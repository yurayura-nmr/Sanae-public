// g++ -std=c++17 -Wall -o test_binding ../binding_utils.cpp test_binding_utils.cpp

#include "../binding_utils.h"
#include <iostream>
#include <cassert>
#include <cmath>

int main() {
    using namespace BindingUtils;

    std::cout << "=== Testing BindingUtils ===\n";

    // ------------------------------------------------------------------
    // 1. Golden-ratio special case: P0 = L0 = Kd = x
    //    sum = 3x, discriminant = 5x^2, PL = x(3 - sqrt(5))/2
    //    fraction = (3 - sqrt(5)) / 2  (independent of x)
    // ------------------------------------------------------------------
    {
        double x = 1.0;
        double frac = CalculateFractionBoundExact(x, x, x);
        double expected = (3.0 - std::sqrt(5.0)) / 2.0;
        std::cout << "P0=L0=Kd=1 -> fraction = " << frac
                   << " (expected " << expected << ")\n";
        assert(std::abs(frac - expected) < 1e-9);
    }

    // Same relation should hold regardless of the common scale x
    {
        double x = 7.5;
        double frac = CalculateFractionBoundExact(x, x, x);
        double expected = (3.0 - std::sqrt(5.0)) / 2.0;
        assert(std::abs(frac - expected) < 1e-9);
        std::cout << "Scale invariance check (x=7.5): PASSED\n";
    }

    // ------------------------------------------------------------------
    // 2. Asymmetry check: swapping P0 and L0 changes the *fraction*,
    //    even though [PL] itself is symmetric in P0/L0 (per the docstring).
    // ------------------------------------------------------------------
    {
        double P0 = 2.0, L0 = 8.0, kd = 1.0;

        double sum = P0 + L0 + kd;
        double disc = sum * sum - 4.0 * P0 * L0;
        double PL = (sum - std::sqrt(disc)) / 2.0; // symmetric in P0/L0

        double frac_normal = CalculateFractionBoundExact(P0, L0, kd);
        double frac_swapped = CalculateFractionBoundExact(L0, P0, kd);

        assert(std::abs(frac_normal - PL / P0) < 1e-9);
        assert(std::abs(frac_swapped - PL / L0) < 1e-9);
        assert(std::abs(frac_normal - frac_swapped) > 1e-6); // must differ

        std::cout << "Asymmetry check: normal=" << frac_normal
                   << " swapped=" << frac_swapped << " -> PASSED\n";
    }

    // ------------------------------------------------------------------
    // 3. Tight-binding limit, ligand in excess: Kd -> 0, L0 >= P0
    //    Nearly all protein should be bound (fraction -> 1).
    // ------------------------------------------------------------------
    {
        double frac = CalculateFractionBoundExact(1.0, 1.0, 1e-9);
        std::cout << "Tight-binding, L0>=P0: fraction = " << frac
                   << " (expect ~1)\n";
        assert(frac > 0.999);
    }

    // ------------------------------------------------------------------
    // 4. Tight-binding limit, ligand sub-stoichiometric: Kd -> 0, L0 < P0
    //    Ligand is the limiting reagent, so fraction -> L0/P0, not 1.
    // ------------------------------------------------------------------
    {
        double P0 = 1.0, L0 = 0.3;
        double frac = CalculateFractionBoundExact(P0, L0, 1e-9);
        std::cout << "Tight-binding, L0<P0: fraction = " << frac
                   << " (expect ~" << (L0 / P0) << ")\n";
        assert(std::abs(frac - L0 / P0) < 1e-6);
    }

    // ------------------------------------------------------------------
    // 5. Large ligand excess: exact solution should converge to the
    //    simple hyperbolic approximation L0/(L0+Kd).
    // ------------------------------------------------------------------
    {
        double P0 = 0.001, L0 = 100.0, kd = 5.0;
        double frac_exact = CalculateFractionBoundExact(P0, L0, kd);
        double frac_approx = L0 / (L0 + kd);
        std::cout << "Excess-ligand check: exact=" << frac_exact
                   << " approx=" << frac_approx << "\n";
        assert(std::abs(frac_exact - frac_approx) < 1e-3);
    }

    // ------------------------------------------------------------------
    // 6. Input validation for CalculateFractionBoundExact
    // ------------------------------------------------------------------
    {
        bool threw;

        threw = false;
        try { CalculateFractionBoundExact(1.0, 1.0, 0.0); }
        catch (const std::invalid_argument&) { threw = true; }
        assert(threw);

        threw = false;
        try { CalculateFractionBoundExact(0.0, 1.0, 1.0); }
        catch (const std::invalid_argument&) { threw = true; }
        assert(threw);

        threw = false;
        try { CalculateFractionBoundExact(1.0, -2.0, 1.0); }
        catch (const std::invalid_argument&) { threw = true; }
        assert(threw);

        std::cout << "CalculateFractionBoundExact validation: PASSED\n";
    }

    // ------------------------------------------------------------------
    // 7. SingleSiteBinding: just verify it runs without throwing and
    //    also enforces the same validation as the exact calculation.
    // ------------------------------------------------------------------
    {
        SingleSiteBinding(1.0, 1.0, 1.0); // prints report, should not throw

        bool threw = false;
        try { SingleSiteBinding(1.0, 1.0, -1.0); }
        catch (const std::invalid_argument&) { threw = true; }
        assert(threw);

        std::cout << "SingleSiteBinding smoke test: PASSED\n";
    }

    // ------------------------------------------------------------------
    // 8. TwoSiteIndependent
    //    NOTE: protein_total does not appear in the formula at all -
    //    only ligand_total, kd1, kd2 affect the result. Confirm that
    //    behavior explicitly so a future refactor doesn't silently
    //    change it without a test noticing.
    // ------------------------------------------------------------------
    {
        double L0 = 5.0, kd1 = 2.0, kd2 = 8.0;
        double theta1 = L0 / (kd1 + L0);
        double theta2 = L0 / (kd2 + L0);
        double expected = 0.5 * (theta1 + theta2);

        double result_P1 = TwoSiteIndependent(1.0, L0, kd1, kd2);
        double result_P100 = TwoSiteIndependent(100.0, L0, kd1, kd2);

        assert(std::abs(result_P1 - expected) < 1e-9);
        assert(std::abs(result_P1 - result_P100) < 1e-12); // protein_total ignored

        std::cout << "TwoSiteIndependent formula check: PASSED (" << expected << ")\n";
    }

    // Degenerate case: kd1 == kd2 should reduce to single-site theta
    {
        double L0 = 3.0, kd = 4.0;
        double result = TwoSiteIndependent(1.0, L0, kd, kd);
        double expected = L0 / (kd + L0);
        assert(std::abs(result - expected) < 1e-9);
        std::cout << "TwoSiteIndependent degenerate (kd1=kd2) check: PASSED\n";
    }

    // Validation
    {
        bool threw = false;
        try { TwoSiteIndependent(1.0, 1.0, 0.0, 1.0); }
        catch (const std::invalid_argument&) { threw = true; }
        assert(threw);

        threw = false;
        try { TwoSiteIndependent(1.0, 0.0, 1.0, 1.0); }
        catch (const std::invalid_argument&) { threw = true; }
        assert(threw);

        std::cout << "TwoSiteIndependent validation: PASSED\n";
    }

    // ------------------------------------------------------------------
    // 9. ThreeSiteIndependent - same pattern, three sites
    // ------------------------------------------------------------------
    {
        double L0 = 6.0, kd1 = 1.0, kd2 = 3.0, kd3 = 9.0;
        double theta1 = L0 / (kd1 + L0);
        double theta2 = L0 / (kd2 + L0);
        double theta3 = L0 / (kd3 + L0);
        double expected = (theta1 + theta2 + theta3) / 3.0;

        double result = ThreeSiteIndependent(1.0, L0, kd1, kd2, kd3);
        assert(std::abs(result - expected) < 1e-9);
        std::cout << "ThreeSiteIndependent formula check: PASSED (" << expected << ")\n";
    }

    // Degenerate: kd1==kd2==kd3 should reduce to single-site theta
    {
        double L0 = 2.0, kd = 5.0;
        double result = ThreeSiteIndependent(1.0, L0, kd, kd, kd);
        double expected = L0 / (kd + L0);
        assert(std::abs(result - expected) < 1e-9);
        std::cout << "ThreeSiteIndependent degenerate check: PASSED\n";
    }

    // Validation
    {
        bool threw = false;
        try { ThreeSiteIndependent(1.0, 1.0, 1.0, -1.0, 1.0); }
        catch (const std::invalid_argument&) { threw = true; }
        assert(threw);
        std::cout << "ThreeSiteIndependent validation: PASSED\n";
    }

    // ------------------------------------------------------------------
    // 10. CalculateOffRateFromKex
    //     Reuse the P0=L0=Kd=1 golden-ratio case to get PL and L_free
    //     analytically, then check koff = kex / (1 + L_free/Kd).
    // ------------------------------------------------------------------
    {
        double P0 = 1.0, L0 = 1.0, kd = 1.0, kex = 50.0;

        double sum = P0 + L0 + kd;
        double disc = sum * sum - 4.0 * P0 * L0;
        double PL = (sum - std::sqrt(disc)) / 2.0;
        double L_free = L0 - PL;
        double expected_koff = kex / (1.0 + L_free / kd);

        double koff = CalculateOffRateFromKex(P0, L0, kd, kex);
        std::cout << "CalculateOffRateFromKex: koff = " << koff
                   << " (expected " << expected_koff << ")\n";
        assert(std::abs(koff - expected_koff) < 1e-9);
    }

    // Sanity bound: koff must be <= kex (since 1 + L_free/Kd >= 1)
    {
        double koff = CalculateOffRateFromKex(2.0, 20.0, 3.0, 100.0);
        assert(koff <= 100.0 + 1e-9);
        std::cout << "CalculateOffRateFromKex upper-bound check: PASSED\n";
    }

    // Validation
    {
        bool threw;

        threw = false;
        try { CalculateOffRateFromKex(1.0, 1.0, 0.0, 10.0); }
        catch (const std::invalid_argument&) { threw = true; }
        assert(threw);

        threw = false;
        try { CalculateOffRateFromKex(1.0, 1.0, 1.0, 0.0); }
        catch (const std::invalid_argument&) { threw = true; }
        assert(threw);

        threw = false;
        try { CalculateOffRateFromKex(-1.0, 1.0, 1.0, 10.0); }
        catch (const std::invalid_argument&) { threw = true; }
        assert(threw);

        std::cout << "CalculateOffRateFromKex validation: PASSED\n";
    }

    std::cout << "\nAll tests passed! BindingUtils is solid.\n";
    return 0;
}
