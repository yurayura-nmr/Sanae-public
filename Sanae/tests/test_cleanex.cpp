// g++ -std=c++17 -Wall -Wextra -I../../include/eigen3/ -o test_cleanex ../cleanex_utils.cpp ../fileio.cpp test_cleanex.cpp

#include "../cleanex_utils.h"
#include <iostream>
#include <fstream>
#include <sstream>
#include <cassert>
#include <cmath>
#include <cstdio>
#include <vector>

namespace
{

    struct ParsedCurve
    {
        std::vector<double> tm_ms;
        std::vector<double> I;
    };

    // Reads back a file written by SimulateCleanex, skipping '#' header lines.
    ParsedCurve ReadCurveFile(const std::string &filename)
    {
        ParsedCurve curve;
        std::ifstream in(filename);
        assert(in.is_open() && "output file should exist and be readable");

        std::string line;
        while (std::getline(in, line))
        {
            if (line.empty() || line[0] == '#')
                continue;
            std::istringstream iss(line);
            double tm, I;
            if (iss >> tm >> I)
            {
                curve.tm_ms.push_back(tm);
                curve.I.push_back(I);
            }
        }
        return curve;
    }

    double AnalyticI(double tm, double k, double R1a, double R1b, double prefactor)
    {
        return prefactor * (-std::exp(-(R1a + k) * tm) + std::exp(-R1b * tm));
    }

} // namespace

int main()
{
    using namespace CleanexUtils;

    std::cout << "=== Testing CleanexUtils::SimulateCleanex ===\n";

    // ------------------------------------------------------------------
    // 1. Basic run: check every point in the output file matches the
    //    analytic formula, and tm is log-spaced from tm_min to tm_max.
    // ------------------------------------------------------------------
    {
        double k = 20.0, R1a = 2.0, R1b = 5.0;
        double tm_min = 0.001, tm_max = 1.0; // seconds
        size_t n_points = 20;
        const std::string fname = "test_cleanex_basic.dat";

        SimulateCleanex(k, R1a, R1b, tm_min, tm_max, n_points,
                        "Basic test case", fname);

        double prefactor = k / (R1a - R1b + k);
        ParsedCurve curve = ReadCurveFile(fname);

        assert(curve.tm_ms.size() == n_points);
        assert(curve.I.size() == n_points);

        // First and last tm should match tm_min/tm_max (in ms)
        assert(std::abs(curve.tm_ms.front() - tm_min * 1000.0) < 1e-6);
        assert(std::abs(curve.tm_ms.back() - tm_max * 1000.0) < 1e-6);

        // Every point should match the analytic curve
        for (size_t i = 0; i < curve.tm_ms.size(); ++i)
        {
            double tm_s = curve.tm_ms[i] / 1000.0;
            double expected = AnalyticI(tm_s, k, R1a, R1b, prefactor);
            assert(std::abs(curve.I[i] - expected) < 1e-6);
        }

        // Log-spacing check: ratio of consecutive tm values should be constant
        double ratio_expected = std::pow(tm_max / tm_min, 1.0 / (n_points - 1));
        for (size_t i = 1; i < curve.tm_ms.size(); ++i)
        {
            double ratio = curve.tm_ms[i] / curve.tm_ms[i - 1];
            assert(std::abs(ratio - ratio_expected) < 1e-6);
        }

        std::cout << "Basic run + analytic match + log-spacing: PASSED\n";
        std::remove(fname.c_str());
    }

    // ------------------------------------------------------------------
    // 2. n_points == 1 edge case: must not produce NaN (previously would,
    //    via division by (n_points - 1) == 0).
    // ------------------------------------------------------------------
    {
        const std::string fname = "test_cleanex_single.dat";
        SimulateCleanex(15.0, 2.0, 5.0, 0.01, 0.2, 1, "Single point", fname);

        ParsedCurve curve = ReadCurveFile(fname);
        assert(curve.tm_ms.size() == 1);
        assert(!std::isnan(curve.tm_ms[0]));
        assert(!std::isnan(curve.I[0]));
        // With frac defined as 0.0 for n_points==1, tm should equal tm_min
        assert(std::abs(curve.tm_ms[0] - 0.01 * 1000.0) < 1e-6);

        std::cout << "n_points == 1 edge case (no NaN): PASSED\n";
        std::remove(fname.c_str());
    }

    // ------------------------------------------------------------------
    // 3. n_points == 0: should simply return without creating a file.
    // ------------------------------------------------------------------
    {
        const std::string fname = "test_cleanex_zero.dat";
        std::remove(fname.c_str()); // ensure clean slate
        SimulateCleanex(15.0, 2.0, 5.0, 0.01, 0.2, 0, "Zero points", fname);

        std::ifstream in(fname);
        assert(!in.is_open() && "no file should be written for n_points == 0");
        std::cout << "n_points == 0 guard: PASSED\n";
    }

    // ------------------------------------------------------------------
    // 4. Empty filename: should not write any file, but should still run
    //    (prints to stdout) without crashing.
    // ------------------------------------------------------------------
    {
        SimulateCleanex(15.0, 2.0, 5.0, 0.01, 0.2, 5, "No file requested", "");
        std::cout << "Empty filename (stdout-only) run: PASSED\n";
    }

    // ------------------------------------------------------------------
    // 5. Degenerate case: R1a - R1b + k ≈ 0. Should not crash, and should
    //    fall back to the TINY-clamped denominator rather than blow up.
    // ------------------------------------------------------------------
    {
        double R1a = 5.0, k = 3.0;
        double R1b = R1a + k; // forces denom == 0 exactly
        const std::string fname = "test_cleanex_degenerate.dat";

        SimulateCleanex(k, R1a, R1b, 0.001, 0.1, 10, "Degenerate case", fname);

        ParsedCurve curve = ReadCurveFile(fname);
        for (double val : curve.I)
        {
            assert(std::isfinite(val));
        }
        std::cout << "Degenerate denominator case (finite output): PASSED\n";
        std::remove(fname.c_str());
    }

    // ------------------------------------------------------------------
    // 6. Prefactor sanity check across a non-degenerate case, read from
    //    the file header comment.
    // ------------------------------------------------------------------
    {
        double k = 10.0, R1a = 1.0, R1b = 4.0;
        const std::string fname = "test_cleanex_prefactor.dat";
        SimulateCleanex(k, R1a, R1b, 0.005, 0.3, 8, "Prefactor check", fname);

        std::ifstream in(fname);
        std::string line;
        double file_prefactor = std::numeric_limits<double>::quiet_NaN();
        while (std::getline(in, line))
        {
            if (line.rfind("# prefactor", 0) == 0)
            {
                std::istringstream iss(line);
                std::string hash, word, eq;
                iss >> hash >> word >> eq >> file_prefactor;
                break;
            }
        }
        double expected_prefactor = k / (R1a - R1b + k);
        assert(std::abs(file_prefactor - expected_prefactor) < 1e-9);
        std::cout << "Prefactor header check: PASSED (" << file_prefactor << ")\n";
        std::remove(fname.c_str());
    }

    std::cout << "\nAll tests passed! CleanexUtils is solid.\n";
    return 0;
}
