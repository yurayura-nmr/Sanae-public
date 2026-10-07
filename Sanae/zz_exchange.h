#pragma once

#include <string>
#include <vector>

/**
 * @file zz_exchange.h
 * @brief ZZ exchange experiment simulation for slow conformational exchange.
 *
 * Simulates peak intensities as a function of mixing time T for the
 * ZZ (longitudinal) exchange experiment (Farrow et al. 1994).
 *
 * Valid exchange regime: kex = 0.1 to 10 s^-1 (slow exchange)
 * Requirement: both peaks visible -> pB >= 0.1
 * Requirement: kex > R1 for accurate determination
 *
 * The longitudinal Bloch-McConnell equations during mixing time T:
 *
 *   d/dt [dMz_A]  =  [-R1A - k_ab    k_ba    ] [dMz_A]
 *        [dMz_B]     [  k_ab      -R1B - k_ba] [dMz_B]
 *
 * Four observable intensities vs T:
 *   a11 = autopeak  A    (wN_A, wH_A)  starts at 1, decays
 *   a22 = autopeak  B    (wN_B, wH_B)  starts at 1, decays
 *   a12 = crosspeak B->A (wN_B, wH_A)  starts at 0, rises then decays
 *   a21 = crosspeak A->B (wN_A, wH_B)  starts at 0, rises then decays
 *
 * When R1A = R1B = R1 (simplified case):
 *   a11 = [pA + pB*exp(-kex * T)] * exp(-R1 * T)
 *   a22 = [pB + pA*exp(-kex * T)] * exp(-R1 * T)
 *   a12 = pA * [1 - exp(-kex * T)] * exp(-R1 * T)
 *   a21 = pB * [1 - exp(-kex * T)] * exp(-R1 * T)
 *
 * @note DIVERGENCE FROM SOURCE LECTURE NOTES:
 *       Pramodh Vallurupalli's lecture notes ("Chemical Exchange", TIFR
 *       Mumbai NMR Workshop, 2009; Baldwin lab), Section 5, Eq 26, give
 *       these same expressions with exp(-2*kex*T) instead of exp(-kex*T).
 *       We believe this is a transcription error in the notes: the
 *       underlying rate matrix they themselves give in Eq 24,
 *           d/dt [DMzA; DMzB] = [[-k1, k-1], [k1, -k-1]] [DMzA; DMzB],
 *       has eigenvalues 0 and -(k1+k-1) = -kex (standard result for
 *       reversible first-order two-state kinetics; consistent with their
 *       own Eq 7 for the analogous concentration-relaxation case). Direct
 *       substitution of the exp(-kex*T) solution back into Eq 24 confirms
 *       it satisfies that equation exactly; exp(-2*kex*T) does not. This
 *       is also confirmed numerically: integrating Eq 24 (as this file
 *       does, via Runge-Kutta, independent of either closed form) agrees
 *       with the exp(-kex*T) expressions above, not exp(-2*kex*T). See
 *       test_zz_exchange.cpp for the numerical cross-check.
 *
 * Reference: Farrow, Zhang, Forman-Kay, Kay (1994) J Biomol NMR 4, 727
 *            Vallurupalli, "Chemical Exchange" lecture notes, Section 5,
 *            Eq 24 (rate matrix; used as-is) and Eq 26 (closed-form
 *            solution; corrected here per the note above).
 *
 * @note This is NOT an FID simulation. Output is peak intensity vs T,
 *       not a time-domain signal. Use nmr2stateExchange for lineshapes.
 * 
 * @todo SimulateZZExchange currently prints its full banner/table output
 *       unconditionally via printf/std::cout, with no way to silence it
 *       (unlike Sanae::Solution, which has a verbose_ flag/SetVerbose()).
 *       Consider adding a similar verbose parameter (default true, to
 *       preserve current behavior) so callers/tests can opt out of the
 *       console output without needing to redirect stdout. Same pattern
 *       as the verbose parameter added to R1_inversion_recovery.
 */

namespace Sanae
{
	/**
	 * @brief Result of a ZZ exchange simulation: mixing times and the
	 *        four observable peak intensities at each mixing time.
	 *
	 * All vectors are the same length (n_points) and index-aligned:
	 * T[i], a11[i], a22[i], a12[i], a21[i] all refer to the same
	 * mixing time point. Empty vectors indicate the simulation could
	 * not be run (e.g. kex = 0 -- see SimulateZZExchange's @note).
	 */
	struct ZZExchangeResult
	{
		std::vector<double> T;
		std::vector<double> a11; // autopeak A
		std::vector<double> a22; // autopeak B
		std::vector<double> a12; // crosspeak B->A
		std::vector<double> a21; // crosspeak A->B
	};

	/**
	 * @brief Simulate ZZ exchange peak intensities vs mixing time.
	 *
	 * Uses Runge-Kutta integration of the longitudinal Bloch-McConnell
	 * equations — consistent with SANAE philosophy of numerical integration
	 * rather than analytical solutions, even though analytical forms exist
	 * for the R1A = R1B case.
	 *
	 * @param k_ab      A -> B exchange rate        (s^-1)
	 * @param k_ba      B -> A exchange rate        (s^-1)
	 * @param R1_A      Longitudinal relaxation, A  (s^-1)
	 * @param R1_B      Longitudinal relaxation, B  (s^-1)
	 * @param T_min     Minimum mixing time         (s)
	 * @param T_max     Maximum mixing time         (s)
	 * @param n_points  Number of mixing time points
	 * @param label     Label for output
	 * @param filename  Output file (T vs a11 a22 a12 a21). Pass an
	 *                   empty string to skip writing a file.
	 *
	 * @return
	 *     The simulated curves (T, a11, a22, a12, a21), index-aligned,
	 *     for direct use by callers/tests without needing to read them
	 *     back from filename. See ZZExchangeResult.
	 *
	 * @note If k_ab + k_ba (kex) is zero or negative, the populations
	 *       pA = k_ba/kex and pB = k_ab/kex are undefined (division by
	 *       zero). In that case the function prints an error and
	 *       returns an empty ZZExchangeResult (all vectors empty)
	 *       without writing a file, rather than propagating NaN.
	 */
	ZZExchangeResult SimulateZZExchange(
		double k_ab,
		double k_ba,
		double R1_A,
		double R1_B,
		double T_min = 0.01,
		double T_max = 1.0,
		int n_points = 50,
		const std::string &label = "ZZ exchange simulation",
		const std::string &filename = "zz_exchange.txt");
}