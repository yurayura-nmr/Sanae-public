#pragma once

#include <string>
#include <vector>

/**
 * @file   cpmg_3state.h
 * @brief  3-state CPMG relaxation dispersion via matrix exponential (Eigen).
 *
 * Implements the general matrix method for 3-state exchange (A <-> B, A <-> C),
 * corresponding to the GLOVE model cpmg_matrix3f:
 *
 *   Nature 430, 586-590 (2004)
 *
 * The Liouvillian L = K + i*Omega is propagated through one CPMG echo
 * (free precession + 180 degree refocusing) using Eigen matrix exponentials.
 * R2eff is recovered from the decay of the summed transverse magnetisation.
 *
 * @par Topology (cpmg_matrix3f):
 *   State A (major) exchanges with B and C independently.
 *   No direct B <-> C exchange (kBC = kCB = 0).
 *
 * @par Sign convention and degeneracy:
 *   dwAB and dwAC enter the Liouvillian as independent signed quantities
 *   (imaginary diagonal elements i*dwAB, i*dwAC), so the dispersion
 *   profile IS sensitive to their relative sign.  However, simultaneously
 *   negating both:
 *
 *     { dwAB, dwAC }  ->  { -dwAB, -dwAC }
 *
 *   is equivalent to complex-conjugating L, leaving all R2eff values
 *   unchanged.  This global sign degeneracy cannot be resolved from CPMG
 *   data alone; only the relative sign of dwAB vs dwAC is observable.
 *
 *   References:
 *     Skrynnikov, Dahlquist & Kay, JACS 124, 12352 (2002)  -- sign degeneracy
 *     Korzhnev et al., Nature 430, 586 (2004)              -- 3-state matrix method
 *
 * @todo The nu_CPMG axis convention used internally (nu = N / Trelax) differs
 *       by a factor of 2 from the standard literature definition
 *       nu_CPMG = 1 / (4*tcp) = N / (2*Trelax) (see e.g. Vallurupalli's
 *       "Chemical Exchange" lecture notes, Eq 27, and the analogous
 *       derivation used for zz_exchange.h). The physics at each computed
 *       point (the actual N, tcp, and resulting R2eff) is unaffected --
 *       only the printed/stored nu_CPMG label would shift by 2x if this
 *       is confirmed to be an error rather than an intentional alternate
 *       convention. Deferred pending review; see BuildNAxis() in the .cpp
 *       for where the fix would go if confirmed.
 *
 * @par Typical usage:
 * @code
 *   // Simulate a profile using parameters from a GLOVE fit
 *   Sanae::Matrix3State(
 *       1.2, -0.8,         // dwAB, dwAC (ppm, signed, from GLOVE)
 *       120.0, 1800.0,     // kAB, kBA (s^-1)
 *       60.0,  900.0,      // kAC, kCA (s^-1)
 *       15.0, 600.0, 0.101,
 *       "dV70-HN distal Ub");
 *
 *   // Confirm global flip degeneracy and that relative flip differs
 *   Sanae::Matrix3StateSignDemo(
 *       1.2, -0.8,
 *       120.0, 1800.0, 60.0, 900.0,
 *       15.0, 600.0, 0.101);
 * @endcode
 */
namespace Sanae
{
	/**
	 * @brief Result of a Matrix3State simulation: the CPMG frequency axis
	 *        and the corresponding R2eff values.
	 *
	 * nu and R2eff are the same length and index-aligned: nu[i], R2eff[i]
	 * refer to the same dispersion point. Empty vectors indicate the
	 * simulation could not be run (see Matrix3State's @note on rate
	 * validity).
	 */
	struct Matrix3StateResult
	{
		std::vector<double> nu;
		std::vector<double> R2eff;
	};

	/**
	 * @brief Result of a Matrix3StateSignDemo run: the CPMG frequency axis
	 *        and the three R2eff curves (original, global flip, relative
	 *        flip), index-aligned with nu.
	 *
	 * Empty vectors indicate the simulation could not be run (see
	 * Matrix3StateSignDemo's @note on rate validity).
	 */
	struct Matrix3StateSignDemoResult
	{
		std::vector<double> nu;
		std::vector<double> R2_orig;
		std::vector<double> R2_gflip;
		std::vector<double> R2_rflip;
	};

	/**
	 * @brief Compute and print a CPMG relaxation dispersion profile
	 *        using the 3-state matrix method (cpmg_matrix3f topology).
	 *
	 * Equivalent to running GLOVE's cpmg_matrix3f model but as a forward
	 * simulator.  Parameters map directly onto GLOVE output, so fitted
	 * values can be plugged in without conversion.
	 *
	 * nu_CPMG axis is log-spaced, matching CarverRichards() convention.
	 *
	 * @param dwAB_ppm  Chemical shift difference Omega_B - Omega_A   (ppm, signed)
	 * @param dwAC_ppm  Chemical shift difference Omega_C - Omega_A   (ppm, signed)
	 * @param kAB       A->B exchange rate                      (s^-1)
	 * @param kBA       B->A exchange rate                      (s^-1)
	 * @param kAC       A->C exchange rate                      (s^-1)
	 * @param kCA       C->A exchange rate                      (s^-1)
	 * @param R20       Intrinsic transverse relaxation rate   (s^-1)
	 * @param B0_MHz    Spectrometer proton frequency          (MHz)
	 * @param gamma     Gyromagnetic ratio relative to 1H
	 *                  (e.g. 0.101 for 15N, 0.251 for 13C, 1.0 for 1H)
	 * @param nu_min    Minimum CPMG frequency 1/tcp           (Hz, default: 25.0)
	 * @param nu_max    Maximum CPMG frequency 1/tcp           (Hz, default: 2000.0)
	 * @param n_points  Number of log-spaced dispersion points (default: 40)
	 * @param label     Optional description printed in header (default: "")
	 * @param filename  If non-empty, write tab-delimited data to this path
	 *
	 * @return
	 *     The simulated dispersion curve (nu, R2eff), index-aligned. See
	 *     Matrix3StateResult.
	 *
	 * @note
	 *     Populations pA, pB, pC are derived from kAB, kBA, kAC, kCA via
	 *     detailed balance, which requires kBA > 0 (if kAB > 0) and
	 *     kCA > 0 (if kAC > 0) to be well-defined. Any rate passed as
	 *     negative, or an irreversible one-way rate (e.g. kAB > 0 with
	 *     kBA <= 0), is invalid: the function prints an error and returns
	 *     an empty Matrix3StateResult (all vectors empty) without writing
	 *     a file, rather than silently computing an incorrect population.
	 */
	Matrix3StateResult Matrix3State(
		double dwAB_ppm,
		double dwAC_ppm,
		double kAB,
		double kBA,
		double kAC,
		double kCA,
		double R20,
		double B0_MHz,
		double gamma,
		double nu_min = 25.0,
		double nu_max = 2000.0,
		int n_points = 40,
		const std::string &label = "",
		const std::string &filename = "");

	/**
	 * @brief Demonstrate the global sign degeneracy of 3-state CPMG data.
	 *
	 * Computes three dispersion profiles and prints them side-by-side:
	 *
	 *   (1) Original:      { +dwAB, +dwAC }
	 *   (2) Global flip:   { -dwAB, -dwAC }  -- numerically identical to (1)
	 *   (3) Relative flip: { -dwAB, +dwAC }  -- differs from (1) if dwAB != dwAC
	 *
	 * The near-zero difference column (orig - gflip) confirms that the global
	 * sign cannot be recovered from CPMG data alone, while the difference
	 * between (1) and (3) confirms that the relative sign IS observable.
	 *
	 * Intended for inclusion as a supplementary figure or rebuttal evidence.
	 *
	 * @param dwAB_ppm  Chemical shift difference Omega_B - Omega_A   (ppm, signed)
	 * @param dwAC_ppm  Chemical shift difference Omega_C - Omega_A   (ppm, signed)
	 * @param kAB       A->B exchange rate                      (s^-1)
	 * @param kBA       B->A exchange rate                      (s^-1)
	 * @param kAC       A->C exchange rate                      (s^-1)
	 * @param kCA       C->A exchange rate                      (s^-1)
	 * @param R20       Intrinsic transverse relaxation rate   (s^-1)
	 * @param B0_MHz    Spectrometer proton frequency          (MHz)
	 * @param gamma     Gyromagnetic ratio relative to 1H
	 * @param nu_min    Minimum CPMG frequency                 (Hz, default: 25.0)
	 * @param nu_max    Maximum CPMG frequency                 (Hz, default: 2000.0)
	 * @param n_points  Number of dispersion points            (default: 40)
	 * @param filename  If non-empty, write columns
	 *                  nu | R2_orig | R2_gflip | R2_rflip to this path
	 *
	 * @return
	 *     The simulated nu axis and all three R2eff curves, index-aligned.
	 *     See Matrix3StateSignDemoResult.
	 *
	 * @note
	 *     Same rate-validity requirements as Matrix3State; see its @note.
	 *     Returns an empty Matrix3StateSignDemoResult if rates are invalid.
	 */
	Matrix3StateSignDemoResult Matrix3StateSignDemo(
		double dwAB_ppm,
		double dwAC_ppm,
		double kAB,
		double kBA,
		double kAC,
		double kCA,
		double R20,
		double B0_MHz,
		double gamma,
		double nu_min = 25.0,
		double nu_max = 2000.0,
		int n_points = 40,
		const std::string &filename = "");

} // namespace Sanae