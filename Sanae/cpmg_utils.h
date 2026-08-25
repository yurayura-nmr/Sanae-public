/**
 * @file   cpmg_utils.h
 * @brief  Utility functions for NMR CPMG relaxation dispersion calculations.
 *
 * Implements the Carver & Richards equation for 2-state conformational exchange:
 *
 *   R2eff = R20 + (kex/2) * { 1 - (2/tcp·kex) * acosh[ D+ cosh(η+) - D- cosh(η-) ] }
 *
 * Reference: Carver & Richards, J Magn Reson 6, 89-105 (1972).
 *
 * @param gamma  Gyromagnetic ratio RELATIVE TO 1H (dimensionless, e.g. ~0.101
 *               for 15N). B0_MHz is the PROTON spectrometer frequency.
 *               NOTE: this differs from Environment.cpp's gamma convention
 *               (absolute MHz/T) — do not confuse the two.
 * @par Typical usage in Sanae.cpp:
 * @code
 *   // 15N at 700 MHz, fast exchange
 *   DispersionUtils::CarverRichards(0.95, 5000.0, 1.0, 700.0, 0.101, 10.0);
 * @endcode
 *
 * @todo Make compliant for any nucleus by simply passing char H, N, C, F instead of gamma relative to 1H.
 */

#pragma once

#include <string>
#include <vector>

namespace DispersionUtils
{
	/**
	 * @brief Compute and print a CPMG relaxation dispersion profile
	 *        using the Carver-Richards equation (2-state exchange).
	 *
	 * Output is printed to stdout as a formatted table. The nu_CPMG axis
	 * is log-spaced for better coverage of the dispersion. Note that
	 * nu_CPMG (Hz) == 1/tcp (s⁻¹) numerically — this is the same x-axis
	 * as used by GLOVE.
	 *
	 * @param pA        Population of major state A            (unitless, e.g. 0.95)
	 * @param kex       Total exchange rate kAB + kBA          (s⁻¹)
	 * @param dw_ppm    Chemical shift difference |δA - δB|    (ppm)
	 * @param B0_MHz    Spectrometer proton frequency          (MHz, e.g. 700.0)
	 * @param gamma     Gyromagnetic ratio relative to 1H
	 *                  (e.g. 0.101 for 15N, 0.251 for 13C, 1.0 for 1H, XX for 19F)
	 * @param R20       Intrinsic transverse relaxation rate   (s⁻¹)
	 * @param nu_min    Minimum CPMG frequency 1/tcp           (Hz, default: 25.0)
	 * @param nu_max    Maximum CPMG frequency 1/tcp           (Hz, default: 2000.0)
	 * @param n_points  Number of log-spaced dispersion points (default: 40)
	 * @param label     Optional description printed in header (default: "")
	 */

	void CarverRichards(double pA, double kex, double dw_ppm, double B0_MHz,
						double gamma, double R20, double nu_min, double nu_max,
						size_t n_points, const std::string &label,
						const std::string &filename = "");

	/**
	 * @brief Simplified Meiboom equation for CPMG relaxation dispersion.
	 *
	 * This is an alternative (less accurate) implementation of the
	 * exchange, valid in the fast exchange limit
	 * or when pA ≈ 1. It is kept for legacy/comparison purposes.
	 *
	 * @param pA        Population of major state A (unitless)
	 * @param kex       Total exchange rate (s⁻¹)
	 * @param dw_ppm    Chemical shift difference (ppm)
	 * @param B0_MHz    Spectrometer proton frequency (MHz)
	 * @param gamma     Gyromagnetic ratio relative to 1H
	 * @param R20       Intrinsic transverse relaxation rate (s⁻¹)
	 * @param nu_min    Minimum CPMG frequency (Hz)
	 * @param nu_max    Maximum CPMG frequency (Hz)
	 * @param n_points  Number of log-spaced dispersion points
	 * @param label     Optional description printed in header
	 */

	void Meiboom(double pA, double kex, double dw_ppm, double B0_MHz,
				 double gamma, double R20, double nu_min, double nu_max,
				 size_t n_points, const std::string &label,
				 const std::string &filename = "");

	/**
	 * @brief Compute Rex and field-dependence parameter alpha at multiple B0 fields.
	 *
	 * Uses the unified Rex expression (valid when pA > 0.7, R2A ≈ R2B):
	 *
	 *   Rex = pA*pB*kex / (1 + (kex/dw)²)
	 *
	 * where dw = dw_ppm * gamma * B0 * 2pi (rad/s), so dw ∝ B0.
	 *
	 * The field-dependence parameter alpha (Palmer et al. 2001):
	 *   alpha = d(ln Rex) / d(ln dw)
	 *
	 * Exchange regime from alpha:
	 *   0 ≤ alpha < 1 : slow exchange   — Rex independent of B0
	 *   alpha = 1     : intermediate    — Rex linear in B0
	 *   1 < alpha ≤ 2 : fast exchange   — Rex ∝ B0²
	 *
	 * Limits:
	 *   Slow  (kex/dw → 0): Rex → pB*kex        (B0-independent)
	 *   Fast  (kex/dw → ∞): Rex → pA*pB*dw²/kex (Rex ∝ B0²)
	 *
	 * Reference: Palmer et al. (2001) Chem Rev 101, 3457
	 *            Lecture notes Section 4.3
	 *
	 * @param kex       Total exchange rate              (s⁻¹)
	 * @param pA        Population of major state        (> 0.7 for approximation)
	 * @param dw_ppm    Chemical shift difference        (ppm)
	 * @param gamma     Gyromagnetic ratio factor        (e.g. 0.101 for 15N)
	 * @param B0_list   Vector of B0 fields to evaluate  (MHz, 1H frequency)
	 * @param label     Label for output
	 * @param filename  Output file
	 */

	void RexVsField(
		double kex,
		double pA,
		double dw_ppm,
		double gamma,
		const std::vector<double> &B0_list,
		const std::string &label = "Rex vs B0 field dependence",
		const std::string &filename = "rex_vs_field.txt");

	/**
	 * @brief Simulate R1rho relaxation dispersion profile.
	 *
	 * Implements Eq 32-33 (Vallurupalli lecture notes / Palmer 2001):
	 *
	 *   R1rho = R1*cos²θ + R2*sin²θ                              [Eq 32]
	 *
	 *   with Rex contribution (fast exchange limit, omega1 >> dw):
	 *
	 *   R1rho = R1rho(weff→∞) + sin²θ * pA*pB*dw²*kex / (kex² + weff²)  [Eq 33]
	 *
	 * where:
	 *   tan θ   = omega1 / dOmega
	 *   dOmega  = pA*wA + pB*wB - wRF   (offset of pop-weighted CS from carrier)
	 *   weff    = sqrt(omega1² + dOmega²)
	 *
	 * Valid regime: kex = 5,000 – 50,000 s⁻¹ (faster than CPMG)
	 * Assumption:  fast exchange limit, weff same for states A and B (omega1 >> dw)
	 *
	 * Two scan modes:
	 *   (1) Vary omega1 at fixed offset dOmega  → R1rho vs spin-lock field strength
	 *   (2) Vary dOmega at fixed omega1         → off-resonance R1rho profile
	 *
	 * Reference: Palmer, Kroenke, Loria (2001) Methods Enzymol 339, 204
	 *            Vallurupalli lecture notes Eq 32-33
	 * @param scan_mode  "omega1"  : vary spin-lock field at fixed dOmega (default)
	 *                   "offset"  : vary dOmega offset at fixed omega1
	 *                               → off-resonance R1rho profile (Figure 13)
	 */

	void R1rhoDispersion(
		double kex,
		double pA,
		double dw_ppm,
		double b0_mhz,
		double gamma,
		double R1,
		double R2,
		double omega1_min_hz,
		double omega1_max_hz,
		double dOmega_hz = 0.0,
		size_t n_points = 40,
		const std::string &scan_mode = "omega1", // "omega1" or "offset"
		double dOmega_min_hz = -500.0,			 // used when scan_mode="offset"
		double dOmega_max_hz = 500.0,			 // used when scan_mode="offset"
		const std::string &label = "R1rho dispersion",
		const std::string &filename = "r1rho_dispersion.txt");

} // namespace DispersionUtils
