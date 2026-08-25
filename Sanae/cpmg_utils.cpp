/**
 * @file   cpmg_utils.cpp
 * @brief  Implementation of CPMG relaxation dispersion utilities.
 *
 * Implements the Carver & Richards (1972) equation for 2-state exchange.
 * The internal calcR2eff() function mirrors the calc_R2() routine from
 * the GLOVE fitting program (Richards et al.).
 *
 * @note  The x-axis convention follows GLOVE: nu_CPMG in Hz == 1/tcp in s⁻¹.
 *        These are numerically identical; no 2π conversion is needed.
 * @todo  Move output from RexVsField and R1rhoDispersion to fileio.cpp
 */

#include "fileio.h"
#include "cpmg_utils.h"

#include <cmath>
#include <cstdio>
#include <string>
#include <vector>
#include <fstream>
#include <iostream>
#include <algorithm>

namespace DispersionUtils
{
	static constexpr double PI = 3.14159265358979323846;

	// ---------------------------------------------------------------------------
	/**
	 * @brief  Compute R2eff for a single CPMG tcp pulse spacing (internal).
	 *
	 * Full Carver-Richards expression:
	 * @code
	 *   psi  = kex² - dw²
	 *   zeta = 2 · dw · kex · (pA - pB)
	 *   norm = sqrt(psi² + zeta²)
	 *   D±   = ½ ( ±1 + (kex² + dw²) / norm )
	 *   η±   = (tcp / 2√2) · sqrt(norm ± psi)        — both real since norm ≥ |psi|
	 *   R2eff = R20 + (kex/2) · { 1 - (2/tcp·kex) · acosh[ D+·cosh(η+) - D-·cosh(η-) ] }
	 * @endcode
	 *
	 * @param tcp   Pulse spacing, 1/nu_CPMG                      (s)
	 * @param kex   Total exchange rate kAB + kBA                 (s⁻¹)
	 * @param dw    Chemical shift difference in angular frequency (rad/s)
	 * @param pA    Population of major state A                   (unitless)
	 * @param R20   Intrinsic transverse relaxation rate          (s⁻¹)
	 * @return      Effective transverse relaxation rate R2eff    (s⁻¹)
	 */

	static double calcR2eff(double tcp, // pulse spacing (s)
							double kex, // total exchange rate (s⁻¹)
							double dw,	// chemical shift difference (rad/s)
							double pA,	// population of major state
							double R20) // intrinsic R2 (s⁻¹)
	{
		const double pB = 1.0 - pA;
		const double dw2 = dw * dw;
		const double kex2 = kex * kex;

		const double pop = 1.0 - 2.0 * pB; // pa - pb
		// or: const double pop = sqrt(1.0 - 4.0 * pA * pB);

		const double psi = kex2 - dw2;
		const double zeta = 2.0 * dw * kex * pop;
		const double zeta2 = zeta * zeta;
		const double norm = std::sqrt(psi * psi + zeta2);

		const double inorm = 1.0 / norm;
		const double pdw2 = (kex2 + dw2) * inorm;
		const double Dp = 0.5 * (1.0 + pdw2);
		const double Dm = 0.5 * (-1.0 + pdw2);

		const double eta_p = tcp * std::sqrt(0.5 * (psi + norm));
		const double eta_m = tcp * std::sqrt(0.5 * (norm - psi));

		double arg = Dp * std::cosh(eta_p) - Dm * std::cos(eta_m);
		if (arg < 1.0)
			arg = 1.0;

		const double i_tcp = 1.0 / tcp;
		return R20 + 0.5 * (kex - i_tcp * std::acosh(arg));
	}

	static double calcR2effMeiboom(double tcp, // pulse spacing (s)
								   double kex, // exchange rate (s^-1)
								   double dw,  // chemical shift difference (rad/s)
								   double pA,  // population of major state
								   double R20) // intrinsic R2
	{
		const double pB = 1.0 - pA;

		// pA*pB*dw^2
		const double pdw = pA * pB * dw * dw;

		const double x = 0.5 * kex * tcp; // kex*tcp/2

		// avoid divide-by-zero for very small x
		double tanh_term;
		if (std::abs(x) < 1e-12)
			tanh_term = 1.0;
		else
			tanh_term = std::tanh(x) / x;

		return R20 + (pdw / kex) * (1.0 - tanh_term);
	}

	// ---------------------------------------------------------------------------
	// Public API
	// ---------------------------------------------------------------------------
	void CarverRichards(double pA,
						double kex,
						double dw_ppm,
						double B0_MHz,
						double gamma,
						double R20,
						double nu_min,
						double nu_max,
						size_t n_points,
						const std::string &label,
						const std::string &filename)
	{
		// Convert ppm → rad/s: dw(rad/s) = dw(ppm) * gamma * B0(MHz) * 2π
		const double dw_rads = dw_ppm * gamma * B0_MHz * 2.0 * PI;

		const double pB = 1.0 - pA;

		// Header
		printf("\n");
		printf("========================================================\n");
		printf("  Carver-Richards CPMG Dispersion Profile\n");
		if (!label.empty())
			printf("  %s\n", label.c_str());
		printf("--------------------------------------------------------\n");
		printf("  pA = %.4f   pB = %.4f\n", pA, pB);
		printf("  kex = %.1f s⁻¹\n", kex);
		printf("  dw = %.3f ppm  (%.2f rad/s)\n", dw_ppm, dw_rads);
		printf("  B0 = %.1f MHz   gamma = %.4f\n", B0_MHz, gamma);
		printf("  R20 = %.2f s⁻¹\n", R20);
		printf("  nu_CPMG range: %.1f – %.1f Hz (%zu points)\n",
			   nu_min, nu_max, n_points);
		printf("--------------------------------------------------------\n");
		// nu_CPMG (Hz) == 1/tcp (s⁻¹) — numerically identical, GLOVE plots the latter
		printf("  %-16s  %s\n", "nu_CPMG (Hz) = 1/tcp (s⁻¹) [GLOVE]", "R2eff (s⁻¹)");
		printf("  %-16s  %s\n", "------------------------------------", "-----------");

		std::vector<double> nu_vec(n_points), R2eff_vec(n_points);

		// Scan over nu_CPMG = 1/tcp
		for (size_t i = 0; i < n_points; ++i)
		{
			double nu = nu_min * std::pow(nu_max / nu_min,
										  static_cast<double>(i) / (n_points - 1));
			double tcp = 1.0 / nu;
			double R2eff = calcR2eff(tcp, kex, dw_rads, pA, R20);

			nu_vec[i] = nu;
			R2eff_vec[i] = R2eff;

			printf("  %-14.2f  %.4f\n", nu, R2eff);
		}

		printf("========================================================\n\n");

		if (!filename.empty() &&
			Sanae::writeCpmgCurveToFile(filename, nu_vec, R2eff_vec, kex, dw_ppm,
										B0_MHz, pA, R20, "Carver-Richards"))
			std::cout << "  [CPMG] Curve written to: " << filename << "\n";
	}

	void Meiboom(double pA,
				 double kex,
				 double dw_ppm,
				 double B0_MHz,
				 double gamma,
				 double R20,
				 double nu_min,
				 double nu_max,
				 size_t n_points,
				 const std::string &label,
				 const std::string &filename)
	{
		const double dw_rads = dw_ppm * gamma * B0_MHz * 2.0 * PI;

		printf("\n");
		printf("========================================================\n");
		printf("  Luz-Meiboom CPMG Dispersion Profile\n");
		if (!label.empty())
			printf("  %s\n", label.c_str());
		printf("--------------------------------------------------------\n");

		std::vector<double> nu_vec(n_points), R2eff_vec(n_points);

		for (size_t i = 0; i < n_points; ++i)
		{
			double nu = nu_min * std::pow(nu_max / nu_min,
										  static_cast<double>(i) / (n_points - 1));

			double tcp = 1.0 / nu;
			double R2eff = calcR2effMeiboom(tcp, kex, dw_rads, pA, R20);

			nu_vec[i] = nu;
			R2eff_vec[i] = R2eff;

			printf("  %-14.2f  %.4f\n", nu, R2eff);
		}

		printf("========================================================\n");

		if (!filename.empty() &&
			Sanae::writeCpmgCurveToFile(filename, nu_vec, R2eff_vec, kex, dw_ppm,
										B0_MHz, pA, R20, "Luz-Meiboom"))
			std::cout << "  [CPMG] Curve written to: " << filename << "\n";
	}

	void RexVsField(
		double kex,
		double pA,
		double dw_ppm,
		double gamma,
		const std::vector<double> &B0_list,
		const std::string &label,
		const std::string &filename)
	{
		const double pB = 1.0 - pA;

		// --- Print header ---
		printf("\n");
		printf("================================================================\n");
		printf("  %s\n", label.c_str());
		printf("  Rex = pA*pB*kex / (1 + (kex/dw)²)  [Palmer et al. 2001]\n");
		printf("  alpha = d(ln Rex)/d(ln dw)  — field dependence parameter\n");
		printf("----------------------------------------------------------------\n");
		printf("  kex    = %.1f s⁻¹\n", kex);
		printf("  pA     = %.3f   pB = %.3f\n", pA, pB);
		printf("  dw     = %.3f ppm\n", dw_ppm);
		printf("  gamma  = %.4f\n", gamma);
		printf("\n");
		printf("  NOTE: approximation valid when pA > 0.7 and R2A ≈ R2B.\n");
		printf("  For pA < 0.7 use Carver-Richards (run_carver_richards).\n");
		printf("----------------------------------------------------------------\n");

		// Regime limits for reference
		// Slow limit: Rex → pB*kex (B0-independent)
		// Fast limit: Rex → pA*pB*dw²/kex (Rex ∝ B0²)
		printf("  Limiting values (B0-independent reference):\n");
		printf("    Slow exchange limit  Rex → pB*kex       = %.2f s⁻¹\n",
			   pB * kex);
		printf("    (Rex approaches this from below as dw decreases / B0 decreases)\n");
		printf("----------------------------------------------------------------\n");

		printf("  %-8s  %-8s  %-10s  %-10s  %-8s  %-10s  %s\n",
			   "B0 (MHz)", "B0 (T)", "dw (rad/s)", "kex/dw",
			   "Rex", "alpha", "Regime");
		printf("  %-8s  %-8s  %-10s  %-10s  %-8s  %-10s  %s\n",
			   "--------", "--------", "----------", "----------",
			   "--------", "----------", "----------");

		// Storage for output and alpha calculation
		std::vector<double> B0_vec, dw_vec, Rex_vec, alpha_vec;

		// Need two adjacent points to compute alpha numerically
		// alpha = d(ln Rex)/d(ln dw) ≈ Δ(ln Rex)/Δ(ln dw)
		// Compute at all fields first, then differentiate

		for (double B0 : B0_list)
		{
			double dw_rads = dw_ppm * gamma * B0 * 2.0 * PI;
			double ratio = kex / dw_rads;
			double Rex = pA * pB * kex / (1.0 + ratio * ratio);

			B0_vec.push_back(B0);
			dw_vec.push_back(dw_rads);
			Rex_vec.push_back(Rex);
		}

		// Compute alpha numerically from adjacent points
		// alpha_i = d(ln Rex)/d(ln dw) at point i
		// Use central differences for interior points, forward/backward at edges
		alpha_vec.resize(B0_vec.size(), 0.0);
		int N = static_cast<int>(B0_vec.size());

		for (int i = 0; i < N; ++i)
		{
			int lo = std::max(0, i - 1);
			int hi = std::min(N - 1, i + 1);
			if (hi == lo)
			{
				alpha_vec[i] = 0.0;
				continue;
			}
			double d_lnRex = std::log(Rex_vec[hi]) - std::log(Rex_vec[lo]);
			double d_lndw = std::log(dw_vec[hi]) - std::log(dw_vec[lo]);
			alpha_vec[i] = (std::abs(d_lndw) > 1e-12) ? d_lnRex / d_lndw : 0.0;
		}

		// Print table
		for (int i = 0; i < N; ++i)
		{
			double B0 = B0_vec[i];
			double B0_T = 2.0 * PI * B0 * 1.0e6 / 267520000.0; // T, from gH
			double dw = dw_vec[i];
			double ratio = kex / dw;
			double Rex = Rex_vec[i];
			double alpha = alpha_vec[i];

			// Regime from alpha
			std::string regime;
			if (alpha < 0.5)
				regime = "Slow";
			else if (alpha < 1.5)
				regime = "Intermediate";
			else
				regime = "Fast";

			printf("  %-8.0f  %-8.2f  %-10.1f  %-10.4f  %-8.2f  %-10.4f  %s\n",
				   B0, B0_T, dw, ratio, Rex, alpha, regime.c_str());
		}

		printf("================================================================\n");
		printf("  Fast limit Rex(B0) ∝ B0²: confirm by alpha → 2 at high field\n");
		printf("  Slow limit Rex(B0) const: confirm by alpha → 0 at low field\n");
		printf("================================================================\n\n");

		// --- Write output file ---
		if (filename.empty())
			return;

		std::ofstream out(filename);
		if (!out.is_open())
		{
			std::cerr << "  [REX] Error: could not open: " << filename << "\n";
			return;
		}

		out << "# SANAE Rex vs B0 field dependence\n";
		out << "# kex   = " << kex << " s⁻¹\n";
		out << "# pA    = " << pA << "\n";
		out << "# dw    = " << dw_ppm << " ppm\n";
		out << "# gamma = " << gamma << "\n";
		out << "# B0_MHz\tdw_rads\tkex_over_dw\tRex\talpha\n";

		for (int i = 0; i < N; ++i)
			out << B0_vec[i] << "\t"
				<< dw_vec[i] << "\t"
				<< kex / dw_vec[i] << "\t"
				<< Rex_vec[i] << "\t"
				<< alpha_vec[i] << "\n";

		std::cout << "  [REX] Field dependence written to: " << filename << "\n";
	}

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
		double dOmega_hz,
		size_t n_points,
		const std::string &scan_mode,
		double dOmega_min_hz,
		double dOmega_max_hz,
		const std::string &label,
		const std::string &filename)
	{
		double pB = 1.0 - pA;

		// Convert chemical shift difference from ppm to rad/s
		double dw_rads = dw_ppm * gamma * b0_mhz * 2.0 * PI;

		// Configuration values (rad/s)
		double dOmega_fixed_rads = dOmega_hz * 2.0 * PI;
		double omega1_fixed_rads = omega1_min_hz * 2.0 * PI;
		bool scan_offset = (scan_mode == "offset");

		// --- Print header ---
		printf("\n");
		printf("================================================================\n");
		printf("  %s\n", label.c_str());
		printf("  Model: Miloushev (General Exchange Regime)\n");
		printf("  Accounts for asymmetric, state-specific tilt angles (theta_A, theta_B)\n");
		printf("----------------------------------------------------------------\n");
		printf("  kex    = %.1f s⁻¹\n", kex);
		printf("  pA     = %.3f   pB = %.3f\n", pA, pB);
		printf("  dw     = %.3f ppm  (%.2f rad/s)\n", dw_ppm, dw_rads);
		printf("  R1     = %.3f s⁻¹\n", R1);
		printf("  R2     = %.3f s⁻¹\n", R2);
		printf("  B0     = %.0f MHz   gamma = %.4f\n", b0_mhz, gamma);
		printf("----------------------------------------------------------------\n");

		if (scan_offset)
		{
			printf("  Scan mode: OFF-RESONANCE — vary dOmega at fixed omega1\n");
			printf("  Fixed omega1  = %.1f Hz\n", omega1_min_hz);
			printf("  dOmega range  = %.1f – %.1f Hz\n", dOmega_min_hz, dOmega_max_hz);
		}
		else
		{
			printf("  Scan mode: SPIN-LOCK FIELD — vary omega1 at fixed dOmega\n");
			printf("  Fixed dOmega  = %.2f Hz\n", dOmega_hz);
			printf("  omega1 range  = %.1f – %.1f Hz\n", omega1_min_hz, omega1_max_hz);
		}

		printf("----------------------------------------------------------------\n");
		printf("  %-12s  %-12s  %-10s  %-8s  %-10s  %-10s  %-10s\n",
			   scan_offset ? "dOmega(Hz)" : "nu1 (Hz)",
			   "weff_avg(Hz)", "theta_avg", "sin²θ_avg",
			   "Rex", "R1rho", "R2_obs");
		printf("  %-12s  %-12s  %-10s  %-8s  %-10s  %-10s  %-10s\n",
			   "------------", "------------", "----------", "--------",
			   "----------", "----------", "----------");

		std::vector<double> x_vec, R1rho_vec, R2obs_vec;

		for (size_t i = 0; i < n_points; ++i)
		{
			double omega1_cur, dOmega_cur;

			if (scan_offset)
			{
				double frac = static_cast<double>(i) / (n_points - 1);
				double dOmega_hz_cur = dOmega_min_hz + (dOmega_max_hz - dOmega_min_hz) * frac;
				dOmega_cur = dOmega_hz_cur * 2.0 * PI;
				omega1_cur = omega1_fixed_rads;
			}
			else
			{
				double frac = static_cast<double>(i) / (n_points - 1);
				double nu1 = omega1_min_hz * std::pow(omega1_max_hz / omega1_min_hz, frac);
				omega1_cur = nu1 * 2.0 * PI;
				dOmega_cur = dOmega_fixed_rads;
			}

			// --- Miloushev-Bushweller Model Equations ---

			// 1. Calculate specific state offsets relative to the carrier position
			// Assumes carrier (dOmega_cur) is referenced to the population-weighted center.
			double wA = dOmega_cur - pB * dw_rads;
			double wB = dOmega_cur + pA * dw_rads;

			// 2. State-specific effective fields and geometric parameters
			double weffA = std::sqrt(omega1_cur * omega1_cur + wA * wA);
			double weffB = std::sqrt(omega1_cur * omega1_cur + wB * wB);

			double thetaA = std::atan2(omega1_cur, wA);
			double thetaB = std::atan2(omega1_cur, wB);

			double sinA = std::sin(thetaA);
			double cosA = std::cos(thetaA);
			double sinB = std::sin(thetaB);
			double cosB = std::cos(thetaB);

			// Average/effective framing for baseline removal handling
			double theta_avg = pA * thetaA + pB * thetaB;
			double cos2_avg = pA * (cosA * cosA) + pB * (cosB * cosB);
			double sin2_avg = pA * (sinA * sinA) + pB * (sinB * sinB);
			double weff_avg = pA * weffA + pB * weffB;

			// 3. Intrinsic baseline relaxation tracking
			double R1rho_0 = R1 * cos2_avg + R2 * sin2_avg;

			// 4. Miloushev and Bushweller analytical exchange term (Rex)
			// Evaluates the delta frequency shifts scaled against effective fields
			double delta_w_eff = weffA - weffB;
			double sin_delta_theta = std::sin(thetaA - thetaB);

			double numerator = pA * pB * kex * (delta_w_eff * delta_w_eff + omega1_cur * omega1_cur * sin_delta_theta * sin_delta_theta);
			double denominator = kex * kex + (pB * weffA + pA * weffB) * (pB * weffA + pA * weffB);

			double Rex_contrib = numerator / denominator;

			// Total R1rho
			double R1rho = R1rho_0 + Rex_contrib;

			// Extract observed R2 using the customized geometric denominator
			double R2_obs = (sin2_avg > 1e-6) ? (R1rho - R1 * cos2_avg) / sin2_avg : R2;

			double x_val = scan_offset ? dOmega_cur / (2.0 * PI) : omega1_cur / (2.0 * PI);

			x_vec.push_back(x_val);
			R1rho_vec.push_back(R1rho);
			R2obs_vec.push_back(R2_obs);

			printf("  %-12.2f  %-12.2f  %-10.2f  %-8.4f  %-10.4f  %-10.4f  %-10.4f\n",
				   x_val,
				   weff_avg / (2.0 * PI),
				   theta_avg * 180.0 / PI,
				   sin2_avg,
				   Rex_contrib,
				   R1rho,
				   R2_obs);
		}

		printf("================================================================\n\n");

		// --- Write output file ---
		if (filename.empty())
			return;

		std::ofstream out(filename);
		if (!out.is_open())
		{
			std::cerr << "  [R1rho] Error: could not open file: " << filename << "\n";
			return;
		}

		out << "# R1rho dispersion simulation — Miloushev & Bushweller Analytical Model\n";
		out << "# scan_mode  = " << scan_mode << "\n";
		out << "# kex        = " << kex << " s⁻¹\n";
		out << "# pA         = " << pA << "\n";
		out << "# dw_ppm     = " << dw_ppm << " ppm\n";
		out << "# R1         = " << R1 << " s⁻¹\n";
		out << "# R2         = " << R2 << " s⁻¹\n";

		if (scan_offset)
		{
			out << "# omega1_fixed = " << omega1_min_hz << " Hz\n";
			out << "# dOmega_range = " << dOmega_min_hz << " to " << dOmega_max_hz << " Hz\n";
			out << "# dOmega_hz\tweff_avg_hz\ttheta_avg_deg\tsin2theta_avg\tRex\tR1rho\tR2obs\n";
		}
		else
		{
			out << "# dOmega_fixed = " << dOmega_hz << " Hz\n";
			out << "# omega1_range = " << omega1_min_hz << " to " << omega1_max_hz << " Hz\n";
			out << "# nu1_hz\tweff_avg_hz\ttheta_avg_deg\tsin2theta_avg\tRex\tR1rho\tR2obs\n";
		}

		for (size_t i = 0; i < x_vec.size(); ++i)
		{
			double o1 = scan_offset ? omega1_fixed_rads : x_vec[i] * 2.0 * PI;
			double dO = scan_offset ? x_vec[i] * 2.0 * PI : dOmega_fixed_rads;

			double wA = dO - pB * dw_rads;
			double wB = dO + pA * dw_rads;
			double weffA = std::sqrt(o1 * o1 + wA * wA);
			double weffB = std::sqrt(o1 * o1 + wB * wB);

			double thetaA = std::atan2(o1, wA);
			double thetaB = std::atan2(o1, wB);
			double sin2_avg = pA * std::sin(thetaA) * std::sin(thetaA) + pB * std::sin(thetaB) * std::sin(thetaB);
			double cos2_avg = 1.0 - sin2_avg;

			out << x_vec[i] << "\t"
				<< (pA * weffA + pB * weffB) / (2.0 * PI) << "\t"
				<< (pA * thetaA + pB * thetaB) * 180.0 / PI << "\t"
				<< sin2_avg << "\t"
				<< R1rho_vec[i] - (R1 * cos2_avg + R2 * sin2_avg) << "\t"
				<< R1rho_vec[i] << "\t"
				<< R2obs_vec[i] << "\n";
		}

		std::cout << "  [R1rho] Dispersion profile written to: " << filename << "\n";
	}

} // namespace DispersionUtils
