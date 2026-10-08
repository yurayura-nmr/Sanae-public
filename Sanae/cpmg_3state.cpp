#include "cpmg_3state.h"
#include <cmath>
#include <cstdio>
#include <fstream>
#include <iostream>
#include <complex>
#include <vector>
#include <algorithm> // std::max, std::min (portability: not relied on transitively)
#include <Eigen/Dense>
#include <unsupported/Eigen/MatrixFunctions> // matrix exponential

namespace Sanae
{

	// ---------------------------------------------------------------------------
	// Internal helpers
	// ---------------------------------------------------------------------------
	namespace
	{

		constexpr double PI = 3.14159265358979323846;

		/**
		 * Build 3x3 complex Liouvillian for cpmg_matrix3f topology (A<->B, A<->C).
		 *
		 *   L = K + i*Omega
		 *
		 * Exchange matrix K (real):
		 *   K = [ -(kAB+kAC)   kBA          kCA         ]
		 *       [   kAB       -kBA          0           ]
		 *       [   kAC        0           -kCA         ]
		 *
		 * Frequency matrix Omega (imaginary diagonal):
		 *   Omega = diag(0, dwAB, dwAC)   in rad/s, relative to state A
		 *
		 * The 180 degree refocusing pulse is modelled by conjugating L,
		 * which negates all imaginary parts (i.e. flips sign of all dw offsets).
		 * This is exact and parallels the sign-flip of dw in the 2-state code.
		 *
		 * @param dwAB   Omega_B - Omega_A  (rad/s, signed)
		 * @param dwAC   Omega_C - Omega_A  (rad/s, signed)
		 * @param kAB    A->B  (s^-1)
		 * @param kBA    B->A  (s^-1)
		 * @param kAC    A->C  (s^-1)
		 * @param kCA    C->A  (s^-1)
		 * @param R20    Intrinsic R2, assumed equal for all states (s^-1)
		 */
		Eigen::Matrix3cd buildL(
			double dwAB, double dwAC,
			double kAB, double kBA,
			double kAC, double kCA,
			double R20)
		{
			using cd = std::complex<double>;
			const cd I(0.0, 1.0);

			Eigen::Matrix3cd L = Eigen::Matrix3cd::Zero();

			// Exchange (real) part -- cpmg_matrix3f topology, no B<->C
			L(0, 0) = -(kAB + kAC);
			L(0, 1) = kBA;
			L(0, 2) = kCA;
			L(1, 0) = kAB;
			L(1, 1) = -kBA;
			L(2, 0) = kAC;
			L(2, 2) = -kCA;

			// Relaxation + frequency offset on diagonal
			// State A: reference frequency 0
			// State B: offset +dwAB  (rad/s)
			// State C: offset +dwAC  (rad/s)
			L(0, 0) -= R20;
			L(1, 1) -= R20 - I * dwAB;
			L(2, 2) -= R20 - I * dwAC;

			return L;
		}

		/**
		 * Raise a 3x3 matrix to integer power via repeated squaring.
		 */
		Eigen::Matrix3cd matpow(const Eigen::Matrix3cd &M, int n)
		{
			Eigen::Matrix3cd result = Eigen::Matrix3cd::Identity();
			Eigen::Matrix3cd base = M;
			while (n > 0)
			{
				if (n & 1)
					result = result * base;
				base = base * base;
				n >>= 1;
			}
			return result;
		}

		/**
		 * Compute R2eff for one (nu_CPMG, Trelax) point given a Liouvillian L.
		 *
		 * Pulse sequence per echo unit -- mirrors cpmg_bloch_mcconnell.cpp Eq 29:
		 *
		 *   block = exp(L* * tcp) * exp(L * tcp) * exp(L * tcp) * exp(L* * tcp)
		 *
		 * where L* = L.conjugate() is the propagator after the 180 pulse.
		 * Applied N/2 times: total time = 4*tcp*(N/2) = 2*N*tcp = Trelax.
		 *
		 * R2eff = -(1/Trelax) * ln( |sum_i M_i(Trelax)| / |sum_i M_i(0)| )
		 *
		 * @param L       Liouvillian (3x3 complex)
		 * @param pA,pB,pC  Equilibrium populations (initial magnetisation)
		 * @param tcp     Half-echo spacing = Trelax / (2*N)   (s)
		 * @param N       Number of pi pulses (even)
		 * @param Trelax  Total CPMG relaxation period          (s)
		 */
		double calcR2eff3(
			const Eigen::Matrix3cd &L,
			double pA, double pB, double pC,
			double tcp, int N, double Trelax)
		{
			Eigen::Matrix3cd Lc = L.conjugate(); // after 180 pulse

			Eigen::Matrix3cd eL = (L * tcp).exp();
			Eigen::Matrix3cd eLc = (Lc * tcp).exp();

			// Four-exponential symmetric block (matches 2-state convention)
			Eigen::Matrix3cd block = eL * eLc * eLc * eL;

			Eigen::Matrix3cd U = matpow(block, N / 2);

			// Initial magnetisation
			Eigen::Vector3cd M0(pA, pB, pC);
			Eigen::Vector3cd Mf = U * M0;

			double mag0 = std::abs(M0(0) + M0(1) + M0(2));
			double magf = std::abs(Mf(0) + Mf(1) + Mf(2));

			if (mag0 < 1e-12 || magf < 1e-12)
				return 0.0;

			return -(1.0 / Trelax) * std::log(magf / mag0);
		}

		/**
		 * Derive equilibrium populations from rate constants (detailed balance).
		 * pB/pA = kAB/kBA,  pC/pA = kAC/kCA,  pA+pB+pC = 1.
		 *
		 * Requires kBA > 0 whenever kAB > 0 (and symmetrically for C); an
		 * irreversible one-way rate (forward > 0, backward <= 0) makes
		 * detailed balance undefined, so this is rejected rather than
		 * silently producing a population that does not reflect true
		 * equilibrium. All rates must also be non-negative.
		 *
		 * Returns false (with error set) if the rates are invalid; pA, pB,
		 * pC are left unmodified in that case.
		 */
		bool populations(double kAB, double kBA, double kAC, double kCA,
						  double &pA, double &pB, double &pC,
						  std::string &error)
		{
			if (kAB < 0.0 || kBA < 0.0 || kAC < 0.0 || kCA < 0.0)
			{
				error = "all exchange rates must be non-negative";
				return false;
			}
			if (kAB > 0.0 && kBA <= 0.0)
			{
				error = "kAB > 0 but kBA <= 0 (irreversible A->B exchange; "
						"detailed balance undefined)";
				return false;
			}
			if (kAC > 0.0 && kCA <= 0.0)
			{
				error = "kAC > 0 but kCA <= 0 (irreversible A->C exchange; "
						"detailed balance undefined)";
				return false;
			}

			double rB = (kBA > 0.0) ? kAB / kBA : 0.0;
			double rC = (kCA > 0.0) ? kAC / kCA : 0.0;
			double norm = 1.0 + rB + rC;
			pA = 1.0 / norm;
			pB = rB / norm;
			pC = rC / norm;
			return true;
		}

		/**
		 * One (N, nu) point of the CPMG frequency axis.
		 */
		struct NAxisPoint
		{
			int N;
			double nu;
		};

		/**
		 * Build the N axis (even integers, log-spaced, deduplicated by nu)
		 * shared by Matrix3State and Matrix3StateSignDemo. Previously this
		 * logic was duplicated in both functions; factored out here so a
		 * fix only needs to be made in one place.
		 *
		 * NOTE: nu uses the convention nu = N / Trelax -- see the @todo in
		 * cpmg_3state.h regarding this convention (deferred).
		 */
		std::vector<NAxisPoint> BuildNAxis(
			double nu_min, double nu_max, int n_points, double Trelax)
		{
			int N_min = 2 * std::max(1,
									  static_cast<int>(std::round(nu_min * Trelax / 2.0)));
			int N_max = 2 * std::max(1,
									  static_cast<int>(std::round(nu_max * Trelax / 2.0)));
			if (N_max <= N_min)
				N_max = N_min + 2;

			std::vector<NAxisPoint> axis;
			for (int i = 0; i < n_points; ++i)
			{
				double frac = static_cast<double>(i) / (n_points - 1);
				double N_raw = N_min * std::pow(
										   static_cast<double>(N_max) / N_min, frac);

				int N = 2 * std::max(1,
									  static_cast<int>(std::round(N_raw / 2.0)));
				N = std::max(N_min, std::min(N_max, N));

				double nu = static_cast<double>(N) / Trelax;
				if (!axis.empty() && std::abs(nu - axis.back().nu) < 0.01)
					continue;

				axis.push_back(NAxisPoint{N, nu});
			}
			return axis;
		}

		/**
		 * Walk a pre-built N axis and print/store one dispersion profile.
		 */
		void runProfile(
			const Eigen::Matrix3cd &L,
			double pA, double pB, double pC,
			double Trelax,
			const std::vector<NAxisPoint> &axis,
			std::vector<double> &nu_vec,
			std::vector<double> &R2_vec)
		{
			for (const auto &pt : axis)
			{
				double tcp = Trelax / (2.0 * pt.N);
				double R2eff = calcR2eff3(L, pA, pB, pC, tcp, pt.N, Trelax);

				nu_vec.push_back(pt.nu);
				R2_vec.push_back(R2eff);

				printf("  %-14.2f  %-8d  %-10.4f  %-10.4f\n",
					   pt.nu, pt.N, tcp * 1000.0, R2eff);
			}
		}

	} // anonymous namespace

	// ---------------------------------------------------------------------------
	// Public API
	// ---------------------------------------------------------------------------

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
		double nu_min,
		double nu_max,
		int n_points,
		const std::string &label,
		const std::string &filename)
	{
		const double Trelax = 0.04; // standard CPMG relaxation delay (s)

		double dwAB = dwAB_ppm * gamma * B0_MHz * 2.0 * PI; // rad/s
		double dwAC = dwAC_ppm * gamma * B0_MHz * 2.0 * PI;

		double pA, pB, pC;
		std::string error;
		if (!populations(kAB, kBA, kAC, kCA, pA, pB, pC, error))
		{
			std::cerr << "  [3STATE] Error: " << error << ". Aborting simulation.\n";
			return Matrix3StateResult{};
		}

		Eigen::Matrix3cd L = buildL(dwAB, dwAC, kAB, kBA, kAC, kCA, R20);

		// --- Header ---
		printf("\n");
		printf("================================================================\n");
		printf("  3-State Matrix CPMG Dispersion [Korzhnev et al. Nature 2004]\n");
		printf("  Topology: A<->B, A<->C  (cpmg_matrix3f)\n");
		if (!label.empty())
			printf("  %s\n", label.c_str());
		printf("----------------------------------------------------------------\n");
		printf("  pA=%.4f  pB=%.4f  pC=%.4f\n", pA, pB, pC);
		printf("  kAB=%.1f  kBA=%.1f  kAC=%.1f  kCA=%.1f s^-1\n",
			   kAB, kBA, kAC, kCA);
		printf("  dwAB = %+.3f ppm (%+.2f rad/s)\n", dwAB_ppm, dwAB);
		printf("  dwAC = %+.3f ppm (%+.2f rad/s)\n", dwAC_ppm, dwAC);
		printf("  R20 = %.2f s^-1   B0 = %.0f MHz   Trelax = %.3f s\n",
			   R20, B0_MHz, Trelax);
		printf("----------------------------------------------------------------\n");
		printf("  %-14s  %-8s  %-10s  %-10s\n",
			   "nu_CPMG (Hz)", "N", "tcp (ms)", "R2eff (s^-1)");
		printf("  %-14s  %-8s  %-10s  %-10s\n",
			   "--------------", "--------", "----------", "----------");

		std::vector<NAxisPoint> axis = BuildNAxis(nu_min, nu_max, n_points, Trelax);

		Matrix3StateResult result;
		runProfile(L, pA, pB, pC, Trelax, axis, result.nu, result.R2eff);

		printf("================================================================\n\n");

		if (filename.empty())
			return result;

		std::ofstream out(filename);
		if (!out.is_open())
		{
			std::cerr << "  [3STATE] Cannot open: " << filename << "\n";
			return result;
		}
		out << "# 3-state matrix CPMG -- " << label << "\n";
		out << "# dwAB_ppm=" << dwAB_ppm << "  dwAC_ppm=" << dwAC_ppm << "\n";
		out << "# kAB=" << kAB << "  kBA=" << kBA
			<< "  kAC=" << kAC << "  kCA=" << kCA << " s^-1\n";
		out << "# nu_CPMG_Hz\tR2eff_s\n";
		for (size_t i = 0; i < result.nu.size(); ++i)
			out << result.nu[i] << "\t" << result.R2eff[i] << "\n";

		std::cout << "  [3STATE] Written to: " << filename << "\n";

		return result;
	}

	// ---------------------------------------------------------------------------

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
		double nu_min,
		double nu_max,
		int n_points,
		const std::string &filename)
	{
		const double Trelax = 0.04;

		double dwAB = dwAB_ppm * gamma * B0_MHz * 2.0 * PI;
		double dwAC = dwAC_ppm * gamma * B0_MHz * 2.0 * PI;

		double pA, pB, pC;
		std::string error;
		if (!populations(kAB, kBA, kAC, kCA, pA, pB, pC, error))
		{
			std::cerr << "  [SIGNDEMO] Error: " << error << ". Aborting simulation.\n";
			return Matrix3StateSignDemoResult{};
		}

		// Three Liouvillians: original, global flip, relative flip
		Eigen::Matrix3cd L_orig = buildL(dwAB, dwAC, kAB, kBA, kAC, kCA, R20);
		Eigen::Matrix3cd L_gflip = buildL(-dwAB, -dwAC, kAB, kBA, kAC, kCA, R20);
		Eigen::Matrix3cd L_rflip = buildL(-dwAB, dwAC, kAB, kBA, kAC, kCA, R20);

		printf("\n");
		printf("================================================================\n");
		printf("  3-State CPMG Global Sign Degeneracy Demo\n");
		printf("  Ref: Skrynnikov, Dahlquist & Kay, JACS 124, 12352 (2002)\n");
		printf("----------------------------------------------------------------\n");
		printf("  Original:      dwAB=%+.3f  dwAC=%+.3f ppm\n",
			   dwAB_ppm, dwAC_ppm);
		printf("  Global flip:   dwAB=%+.3f  dwAC=%+.3f ppm  (should be identical)\n",
			   -dwAB_ppm, -dwAC_ppm);
		printf("  Relative flip: dwAB=%+.3f  dwAC=%+.3f ppm  (should differ)\n",
			   -dwAB_ppm, dwAC_ppm);
		printf("----------------------------------------------------------------\n");
		printf("  %-14s  %-12s  %-12s  %-12s  %s\n",
			   "nu_CPMG (Hz)", "R2_orig", "R2_gflip", "R2_rflip", "orig-gflip");
		printf("  %-14s  %-12s  %-12s  %-12s  %s\n",
			   "--------------", "------------", "------------", "------------", "----------");

		std::vector<NAxisPoint> axis = BuildNAxis(nu_min, nu_max, n_points, Trelax);

		Matrix3StateSignDemoResult result;

		for (const auto &pt : axis)
		{
			double tcp = Trelax / (2.0 * pt.N);

			double r_orig = calcR2eff3(L_orig, pA, pB, pC, tcp, pt.N, Trelax);
			double r_gflip = calcR2eff3(L_gflip, pA, pB, pC, tcp, pt.N, Trelax);
			double r_rflip = calcR2eff3(L_rflip, pA, pB, pC, tcp, pt.N, Trelax);
			double diff = r_orig - r_gflip;

			printf("  %-14.2f  %-12.4f  %-12.4f  %-12.4f  %+.2e\n",
				   pt.nu, r_orig, r_gflip, r_rflip, diff);

			result.nu.push_back(pt.nu);
			result.R2_orig.push_back(r_orig);
			result.R2_gflip.push_back(r_gflip);
			result.R2_rflip.push_back(r_rflip);
		}

		printf("================================================================\n");
		printf("  orig-gflip ~ 0 everywhere  ->  global flip is unobservable\n");
		printf("  R2_rflip != R2_orig        ->  relative sign IS observable\n");
		printf("================================================================\n\n");

		if (filename.empty())
			return result;

		std::ofstream out(filename);
		if (!out.is_open())
		{
			std::cerr << "  [SIGNDEMO] Cannot open: " << filename << "\n";
			return result;
		}
		out << "# 3-state CPMG global sign degeneracy demo\n";
		out << "# dwAB_orig=" << dwAB_ppm << " dwAC_orig=" << dwAC_ppm << " ppm\n";
		out << "# kAB=" << kAB << " kBA=" << kBA
			<< " kAC=" << kAC << " kCA=" << kCA << " s^-1\n";
		out << "# nu_CPMG_Hz\tR2_orig\tR2_global_flip\tR2_relative_flip\n";
		for (size_t i = 0; i < result.nu.size(); ++i)
			out << result.nu[i] << "\t"
				<< result.R2_orig[i] << "\t"
				<< result.R2_gflip[i] << "\t"
				<< result.R2_rflip[i] << "\n";

		std::cout << "  [SIGNDEMO] Written to: " << filename << "\n";

		return result;
	}

} // namespace Sanae