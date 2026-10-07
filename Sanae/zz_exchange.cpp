#include "zz_exchange.h"
#include "RungeKutta.h"
#include <Eigen/Dense>
#include <cmath>
#include <cstdio>
#include <fstream>
#include <iostream>
#include <vector>

namespace Sanae
{

	// ---------------------------------------------------------------------------
	// Propagate from T=0 to T=T_mix using n_steps RK4 steps (via the shared
	// Sanae RungeKutta integrator) for the longitudinal Bloch-McConnell system:
	//
	//   d/dt [dMz_A]  =  [-R1A - k_ab    k_ba    ] [dMz_A]
	//        [dMz_B]     [  k_ab      -R1B - k_ba] [dMz_B]
	//
	// This system is purely real (no chemical shift term, unlike the
	// transverse case), so A is built as a complex matrix with zero
	// imaginary parts purely to satisfy RungeKutta's MatrixXcd interface.
	//
	// Returns [a11, a22, a12, a21] at T_mix. Each element aij is computed
	// from a specific initial condition:
	//   a11, a21 : start with dMz_A=1, dMz_B=0  (A magnetization)
	//   a12, a22 : start with dMz_A=0, dMz_B=1  (B magnetization)
	// ---------------------------------------------------------------------------
	static void propagate(
		double k_ab, double k_ba,
		double R1_A, double R1_B,
		double T_mix, int n_steps,
		double& a11, double& a12,
		double& a21, double& a22)
	{
		double dt = T_mix / n_steps;

		Eigen::MatrixXcd A(2, 2);
		A(0, 0) = std::complex<double>(-(R1_A + k_ab), 0.0);
		A(0, 1) = std::complex<double>(k_ba, 0.0);
		A(1, 0) = std::complex<double>(k_ab, 0.0);
		A(1, 1) = std::complex<double>(-(R1_B + k_ba), 0.0);

		// RungeKutta expects M pre-sized to the total number of columns
		// (including the initial condition at col 0) and fills columns
		// 1..num_steps-1 in place; the final state after n_steps RK4
		// steps therefore lands at column n_steps, so M needs n_steps+1
		// columns and we pass num_steps = n_steps+1.
		const int num_cols = n_steps + 1;

		// Column 1: initial A magnetization -> gives a11, a21
		Eigen::MatrixXcd M1(2, num_cols);
		M1.col(0) = Eigen::Vector2cd(1.0, 0.0);
		RungeKutta(A, M1, dt, num_cols);

		// Column 2: initial B magnetization -> gives a12, a22
		Eigen::MatrixXcd M2(2, num_cols);
		M2.col(0) = Eigen::Vector2cd(0.0, 1.0);
		RungeKutta(A, M2, dt, num_cols);

		a11 = M1(0, n_steps).real(); // A stays in A
		a21 = M1(1, n_steps).real(); // A transfers to B
		a12 = M2(0, n_steps).real(); // B transfers to A
		a22 = M2(1, n_steps).real(); // B stays in B
	}

	ZZExchangeResult SimulateZZExchange(
		double k_ab,
		double k_ba,
		double R1_A,
		double R1_B,
		double T_min,
		double T_max,
		int n_points,
		const std::string& label,
		const std::string& filename)
	{
		double kex = k_ab + k_ba;

		// Guard against division by zero (and negative kex, which is
		// unphysical): pA = k_ba/kex and pB = k_ab/kex are undefined
		// when kex <= 0. Documented range is kex = 0.1 to 10 s^-1, so
		// this is already out of spec, but fail loudly rather than
		// silently propagating NaN through the rest of the function.
		if (kex <= 0.0)
		{
			std::cerr << "  [ZZ] Error: k_ab + k_ba (kex) = " << kex
				<< " <= 0; populations pA, pB are undefined. "
				<< "Aborting simulation.\n";
			return ZZExchangeResult{};
		}

		double pA = k_ba / kex;
		double pB = k_ab / kex;

		// --- Sanity checks for ZZ experiment validity ---
		printf("\n");
		printf("================================================================\n");
		printf("  %s\n", label.c_str());
		printf("  Farrow et al. (1994) J Biomol NMR 4, 727\n");
		printf("----------------------------------------------------------------\n");
		printf("  k_ab  = %.3f s^-1   k_ba = %.3f s^-1\n", k_ab, k_ba);
		printf("  kex   = %.3f s^-1\n", kex);
		printf("  pA    = %.3f       pB   = %.3f\n", pA, pB);
		printf("  R1_A  = %.3f s^-1   R1_B = %.3f s^-1\n", R1_A, R1_B);
		printf("  T range: %.3f - %.3f s (%d points)\n",
			T_min, T_max, n_points);
		printf("----------------------------------------------------------------\n");

		// Validity warnings
		if (pB < 0.1)
			printf("  WARNING: pB = %.3f < 0.1 -- minor peak may not be visible\n", pB);
		if (kex > 10.0)
			printf("  WARNING: kex = %.1f > 10 s^-1 -- approaching intermediate exchange\n"
				"           consider nmr2stateExchange for lineshape instead\n",
				kex);
		if (kex < 0.1)
			printf("  WARNING: kex = %.3f < 0.1 s^-1 -- exchange very slow,\n"
				"           crosspeak may not be visible within T_max\n",
				kex);
		double R1_avg = 0.5 * (R1_A + R1_B);
		if (kex < R1_avg)
			printf("  WARNING: kex (%.3f) < R1 (%.3f) -- magnetization lost before\n"
				"           exchange occurs; consider longer T1 or larger kex\n",
				kex, R1_avg);

		// Analytical check (R1A = R1B case) for validation.
		// NOTE: uses exp(-kex*T), not exp(-2*kex*T) -- see @note in
		// zz_exchange.h for why this differs from the source lecture
		// notes' Eq 26 (believed to be a transcription error there).
		bool equal_R1 = std::abs(R1_A - R1_B) < 0.01;
		if (equal_R1)
			printf("  R1A ~ R1B: analytical expressions apply -- numerical\n"
				"  result should match closely (use as sanity check)\n");

		printf("----------------------------------------------------------------\n");
		printf("  %-8s  %-10s  %-10s  %-10s  %-10s",
			"T (s)", "a11(AA)", "a22(BB)", "a12(BA)", "a21(AB)");
		if (equal_R1)
			printf("  %-10s  %-10s  %-10s  %-10s",
				"a11_anal", "a22_anal", "a12_anal", "a21_anal");
		printf("\n");

		// Storage
		ZZExchangeResult result;
		result.T.resize(n_points);
		result.a11.resize(n_points);
		result.a22.resize(n_points);
		result.a12.resize(n_points);
		result.a21.resize(n_points);

		// Number of RK4 steps per mixing time -- 1000 per second is accurate
		constexpr int STEPS_PER_SECOND = 1000;

		for (int i = 0; i < n_points; ++i)
		{
			// Log-spaced mixing times
			double T = T_min * std::pow(T_max / T_min,
				static_cast<double>(i) / (n_points - 1));

			int n_steps = std::max(10, static_cast<int>(T * STEPS_PER_SECOND));

			double a11, a12, a21, a22;
			propagate(k_ab, k_ba, R1_A, R1_B, T, n_steps,
				a11, a12, a21, a22);

			result.T[i] = T;
			result.a11[i] = a11;
			result.a22[i] = a22;
			result.a12[i] = a12;
			result.a21[i] = a21;

			printf("  %-8.4f  %-10.6f  %-10.6f  %-10.6f  %-10.6f",
				T, a11, a22, a12, a21);

			// Analytical comparison (equal R1 case only)
			if (equal_R1)
			{
				double R1 = R1_A;
				double e1 = std::exp(-kex * T);
				double e2 = std::exp(-R1 * T);
				double a11a = (pA + pB * e1) * e2;
				double a22a = (pB + pA * e1) * e2;
				double a12a = pA * (1.0 - e1) * e2;
				double a21a = pB * (1.0 - e1) * e2;
				printf("  %-10.6f  %-10.6f  %-10.6f  %-10.6f",
					a11a, a22a, a12a, a21a);
			}
			printf("\n");
		}

		printf("================================================================\n\n");

		// --- Write output file ---
		if (filename.empty())
			return result;

		std::ofstream out(filename);
		if (!out.is_open())
		{
			std::cerr << "  [ZZ] Error: could not open: " << filename << "\n";
			return result;
		}

		out << "# SANAE ZZ exchange simulation\n";
		out << "# Farrow et al. (1994) J Biomol NMR 4, 727\n";
		out << "# k_ab = " << k_ab << " s^-1\n";
		out << "# k_ba = " << k_ba << " s^-1\n";
		out << "# kex  = " << kex << " s^-1\n";
		out << "# pA   = " << pA << "\n";
		out << "# pB   = " << pB << "\n";
		out << "# R1_A = " << R1_A << " s^-1\n";
		out << "# R1_B = " << R1_B << " s^-1\n";
		out << "# T_s\ta11_AA\ta22_BB\ta12_BA\ta21_AB\n";

		for (int i = 0; i < n_points; ++i)
			out << result.T[i] << "\t"
			<< result.a11[i] << "\t"
			<< result.a22[i] << "\t"
			<< result.a12[i] << "\t"
			<< result.a21[i] << "\n";

		std::cout << "  [ZZ] Exchange curves written to: " << filename << "\n";

		return result;
	}

}