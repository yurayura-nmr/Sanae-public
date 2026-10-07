//  cl /EHsc /std:c++17 /I"..\..\include\eigen3" /Fe:test_zz_exchange.exe ^
//     ../zz_exchange.cpp ../RungeKutta.cpp test_zz_exchange.cpp

#include "../zz_exchange.h"

#include <Eigen/Dense>
#include <iostream>
#include <cassert>
#include <cmath>
#include <cstdio>

int main()
{
	std::cout << "=== Testing Sanae::SimulateZZExchange ===\n";

	// -------------------------------------------------------------
	std::cout << "\n--- kex = 0 guard (division-by-zero protection) ---\n";
	{
		// k_ab = k_ba = 0 => kex = 0 => pA, pB undefined
		Sanae::ZZExchangeResult result = Sanae::SimulateZZExchange(
			0.0, 0.0, 1.0, 1.0, 0.01, 1.0, 10,
			"kex=0 guard test", "");

		bool all_empty = result.T.empty() && result.a11.empty() &&
						  result.a22.empty() && result.a12.empty() &&
						  result.a21.empty();

		std::cout << "T.size() = " << result.T.size()
				  << " (expected 0)\n";

		assert(all_empty);
		std::cout << "Zero Kex Guard Check: "
				  << (all_empty ? "PASSED" : "FAILED")
				  << " (empty result, no crash, no NaN)\n";
	}

	// -------------------------------------------------------------
	std::cout << "\n--- Negative kex guard ---\n";
	{
		// k_ab = -1 is unphysical; k_ab + k_ba < 0
		Sanae::ZZExchangeResult result = Sanae::SimulateZZExchange(
			-1.0, 0.5, 1.0, 1.0, 0.01, 1.0, 10,
			"negative kex guard test", "");

		bool all_empty = result.T.empty();
		assert(all_empty);
		std::cout << "Negative Kex Guard Check: "
				  << (all_empty ? "PASSED" : "FAILED")
				  << " (empty result)\n";
	}

	// -------------------------------------------------------------
	std::cout << "\n--- Basic run produces well-formed, index-aligned curves ---\n";
	{
		int n_points = 20;
		Sanae::ZZExchangeResult result = Sanae::SimulateZZExchange(
			1.0, 1.0, 2.0, 2.0, 0.01, 1.0, n_points,
			"basic run test", "");

		bool sizes_ok =
			static_cast<int>(result.T.size()) == n_points &&
			static_cast<int>(result.a11.size()) == n_points &&
			static_cast<int>(result.a22.size()) == n_points &&
			static_cast<int>(result.a12.size()) == n_points &&
			static_cast<int>(result.a21.size()) == n_points;

		assert(sizes_ok);
		std::cout << "Result Size Check: "
				  << (sizes_ok ? "PASSED" : "FAILED")
				  << " (all vectors have n_points=" << n_points << " entries)\n";

		// Mixing times should be increasing (log-spaced ascending)
		bool T_increasing = true;
		for (int i = 1; i < n_points; ++i)
			if (result.T[i] <= result.T[i - 1])
				T_increasing = false;

		assert(T_increasing);
		std::cout << "T Monotonic Increasing Check: "
				  << (T_increasing ? "PASSED" : "FAILED") << "\n";
	}

	// -------------------------------------------------------------
	std::cout << "\n--- Numerical result matches corrected analytical formula (equal R1) ---\n";
	{
		// Symmetric case: k_ab = k_ba = k => kex = 2k, pA = pB = 0.5
		double k_ab = 1.0, k_ba = 1.0;
		double R1 = 2.0;
		double kex = k_ab + k_ba;
		double pA = k_ba / kex;
		double pB = k_ab / kex;

		int n_points = 15;
		Sanae::ZZExchangeResult result = Sanae::SimulateZZExchange(
			k_ab, k_ba, R1, R1, 0.01, 1.0, n_points,
			"analytical cross-check", "");

		double max_abs_err_a11 = 0.0;
		double max_abs_err_a21 = 0.0;

		for (int i = 0; i < n_points; ++i)
		{
			double T = result.T[i];
			double e1 = std::exp(-kex * T);
			double e2 = std::exp(-R1 * T);

			double a11_expected = (pA + pB * e1) * e2;
			double a21_expected = pB * (1.0 - e1) * e2;

			double err_a11 = std::abs(result.a11[i] - a11_expected);
			double err_a21 = std::abs(result.a21[i] - a21_expected);

			max_abs_err_a11 = std::max(max_abs_err_a11, err_a11);
			max_abs_err_a21 = std::max(max_abs_err_a21, err_a21);
		}

		std::cout << "Max |a11 - a11_expected| = " << max_abs_err_a11 << "\n";
		std::cout << "Max |a21 - a21_expected| = " << max_abs_err_a21 << "\n";

		bool a11_ok = max_abs_err_a11 < 1e-4;
		bool a21_ok = max_abs_err_a21 < 1e-4;

		assert(a11_ok);
		assert(a21_ok);
		std::cout << "Analytical Match Check (exp(-kex*T) formula): "
				  << ((a11_ok && a21_ok) ? "PASSED" : "FAILED") << "\n";
	}

	// -------------------------------------------------------------
	std::cout << "\n--- Confirms the formula fix: exp(-2*kex*T) would NOT match ---\n";
	{
		// Same setup as above, but deliberately check against the
		// OLD (incorrect) exp(-2*kex*T) formula from the original
		// lecture-notes transcription, to make sure the numerical
		// result has actually diverged from it (i.e. we are not
		// accidentally still matching the old, wrong formula).
		double k_ab = 1.0, k_ba = 1.0;
		double R1 = 2.0;
		double kex = k_ab + k_ba;
		double pA = k_ba / kex;
		double pB = k_ab / kex;

		int n_points = 15;
		Sanae::ZZExchangeResult result = Sanae::SimulateZZExchange(
			k_ab, k_ba, R1, R1, 0.01, 1.0, n_points,
			"old formula divergence check", "");

		// Use the largest T point, where the two formulas diverge most
		double T = result.T.back();
		double e1_old = std::exp(-2.0 * kex * T); // the old, incorrect formula
		double e2 = std::exp(-R1 * T);
		double a11_old_formula = (pA + pB * e1_old) * e2;

		double diff = std::abs(result.a11.back() - a11_old_formula);

		std::cout << "At T = " << T << ":\n";
		std::cout << "  Numerical a11        = " << result.a11.back() << "\n";
		std::cout << "  Old formula (2*kex)  = " << a11_old_formula << "\n";
		std::cout << "  |difference|          = " << diff << "\n";

		bool diverges_from_old = diff > 1e-3;
		assert(diverges_from_old);
		std::cout << "Old Formula Divergence Check: "
				  << (diverges_from_old ? "PASSED" : "FAILED")
				  << " (confirms exp(-kex*T), not exp(-2*kex*T), is correct)\n";
	}

	// -------------------------------------------------------------
	std::cout << "\n--- Symmetry at long mixing time (equal rates => equal populations) ---\n";
	{
		// Fully symmetric case (k_ab = k_ba, R1_A = R1_B): at long T,
		// a11 and a22 should converge to the same value (pA = pB = 0.5),
		// and a12, a21 should likewise converge to each other.
		Sanae::ZZExchangeResult result = Sanae::SimulateZZExchange(
			2.0, 2.0, 1.0, 1.0, 5.0, 5.0, 1,
			"long-T symmetry check", "");

		double diff_auto = std::abs(result.a11[0] - result.a22[0]);
		double diff_cross = std::abs(result.a12[0] - result.a21[0]);

		std::cout << "a11 = " << result.a11[0] << ", a22 = " << result.a22[0]
				  << ", |diff| = " << diff_auto << "\n";
		std::cout << "a12 = " << result.a12[0] << ", a21 = " << result.a21[0]
				  << ", |diff| = " << diff_cross << "\n";

		bool symmetric = diff_auto < 1e-6 && diff_cross < 1e-6;
		assert(symmetric);
		std::cout << "Symmetric Populations Check: "
				  << (symmetric ? "PASSED" : "FAILED") << "\n";
	}

	// -------------------------------------------------------------
	std::cout << "\n--- Empty filename skips file write but still returns data ---\n";
	{
		Sanae::ZZExchangeResult result = Sanae::SimulateZZExchange(
			1.0, 1.0, 2.0, 2.0, 0.01, 1.0, 5,
			"no file write test", ""); // empty filename

		bool nonempty = !result.T.empty();
		assert(nonempty);
		std::cout << "No-File-Write Still Returns Data Check: "
				  << (nonempty ? "PASSED" : "FAILED") << "\n";
	}

	std::cout << "\nAll tests passed! zz_exchange.cpp (SimulateZZExchange) is solid.\n";
	return 0;
}