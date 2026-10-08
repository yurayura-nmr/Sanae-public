//  cl /EHsc /std:c++17 /I"..\..\include\eigen3" /Fe:test_cpmg_3state.exe ^
//     ../cpmg_3state.cpp test_cpmg_3state.cpp

#include "../cpmg_3state.h"

#include <iostream>
#include <cassert>
#include <cmath>
#include <algorithm>

int main()
{
	std::cout << "=== Testing Sanae::Matrix3State / Matrix3StateSignDemo ===\n";

	// -------------------------------------------------------------
	std::cout << "\n--- Rate validity guards ---\n";
	{
		// Negative rate
		Sanae::Matrix3StateResult r1 = Sanae::Matrix3State(
			1.0, -1.0, -5.0, 100.0, 50.0, 200.0, 15.0, 600.0, 0.101,
			25.0, 2000.0, 20, "negative rate test", "");
		bool empty1 = r1.nu.empty() && r1.R2eff.empty();
		assert(empty1);
		std::cout << "Negative Rate Guard Check: "
				  << (empty1 ? "PASSED" : "FAILED") << "\n";

		// Irreversible A->B (kAB > 0, kBA <= 0)
		Sanae::Matrix3StateResult r2 = Sanae::Matrix3State(
			1.0, -1.0, 100.0, 0.0, 50.0, 200.0, 15.0, 600.0, 0.101,
			25.0, 2000.0, 20, "irreversible AB test", "");
		bool empty2 = r2.nu.empty();
		assert(empty2);
		std::cout << "Irreversible A->B Guard Check: "
				  << (empty2 ? "PASSED" : "FAILED") << "\n";

		// Irreversible A->C (kAC > 0, kCA <= 0)
		Sanae::Matrix3StateResult r3 = Sanae::Matrix3State(
			1.0, -1.0, 100.0, 400.0, 50.0, 0.0, 15.0, 600.0, 0.101,
			25.0, 2000.0, 20, "irreversible AC test", "");
		bool empty3 = r3.nu.empty();
		assert(empty3);
		std::cout << "Irreversible A->C Guard Check: "
				  << (empty3 ? "PASSED" : "FAILED") << "\n";

		// Same guards should also empty out Matrix3StateSignDemo
		Sanae::Matrix3StateSignDemoResult r4 = Sanae::Matrix3StateSignDemo(
			1.0, -1.0, 100.0, 0.0, 50.0, 200.0, 15.0, 600.0, 0.101,
			25.0, 2000.0, 20, "");
		bool empty4 = r4.nu.empty() && r4.R2_orig.empty() &&
					  r4.R2_gflip.empty() && r4.R2_rflip.empty();
		assert(empty4);
		std::cout << "SignDemo Guard Check: "
				  << (empty4 ? "PASSED" : "FAILED") << "\n";
	}

	// -------------------------------------------------------------
	std::cout << "\n--- No exchange: R2eff must equal R20 exactly, for every nu ---\n";
	{
		// With all rates zero, there is no exchange contribution (Rex=0)
		// and no coupling between states, so the (purely-A) magnetization
		// simply decays at R20 regardless of the pulse train. This is a
		// decisive, closed-form correctness check of the whole pipeline
		// (buildL, matpow, calcR2eff3, BuildNAxis, runProfile) without
		// needing access to any internal (anonymous-namespace) function.
		double R20 = 12.0;
		Sanae::Matrix3StateResult result = Sanae::Matrix3State(
			0.0, 0.0,      // dw values irrelevant: populations are 100% in A
			0.0, 0.0,      // kAB, kBA = 0 (no A<->B exchange)
			0.0, 0.0,      // kAC, kCA = 0 (no A<->C exchange)
			R20, 600.0, 0.101,
			25.0, 2000.0, 15,
			"no-exchange sanity check", "");

		bool nonempty = !result.R2eff.empty();
		assert(nonempty);

		double max_err = 0.0;
		for (double r2eff : result.R2eff)
			max_err = std::max(max_err, std::abs(r2eff - R20));

		std::cout << "R20 = " << R20 << ", max|R2eff - R20| over "
				  << result.R2eff.size() << " points = " << max_err << "\n";

		bool matches = max_err < 1e-9;
		assert(matches);
		std::cout << "No-Exchange R2eff=R20 Check: "
				  << (matches ? "PASSED" : "FAILED") << "\n";
	}

	// -------------------------------------------------------------
	std::cout << "\n--- Valid exchange case produces well-formed, ordered curve ---\n";
	{
		Sanae::Matrix3StateResult result = Sanae::Matrix3State(
			2.0, -1.5,
			100.0, 400.0,   // kAB, kBA
			50.0, 200.0,    // kAC, kCA
			15.0, 600.0, 0.101,
			25.0, 2000.0, 20,
			"valid exchange test", "");

		bool sizes_match = result.nu.size() == result.R2eff.size() && !result.nu.empty();
		assert(sizes_match);
		std::cout << "Result Size Match Check: "
				  << (sizes_match ? "PASSED" : "FAILED")
				  << " (" << result.nu.size() << " points)\n";

		bool nu_increasing = true;
		for (std::size_t i = 1; i < result.nu.size(); ++i)
			if (result.nu[i] <= result.nu[i - 1])
				nu_increasing = false;
		assert(nu_increasing);
		std::cout << "Nu Axis Strictly Increasing Check: "
				  << (nu_increasing ? "PASSED" : "FAILED") << "\n";

		// Qualitative CPMG behavior: Rex should be quenched at high
		// pulsing frequency, so R2eff at the lowest nu should exceed
		// R2eff at the highest nu for a case with real exchange present.
		bool rex_quenches = result.R2eff.front() > result.R2eff.back();
		std::cout << "R2eff[first nu] = " << result.R2eff.front()
				  << ", R2eff[last nu] = " << result.R2eff.back() << "\n";
		assert(rex_quenches);
		std::cout << "Rex Quenching Trend Check: "
				  << (rex_quenches ? "PASSED" : "FAILED") << "\n";
	}

	// -------------------------------------------------------------
	std::cout << "\n--- Global sign degeneracy: {dwAB,dwAC} -> {-dwAB,-dwAC} is unobservable ---\n";
	{
		Sanae::Matrix3StateSignDemoResult result = Sanae::Matrix3StateSignDemo(
			2.0, -1.5,
			100.0, 400.0,
			50.0, 200.0,
			15.0, 600.0, 0.101,
			25.0, 2000.0, 20,
			"");

		bool nonempty = !result.nu.empty();
		assert(nonempty);

		double max_diff_orig_gflip = 0.0;
		for (std::size_t i = 0; i < result.nu.size(); ++i)
			max_diff_orig_gflip = std::max(max_diff_orig_gflip,
				std::abs(result.R2_orig[i] - result.R2_gflip[i]));

		std::cout << "max|R2_orig - R2_gflip| over " << result.nu.size()
				  << " points = " << max_diff_orig_gflip << "\n";

		bool globally_degenerate = max_diff_orig_gflip < 1e-9;
		assert(globally_degenerate);
		std::cout << "Global Sign Degeneracy Check: "
				  << (globally_degenerate ? "PASSED" : "FAILED") << "\n";
	}

	// -------------------------------------------------------------
	std::cout << "\n--- Relative sign IS observable: {-dwAB,+dwAC} differs from original ---\n";
	{
		Sanae::Matrix3StateSignDemoResult result = Sanae::Matrix3StateSignDemo(
			2.0, -1.5,
			100.0, 400.0,
			50.0, 200.0,
			15.0, 600.0, 0.101,
			25.0, 2000.0, 20,
			"");

		double max_diff_orig_rflip = 0.0;
		for (std::size_t i = 0; i < result.nu.size(); ++i)
			max_diff_orig_rflip = std::max(max_diff_orig_rflip,
				std::abs(result.R2_orig[i] - result.R2_rflip[i]));

		std::cout << "max|R2_orig - R2_rflip| over " << result.nu.size()
				  << " points = " << max_diff_orig_rflip << "\n";

		// Should be clearly nonzero (unlike the global-flip case above),
		// since dwAB != dwAC here.
		bool relative_sign_observable = max_diff_orig_rflip > 1e-3;
		assert(relative_sign_observable);
		std::cout << "Relative Sign Observability Check: "
				  << (relative_sign_observable ? "PASSED" : "FAILED") << "\n";
	}

	// -------------------------------------------------------------
	std::cout << "\n--- Shared N-axis consistency between Matrix3State and SignDemo ---\n";
	{
		// Both functions build their nu axis via the same shared
		// BuildNAxis() helper with the same fixed internal Trelax, so
		// for identical nu_min/nu_max/n_points they must produce
		// identical nu vectors -- this exercises the de-duplication
		// refactor that replaced two independently-maintained copies.
		Sanae::Matrix3StateResult r1 = Sanae::Matrix3State(
			1.0, -1.0, 100.0, 400.0, 50.0, 200.0, 15.0, 600.0, 0.101,
			25.0, 2000.0, 20, "", "");
		Sanae::Matrix3StateSignDemoResult r2 = Sanae::Matrix3StateSignDemo(
			1.0, -1.0, 100.0, 400.0, 50.0, 200.0, 15.0, 600.0, 0.101,
			25.0, 2000.0, 20, "");

		bool same_size = r1.nu.size() == r2.nu.size();
		assert(same_size);

		bool same_values = true;
		for (std::size_t i = 0; i < r1.nu.size() && same_values; ++i)
			if (std::abs(r1.nu[i] - r2.nu[i]) > 1e-9)
				same_values = false;

		std::cout << "Matrix3State nu.size() = " << r1.nu.size()
				  << ", Matrix3StateSignDemo nu.size() = " << r2.nu.size() << "\n";

		assert(same_values);
		std::cout << "Shared N-Axis Consistency Check: "
				  << ((same_size && same_values) ? "PASSED" : "FAILED") << "\n";
	}

	// -------------------------------------------------------------
	std::cout << "\n--- Empty filename skips file write but still returns data ---\n";
	{
		Sanae::Matrix3StateResult result = Sanae::Matrix3State(
			1.0, -1.0, 100.0, 400.0, 50.0, 200.0, 15.0, 600.0, 0.101,
			25.0, 2000.0, 10, "no file write test", "");
		bool nonempty = !result.nu.empty();
		assert(nonempty);
		std::cout << "No-File-Write Still Returns Data Check: "
				  << (nonempty ? "PASSED" : "FAILED") << "\n";
	}

	std::cout << "\nAll tests passed! cpmg_3state.cpp (Matrix3State / Matrix3StateSignDemo) is solid.\n";
	return 0;
}