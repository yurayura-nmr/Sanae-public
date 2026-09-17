// g++ -std=c++17 -Wall -I../../include/eigen3/ -o test_r1_inversion_recovery ../r1_inversion_recovery.cpp ../Solution.cpp ../State.cpp ../Environment.cpp ../RungeKutta.cpp ../fileio.cpp test_r1_inversion_recovery.cpp

#include "../r1_inversion_recovery.h"

#include <Eigen/Dense>
#include <iostream>
#include <cassert>
#include <cmath>
#include <cstdio>
#include <string>

int main()
{
	std::cout << "=== Testing Sanae::R1_inversion_recovery ===\n";

	const std::string test_output_file = "test_r1_inversion_recovery_tmp.txt";

	// -------------------------------------------------------------
	std::cout << "\n--- Basic run: default parameters ---\n";
	{
		Eigen::MatrixXcd M = Sanae::R1_inversion_recovery(
			5.0, 7.0, 1.0, 0.1, 'N', test_output_file, false);

		bool nonempty = (M.rows() > 0 && M.cols() > 0);
		assert(nonempty);
		std::cout << "Trace Shape Check: "
				  << (nonempty ? "PASSED" : "FAILED")
				  << " (rows=" << M.rows() << ", cols=" << M.cols() << ")\n";
	}

	// -------------------------------------------------------------
	std::cout << "\n--- Initial magnetization matches Mx_iR = exp(-tau*R1) (Meq=1) ---\n";
	{
		double R1 = 5.0;
		double tau = 0.1;
		double expected_Mx_iR = std::exp(-tau * R1);

		Eigen::MatrixXcd M = Sanae::R1_inversion_recovery(
			R1, 7.0, 1.0, tau, 'N', test_output_file, false);

		double actual_M0 = M(0, 0).real();

		std::cout << "Expected M(0) = " << expected_Mx_iR << "\n";
		std::cout << "Actual   M(0) = " << actual_M0 << "\n";

		bool matches = std::abs(actual_M0 - expected_Mx_iR) < 1e-9;
		assert(matches);
		std::cout << "Initial Magnetization Check: "
				  << (matches ? "PASSED" : "FAILED")
				  << " (matches documented Meq=1 formula)\n";
	}

	// -------------------------------------------------------------
	std::cout << "\n--- tau = 0 gives no relaxation decay yet ---\n";
	{
		double R1 = 5.0;
		double tau = 0.0;
		double expected_Mx_iR = std::exp(-tau * R1); // == 1.0

		Eigen::MatrixXcd M = Sanae::R1_inversion_recovery(
			R1, 7.0, 1.0, tau, 'N', test_output_file, false);

		double actual_M0 = M(0, 0).real();

		bool matches = std::abs(actual_M0 - expected_Mx_iR) < 1e-9;
		assert(matches);
		std::cout << "Tau=0 Boundary Check: "
				  << (matches ? "PASSED" : "FAILED")
				  << " (M(0) = " << actual_M0 << ", expected 1.0)\n";
	}

	// -------------------------------------------------------------
	std::cout << "\n--- Larger R1 decays faster (monotonic sanity check) ---\n";
	{
		double tau = 0.1;

		Eigen::MatrixXcd M_slow = Sanae::R1_inversion_recovery(
			2.0, 7.0, 1.0, tau, 'N', test_output_file, false);
		Eigen::MatrixXcd M_fast = Sanae::R1_inversion_recovery(
			20.0, 7.0, 1.0, tau, 'N', test_output_file, false);

		double M0_slow = M_slow(0, 0).real();
		double M0_fast = M_fast(0, 0).real();

		std::cout << "M(0) with R1=2.0  = " << M0_slow << "\n";
		std::cout << "M(0) with R1=20.0 = " << M0_fast << "\n";

		bool ordering_ok = M0_fast < M0_slow;
		assert(ordering_ok);
		std::cout << "R1 Ordering Check: "
				  << (ordering_ok ? "PASSED" : "FAILED")
				  << " (higher R1 -> smaller initial Mx)\n";
	}

	// -------------------------------------------------------------
	std::cout << "\n--- Runs across supported nuclei without crashing ---\n";
	{
		bool all_ok = true;
		for (char nucleus : {'H', 'N', 'C', 'F'})
		{
			Eigen::MatrixXcd M = Sanae::R1_inversion_recovery(
				5.0, 7.0, 1.0, 0.1, nucleus, test_output_file, false);
			if (M.rows() == 0 || M.cols() == 0)
				all_ok = false;
		}
		assert(all_ok);
		std::cout << "Nucleus Sweep Check: "
				  << (all_ok ? "PASSED" : "FAILED")
				  << " (H, N, C, F all produced a trace)\n";
	}

	// -------------------------------------------------------------
	// Cleanup: remove the temp output file created by Solution's
	// internal writeMagnetizationToFile() call.
	std::remove(test_output_file.c_str());

	std::cout << "\nAll tests passed! r1_inversion_recovery.cpp (R1_inversion_recovery) is solid.\n";
	return 0;
}