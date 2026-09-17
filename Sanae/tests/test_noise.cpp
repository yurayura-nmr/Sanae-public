//  g++ -std=c++17 -Wall -I../../include/eigen3/ -o test_noise ../noise.cpp test_noise.cpp

#include "../noise.h"

#include <Eigen/Dense>
#include <iostream>
#include <cassert>
#include <cmath>
#include <complex>
#include <vector>

// ---------------------------------------------------------------------
// Helpers
// ---------------------------------------------------------------------

// Builds a deterministic non-trivial FID matrix for testing.
// Values are distinct per-cell so any accidental mutation is detectable.
Eigen::MatrixXcd MakeTestMatrix(int rows, int cols)
{
	Eigen::MatrixXcd M(rows, cols);
	for (int r = 0; r < rows; ++r)
		for (int c = 0; c < cols; ++c)
			M(r, c) = std::complex<double>(
				1.0 + 0.1 * r + c,
				0.5 - 0.1 * r + 0.01 * c);
	return M;
}

bool MatricesEqual(const Eigen::MatrixXcd &A, const Eigen::MatrixXcd &B)
{
	if (A.rows() != B.rows() || A.cols() != B.cols())
		return false;
	return A == B;
}

// ---------------------------------------------------------------------
// Main
// ---------------------------------------------------------------------

int main()
{
	std::cout << "=== Testing AddGaussianNoise (noise.cpp) ===\n";

	// -------------------------------------------------------------
	std::cout << "\n--- sigma = 0.0 (disabled) ---\n";
	{
		Eigen::MatrixXcd M = MakeTestMatrix(3, 5);
		Eigen::MatrixXcd M_ref = M;

		AddGaussianNoise(M, 0.0, 42);

		bool unchanged = MatricesEqual(M, M_ref);
		assert(unchanged);
		std::cout << "Zero Sigma Check: "
				  << (unchanged ? "PASSED" : "FAILED")
				  << " (matrix unchanged)\n";
	}

	// -------------------------------------------------------------
	std::cout << "\n--- sigma < 0.0 (treated as disabled) ---\n";
	{
		Eigen::MatrixXcd M = MakeTestMatrix(3, 5);
		Eigen::MatrixXcd M_ref = M;

		AddGaussianNoise(M, -0.05, 42);

		bool unchanged = MatricesEqual(M, M_ref);
		assert(unchanged);
		std::cout << "Negative Sigma Check: "
				  << (unchanged ? "PASSED" : "FAILED")
				  << " (matrix unchanged)\n";
	}

	// -------------------------------------------------------------
	std::cout << "\n--- Column 0 (t=0) protection ---\n";
	{
		Eigen::MatrixXcd M = MakeTestMatrix(4, 10);
		Eigen::MatrixXcd M_ref = M;

		// Deliberately large sigma to make any accidental mutation obvious
		AddGaussianNoise(M, 1.0, 42);

		bool col0_intact = true;
		for (Eigen::Index row = 0; row < M.rows(); ++row)
			if (M(row, 0) != M_ref(row, 0))
				col0_intact = false;

		assert(col0_intact);
		std::cout << "Column 0 Integrity Check: "
				  << (col0_intact ? "PASSED" : "FAILED")
				  << " (t=0 left exact)\n";

		bool any_changed = false;
		for (Eigen::Index col = 1; col < M.cols(); ++col)
			for (Eigen::Index row = 0; row < M.rows(); ++row)
				if (M(row, col) != M_ref(row, col))
					any_changed = true;

		assert(any_changed);
		std::cout << "Noise Applied Check: "
				  << (any_changed ? "PASSED" : "FAILED")
				  << " (cols >= 1 perturbed)\n";
	}

	// -------------------------------------------------------------
	std::cout << "\n--- Seed reproducibility ---\n";
	{
		Eigen::MatrixXcd M1 = MakeTestMatrix(5, 20);
		Eigen::MatrixXcd M2 = M1;

		AddGaussianNoise(M1, 0.02, 123);
		AddGaussianNoise(M2, 0.02, 123);

		bool identical = MatricesEqual(M1, M2);
		assert(identical);
		std::cout << "Same Seed Reproducibility Check: "
				  << (identical ? "PASSED" : "FAILED")
				  << " (seed=123 == seed=123)\n";

		Eigen::MatrixXcd M3 = MakeTestMatrix(5, 20);
		Eigen::MatrixXcd M4 = M3;

		AddGaussianNoise(M3, 0.02, 1);
		AddGaussianNoise(M4, 0.02, 2);

		bool different = !MatricesEqual(M3, M4);
		assert(different);
		std::cout << "Different Seed Divergence Check: "
				  << (different ? "PASSED" : "FAILED")
				  << " (seed=1 != seed=2)\n";
	}

	// -------------------------------------------------------------
	std::cout << "\n--- Noise statistics (N(0, sigma)) ---\n";
	{
		const int rows = 20;
		const int cols = 2000;
		const double sigma = 0.05;

		Eigen::MatrixXcd M = Eigen::MatrixXcd::Zero(rows, cols);
		AddGaussianNoise(M, sigma, 7);

		std::vector<double> samples;
		samples.reserve(rows * (cols - 1) * 2);
		for (Eigen::Index col = 1; col < M.cols(); ++col)
			for (Eigen::Index row = 0; row < M.rows(); ++row)
			{
				samples.push_back(M(row, col).real());
				samples.push_back(M(row, col).imag());
			}

		double mean = 0.0;
		for (double s : samples)
			mean += s;
		mean /= samples.size();

		double var = 0.0;
		for (double s : samples)
			var += (s - mean) * (s - mean);
		var /= samples.size();
		double stddev = std::sqrt(var);

		std::cout << "Sample mean   = " << mean << "  (expected ~0.0)\n";
		std::cout << "Sample stddev = " << stddev << "  (expected ~" << sigma << ")\n";

		bool mean_ok = std::abs(mean) < 0.005;
		bool stddev_ok = std::abs(stddev - sigma) < 0.005;
		assert(mean_ok);
		assert(stddev_ok);
		std::cout << "Noise Distribution Check: "
				  << ((mean_ok && stddev_ok) ? "PASSED" : "FAILED")
				  << " (matches N(0, sigma))\n";
	}

	// -------------------------------------------------------------
	std::cout << "\n--- Degenerate single-column matrix ---\n";
	{
		// Only col=0 (t=0) exists — loop body never executes
		Eigen::MatrixXcd M = MakeTestMatrix(3, 1);
		Eigen::MatrixXcd M_ref = M;

		AddGaussianNoise(M, 0.05, 42);

		bool unchanged = MatricesEqual(M, M_ref);
		assert(unchanged);
		std::cout << "Single-Column Safety Check: "
				  << (unchanged ? "PASSED" : "FAILED")
				  << " (no crash, no-op on t=0 only)\n";
	}

	std::cout << "\nAll tests passed! noise.cpp (AddGaussianNoise) is solid.\n";
	return 0;
}