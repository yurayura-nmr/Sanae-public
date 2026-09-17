#include "noise.h"
#include <random>

void AddGaussianNoise(Eigen::MatrixXcd &M, double sigma, int seed)
{
	// No-op if noise is disabled — avoids any random number overhead
	if (sigma <= 0.0)
		return;

	// Seed strategy: use random_device for non-reproducible noise (default),
	// or a fixed seed for reproducible simulations (e.g. for unit tests)
	std::mt19937 gen;
	if (seed < 0)
	{
		std::random_device rd;
		gen.seed(rd());
	}
	else
	{
		gen.seed(static_cast<unsigned int>(seed));
	}

	std::normal_distribution<double> dist(0.0, sigma);

	// Add independent Gaussian noise to real (Mx) and imaginary (My) channels
	// Skip col=0 — this is M(t=0), the initial condition, which must remain exact
	for (Eigen::Index col = 1; col < M.cols(); ++col)
		for (Eigen::Index row = 0; row < M.rows(); ++row)
			M(row, col) += std::complex<double>(dist(gen), dist(gen));
}