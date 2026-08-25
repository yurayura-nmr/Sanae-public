#include "cleanex_utils.h"
#include "fileio.h"
#include <cmath>
#include <cstdio>
#include <iostream>
#include <limits>
#include <vector>

namespace CleanexUtils
{

	void SimulateCleanex(
		double k,
		double R1a,
		double R1b,
		double tm_min,
		double tm_max,
		size_t n_points,
		const std::string &label,
		const std::string &filename)
	{
		if (n_points == 0)
		{
			std::cerr << "  [CLEANEX] Warning: n_points == 0, nothing to simulate.\n";
			return;
		}

		// --- Derived quantities ---
		constexpr double TINY = 1.0e-10;
		double denom = R1a - R1b + k;
		if (std::abs(denom) < TINY)
		{
			std::cerr << "  [CLEANEX] Warning: R1a - R1b + k \xE2\x89\x88 0 \xE2\x80\x94 degenerate case.\n"
					  << "            Perturb R1a or R1b slightly.\n";
			denom = TINY;
		}
		double prefactor = k / denom;

		// Find analytical maximum of I(tm) for reference
		double tm_peak = std::numeric_limits<double>::quiet_NaN();
		if ((R1a + k) > R1b && R1b > 0.0)
			tm_peak = std::log((R1a + k) / R1b) / ((R1a + k) - R1b);

		// --- Print header ---
		printf("\n========================================================\n");
		printf("  %s\n", label.c_str());
		printf("  CLEANEX-PM: I(tm) = k/(R1a-R1b+k) * (-exp(-(R1a+k)*tm) + exp(-R1b*tm))\n");
		printf("--------------------------------------------------------\n");
		printf("  k    = %.2f s\xE2\x81\xBB\xC2\xB9   (tau_ex = %.2f ms)\n", k, 1000.0 / k);
		printf("  R1a  = %.3f s\xE2\x81\xBB\xC2\xB9  (water R1)\n", R1a);
		printf("  R1b  = %.3f s\xE2\x81\xBB\xC2\xB9  (NH R1)\n", R1b);
		printf("  prefactor k/(R1a-R1b+k) = %.4f\n", prefactor);
		if (!std::isnan(tm_peak))
			printf("  Predicted peak at tm = %.1f ms\n", tm_peak * 1000.0);
		printf("  tm range: %.1f \xE2\x80\x93 %.1f ms  (%zu points)\n",
			   tm_min * 1000.0, tm_max * 1000.0, n_points);
		printf("--------------------------------------------------------\n");
		printf("  %-14s  %-12s\n", "tm (ms)", "I(tm)");
		printf("  %-14s  %-12s\n", "--------------", "------------");

		// --- Simulate curve ---
		std::vector<double> tm_vec(n_points);
		std::vector<double> I_vec(n_points);

		for (size_t i = 0; i < n_points; ++i)
		{
			double frac = (n_points == 1)
							  ? 0.0
							  : static_cast<double>(i) / static_cast<double>(n_points - 1);

			// Log-spaced mixing times give better coverage of the rising edge
			double tm = tm_min * std::pow(tm_max / tm_min, frac);

			double expR1a = std::exp(-(R1a + k) * tm);
			double expR1b = std::exp(-R1b * tm);
			double I = prefactor * (-expR1a + expR1b);

			tm_vec[i] = tm;
			I_vec[i] = I;

			printf("  %-14.4f  %-12.6f\n", tm * 1000.0, I);
		}

		printf("========================================================\n\n");

		// --- Write output file (delegated to fileio) ---
		if (filename.empty())
			return;

		if (Sanae::writeCleanexCurveToFile(filename, tm_vec, I_vec, k, R1a, R1b, prefactor))
			std::cout << "  [CLEANEX] Curve written to: " << filename << "\n";
	}

} // namespace CleanexUtils