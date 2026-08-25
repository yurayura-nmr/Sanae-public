/**
 * @file cleanex_utils.h
 * @brief CLEANEX-PM hydrogen exchange rate simulation.
 *
 * Simulates the CLEANEX-PM NMR experiment for measuring amide proton
 * hydrogen exchange rates with bulk water.
 *
 * The CLEANEX-PM signal intensity as a function of mixing time tm is:
 *
 *   I(tm) = k / (R1a - R1b + k) * (-exp(-(R1a+k)*tm) + exp(-R1b*tm))
 *
 * where:
 *   k   = NH hydrogen exchange rate with water  (s⁻¹)
 *   R1a = R1 of bulk water                      (s⁻¹, typically ~0.3 s⁻¹)
 *   R1b = R1 of NH proton                       (s⁻¹, typically 1-3 s⁻¹)
 *   tm  = CLEANEX mixing time                   (s)
 *
 * The curve rises from 0, reaches a maximum, then decays.
 * The maximum position and height are sensitive to k.
 * Fast exchanging residues (loops, termini) show earlier, higher peaks.
 * Slow exchanging residues (buried beta-sheet) show later, lower peaks.
 *
 * Reference: Hwang et al. JACS 1998, Pellecchia et al. JBNMR 1999.
 * Fitting code reference: GLOVE func_cleanex.cpp
 *
 * @note R1a is the water R1 — measure separately or use ~0.3 s⁻¹ at 298K.
 *       R1b is the NH proton R1 — can be measured from standard R1 experiment
 *       or estimated from RelaxationUtils::EstimateR2_N scaled appropriately.
 */

#pragma once

#include <string>
#include <vector>

namespace CleanexUtils
{

	/**
	 * @brief Simulate CLEANEX-PM intensity curve I(tm).
	 *
	 * @param k          NH exchange rate              (s⁻¹)
	 * @param R1a        R1 of bulk water              (s⁻¹, default 0.3)
	 * @param R1b        R1 of NH proton               (s⁻¹)
	 * @param tm_min     Minimum mixing time           (s)
	 * @param tm_max     Maximum mixing time           (s)
	 * @param n_points   Number of mixing time points
	 * @param label      Label for output
	 * @param filename   Output file (tm vs I(tm))
	 */

	void SimulateCleanex(
		double k,
		double R1a,
		double R1b,
		double tm_min = 0.001,
		double tm_max = 0.5,
		size_t n_points = 50,
		const std::string &label = "CLEANEX-PM simulation",
		const std::string &filename = "cleanex_curve.txt");

} // namespace CleanexUtils
