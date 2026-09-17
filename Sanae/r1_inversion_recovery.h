#pragma once

#include <Eigen/Dense>
#include <string>

/**
 * @file r1_inversion_recovery.h
 * @brief R1 (longitudinal) inversion-recovery simulation for a single spin.
 *
 * This function simulates the evolution of longitudinal magnetization following
 * a 180 degree inversion pulse. The magnetization obeys the standard
 * inversion-recovery relation:
 *
 *     Mz(tau) = Meq * (1 - 2 * exp(-R1 * tau))
 *
 * @note This implementation assumes Meq = 1 (normalized equilibrium
 *       magnetization). The current code computes only the bare decay
 *       factor exp(-R1 * tau), which is equivalent to the full formula
 *       above only under that Meq = 1 assumption. If Meq != 1 is ever
 *       needed, the full expression must be used instead.
 *
 * A numerical integrator (1-state Bloch equation solver) is used internally to
 * propagate the magnetization under longitudinal relaxation (R1), transverse
 * relaxation (R2), and chemical shift evolution.
 *
 * The function sets up:
 *  - A spectrometer environment (e.g., 700 MHz, selected nucleus)
 *  - A single spin state with user-specified parameters (R1, R2, shift)
 *  - An inversion-recovery delay `tau`
 *  - Numerical propagation via `Sanae::Solution(1)`
 *
 * @note
 * This is an educational, minimal example meant to illustrate T1 dynamics.
 * It is analogous in spirit to `bloch_relax()` and may be used in the same style.
 *
 * @todo Currently returns void; result is only observable via whatever
 *       Sanae::Solution::SetupOneState() does internally (e.g. file output).
 *       Consider returning the simulated magnetization trace (or a scalar
 *       Mz(tau) value) so this function is directly testable and so callers
 *       can use the trace programmatically (e.g. for ML training data
 *       generation from simulated R1 experiments/spectra).
 * @todo Sanae::Solution::SetupOneState() is a misnomer (it solves and writes
 *       output, not just "sets up") and itself returns void. A future
 *       revision should have it return the Eigen::MatrixXcd trace, which
 *       will also let this function return that trace. Renaming and
 *       signature changes are deferred to avoid breaking existing tests
 *       (test_solution.cpp) until that file is revisited deliberately.
 * @todo Extend to simulate a full inversion-recovery experiment: a series
 *       of tau values producing a full simulated spectral series (optionally
 *       with noise via AddGaussianNoise), rather than a single-tau scalar.
 *
 * @see bloch_relax(), nmr2stateExchange()
 *
 * @param R1
 *     Longitudinal relaxation rate in s^-1.
 *     Typical values:
 *     - Proteins at high field: 0.5-5 s^-1
 *     - Small molecules: 1-10 s^-1
 *
 * @param R2
 *     Transverse relaxation rate in s^-1.
 *     Affects transverse magnetization (spectral linewidth) after the delay.
 *
 * @param CSppm
 *     Chemical shift of the spin in parts per million (ppm).
 *     Position of peak in simulated spectrum.
 *
 * @param tau
 *     Inversion-recovery delay time in seconds.
 *     This is the time between the 180 degree inversion pulse and acquisition
 *     (or the subsequent 90 degree pulse). Determines how much magnetization
 *     recovers toward equilibrium.
 *
 * @param nucleus
 *     NMR nucleus type, e.g. `'H'`, `'N'`, `'C'`.
 *     Determines the gyromagnetic ratio and the chemical shift frequency.
 *
 * @param output_filename
 *     Path to write the simulated magnetization trace to. Defaults to
 *     "r1_inversion_recovery_output.txt". Pass a custom path (e.g. a temp
 *     file) if you don't want to overwrite a shared default output file —
 *     useful for tests.
 *
 * @return
 *     The simulated magnetization trace M(t), as returned by
 *     Sanae::Solution::GetLastMagnetization() after propagation.
 */
namespace Sanae
{
	Eigen::MatrixXcd R1_inversion_recovery(
		double R1 = 5.0,
		double R2 = 7.0,
		double CSppm = 1.0,
		double tau = 0.1,
		char nucleus = 'N',
		const std::string &output_filename = "r1_inversion_recovery_output.txt",
		bool verbose = false);
}