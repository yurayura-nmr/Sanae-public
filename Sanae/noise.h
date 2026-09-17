#pragma once

#include <Eigen/Dense>

/**
 * @file noise.h
 * @brief Adds Gaussian white noise to a simulated NMR FID matrix.
 *
 * Models receiver noise in a real NMR experiment.
 * Here, noise is applied after numerical integration of the
 * Bloch-McConnell equations — it is purely an instrumental/observational
 * effect, not part of the actual spin dynamics.
 *
 * Noise is drawn from N(0, sigma) independently for real and
 * imaginary components of each point
 * (i.e, uncorrelated receiver noise on Mx and My channels).
 *
 * @note col=0 is never modified — it holds M(t=0), the initial
 *       condition set before integration, which should be exact.
 */

/**
 * @brief Add Gaussian white noise to an NMR FID matrix in-place.
 *
 * @param M      FID matrix (rows = states, cols = time points).
 *               Modified in-place. Pass as reference from Solution.cpp.
 * @param sigma  Standard deviation of noise (same units as magnetization).
 *               Typical values relative to M0 = 1.0:
 *               0.001 = excellent S/N
 *               0.01  = good S/N  (realistic high-field NMR)
 *               0.05  = poor S/N
 *               0.0   = no noise added (returns immediately)
 *               Negative values are also treated as "no noise" and
 *               return immediately without modifying M.
 * @param seed   Random seed for reproducibility. Default -1 uses
 *               std::random_device (non-reproducible).
 *
 * @example
 *   // After RungeKutta, before file write:
 *   AddGaussianNoise(M, noise_level_, -1);
 */
void AddGaussianNoise(
    Eigen::MatrixXcd &M,
    double sigma,
    int seed = -1);
