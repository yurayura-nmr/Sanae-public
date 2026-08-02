/**
 * @file RungeKutta.h
 * @brief Fourth-order Runge–Kutta integrator for complex matrices.
 *
 * Performs numerical integration of the system of linear differential equations:
 * \f$ \frac{dM}{dt} = A \cdot M \f$
 * using the classic 4th-order Runge–Kutta method.
 *
 * @param A         System matrix defining the evolution (Eigen::MatrixXcd).
 * @param M         State vector/matrix to be updated (Eigen::MatrixXcd).
 * @param dt        Time step for integration.
 * @param num_steps Number of integration steps to perform.
 *
 * @note This is used by Sanae for creating the actual FID, i.e., we are using
 *       numerical solutions, never setting up any analytical approximation.
 *       It is not elegant, it is brute force, but it helps us check all
 *       scenarios without keeping track of which approximation is valid to what
 *       degree under what condition.
 */

#pragma once

#include <Eigen/Core>

void RungeKutta(const Eigen::MatrixXcd &A, Eigen::MatrixXcd &M, double dt, int num_steps);