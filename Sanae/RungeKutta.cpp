#include "RungeKutta.h"

// @todo Consider adaptive step size if stability issues arise.

void RungeKutta(const Eigen::MatrixXcd &A, Eigen::MatrixXcd &M, double dt, int num_steps)
{
	// Loop through timesteps and solve using 4th-order Runge-Kutta method
	for (int i = 1; i < num_steps; ++i)
	{
		Eigen::MatrixXcd k1 = A * M.col(i - 1);					  // k1 = A *  M(i-1)
		Eigen::MatrixXcd k2 = A * (M.col(i - 1) + 0.5 * dt * k1); // k2 = A * (M(i-1) + dt/2 * k1)
		Eigen::MatrixXcd k3 = A * (M.col(i - 1) + 0.5 * dt * k2); // k3 = A * (M(i-1) + dt/2 * k2)
		Eigen::MatrixXcd k4 = A * (M.col(i - 1) + dt * k3);		  // k4 = A * (M(i-1) + dt * k3)

		// Update M using the weighted sum of the Runge-Kutta terms
		M.col(i) = M.col(i - 1) + (dt / 6.0) * (k1 + 2.0 * k2 + 2.0 * k3 + k4);
	}
}
