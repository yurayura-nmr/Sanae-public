#pragma once

#include <Eigen/Dense>

/**
 * @file numerical_solvers.h
 * @brief Generic 3-variable Newton-Raphson solver with a
 *        positivity-preserving backtracking line search.
 *
 * Intended for systems of 3 nonlinear equations in 3 unknowns where the
 * unknowns represent physically positive quantities (e.g. populations,
 * concentrations, fractional state occupancies) that must never be
 * allowed to go negative or zero during the iterative solve, even as
 * intermediate (non-converged) values.
 *
 * The caller supplies a function that computes both the residual vector
 * F(x) and its Jacobian J(x) at a given point x. The solver iterates
 * x_new = x + alpha * dx, where dx solves J * dx = -F, and alpha is
 * halved (backtracking line search) until x_new has strictly positive
 * components, or the line search budget is exhausted.
 */

namespace NumericalUtils
{

	/**
	 * @brief Solve F(x) = 0 for a 3-variable system using Newton-Raphson
	 *        with a positivity-preserving backtracking line search.
	 *
	 * @tparam ResidualFunc
	 *     Callable with signature:
	 *         void f(const Eigen::Vector3d& x, Eigen::Vector3d& F, Eigen::Matrix3d& J)
	 *     F and J are output parameters: on each call, the caller must
	 *     fill in the residual vector F(x) and the Jacobian J(x) = dF/dx
	 *     evaluated at the current x.
	 *
	 * @param residualAndJacobian
	 *     User-supplied function computing F(x) and J(x). See ResidualFunc.
	 * @param x
	 *     Initial guess on input; updated in-place with the solution
	 *     (converged or not) on return.
	 * @param tol
	 *     Convergence tolerance on the residual norm ||F(x)||. The solver
	 *     returns true as soon as ||F(x)|| < tol.
	 * @param maxIter
	 *     Maximum number of Newton iterations before giving up.
	 *
	 * @return
	 *     true if ||F(x)|| < tol was reached within maxIter iterations;
	 *     false otherwise. On false, x holds the last iterate reached,
	 *     which may still be a reasonable (non-converged) estimate.
	 *
	 * @note
	 *     The line search only accepts a step if every component of
	 *     x_new is strictly positive (x_new.minCoeff() > 0.0). If no
	 *     such step is found within maxLS halvings of alpha, x is left
	 *     unchanged for that iteration.
	 */
	template <typename ResidualFunc>
	bool newtonRaphson(ResidualFunc residualAndJacobian,
		Eigen::Vector3d& x,
		double tol = 1e-12,
		int maxIter = 50)
	{
		using Eigen::Matrix3d;
		using Eigen::Vector3d;

		constexpr int maxLS = 10; // Line search max iterations

		Vector3d F;
		Matrix3d J;

		for (int iter = 0; iter < maxIter; ++iter)
		{
			// user must provide residual and Jacobian in this signature:
			//     void f(const Vector3d& x, Vector3d& F, Matrix3d& J)
			residualAndJacobian(x, F, J);

			if (F.norm() < tol)
				return true;

			Vector3d dx = J.partialPivLu().solve(-F);

			// Simple positivity-preserving backtracking line search
			double alpha = 1.0;
			for (int ls = 0; ls < maxLS; ++ls)
			{
				Vector3d x_new = x + alpha * dx;
				if (x_new.minCoeff() > 0.0)
				{
					x = x_new;
					break;
				}
				alpha *= 0.5;
			}
		}

		// One final check: the last update inside the loop is never
		// re-evaluated against tol before the loop exits via the
		// iteration count, so check it here before giving up.
		residualAndJacobian(x, F, J);
		return F.norm() < tol;
	}

} // namespace NumericalUtils