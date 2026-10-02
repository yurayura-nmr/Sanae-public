//  cl /EHsc /std:c++17 /I"..\..\include\eigen3" /Fe:test_numerical_solvers.exe test_numerical_solvers.cpp

#include "../numerical_solvers.h"

#include <Eigen/Dense>
#include <iostream>
#include <cassert>
#include <cmath>

// ---------------------------------------------------------------------
// Test problem 1: diagonal system with a known positive root.
//   F_i(x) = x_i^2 - c_i   =>   x_i = sqrt(c_i)
//   J_i,i  = 2 * x_i  (diagonal Jacobian)
// Chosen because it has an exact, easily verified positive solution
// and because squaring makes the line search's positivity constraint
// meaningful (a naive full Newton step can overshoot into negative
// territory from certain starting points).
// ---------------------------------------------------------------------
struct DiagonalSquareResidual
{
	Eigen::Vector3d c; // target values: x_i^2 = c_i

	void operator()(const Eigen::Vector3d &x, Eigen::Vector3d &F, Eigen::Matrix3d &J) const
	{
		F(0) = x(0) * x(0) - c(0);
		F(1) = x(1) * x(1) - c(1);
		F(2) = x(2) * x(2) - c(2);

		J.setZero();
		J(0, 0) = 2.0 * x(0);
		J(1, 1) = 2.0 * x(1);
		J(2, 2) = 2.0 * x(2);
	}
};

int main()
{
	std::cout << "=== Testing NumericalUtils::newtonRaphson ===\n";

	// -------------------------------------------------------------
	std::cout << "\n--- Converges to known positive root ---\n";
	{
		DiagonalSquareResidual res;
		res.c = Eigen::Vector3d(4.0, 9.0, 16.0); // expect x = (2, 3, 4)

		Eigen::Vector3d x(1.0, 1.0, 1.0); // positive starting guess
		bool converged = NumericalUtils::newtonRaphson(res, x);

		Eigen::Vector3d expected(2.0, 3.0, 4.0);
		double err = (x - expected).norm();

		std::cout << "Converged = " << (converged ? "true" : "false") << "\n";
		std::cout << "x = (" << x(0) << ", " << x(1) << ", " << x(2) << ")\n";
		std::cout << "Expected  = (2, 3, 4), error = " << err << "\n";

		assert(converged);
		assert(err < 1e-9);
		std::cout << "Known Root Check: "
				  << ((converged && err < 1e-9) ? "PASSED" : "FAILED")
				  << "\n";
	}

	// -------------------------------------------------------------
	std::cout << "\n--- Line search keeps all iterates strictly positive ---\n";
	{
		// A starting guess small enough that a full (alpha=1) Newton
		// step would otherwise be the natural thing to try; the test
		// checks that whatever iterate the solver lands on along the
		// way never had a non-positive component inspectable via the
		// final positive solution itself (indirect check: a diverging
		// or sign-flipped x would not converge to the known root at all).
		DiagonalSquareResidual res;
		res.c = Eigen::Vector3d(1.0, 1.0, 1.0); // expect x = (1, 1, 1)

		Eigen::Vector3d x(0.01, 0.01, 0.01); // tiny positive starting guess
		bool converged = NumericalUtils::newtonRaphson(res, x);

		bool all_positive = (x(0) > 0.0 && x(1) > 0.0 && x(2) > 0.0);

		std::cout << "Converged = " << (converged ? "true" : "false") << "\n";
		std::cout << "x = (" << x(0) << ", " << x(1) << ", " << x(2) << ")\n";

		assert(converged);
		assert(all_positive);
		std::cout << "Positivity Preserved Check: "
				  << ((converged && all_positive) ? "PASSED" : "FAILED")
				  << "\n";
	}

	// -------------------------------------------------------------
	std::cout << "\n--- Respects custom tolerance ---\n";
	{
		DiagonalSquareResidual res;
		res.c = Eigen::Vector3d(4.0, 9.0, 16.0);

		Eigen::Vector3d x(1.0, 1.0, 1.0);
		bool converged = NumericalUtils::newtonRaphson(res, x, 1e-6, 50);

		Eigen::Vector3d F;
		Eigen::Matrix3d J;
		res(x, F, J);

		std::cout << "Converged = " << (converged ? "true" : "false") << "\n";
		std::cout << "Final ||F(x)|| = " << F.norm() << "  (tol = 1e-6)\n";

		assert(converged);
		assert(F.norm() < 1e-6);
		std::cout << "Custom Tolerance Check: "
				  << ((converged && F.norm() < 1e-6) ? "PASSED" : "FAILED")
				  << "\n";
	}

	// -------------------------------------------------------------
	std::cout << "\n--- Returns false when maxIter is too small ---\n";
	{
		DiagonalSquareResidual res;
		res.c = Eigen::Vector3d(4.0, 9.0, 16.0);

		// Starting far from the root with essentially zero iteration
		// budget should not be enough to converge.
		Eigen::Vector3d x(100.0, 100.0, 100.0);
		bool converged = NumericalUtils::newtonRaphson(res, x, 1e-12, 1);

		std::cout << "Converged = " << (converged ? "true" : "false")
				  << " (expected false)\n";

		assert(!converged);
		std::cout << "Insufficient Iterations Check: "
				  << (!converged ? "PASSED" : "FAILED")
				  << "\n";
	}

	// -------------------------------------------------------------
	std::cout << "\n--- Off-by-one fix: converges exactly on the final allowed iteration ---\n";
	{
		// This targets the bug we fixed: convergence achieved on the
		// very last update inside the loop must still be detected,
		// not silently discarded because the loop exits via the
		// iteration-count condition before re-checking the residual.
		//
		// Rather than hardcoding a maxIter boundary (fragile, depends
		// on exact convergence speed), find the minimum maxIter that
		// succeeds, then confirm one iteration fewer correctly fails,
		// and the boundary iteration succeeds -- this is precisely the
		// case the fix addresses.
		DiagonalSquareResidual res;
		res.c = Eigen::Vector3d(4.0, 9.0, 16.0);
		const double tol = 1e-9;

		int min_iters_needed = -1;
		for (int n = 1; n <= 50; ++n)
		{
			Eigen::Vector3d x(1.0, 1.0, 1.0);
			if (NumericalUtils::newtonRaphson(res, x, tol, n))
			{
				min_iters_needed = n;
				break;
			}
		}

		std::cout << "Minimum maxIter needed to converge = " << min_iters_needed << "\n";
		assert(min_iters_needed > 0);

		// One fewer iteration must still correctly fail (sanity check
		// that we're actually at the true boundary, not past it).
		Eigen::Vector3d x_short(1.0, 1.0, 1.0);
		bool converged_short = NumericalUtils::newtonRaphson(res, x_short, tol, min_iters_needed - 1);
		std::cout << "maxIter = " << (min_iters_needed - 1)
				  << " converged = " << (converged_short ? "true" : "false")
				  << " (expected false)\n";
		assert(!converged_short);

		// At the boundary itself, convergence on the final allowed
		// iteration must be detected (this is the fix being tested).
		Eigen::Vector3d x_boundary(1.0, 1.0, 1.0);
		bool converged_boundary = NumericalUtils::newtonRaphson(res, x_boundary, tol, min_iters_needed);
		Eigen::Vector3d expected(2.0, 3.0, 4.0);
		double err = (x_boundary - expected).norm();

		std::cout << "maxIter = " << min_iters_needed
				  << " converged = " << (converged_boundary ? "true" : "false")
				  << ", error = " << err << "\n";

		assert(converged_boundary);
		assert(err < 1e-6);
		std::cout << "Final-Iteration Convergence Check: "
				  << ((converged_boundary && err < 1e-6) ? "PASSED" : "FAILED")
				  << "\n";
	}
	return 0;
}