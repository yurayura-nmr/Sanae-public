#include "../RungeKutta.h"
#include <iostream>
#include <Eigen/Core>
#include <complex>
#include <cmath>
#include <cassert>

using namespace Eigen;

// Helper to compare matrices within tolerance
bool is_close(const MatrixXcd& a, const MatrixXcd& b, double tol = 1e-9) {
    return (a - b).cwiseAbs().maxCoeff() < tol;
}

int main() {
    std::cout << "=== Testing RungeKutta ===\n";

    // Parameters
    const int n = 3;               // state dimension
    const int num_steps = 1000;
    const double dt = 0.001;
    const double total_time = num_steps * dt;

    // ---------- Test 1: Zero matrix (M should stay constant) ----------
    {
        MatrixXcd A_zero = MatrixXcd::Zero(n, n);
        MatrixXcd M(n, num_steps + 1);
        M.col(0) = VectorXcd::Random(n);  // random initial state

        RungeKutta(A_zero, M, dt, num_steps + 1);  // includes initial col

        // All columns should equal initial
        bool pass = true;
        for (int i = 1; i <= num_steps; ++i) {
            if (!is_close(M.col(i), M.col(0))) {
                pass = false;
                break;
            }
        }
        std::cout << "Test 1 (Zero A): " << (pass ? "PASSED" : "FAILED") << "\n";
        assert(pass);
    }

    // ---------- Test 2: Decay A = -λ I (M should decay) ----------
    {
        const double lambda = 2.0;
        MatrixXcd A_decay = -lambda * MatrixXcd::Identity(n, n);
        MatrixXcd M(n, num_steps + 1);
        M.col(0) = VectorXcd::Random(n);

        RungeKutta(A_decay, M, dt, num_steps + 1);

        // Analytical: M(t) = exp(-λ t) * M(0)
        double t = 0.0;
        bool pass = true;
        for (int i = 0; i <= num_steps; ++i) {
            t = i * dt;
            VectorXcd expected = std::exp(-lambda * t) * M.col(0);
            if (!is_close(M.col(i), expected, 1e-6)) {
                pass = false;
                break;
            }
        }
        std::cout << "Test 2 (Decay): " << (pass ? "PASSED" : "FAILED") << "\n";
        assert(pass);
    }

    // ---------- Test 3: Oscillation A = i * ω (skew-hermitian) ----------
    {
        const double omega = 1.5;
        MatrixXcd A_osc(n, n);
        A_osc.setIdentity();
        A_osc *= std::complex<double>(0.0, omega);  // i*ω

        MatrixXcd M(n, num_steps + 1);
        M.col(0) = VectorXcd::Random(n);

        RungeKutta(A_osc, M, dt, num_steps + 1);

        // Analytical: M(t) = exp(i ω t) * M(0)
        double t = 0.0;
        bool pass = true;
        for (int i = 0; i <= num_steps; ++i) {
            t = i * dt;
            VectorXcd expected = std::exp(std::complex<double>(0.0, omega * t)) * M.col(0);
            if (!is_close(M.col(i), expected, 1e-6)) {
                pass = false;
                break;
            }
        }
        std::cout << "Test 3 (Oscillation): " << (pass ? "PASSED" : "FAILED") << "\n";
        assert(pass);
    }

    std::cout << "\n All tests passed!\n";
    return 0;
}
