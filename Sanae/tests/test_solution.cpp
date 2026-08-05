// g++ -std=c++17 -Wall -I/path/to/eigen -o test_solution ../Solution.cpp ../State.cpp ../Environment.cpp ../RungeKutta.cpp ../fileio.cpp test_solution.cpp

#include "../Solution.h"
#include "../State.h"
#include "../Environment.h"

#include <iostream>
#include <fstream>
#include <cassert>
#include <cmath>
#include <string>

using namespace Sanae;

// Helper to check if a file exists and is non-empty
bool fileExistsAndNonEmpty(const std::string& filename) {
    std::ifstream f(filename);
    return f.good() && f.peek() != std::ifstream::traits_type::eof();
}

int main() {
    std::cout << "=== Integration Test: Solution + State + Environment ===\n";

    // --- 1. Create Environment (700 MHz proton) ---
    Environment env;
    env.SetGamma('H');  // proton, B0=16.4 T by default

    // --- 2. Create two states (A and B) ---
    State A, B;
    A.SetID(0);
    A.SetName("State A");
    A.SetInitMag(1.0);      // initial Mz = 1
    A.SetR2(5.0);           // R2 = 5 s^-1
    A.SetR1(3.0);           // not used in transverse simulation but set anyway
    A.SetCS_ppm(0.0, env);  // no chemical shift offset

    B.SetID(1);
    B.SetName("State B");
    B.SetInitMag(0.5);
    B.SetR2(10.0);
    B.SetR1(4.0);
    B.SetCS_ppm(2.0, env);  // 2 ppm offset

    // --- 3. Create Solution (2‑state) and set numerical parameters ---
    Solution sol(2);
    sol.SetDT(0.0001);          // 100 µs per step
    sol.SetSteps(1000);         // total time = 0.1 s
    sol.SetPF(10);              // print every 10 steps
    sol.SetMethod(1);           // Runge‑Kutta (1) – Euler (0) not implemented

    // --- 4. Set exchange matrix: no exchange (only relaxation) ---
    sol.SetNoExchange(A, B);

    // --- 5. Set output filename ---
    const std::string outFile = "integration_test_fid.txt";
    sol.SetOutputFilename(outFile);

    // --- 6. Run the simulation ---
    sol.SetupTwoState(A, B);

    // --- 7. Verify GetSteps() and GetK() ---
    assert(sol.GetSteps() == 1000);
    Eigen::MatrixXcd K = sol.GetK();
    assert(K.rows() == 2 && K.cols() == 2);
    // No exchange => K should be zero (but diagonal was set to -0 and off-diagonals = 0)
    // Actually SetNoExchange sets cross terms to 0 and diagonal to -0, so K is zero matrix.
    assert(K.isZero(1e-12));

    std::cout << "Getter checks passed.\n";

    // --- 8. Check output file ---
    assert(fileExistsAndNonEmpty(outFile));

    std::ifstream file(outFile);
    std::string line;

    // 8a. Verify header lines
    std::getline(file, line);
    assert(line.find("# SANAE simulation output") != std::string::npos);
    std::getline(file, line);
    assert(line.find("# dt = 0.0001") != std::string::npos);
    std::getline(file, line);
    assert(line.find("# steps = 1000") != std::string::npos);
    std::getline(file, line);
    assert(line.find("# print_every = 10") != std::string::npos);
    std::cout << "Header verified.\n";

    // 8b. Read first data line (should be t=0)
    std::getline(file, line);
    double real0, imag0;
    std::stringstream ss(line);
    ss >> real0 >> imag0;
    // At t=0, M+ = M_A(0) + M_B(0) = 1.0 + 0.5 = 1.5
    // Since CS = 0 for A, 2 ppm for B, at t=0 both are real, so imag=0
    const double tolerance = 1e-6;
    assert(std::abs(real0 - 1.5) < tolerance);
    assert(std::abs(imag0) < tolerance);
    std::cout << "Initial value verified (M+ = 1.5).\n";

    // 8c. Read the last data line (should have decayed)
    std::string lastLine;
    while (std::getline(file, lastLine)) {
        // just read to end
    }
    double realLast, imagLast;
    std::stringstream ssLast(lastLine);
    ssLast >> realLast >> imagLast;
    // Because of relaxation, magnitude should be smaller than 1.5
    double magLast = std::sqrt(realLast*realLast + imagLast*imagLast);
    assert(magLast < 1.5 - 0.1);  // should have decayed significantly
    std::cout << "Decay verified (|M+| at end = " << magLast << ").\n";

    // --- 9. Test SetupOneState as well (optional) ---
    std::cout << "\nTesting SetupOneState...\n";
    State single;
    single.SetID(0);
    single.SetName("Single spin");
    single.SetInitMag(1.0);
    single.SetR2(5.0);
    single.SetCS_ppm(0.0, env);

    Solution sol1(1);
    sol1.SetDT(0.0001);
    sol1.SetSteps(500);
    sol1.SetPF(5);
    sol1.SetMethod(1);
    sol1.SetOutputFilename("single_spin_test.txt");
    sol1.SetupOneState(single);

    assert(fileExistsAndNonEmpty("single_spin_test.txt"));
    // Check first value
    std::ifstream fSingle("single_spin_test.txt");
    // skip header (4 lines)
    for (int i = 0; i < 4; ++i) std::getline(fSingle, line);
    std::getline(fSingle, line);
    double r0, i0;
    std::stringstream ssSingle(line);
    ssSingle >> r0 >> i0;
    assert(std::abs(r0 - 1.0) < tolerance);
    assert(std::abs(i0) < tolerance);
    std::cout << "SetupOneState passed.\n";

    std::cout << "\nAll integration tests passed!\n";
    return 0;
}