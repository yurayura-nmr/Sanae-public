// g++ -std=c++17 -Wall -I../../include/eigen3/ -o test_bloch_relax ../Solution.cpp ../State.cpp ../Environment.cpp ../RungeKutta.cpp ../fileio.cpp ../bloch_relax.cpp test_bloch_relax.cpp

#include "../bloch_relax.h"
#include <iostream>
#include <fstream>
#include <cassert>
#include <sstream>   // for std::stringstream
#include <cmath>     // for std::sqrt, std::abs
#include <string>

// Helper to check if file exists and is non‑empty
bool fileExistsAndNonEmpty(const std::string& filename) {
    std::ifstream f(filename);
    return f.good() && f.peek() != std::ifstream::traits_type::eof();
}

int main() {
    std::cout << "=== Testing bloch_relax ===\n";

    // 1. Run the demo with default parameters
    // This will generate "magnetization_evolution.txt" (default filename)
    bloch_relax();

    // 2. Verify the output file was created
    std::string filename = "magnetization_evolution.txt";
    assert(fileExistsAndNonEmpty(filename));
    std::cout << "File created successfully.\n";

    // 3. Verify the header
    std::ifstream file(filename);
    std::string line;

    std::getline(file, line);
    assert(line.find("# SANAE simulation output") != std::string::npos);

    std::getline(file, line);
    assert(line.find("# dt = ") != std::string::npos);

    std::getline(file, line);
    assert(line.find("# steps = ") != std::string::npos);

    std::getline(file, line);
    assert(line.find("# print_every = ") != std::string::npos);

    std::cout << "Header verified.\n";

    // 4. Read the first data point (t=0)
    std::getline(file, line);
    double real0, imag0;
    std::stringstream ss(line);
    ss >> real0 >> imag0;

    // At t=0, M+ = M0_A + M0_B = 0.5 + 0.5 = 1.0
    const double tolerance = 1e-6;
    assert(std::abs(real0 - 1.0) < tolerance);
    assert(std::abs(imag0) < tolerance);
    std::cout << "Initial value verified: M+ = 1.0\n";

    // 5. Check that decay occurs (read last line)
    std::string lastLine;
    while (std::getline(file, lastLine)) {
        // just loop to the end
    }
    double realLast, imagLast;
    std::stringstream ssLast(lastLine);
    ssLast >> realLast >> imagLast;
    double magLast = std::sqrt(realLast*realLast + imagLast*imagLast);
    // R2 values are 20 and 10, so after 0.2 s (20000 * 10µs) magnetization should be tiny
    assert(magLast < 0.1);
    std::cout << "Decay verified: |M+| at end = " << magLast << "\n";

    std::cout << "\nAll tests passed. bloch_relax works correctly.\n";
    return 0;
}
