// g++ -std=c++17 -Wall -o test_env ../Environment.cpp test_env.cpp

#include "../Environment.h"
#include <iostream>
#include <cassert>
#include <cmath>   // for std::abs

int main() {
    using namespace Sanae;

    std::cout << "=== Testing Sanae::Environment ===\n";

    // 1. Default constructor
    Environment env;
    std::cout << "Default SFO1: " << env.GetSFO1() << " Hz (Expected: 0)\n";
    assert(env.GetSFO1() == 0.0);

    // 2. Set Hydrogen (gamma = 42.58, default B0 = 16.4)
    env.SetGamma('H');
    std::cout << "\n--- After SetGamma('H') ---\n";
    env.DisplayInfo();
    
    // Mathematically: w_mhz = -gamma * B0 = -42.58 * 16.4 = -698.312 MHz
    // SFO1 = w_mhz * 1e6 = -698,312,000 Hz
    double expected_sfo1_h = -42.58 * 16.4 * 1e6;
    assert(std::abs(env.GetSFO1() - expected_sfo1_h) < 1e-6);
    std::cout << "SFO1 Math Check: PASSED (" << expected_sfo1_h << " Hz)\n";

    // 3. Change B0 to 20.0 Tesla
    env.SetB0(20.0);
    std::cout << "\n--- After SetB0(20.0) ---\n";
    env.DisplayInfo();

    // Now w_mhz = -42.58 * 20.0 = -851.6 MHz
    double expected_sfo1_b0 = -42.58 * 20.0 * 1e6;
    assert(std::abs(env.GetSFO1() - expected_sfo1_b0) < 1e-6);
    std::cout << "B0 Update Check: PASSED (" << expected_sfo1_b0 << " Hz)\n";

    // 4. Switch to Nitrogen (gamma = 4.31)
    env.SetGamma('N');
    std::cout << "\n--- After SetGamma('N') with B0=20.0 ---\n";
    env.DisplayInfo();

    // Now w_mhz = -4.31 * 20.0 = -86.2 MHz
    double expected_sfo1_n = -4.31 * 20.0 * 1e6;
    assert(std::abs(env.GetSFO1() - expected_sfo1_n) < 1e-6);
    std::cout << "Nitrogen Switch Check: PASSED (" << expected_sfo1_n << " Hz)\n";

    // 5. Test invalid nuclide exception
    std::cout << "\n--- Testing Invalid Nuclide ---\n";
    try {
        env.SetGamma('X');
        std::cout << "ERROR: SetGamma('X') should have thrown!\n";
        return 1;
    } catch (const std::invalid_argument& e) {
        std::cout << "Caught expected exception: " << e.what() << "\n";
    }

    std::cout << "\nAll tests passed! Environment class is solid.\n";
    return 0;
}
