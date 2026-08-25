// g++ -std=c++17 -Wall -o test_env ../Environment.cpp test_env.cpp

#include "../Environment.h"
#include <iostream>
#include <cassert>
#include <cmath> // for std::abs

int main()
{
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

    // w_mhz = -gamma * B0 = -42.58 * 16.4 (negative, since gamma_H > 0)
    double expected_sfo1_h = -42.58 * 16.4 * 1e6;
    assert(std::abs(env.GetSFO1() - expected_sfo1_h) < 1e-6);
    std::cout << "SFO1 Math Check: PASSED (" << expected_sfo1_h << " Hz)\n";

    // 3. Change B0 to 20.0 Tesla
    env.SetB0(20.0);
    std::cout << "\n--- After SetB0(20.0) ---\n";
    env.DisplayInfo();

    double expected_sfo1_b0 = -42.58 * 20.0 * 1e6;
    assert(std::abs(env.GetSFO1() - expected_sfo1_b0) < 1e-6);
    std::cout << "B0 Update Check: PASSED (" << expected_sfo1_b0 << " Hz)\n";

    // 4. Switch to Nitrogen (gamma = -4.316, now NEGATIVE per updated code)
    env.SetGamma('N');
    std::cout << "\n--- After SetGamma('N') with B0=20.0 ---\n";
    env.DisplayInfo();

    // gamma_N is now -4.316, so w_mhz = -gamma_N * B0 = -(-4.316)*20.0 = +4.316*20.0
    // NOTE: sign flips positive here, unlike H/C/F, because gamma_N itself is negative.
    double expected_sfo1_n = -(-4.316) * 20.0 * 1e6;
    assert(std::abs(env.GetSFO1() - expected_sfo1_n) < 1e-6);
    std::cout << "Nitrogen Switch Check: PASSED (" << expected_sfo1_n << " Hz)\n";
    assert(env.GetSFO1() > 0.0 && "15N SFO1 should be positive given negative gamma");

    // 5. Switch to Carbon (gamma = 10.705, positive like H)
    env.SetGamma('C');
    std::cout << "\n--- After SetGamma('C') with B0=20.0 ---\n";
    env.DisplayInfo();

    double expected_sfo1_c = -10.705 * 20.0 * 1e6;
    assert(std::abs(env.GetSFO1() - expected_sfo1_c) < 1e-6);
    std::cout << "Carbon Switch Check: PASSED (" << expected_sfo1_c << " Hz)\n";

    // 6. Switch to Fluorine (gamma = 40.05, positive)
    env.SetGamma('F');
    std::cout << "\n--- After SetGamma('F') with B0=20.0 ---\n";
    env.DisplayInfo();

    double expected_sfo1_f = -40.05 * 20.0 * 1e6;
    assert(std::abs(env.GetSFO1() - expected_sfo1_f) < 1e-6);
    std::cout << "Fluorine Switch Check: PASSED (" << expected_sfo1_f << " Hz)\n";

    // 7. Test invalid nuclide exception
    std::cout << "\n--- Testing Invalid Nuclide ---\n";
    try
    {
        env.SetGamma('X');
        std::cout << "ERROR: SetGamma('X') should have thrown!\n";
        return 1;
    }
    catch (const std::invalid_argument &e)
    {
        std::cout << "Caught expected exception: " << e.what() << "\n";
    }

    std::cout << "\nAll tests passed! Environment class is solid.\n";
    return 0;
}