// g++ -std=c++17 -Wall -o test_state ../State.cpp ../Environment.cpp test_state.cpp

#include "../State.h"
#include "../Environment.h"
#include <iostream>
#include <cassert>
#include <cmath>

using namespace Sanae;

int main()
{
    std::cout << "=== Testing State ===\n";

    // 1. Create an Environment (default: B0=16.4 T, gamma_H=42.58 MHz/T)
    Environment env;
    env.SetGamma('H');
    std::cout << "SFO1 from env: " << env.GetSFO1() << std::endl;

    // env.GetSFO1() returns -698.312 MHz = -6.98312e8 Hz (since gamma is positive and B0 positive)
    double expected_SFO1 = -42.58 * 16.4 * 1e6; // -6.98312e8 Hz

    // 2. Create a State with default values
    State state;

    // 3. Set properties
    state.SetID(0);
    state.SetName("Proton");
    state.SetInitMag(1.0);
    state.SetR1(3.0);
    state.SetR2(7.0);

    // 4. Set chemical shift: 1 ppm
    const double cs_ppm = 1.0;
    state.SetCS_ppm(cs_ppm, env);

    // 5. Verify getters
    assert(state.GetID() == 0);
    assert(state.GetInitMag() == 1.0);
    assert(state.GetR1() == 3.0);
    assert(state.GetR2() == 7.0);

    double expected_rad_s = cs_ppm * expected_SFO1 / 1e6 * 2.0 * 3.14159265358979323846;
    double actual_rad_s = state.GetCS_rad_s();
    std::cout << "expected_rad_s = " << expected_rad_s << std::endl;
    std::cout << "actual_rad_s   = " << actual_rad_s << std::endl;
    std::cout << "diff = " << std::abs(actual_rad_s - expected_rad_s) << std::endl;
    assert(std::abs(actual_rad_s - expected_rad_s) < 1e-6);

    std::cout << "All getters passed.\n";

    // 6. Display info (manual inspection)
    state.DisplayInfo();

    // 7. Test another state with different parameters
    State state2;
    state2.SetID(1);
    state2.SetName("Nitrogen");
    state2.SetInitMag(0.5);
    state2.SetR2(10.0);
    state2.SetCS_ppm(2.5, env);

    state2.DisplayInfo();

    // 8. Test DisplayInfo_R1 if you uncommented it
    // state2.DisplayInfo_R1(0.01);

    std::cout << "\nAll tests passed. Check the printed output.\n";
    return 0;
}
