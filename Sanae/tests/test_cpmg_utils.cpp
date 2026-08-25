// g++ -std=c++17 -Wall -Wextra -I../../include/eigen3/ -o test_cpmg ../cpmg_utils.cpp ../fileio.cpp test_cpmg_utils.cpp

#include "../cpmg_utils.h"
#include <iostream>
#include <sstream>
#include <fstream>
#include <cassert>
#include <cmath>
#include <cstdio>
#include <vector>
#include <string>

namespace
{

    struct CpmgCurve
    {
        std::vector<double> nu;
        std::vector<double> R2eff;
    };

    CpmgCurve ReadCpmgFile(const std::string &filename)
    {
        CpmgCurve curve;
        std::ifstream in(filename);
        assert(in.is_open() && "output file should exist and be readable");

        std::string line;
        while (std::getline(in, line))
        {
            if (line.empty() || line[0] == '#')
                continue;
            std::istringstream iss(line);
            double nu, r2;
            if (iss >> nu >> r2)
            {
                curve.nu.push_back(nu);
                curve.R2eff.push_back(r2);
            }
        }
        return curve;
    }

} // namespace

int main()
{
    using namespace DispersionUtils;

    std::cout << "=== Testing DispersionUtils (cpmg_utils.cpp) - 15N ===\n";

    // NOTE: "gamma" here is the gyromagnetic ratio RELATIVE TO 1H
    // (dimensionless, ~0.101 for 15N), and B0_MHz is the PROTON
    // spectrometer frequency — not Environment.cpp's MHz/T convention.
    const double gamma_15N_relative = 0.1010;
    const double pA = 0.90;
    const double kex = 2000.0;
    const double dw_ppm = 0.5;
    const double B0_MHz = 700.0;
    const double R20 = 10.0;
    const double nu_min = 40.0;
    const double nu_max = 2000.0;
    const size_t n_points = 15;

    // ------------------------------------------------------------------
    // Carver-Richards
    // ------------------------------------------------------------------
    {
        const std::string fname = "test_cr_15N.dat";
        CarverRichards(pA, kex, dw_ppm, B0_MHz, gamma_15N_relative, R20,
                       nu_min, nu_max, n_points, "CPMG test: 15N", fname);

        CpmgCurve curve = ReadCpmgFile(fname);
        assert(curve.R2eff.size() == n_points);

        // Realistic range check: should stay near R20, not blow up
        for (double r2 : curve.R2eff)
            assert(r2 > R20 - 1e-6 && r2 < R20 + 5.0);

        assert(curve.R2eff.front() > curve.R2eff.back());
        assert(curve.R2eff.back() >= R20 - 1e-6);

        std::cout << "CarverRichards[15N]: R2eff " << curve.R2eff.front()
                  << " -> " << curve.R2eff.back() << " over " << nu_min
                  << "-" << nu_max << " Hz: PASSED\n";
        std::remove(fname.c_str());
    }

    // No-dispersion control: dw=0 should hold R2eff == R20 everywhere
    {
        const std::string fname = "test_cr_15N_nodw.dat";
        CarverRichards(pA, kex, /*dw_ppm=*/0.0, B0_MHz, gamma_15N_relative, R20,
                       nu_min, nu_max, 5, "no dispersion", fname);

        CpmgCurve curve = ReadCpmgFile(fname);
        for (double r2 : curve.R2eff)
            assert(std::abs(r2 - R20) < 1e-6); // full precision now, from fileio

        std::cout << "CarverRichards[15N]: dw=0 recovers R20 everywhere: PASSED\n";
        std::remove(fname.c_str());
    }

    // ------------------------------------------------------------------
    // Meiboom (fast-exchange approximation)
    // ------------------------------------------------------------------
    {
        const std::string fname = "test_meiboom_15N.dat";
        Meiboom(pA, kex, dw_ppm, B0_MHz, gamma_15N_relative, R20,
                nu_min, nu_max, n_points, "Meiboom test: 15N", fname);

        CpmgCurve curve = ReadCpmgFile(fname);
        assert(curve.R2eff.size() == n_points);
        for (double r2 : curve.R2eff)
        {
            assert(!std::isnan(r2));
            assert(r2 > R20 - 1e-6 && r2 < R20 + 5.0); // realistic range
        }
        assert(curve.R2eff.front() >= curve.R2eff.back() - 1e-9);

        std::cout << "Meiboom[15N]: R2eff " << curve.R2eff.front()
                  << " -> " << curve.R2eff.back() << " over " << nu_min
                  << "-" << nu_max << " Hz: PASSED\n";
        std::remove(fname.c_str());
    }

    std::cout << "\nAll tests passed!\n";
    return 0;
}