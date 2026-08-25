#include "fileio.h"

#include <fstream>
#include <iostream>
#include <limits>
#include <vector>

namespace Sanae
{

    bool writeMagnetizationToFile(const Eigen::MatrixXcd &M_plus,
                                  const std::string &filename,
                                  double dt,
                                  int num_steps,
                                  int printFrequency)
    {
        // Transpose so each row corresponds to one time point
        Eigen::MatrixXcd M_plus_transposed = M_plus.transpose();

        std::ofstream outputFile(filename);
        if (!outputFile.is_open())
        {
            std::cerr << "Error opening output file: " << filename << "\n";
            return false;
        }

        // Write simulation header
        outputFile << "# SANAE simulation output\n";
        outputFile << "# dt = " << dt << "\n";
        outputFile << "# steps = " << num_steps << "\n";
        outputFile << "# print_every = " << printFrequency << "\n";

        // Write data (downsampled according to printFrequency)
        for (int i = 0; i < M_plus_transposed.rows(); ++i)
        {
            if (i % printFrequency == 0)
            {
                outputFile << M_plus_transposed(i).real() << " "
                           << M_plus_transposed(i).imag() << "\n";
            }
        }

        return true;
    }

    bool writeCleanexCurveToFile(const std::string &filename,
                                 const std::vector<double> &tm_vec,
                                 const std::vector<double> &I_vec,
                                 double k,
                                 double R1a,
                                 double R1b,
                                 double prefactor)
    {
        if (tm_vec.size() != I_vec.size())
        {
            std::cerr << "Error: tm_vec and I_vec size mismatch\n";
            return false;
        }

        std::ofstream out(filename);
        if (!out.is_open())
        {
            std::cerr << "Error opening output file: " << filename << "\n";
            return false;
        }

        out.precision(std::numeric_limits<double>::max_digits10); // full round-trip precision

        out << "# SANAE CLEANEX-PM simulation\n";
        out << "# k    = " << k << " s\xE2\x81\xBB\xC2\xB9\n";
        out << "# R1a  = " << R1a << " s\xE2\x81\xBB\xC2\xB9  (water R1)\n";
        out << "# R1b  = " << R1b << " s\xE2\x81\xBB\xC2\xB9  (NH R1)\n";
        out << "# prefactor = " << prefactor << "\n";
        out << "# tm_ms\tI_tm\n";

        for (size_t i = 0; i < tm_vec.size(); ++i)
            out << tm_vec[i] * 1000.0 << "\t" << I_vec[i] << "\n";

        return true;
    }

    bool writeCpmgCurveToFile(const std::string &filename,
                              const std::vector<double> &nu_vec,
                              const std::vector<double> &R2eff_vec,
                              double kex,
                              double dw_ppm,
                              double B0_MHz,
                              double pA,
                              double R20,
                              const std::string &model_label)
    {
        if (nu_vec.size() != R2eff_vec.size())
        {
            std::cerr << "Error: nu_vec and R2eff_vec size mismatch\n";
            return false;
        }

        std::ofstream out(filename);
        if (!out.is_open())
        {
            std::cerr << "Error opening output file: " << filename << "\n";
            return false;
        }

        out.precision(std::numeric_limits<double>::max_digits10);

        out << "# SANAE CPMG dispersion simulation (" << model_label << ")\n";
        out << "# kex    = " << kex << " s\xE2\x81\xBB\xC2\xB9\n";
        out << "# dw_ppm = " << dw_ppm << " ppm\n";
        out << "# B0_MHz = " << B0_MHz << " MHz\n";
        out << "# pA     = " << pA << "\n";
        out << "# R20    = " << R20 << " s\xE2\x81\xBB\xC2\xB9\n";
        out << "# nu_CPMG_Hz\tR2eff\n";

        for (size_t i = 0; i < nu_vec.size(); ++i)
            out << nu_vec[i] << "\t" << R2eff_vec[i] << "\n";

        return true;
    }

    void appendToLog(const std::string &message, const Eigen::MatrixXcd &matrix)
    {
        std::ofstream logFile("sanae_simulation.log", std::ios_base::app);
        if (!logFile)
        {
            std::cerr << "Failed to open log file!" << std::endl;
            return;
        }
        logFile << message << matrix << "\n";
    }

    void appendToLog(const std::string &message, const Eigen::VectorXcd &vector)
    {
        std::ofstream logFile("sanae_simulation.log", std::ios_base::app);
        if (!logFile)
        {
            std::cerr << "Failed to open log file!" << std::endl;
            return;
        }
        logFile << message << vector << "\n";
    }

} // namespace Sanae