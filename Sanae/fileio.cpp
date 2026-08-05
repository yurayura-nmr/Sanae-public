#include "fileio.h"
#include <fstream>
#include <iostream>

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