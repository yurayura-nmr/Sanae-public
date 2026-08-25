#pragma once

#include <Eigen/Core>
#include <string>
#include <vector>

namespace Sanae
{

    /**
     * @brief Write the time-domain magnetization (FID) to a text file.
     *
     * The file format includes a header with simulation parameters (dt, steps,
     * print_every) followed by columns: real(M+)  imag(M+).
     *
     * @param M_plus       Eigen matrix of size (num_states x num_steps) containing
     *                     the total magnetization M+ for each time point.
     * @param filename     Output file path.
     * @param dt           Time step (seconds), written to header.
     * @param num_steps    Total number of steps, written to header.
     * @param printFrequency Number of steps between output rows (downsampling).
     * @return true if file was written successfully, false otherwise.
     */
    bool writeMagnetizationToFile(const Eigen::MatrixXcd &M_plus,
                                  const std::string &filename,
                                  double dt,
                                  int num_steps,
                                  int printFrequency);

    // Write a CLEANEX-PM simulated curve to a tab-separated file with a
    // descriptive header. Returns false (and logs to stderr) on failure to open.
    bool writeCleanexCurveToFile(const std::string &filename,
                                 const std::vector<double> &tm_vec,
                                 const std::vector<double> &I_vec,
                                 double k,
                                 double R1a,
                                 double R1b,
                                 double prefactor);

    /**
     * @brief Append a matrix to the simulation log file (sanae_simulation.log).
     *
     * @param message    Descriptive text to prepend.
     * @param matrix     Matrix to log (read‑only).
     */
    void appendToLog(const std::string &message, const Eigen::MatrixXcd &matrix);

    /**
     * @brief Append a vector to the simulation log file (sanae_simulation.log).
     *
     * @param message    Descriptive text to prepend.
     * @param vector     Vector to log (read‑only).
     */
    void appendToLog(const std::string &message, const Eigen::VectorXcd &vector);

} // namespace Sanae