#pragma once

/**
 * @file Solution.h
 * @brief Simulation of NMR exchange dynamics.
 *
 * The `Solution` class manages the numerical integration of state populations
 * or magnetizations under Bloch–McConnell dynamics. For now, it supports two-state
 * and three-state exchange models:
 *
 *  - Construction of exchange matrices for 2- or 3-state systems
 *  - Configurable numerical integration (Euler or Runge–Kutta)
 *  - Logging and output of simulation matrices and vectors
 *  - Helper methods for setting diagonal and cross (exchange) terms
 *
 * @note Designed for modular use with `State` objects representing spin states.
 */

#include <string>
#include "State.h"
#include "RungeKutta.h"

// Using standard literals for easier handling of complex numbers
// Removed using namespace std::literals
// Remember to use literals:: then in Solution.cpp

namespace Sanae
{

	class Solution
	{
	private:
		Eigen::MatrixXcd K_; // Matrix representing exchange terms between states

		// Numerical parameters for the simulation
		double dt_ = 0.001;		 // Time step for the simulation (in seconds)
		int num_steps_ = 2000;	 // Total number of time steps in the simulation
		int printFrequency_ = 1; // Frequency of printing results
		int method_ = 1;		 // Integration method: 0 = Euler, 1 = Runge-Kutta

		std::string output_filename_ = "magnetization_evolution.txt";

	public:
		// Constructor initializes the exchange matrix K_ based on the number of states
		Solution(int size);

		// Setters for customizing simulation parameters
		void SetMethod(const int method);
		void SetPF(const int customPF);
		void SetDT(const double customDT);
		void SetSteps(const int customSteps);
		void SetNumericPars();

		void SetCrossterm(const double k_ij, Sanae::State FirstState, Sanae::State SecondState);
		void SetDiagonaltermsTwoState(const double k_11, const double k_22, Sanae::State FirstState, Sanae::State SecondState);
		void SetDiagonaltermsThreeState(const double k_11, const double k_22, const double k_33, Sanae::State FirstState, Sanae::State SecondState, Sanae::State ThirdState);
		void SetNoExchange(Sanae::State &A, Sanae::State &B);

		void SetupOneState(Sanae::State &StateA);
		void SetupTwoState(Sanae::State &StateA, Sanae::State &StateB);
		void SetupThreeState(Sanae::State &StateA, Sanae::State &StateB, Sanae::State &StateC);

		int GetSteps() const;
		Eigen::MatrixXcd GetK() const;
		void DisplayInfo() const;

		void SetOutputFilename(const std::string &filename);
		void writeLog(const std::string message, Eigen::MatrixXcd &A) const;
		void writeLog(const std::string message, Eigen::VectorXcd &vector) const;
	};
}
