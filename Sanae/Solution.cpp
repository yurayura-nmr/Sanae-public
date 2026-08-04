#include "Solution.h"
#include "State.h"
#include "RungeKutta.h"
#include "noise.h"

#include <fstream>
#include <vector>
#include <iostream>

namespace Sanae
{

	Solution::Solution(int size)
	{
		if (size == 1)
		{
			// Single spin (no exchange, e.g., 1 spin in R1 experiment)
			K_ = Eigen::MatrixXcd::Zero(1, 1);
		}
		else if (size == 2)
		{
			// 2-state exchange (A <-> B)
			K_ = Eigen::MatrixXcd::Zero(2, 2);
		}
		else if (size == 3)
		{
			// 3-state exchange (F <-> I <-> B)
			K_ = Eigen::MatrixXcd::Zero(3, 3);
		}
		else
		{
			std::cout << "Invalid exchange regime (" << size
					  << "). Defaulting to 2-state.\n";
			K_ = Eigen::MatrixXcd::Zero(2, 2);
		};
	}

	// Setters for customizing simulation parameters
	void Solution::SetPF(const int customPF)
	{
		printFrequency_ = customPF;
	} // Set the print frequency

	void Solution::SetDT(const double customDT)
	{
		dt_ = customDT;
	} // Set the time step duration
	void Solution::SetSteps(const int customSteps)
	{
		num_steps_ = customSteps;
	} // Set the number of steps

	void Solution::SetNumericPars()
	{
		// Set default numerical parameters for the simulation
		// Note: Modifying these parameters may affect the time axis in sanaeplot.py; verify any changes
		this->SetMethod(1);	   // Use Runge-Kutta (R-K) method
		this->SetPF(10);	   // Print results every 10 steps
		this->SetDT(0.00001);  // Set timestep to 10 microseconds
		this->SetSteps(20000); // 20000 steps at 10 ƒÊs = 0.2 seconds total simulation time; 2000 print intervals
		this->DisplayInfo();   // Display current settings to the user
	}

	void Solution::SetNoExchange(Sanae::State &A, Sanae::State &B)
	{
		// Configure system for no exchange between states (Bloch relaxation only)
		double k_ab = 0; // Rate of transition from state A to B is zero
		double k_ba = 0; // Rate of transition from state B to A is zero

		this->SetCrossterm(k_ab, A, B); // Set off-diagonal (cross) terms for A <> B transitions
		this->SetCrossterm(k_ba, B, A); // Set off-diagonal (cross) terms for B <> A transitions

		// Set diagonal terms in the exchange matrix
		// Diagonal terms are set in the order upper left to lower right : 11, 22, ...
		this->SetDiagonaltermsTwoState(k_ab, k_ba, A, B); // Diagonal terms: k_11 for A, k_22 for B
	}

	void Solution::writeLog(const std::string message, Eigen::MatrixXcd &A) const
	{
		// Open the log file
		std::ofstream logFile("sanae_simulation.log", std::ios_base::app); // Use append mode

		// Check if log file opened successfully
		if (!logFile)
		{
			std::cerr << "Failed to open log file!" << std::endl;
			return;
		}

		// Write to the log file
		logFile << message << A << "\n";
	}

	// Overload for vectors
	void Solution::writeLog(const std::string message, Eigen::VectorXcd &vector) const
	{
		std::ofstream logFile("sanae_simulation.log", std::ios_base::app);

		if (!logFile)
		{
			std::cerr << "Failed to open log file!" << std::endl;
			return;
		}

		logFile << message << vector << "\n";
	}

	int Solution::GetSteps() const
	{
		return num_steps_;
	}

	Eigen::MatrixXcd Solution::GetK() const
	{
		return K_;
	}

	void Solution::SetMethod(const int method)
	{
		method_ = method;
	} // Set the numerical integration method (Euler or Runge-Kutta)

	void Solution::DisplayInfo() const
	{
		// Display key simulation parameters to the console
		std::cout << "\nInformation on Solution (FID)\n"
				  << "Integration timestep (= FID point spacing)   [s] = " << this->dt_ << "\n"
				  << "Number of timesteps (= total FID points)         = " << this->num_steps_ << "\n"
				  << "Output frequency                                 = " << this->printFrequency_ << "\n"
				  << "Integration method (0: Euler | 1: Runge-Kutta)   = " << this->method_ << "\n"
				  << std::endl;
	}

	void Solution::SetCrossterm(const double k_ij, Sanae::State FirstState, Sanae::State SecondState)
	{
		// Set off-diagonal cross-term in exchange matrix K_ for the exchange between two states
		const double first_id_ = FirstState.GetID();   // Get the ID of the first state
		const double second_id_ = SecondState.GetID(); // Get the ID of the second state
		this->K_(second_id_, first_id_) = k_ij;		   // Update the exchange matrix with the cross-term value

		// Print the updated exchange matrix and the involved states
		std::cout << "\n[ ... Adding cross-term to exchange matrix ... ]\n"
				  << K_ << "\nfor exchange between: " << first_id_ << " and " << second_id_ << "\n";
	}

	void Solution::SetDiagonaltermsTwoState(const double k_11, const double k_22, Sanae::State FirstState, Sanae::State SecondState)
	{
		// Set diagonal terms for two-state exchange in the exchange matrix
		// Example: k_11 represents the decay of state A due to exchange with state B
		// Example: k_22 represents the decay of state B due to exchange with state A
		// In a three-state system, k_22 would represent decay of state B due to exchanges with both A and C
		const double first_id_ = FirstState.GetID();   // Get the ID of the first state
		const double second_id_ = SecondState.GetID(); // Get the ID of the second state

		// Update the diagonal terms: for a two-state system, this is straightforward
		this->K_(first_id_, first_id_) = -k_11;
		this->K_(second_id_, second_id_) = -k_22;

		// Print the updated matrix and the states involved in the decay process
		std::cout << "[ ... Adding diagonal-term to exchange matrix ... ]\n"
				  << K_ << "\nfor decay of: " << first_id_ << " due to the presence of State " << second_id_ << " \n";
	}

	void Solution::SetDiagonaltermsThreeState(const double k_11, const double k_22, const double k_33, Sanae::State FirstState, Sanae::State SecondState, Sanae::State ThirdState)
	{
		// Set diagonal terms for linear three-state exchange in the exchange matrix
		// Example: k_11 represents the decay of state A due to exchange with state B
		// Example: k_22 represents the decay of state B due to exchanges with both A and C
		// Example: k_33 represents the decay of state C due to exchange with B
		const double first_id_ = FirstState.GetID();   // Get the ID of the first state
		const double second_id_ = SecondState.GetID(); // Get the ID of the second state
		const double third_id_ = ThirdState.GetID();   // Get the ID of the third state

		// Update diagonal terms for each state, accounting for contributions from one (A,C) or multiple (B) exchanges
		this->K_(first_id_, first_id_) = -k_11;	  // Decay for first state (A)
		this->K_(second_id_, second_id_) = -k_22; // Decay for second state (B), combines contributions from A and C
		this->K_(third_id_, third_id_) = -k_33;	  // Decay for third state (C)

		// Print the updated matrix and the states involved in the decay process
		std::cout << "[ ... Adding diagonal-term to exchange matrix ... ]\n"
				  << K_ << "\nfor decay of: " << first_id_ << " due to the presence of State " << second_id_ << " \n";
	}

	void Solution::SetupThreeState(Sanae::State &StateA, Sanae::State &StateB, Sanae::State &StateC)
	{
		// Create a matrix to store the time evolution of magnetization M(t) for the 3-state system
		Eigen::MatrixXcd M(3, num_steps_);

		// Define the initial magnetization state at t=0 (M_zero) for the 3-state exchange in the transverse plane
		Eigen::VectorXcd M_zero(3);												 // Magnetization vector for the 3 states
		M_zero << StateA.GetInitMag(), StateB.GetInitMag(), StateC.GetInitMag(); // Initial magnetization for each state
		M.col(0) = M_zero;														 // Set the first column to M(t=0), i.e., the equilibrium magnetization

		// Initialize matrices for solving the ODEs based on the Bloch equations
		Eigen::MatrixXcd A(3, 3);		// Bloch-McConnell matrix (total system evolution)
		Eigen::MatrixXcd R(3, 3);		// Relaxation and chemical shift matrix
		Eigen::MatrixXcd K_total(3, 3); // Exchange matrix (K)
		Eigen::MatrixXcd debug(3, 3);	// Debugging matrix (for intermediate checks)

		// Define relaxation and chemical shift terms for better readability
		double mR2_A = -StateA.GetR2(); // Negative R2 relaxation rate for State A
		double mR2_B = -StateB.GetR2(); // Negative R2 relaxation rate for State B
		double mR2_C = -StateC.GetR2(); // Negative R2 relaxation rate for State C

		double CS_A = StateA.GetCS_rad_s(); // Chemical shift in rad/s for State A
		double CS_B = StateB.GetCS_rad_s(); // Chemical shift in rad/s for State B
		double CS_C = StateC.GetCS_rad_s(); // Chemical shift in rad/s for State C

		// *** Relaxation and chemical shift contributions ***
		// Matrix R: contains relaxation rates and chemical shifts for each state
		//   M_A                    M_B
		R << (mR2_A + 1.0i * CS_A), 0.0, 0.0,
			0.0, (mR2_B + 1.0i * CS_B), 0.0,
			0.0, 0.0, (mR2_C + 1.0i * CS_C);

		// Combine relaxation/chemical shift matrix R with the exchange matrix K to form the full Bloch-McConnell matrix
		A = R + this->GetK();

		// Output the final matrix for inspection
		std::cout << "[... Final matrix is set to ...] :\n"
				  << A << "\n";

		// Output log to enable reproducibility
		writeLog("[... Initial magnetization is set to ...] :\n", M_zero);
		writeLog("[... Final matrix is set to ...] :\n", A);

		if (this->method_ == 1)
		{
			// Solving dM/dt = A * M(0) using the 4th-order Runge-Kutta method
			std::cout << "\n[... Solving d/dt M(t) = A * M(0) using Runge-Kutta method ...]\n"
					  << std::endl;
			RungeKutta(A, M, this->dt_, this->num_steps_);

			// Optional add noise
			//AddGaussianNoise(M, 0.0, -1);
		}
		else if (this->method_ == 0)
		{
			// Euler method is not implemented
			std::cout << "\n[... Euler method is not yet implemented. Use Runge-Kutta. ...]\n"
					  << std::endl;
		}

		// Extract solved magnetization for each state: A, B, and C
		Eigen::MatrixXcd M_A_solved = M.row(0); // Magnetization for State A
		Eigen::MatrixXcd M_B_solved = M.row(1); // Magnetization for State B
		Eigen::MatrixXcd M_C_solved = M.row(2); // Magnetization for State C

		// Combine magnetization for all states (A + B + C), as only M+ is observed in NMR
		Eigen::MatrixXcd M_plus_ABC = M_A_solved + M_B_solved + M_C_solved;

		// Output the initial magnetization and the relaxation matrix (Bloch-McConnell matrix)
		std::cout << "\n[... Magnetization at the start of the experiment ...]\nM(t=0) = \n"
				  << M_zero.real() << std::endl;
		std::cout << "\n[... Relaxation matrix determining the physics of the 3-state exchange process ...]\nA = \n"
				  << A << std::endl;

		// Debug
		// std::cout << "\n[... Time Evolution of total Magnetization A+B throughout the experiment ...]\n[M+(t)] = \n" << M_plus_AB << std::endl;

		// Write the time evolution of the magnetization (M+) to a file
		// Transpose M+ matrix to organize data for output
		Eigen::MatrixXcd M_plus_transposed = M_plus_ABC.transpose();

		// Open output file to store magnetization data
		std::ofstream outputFile(output_filename_);
		if (!outputFile.is_open())
		{
			std::cerr << "Error opening output file!\n";
		}
		else
		{
			// --- Write simulation header ---
			// sanaeplot.py reads dt and steps from here, eliminating manual synchronisation
			// Format: # KEY = VALUE  (lines starting with # are skipped by numpy.loadtxt)
			outputFile << "# SANAE simulation output\n";
			outputFile << "# dt = " << this->dt_ << "\n";
			outputFile << "# steps = " << this->num_steps_ << "\n";
			outputFile << "# print_every = " << this->printFrequency_ << "\n";

			// Loop through the transposed M+ matrix and write the data to the file
			// Only print data at specified intervals (defined by printFrequency_)
			for (int i = 0; i < M_plus_transposed.rows(); ++i)
			{
				if (i % this->printFrequency_ == 0)
				{
					// Write real and imaginary parts (Mx, My) of the magnetization to the file
					outputFile << M_plus_transposed(i).real() << " " << M_plus_transposed(i).imag() << std::endl;
				}
			}
		}
	}

	void Sanae::Solution::SetOutputFilename(const std::string &filename)
	{
		output_filename_ = filename;
	}

	void Solution::SetupTwoState(Sanae::State &StateA, Sanae::State &StateB)
	{
		// Create a matrix to store the time evolution of magnetization in the transverse plane: M(t) for the 2-state system
		Eigen::MatrixXcd M(2, num_steps_);

		// Define the initial state of the magnetization: M(t=0)
		Eigen::VectorXcd M_zero(2);							// Magnetization vector for 2-state system [2x2 matrix]
		M_zero << StateA.GetInitMag(), StateB.GetInitMag(); // Initial magnetization: M_A(0) = M_eq_A, M_B(0) = M_eq_B
		M.col(0) = M_zero;									// Set the initial magnetization for t=0

		// Define the matrices to solve the ODEs (Bloch equations) for magnetization evolution
		Eigen::MatrixXcd A(2, 2);		// Total matrix (Bloch-McConnell matrix)
		Eigen::MatrixXcd R(2, 2);		// Relaxation and chemical shift matrix
		Eigen::MatrixXcd K_total(2, 2); // Exchange matrix (used later)
		Eigen::MatrixXcd debug(2, 2);

		// Define variables to represent relaxation rates and chemical shifts for States A and B
		double mR2_A = -StateA.GetR2();		// Relaxation rate for State A (negative of R2)
		double mR2_B = -StateB.GetR2();		// Relaxation rate for State B
		double CS_A = StateA.GetCS_rad_s(); // Chemical shift for State A (in radians/second)
		double CS_B = StateB.GetCS_rad_s(); // Chemical shift for State B

		// *** Define relaxation and chemical shift matrix ***
		// Diagonal matrix with relaxation and chemical shifts for State A and State B
		//   M_A                    M_B
		R << (mR2_A + 1.0i * CS_A), 0.0,
			0.0, (mR2_B + 1.0i * CS_B);

		// Add exchange terms (from matrix K) to the relaxation matrix to form the full Bloch-McConnell matrix
		A = R + this->GetK();

		// Output the final matrix to be used in time evolution
		std::cout << "[... Final matrix is set to ...] :\n"
				  << A << "\n";

		// Output log to enable reproducibility
		writeLog("[... Initial magnetization is set to ...] :\n", M_zero);
		writeLog("[... Final matrix is set to ...] :\n", A);

		// Solve the ODE using the specified integration method
		if (this->method_ == 1)
		{
			// Solve using 4th-order Runge-Kutta method
			std::cout << "\n[... Solving d/dt M(t) = A * M(0) using Runge-Kutta method ...]\n"
					  << std::endl;
			RungeKutta(A, M, this->dt_, this->num_steps_);
		}
		else if (this->method_ == 0)
		{
			// Placeholder for Euler method (if needed)
		}

		// Combine magnetization components for State A and B into a single quantity M+ (only M+ observed in NMR)
		Eigen::MatrixXcd M_A_solved = M.row(0); // Magnetization for State A
		Eigen::MatrixXcd M_B_solved = M.row(1); // Magnetization for State B

		// Sum of M+_A and M+_B magnetization components (M+ in NMR)
		Eigen::MatrixXcd M_plus_AB = M_A_solved + M_B_solved;

		// Output the initial magnetization and relaxation matrix
		std::cout << "\n[... Magnetization at the start of the experiment ...]\nM(t=0) = \n"
				  << M_zero.real() << std::endl;
		std::cout << "\n[... Relaxation matrix determining the physics of the 2-state exchange process ...]\nA = \n"
				  << A << std::endl;

		// Debug
		// std::cout << "\n[... Time Evolution of total Magnetization A+B throughout the experiment ...]\n[M+(t)] = \n" << M_plus_AB << std::endl;

		// Write the time evolution of the magnetization (M+) to a file
		Eigen::MatrixXcd M_plus_transposed = M_plus_AB.transpose(); // Transpose the matrix for easier file output

		// Open a file to write the magnetization data
		std::ofstream outputFile(output_filename_);
		if (!outputFile.is_open())
		{
			std::cerr << "Error opening output file!\n";
		}
		else
		{
			// --- Write simulation header ---
			// sanaeplot.py reads dt and steps from here, eliminating manual synchronisation
			// Format: # KEY = VALUE  (lines starting with # are skipped by numpy.loadtxt)
			outputFile << "# SANAE simulation output\n";
			outputFile << "# dt = " << this->dt_ << "\n";
			outputFile << "# steps = " << this->num_steps_ << "\n";
			outputFile << "# print_every = " << this->printFrequency_ << "\n";

			// Loop through and write the time evolution data (only at intervals specified by printFrequency_)
			for (int i = 0; i < M_plus_transposed.rows(); ++i)
			{
				if (i % this->printFrequency_ == 0)
				{
					// Write the real and imaginary components (Mx and My) of the magnetization
					outputFile << M_plus_transposed(i).real() << " " << M_plus_transposed(i).imag() << std::endl;
				}
			}
		}
	}

	void Solution::SetupOneState(Sanae::State &StateA)
	{
		// --- Allocate magnetization container ---------------------------------
		Eigen::MatrixXcd M(1, num_steps_);

		// Initial magnetization M(t=0)
		Eigen::VectorXcd M_zero(1);
		M_zero(0) = StateA.GetInitMag();
		M.col(0) = M_zero;

		// --- Construct the 1x1 Bloch matrix A ----------------------------------
		// Relaxation: -R2
		// Chemical shift: + i * ƒ¢ƒÖ
		double mR2 = -StateA.GetR2();
		double w = StateA.GetCS_rad_s(); // rad/s offset

		Eigen::MatrixXcd A(1, 1);
		A(0, 0) = std::complex<double>(mR2, w);

		// Log for reproducibility
		writeLog("[... Initial magnetization is set to ...] :\n", M_zero);
		writeLog("[... Final matrix is set to ...] :\n", A);

		std::cout << "[... Final matrix is set to ...] :\n"
				  << A << "\n";

		// --- Solve ODE dM/dt = A * M ------------------------------------------
		if (this->method_ == 1)
		{
			std::cout << "\n[... Solving d/dt M(t) = A * M(0) using Runge-Kutta method ...]\n";
			RungeKutta(A, M, this->dt_, this->num_steps_);
		}
		else if (this->method_ == 0)
		{
			// Optional Euler implementation
		}

		// --- Extract the only magnetization component (M+) ---------------------
		Eigen::MatrixXcd M_plus = M; // For 1 state, M+ = M
		Eigen::MatrixXcd M_plus_T = M_plus.transpose();

		// --- Output M(t) to file ----------------------------------------------
		std::ofstream outputFile(output_filename_);

		if (!outputFile.is_open())
		{
			std::cerr << "Error opening output file!\n";
			return;
		}
		else
		{
			// TODO - Test this. So far not tested.

			// --- Write simulation header ---
			// sanaeplot.py reads dt and steps from here, eliminating manual synchronisation
			// Format: # KEY = VALUE  (lines starting with # are skipped by numpy.loadtxt)
			outputFile << "# SANAE simulation output\n";
			outputFile << "# dt = " << this->dt_ << "\n";
			outputFile << "# steps = " << this->num_steps_ << "\n";
			outputFile << "# print_every = " << this->printFrequency_ << "\n";

			for (int i = 0; i < M_plus_T.rows(); ++i)
			{
				if (i % this->printFrequency_ == 0)
				{
					outputFile << M_plus_T(i).real() << " "
							   << M_plus_T(i).imag() << "\n";
				}
			}

			std::cout << "\n[... Magnetization at the start ...]\nM(0) = "
					  << M_zero.real() << "\n";
			std::cout << "\n[... Relaxation matrix for 1-state ...]\nA = "
					  << A << "\n";
		}
	}
}
