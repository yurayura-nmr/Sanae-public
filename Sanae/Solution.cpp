#include "Solution.h"
#include "State.h"
#include "RungeKutta.h"
#include "fileio.h"
// #include "noise.h"

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
			// 3-state exchange linear (F <-> I <-> B)
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
		SetMethod(1);	 // Use Runge-Kutta (R-K) method
		SetPF(10);		 // Print results every 10 steps
		SetDT(0.00001);	 // Set timestep to 10 microseconds
		SetSteps(20000); // 20000 steps at 10 microseconds = 0.2 seconds total simulation time; 2000 print intervals
		DisplayInfo();	 // Display current settings to the user
	}

	void Solution::SetNoExchange(Sanae::State &A, Sanae::State &B)
	{
		// Configure system for no exchange between states (Bloch relaxation only)
		double k_ab = 0; // Rate of transition from state A to B is zero
		double k_ba = 0; // Rate of transition from state B to A is zero

		SetCrossterm(k_ab, A, B); // Set off-diagonal (cross) terms for A <> B transitions
		SetCrossterm(k_ba, B, A); // Set off-diagonal (cross) terms for B <> A transitions

		// Set diagonal terms in the exchange matrix
		// Diagonal terms are set in the order upper left to lower right : 11, 22, ...
		SetDiagonaltermsTwoState(k_ab, k_ba, A, B); // Diagonal terms: k_11 for A, k_22 for B
	}

	int Solution::GetSteps() const
	{
		return num_steps_;
	}

	Eigen::MatrixXcd Solution::GetK() const
	{
		return K_;
	}

	void Solution::SetMethod(const int /*method*/)
	{
		// Deprecated: only Runge-Kutta is supported.
		// Intentionally left empty to avoid breaking existing calls.
	}

	void Solution::DisplayInfo() const
	{
		// Display key simulation parameters to the console
		std::cout << "\nInformation on Solution (FID)\n"
				  << "Integration timestep (= FID point spacing)   [s] = " << dt_ << "\n"
				  << "Number of timesteps (= total FID points)         = " << num_steps_ << "\n"
				  << "Output frequency                                 = " << printFrequency_ << "\n"
				  << "Integration method (0: Euler | 1: Runge-Kutta)   = " << method_ << "\n"
				  << std::endl;
	}

	void Solution::SetCrossterm(const double k_ij, Sanae::State FirstState, Sanae::State SecondState)
	{
		// Set off-diagonal cross-term in exchange matrix K_ for the exchange between two states
		const double first_id_ = FirstState.GetID();   // Get the ID of the first state
		const double second_id_ = SecondState.GetID(); // Get the ID of the second state
		K_(second_id_, first_id_) = k_ij;			   // Update the exchange matrix with the cross-term value

		// Debug: Print the updated exchange matrix and the involved states
		// std::cout << "\n[ ... Adding cross-term to exchange matrix ... ]\n"
		//		  << K_ << "\nfor exchange between: " << first_id_ << " and " << second_id_ << "\n";
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
		K_(first_id_, first_id_) = -k_11;
		K_(second_id_, second_id_) = -k_22;

		// Debug: Print the updated matrix and the states involved in the decay process
		// std::cout << "[ ... Adding diagonal-term to exchange matrix ... ]\n"
		//		  << K_ << "\nfor decay of: " << first_id_ << " due to the presence of State " << second_id_ << " \n";
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
		K_(first_id_, first_id_) = -k_11;	// Decay for first state (A)
		K_(second_id_, second_id_) = -k_22; // Decay for second state (B), combines contributions from A and C
		K_(third_id_, third_id_) = -k_33;	// Decay for third state (C)

		// Debug: Print the updated matrix and the states involved in the decay process
		// std::cout << "[ ... Adding diagonal-term to exchange matrix ... ]\n"
		//		  << K_ << "\nfor decay of: " << first_id_ << " due to the presence of State " << second_id_ << " \n";
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
		Eigen::MatrixXcd A(3, 3); // Bloch-McConnell matrix (total system evolution)
		Eigen::MatrixXcd R(3, 3); // Relaxation and chemical shift matrix

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
		R << std::complex<double>(mR2_A, CS_A), 0.0, 0.0,
			0.0, std::complex<double>(mR2_B, CS_B), 0.0,
			0.0, 0.0, std::complex<double>(mR2_C, CS_C);

		// Combine relaxation/chemical shift matrix R with the exchange matrix K to form the full Bloch-McConnell matrix
		A = R + GetK();

		// Output the final matrix for inspection
		std::cout << "[... Final matrix is set to ...] :\n"
				  << A << "\n";

		// Output log to enable reproducibility
		appendToLog("[... Initial magnetization is set to ...] :\n", M_zero);
		appendToLog("[... Final matrix is set to ...] :\n", A);

		if (method_ == 1)
		{
			// Solving dM/dt = A * M(0) using the 4th-order Runge-Kutta method
			std::cout << "\n[... Solving d/dt M(t) = A * M(0) using Runge-Kutta method ...]\n"
					  << std::endl;
			RungeKutta(A, M, dt_, num_steps_);

			// Optional add noise
			// AddGaussianNoise(M, 0.0, -1);
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

		// Open output file to store magnetization data
		if (!writeMagnetizationToFile(M_plus_ABC, output_filename_, dt_, num_steps_, printFrequency_))
		{
			return;
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
		Eigen::MatrixXcd A(2, 2); // Total matrix (Bloch-McConnell matrix)
		Eigen::MatrixXcd R(2, 2); // Relaxation and chemical shift matrix

		// Define variables to represent relaxation rates and chemical shifts for States A and B
		double mR2_A = -StateA.GetR2();		// Relaxation rate for State A (negative of R2)
		double mR2_B = -StateB.GetR2();		// Relaxation rate for State B
		double CS_A = StateA.GetCS_rad_s(); // Chemical shift for State A (in radians/second)
		double CS_B = StateB.GetCS_rad_s(); // Chemical shift for State B

		// *** Define relaxation and chemical shift matrix ***
		// Diagonal matrix with relaxation and chemical shifts for State A and State B
		//   M_A                    M_B
		R << std::complex<double>(mR2_A, CS_A), 0.0,
			0.0, std::complex<double>(mR2_B, CS_B);

		// Add exchange terms (from matrix K) to the relaxation matrix to form the full Bloch-McConnell matrix
		A = R + GetK();

		// Output the final matrix to be used in time evolution
		std::cout << "[... Final matrix is set to ...] :\n"
				  << A << "\n";

		// Output log to enable reproducibility
		appendToLog("[... Initial magnetization is set to ...] :\n", M_zero);
		appendToLog("[... Final matrix is set to ...] :\n", A);

		// Solve the ODE using the specified integration method
		if (method_ == 1)
		{
			// Solve using 4th-order Runge-Kutta method
			std::cout << "\n[... Solving d/dt M(t) = A * M(0) using Runge-Kutta method ...]\n"
					  << std::endl;
			RungeKutta(A, M, dt_, num_steps_);
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

		// --- Write M(t) to file -----------------------------------------------
		if (!writeMagnetizationToFile(M_plus_AB, output_filename_, dt_, num_steps_, printFrequency_))
		{
			return;
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
		// Chemical shift: + i * w
		double mR2 = -StateA.GetR2();
		double w = StateA.GetCS_rad_s(); // rad/s offset

		Eigen::MatrixXcd A(1, 1);
		A(0, 0) = std::complex<double>(mR2, w);

		appendToLog("[... Initial magnetization is set to ...] :\n", M_zero);
		appendToLog("[... Final matrix is set to ...] :\n", A);

		std::cout << "[... Final matrix is set to ...] :\n"
				  << A << "\n";

		// --- Solve ODE dM/dt = A * M ------------------------------------------
		if (method_ == 1)
		{
			std::cout << "\n[... Solving d/dt M(t) = A * M(0) using Runge-Kutta method ...]\n";
			RungeKutta(A, M, dt_, num_steps_);
		}

		// --- Extract the only magnetization component (M+) ---------------------
		Eigen::MatrixXcd M_plus = M; // For 1 state, M+ = M
		Eigen::MatrixXcd M_plus_T = M_plus.transpose();

		// --- Output M(t) to file ----------------------------------------------
		if (!writeMagnetizationToFile(M_plus, output_filename_, dt_, num_steps_, printFrequency_))
		{
			return;
		}

		// (Keep any debug prints you still want – these can also be removed for quiet output)
		std::cout << "\n[... Magnetization at the start ...]\nM(0) = "
				  << M_zero.real() << "\n";
		std::cout << "\n[... Relaxation matrix for 1-state ...]\nA = "
				  << A << "\n";
	}
}
