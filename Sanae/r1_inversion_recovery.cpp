#include "r1_inversion_recovery.h"
//#include "Sanae.h"
#include "State.h"
#include "Solution.h"

namespace Sanae
{
	Eigen::MatrixXcd R1_inversion_recovery(double R1, double R2, double CSppm, double tau, char nucleus, const std::string &output_filename, bool verbose)
	{
		/*
		Simulation of a single spin undergoing longitudinal relaxation with rate R1
		as part of an inversion recovery experiment:		

		Initial condition:  M_z
		Complete inversion: Magnetization (M_z -> -M_z)
		Pulse sequence:     No delay after inversion before detection
							(i.e., the 90 degree pulse is applied immediately)

		Calculate the surviving magnetization after a delay time tau, due to longitudinal relaxation.

		Inversion recovery magnetization: Mz(tau) = Meq * (1 - 2 * exp(-tau / T1))

		NOTE: Meq = 1 is assumed here, so the full formula reduces to
			  a form proportional to exp(-tau * R1).
			  See @note in r1_inversion_recovery.h.
		TODO: if Meq != 1 is ever needed, replace with the full expression:
			   Meq * (1.0 - 2.0 * std::exp(-tau * R1))		
		*/

		// --- NMR environment ----------------------------------------------------
		Sanae::Environment env;
		env.SetGamma(nucleus); // User chooses nucleus (H, N, C, etc.)
		env.DisplayInfo();

		// --- Define the single spin state --------------------------------------
		Sanae::State A;

		A.SetID(0);
		A.SetName("Single spin prepared by inversion recovery");
		A.SetR1(R1);
		A.SetR2(R2);
		A.SetCS_ppm(CSppm, env);

		double Mx_iR = std::exp(-tau * R1);

		// Set magnetization at end of inversion recovery for simulation in transverse
		A.SetInitMag(Mx_iR);
		A.DisplayInfo_R1(tau);

		// --- Numerical simulation setup ----------------------------------------
		Sanae::Solution Solution(1); // Use 1-state integrator
		Solution.SetVerbose(verbose);
		Solution.SetNumericPars();	 // Integration parameters		
		Solution.SetupOneState(A);	 // TODO: rename / change return type — see header @todo

		return Solution.GetLastMagnetization();
	}
}