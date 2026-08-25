//#include "Sanae.h"
#include "Environment.h"
#include "State.h"
#include "Solution.h"

void bloch_relax(
	double R2_A,    // R2 for state A
	double R2_B,    // R2 for state B
	double CSppm_A, // Chemical shift in ppm for A
	double CSppm_B, // Chemical shift in ppm for B
	double M0_A,    // Initial magnetization for A
	double M0_B,    // Initial magnetization for B
	char nucleus)
{
	// Define the NMR magnet
	Sanae::Environment env;
	env.SetGamma(nucleus);
	env.DisplayInfo();

	// Define states
	Sanae::State A, B;

	// State A
	A.SetID(0);
	A.SetName("Freely precessing spin A with only relaxation (no exchange)");
	A.SetR2(R2_A);
	A.SetCS_ppm(CSppm_A, env);
	A.SetInitMag(M0_A);

	// State B
	B.SetID(1);
	B.SetName("Freely precessing spin B with only relaxation (no exchange)");
	B.SetR2(R2_B);
	B.SetCS_ppm(CSppm_B, env);
	B.SetInitMag(M0_B);

	A.DisplayInfo();
	B.DisplayInfo();

	// *** Set the exchange rates *** (no exchange)
	Sanae::Solution sol(2);
	//sol.SetNumericPars();  // Uses default: dt=10µs, steps=20000, print=10
	sol.SetupTwoState(A, B);
}
