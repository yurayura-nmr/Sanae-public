#include "Environment.h"
//#include "Sanae.h"
#include "State.h"

// Unused but possible future idea
// In private function, also calculate R1rho - basically just for my memory aid
// Same for R2 in CPMG - Ulrich: basically monoexponential approximation to get exponential decay rate.
// R1rho_eff = (- 1 / T_relax) * ln (I / I_0)

// CPMG: R2 dispersion can detect exchange if k_ex ~~ 1/tau_CP
//       tau_CP = 2 * tau (delay between 180-pulses in CPMG pulse sequence)
// RC element averages contributions of in-phase and anti-phase coherences to R2
// R2_eff = (- 1 / T_CPMG) * ln (I / I_0)

namespace Sanae
{
	static constexpr double PI = 3.14159265358979323846;

	void State::SetID(const int ID)
	{
		ID_ = ID;
	}

	void State::SetName(const std::string &name)
	{
		name_ = name;
	}

	void State::SetInitMag(const double M0)
	{
		M0_ = M0;
	}

	void State::SetR2(const double R2)
	{
		R2_ = R2;
	}

	void State::SetR1(const double R1)
	{
		R1_ = R1;
	}

	void State::SetCS_ppm(const double cs_ppm, const Sanae::Environment &env)
	{
		cs_ppm_ = cs_ppm;
		cs_Hz_ = cs_ppm * env.GetSFO1() / (1E6); // Convert ppm to Hz
		cs_rad_s_ = cs_Hz_ * 2.0 * PI;			 // Convert chemical_shift_hz from Hz to rad/s
	}

	double State::GetInitMag() const { return M0_; }
	double State::GetR1() const { return R1_; }
	double State::GetR2() const { return R2_; }
	double State::GetCS_rad_s() const { return cs_rad_s_; }
	int State::GetID() const { return ID_; }

	void State::DisplayInfo() const
	{
		std::cout << "State Information on state [" << this->ID_ << "] named <" << this->name_ << "> \n"
				  << "Chemical shift   [ppm]   = " << this->cs_ppm_ << "\n"
				  << "Chemical shift   [Hz]    = " << this->cs_Hz_ << "\n"
				  << "Chemical shift   [rad/s] = " << this->cs_rad_s_ << "\n"
				  << "R2               [/s]    = " << this->R2_ << "\n"
				  << "Magnetization M0 (at t0) = " << this->M0_ << "\n"
				  << std::endl;
	}

	void State::DisplayInfo_R1(const double tau_) const
	{
		std::cout << "State Information on state [" << this->ID_ << "] named <" << this->name_ << "> \n"
				  << "Chemical shift   [ppm]   = " << this->cs_ppm_ << "\n"
				  << "Chemical shift   [Hz]    = " << this->cs_Hz_ << "\n"
				  << "Chemical shift   [rad/s] = " << this->cs_rad_s_ << "\n"
				  << "R1               [/s]    = " << this->R1_ << "\n"
				  << "R2               [/s]    = " << this->R2_ << "\n"
				  << "Waiting time tau [ s]    = " << tau_ << "\n"
				  << "Magnetization M0 (at t0) = " << this->M0_ << "\n"
				  << std::endl;
	}
}
