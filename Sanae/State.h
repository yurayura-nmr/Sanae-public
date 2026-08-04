#pragma once

/**
 * @file State.h
 * @brief Represents a single spin or state in NMR simulations.
 *
 * The `State` class encapsulates all relevant properties of a spin/state
 * used in Bloch–McConnell or relaxation simulations, including:
 *  - Initial magnetization
 *  - Longitudinal (R1) and transverse (R2) relaxation rates
 *  - Chemical shift (ppm, Hz, rad/s)
 *  - Identification and name for bookkeeping
 *
 * @note Designed to integrate seamlessly with `Environment` and `Solution` classes.
 */

#include <string>
#include <iostream>
#include "Environment.h"

namespace Sanae
{

	class State
	{

	public:
		// Methods
		void SetID(const int ID);
		void SetName(const std::string &name);

		void SetInitMag(const double M0);
		void SetR1(const double R1);
		void SetR2(const double R2);
		void SetCS_ppm(const double cs_ppm, const Sanae::Environment &env);

		double GetInitMag() const;
		double GetR1() const;
		double GetR2() const;
		double GetCS_rad_s() const;
		int GetID() const;

		void DisplayInfo() const;
		void DisplayInfo_R1(const double) const;

	private:
		int ID_ = 0; // Position in Bloch-McConnell matrix
		std::string name_ = "NoNameSpin";

		double R2_ = 7.0; // Relaxation rate
		double R1_ = 3.0;
		double cs_ppm_ = 0.0;	// Chemical shift in ppm
		double cs_Hz_ = 0.0;	// Chemical shift in Hz
		double cs_rad_s_ = 0.0; // Chemical shift in rad/s
		double M0_ = 1.0;		// Initial magnetization for this state
	};

}
