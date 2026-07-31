/**
* @file Environment.h
* @brief Defines the NMR environment parameters for Sanae simulations.
*
* Represents spectrometer and spin-environment settings such as magnetic
* field strength, gyromagnetic ratio, and derived Larmor frequencies.
*
* @note Only a minimal subset of fields is currently used. Additional
*       attributes (e.g., B1 fields, offsets, tilt angles) are reserved
*       for future extensions.
*/

#pragma once

namespace Sanae
{
	class Environment
	{
	private:
		// Attributes for NMR settings
		double gamma_ = 42.58;	// [MHz / T]      Default: gamma of 1H spin
		double gamma_rads_ = 0; // (rad / s) / T  Calculated in SetGamma()
		double B0_ = 16.4;		// [T]            Static magnetic field; Default: 700 MHz spectrometer.
		double SFO1_ = 0;		// [Hz]
		double w_mhz_ = 0;		// [MHz]          Larmor frequency

		// Private method to calculate SFO1 from Larmor Frequency
		void SetSFO1()
		{
			SFO1_ = 1.0 * w_mhz_ * 1E6;
		}

	public:
		// Constructor
		Environment() = default;

		// Set magnetic field strength
		void SetB0(const double customB0);

		// Set gamma for different nuclides ('H' for Hydrogen, 'N' for Nitrogen)
		void SetGamma(const char nuclide);

		// Display information about the environment
		void DisplayInfo() const;

		// Return the Larmor frequency (SFO1)
		double GetSFO1() const;
	};
}
