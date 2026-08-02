#include "Sanae.h"
#include "Environment.h"

#include <iostream>
#include <stdexcept>

namespace Sanae
{
	static constexpr double PI = 3.14159265358979323846;

	// ---------------------------------------------------------
	//  Set magnetic field B0 (Tesla)
	// ---------------------------------------------------------
	void Environment::SetB0(const double customB0)
	{
		B0_ = customB0;

		// Recalculate Larmor frequency based on current gamma
        w_mhz_ = (-gamma_rads_ * B0_ / (2.0 * PI)) * 1e-6;

        // Update cached SFO1_ from the new w_mhz_
        SetSFO1();
	}

	// ---------------------------------------------------------
	//  Set gyromagnetic ratio based on nuclide
	// ---------------------------------------------------------
	void Environment::SetGamma(const char nuclide)
	{
		switch (nuclide)
		{
		case 'H':
			gamma_ = 42.58; // MHz/T  (1H)
			break;

		case 'N':
			gamma_ = 4.31; // MHz/T  (15N)
			break;

		case 'C':
			gamma_ = 10.705; // MHz/T (13C)
			break;

		default:
			throw std::invalid_argument(
				"Invalid nuclide. Use 'H' (1H), 'N' (15N), or 'C' (13C).");
		}

		// Convert to (rad/s)/T
		gamma_rads_ = gamma_ * 2.0 * PI * 1e6;

		// Larmor frequency in MHz (signed)
		w_mhz_ = (-gamma_rads_ * B0_ / (2.0 * PI)) * 1e-6;

		SetSFO1();
	}

	// ---------------------------------------------------------
	//  Display current environment values
	// ---------------------------------------------------------
	void Environment::DisplayInfo() const
	{
		std::cout << "\n[ NMR Environment Parameters ]\n"
			<< "Gamma           [MHz/T]      = " << gamma_ << "\n"
			<< "Gamma   [(rad/s)/T]          = " << gamma_rads_ << "\n"
			<< "B0              [T]          = " << B0_ << "\n"
			<< "Larmor (w/2pi)  [MHz]        = " << w_mhz_ << "\n"
			<< "SFO1            [Hz]         = " << SFO1_ << "\n"
			<< "SFO1            [MHz]        = " << SFO1_ / 1e6
			<< "\n"
			<< std::endl;
	}

	// ---------------------------------------------------------
	//  Return Larmor frequency in Hz
	// ---------------------------------------------------------
	double Environment::GetSFO1() const
	{
		return SFO1_;
	}

} // namespace Sanae
