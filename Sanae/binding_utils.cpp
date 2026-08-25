#include <iostream>
#include <cmath>
#include <algorithm>

#include "binding_utils.h"

namespace BindingUtils
{

	/**
	 * Calculate fraction of protein binding sites occupied using the exact quadratic solution.
	 *
	 * For equilibrium: P + L <-> PL
	 * Kd = [P][L]/[PL]
	 *
	 * Solves: [PL] = (P₀ + L₀ + Kd - √((P₀ + L₀ + Kd)² - 4P₀L₀)) / 2
	 *
	 * NOTE:
	 * This implementation returns the fraction of protein binding sites occupied,
	 * i.e., [PL] / P₀. While the quadratic expression for [PL] is symmetric in P₀ and L₀,
	 * the final normalization is not.
	 * Interchanging P₀ and L₀ will therefore change the result.
	 *
	 * A common mistake is to swap protein and ligand concentrations when calling this method,
	 * which leads to physically inconsistent or misleading occupancy values. Please ensure that:
	 *   - protein_total corresponds to P₀ (binding sites)
	 *   - ligand_total corresponds to L₀ (binding partner)
	 *
	 * @param protein_total Total protein concentration (μM)
	 * @param ligand_total Total ligand concentration (μM)
	 * @param kd Equilibrium dissociation constant (μM)
	 * @return Fraction of protein binding sites occupied (0.0 to 1.0)
	 */

	double CalculateFractionBoundExact(double protein_total, double ligand_total, double kd)
	{
		if (kd <= 0)
		{
			throw std::invalid_argument("Kd must be positive");
		}
		if (protein_total <= 0 || ligand_total <= 0)
		{
			throw std::invalid_argument("Protein and ligand concentrations must be positive");
		}

		// Quadratic formula solution for [PL] (protein-ligand complex concentration)
		double sum = protein_total + ligand_total + kd;
		double discriminant = sum * sum - 4.0 * protein_total * ligand_total;

		double complex_concentration = (sum - std::sqrt(discriminant)) / 2.0;

		// Fraction bound = [PL] / [P_total]
		return complex_concentration / protein_total;
	}

	/**
	 * Comprehensive single-site binding analysis with detailed output.
	 * Uses exact quadratic solution - valid for all concentration ranges.
	 *
	 * @param protein_total Total protein concentration (μM)
	 * @param ligand_total Total ligand concentration (μM)
	 * @param kd Equilibrium dissociation constant (μM)
	 */

	void SingleSiteBinding(double protein_total, double ligand_total, double kd)
	{
		if (kd <= 0)
		{
			throw std::invalid_argument("Kd must be positive");
		}
		if (protein_total <= 0 || ligand_total <= 0)
		{
			throw std::invalid_argument("Protein and ligand concentrations must be positive");
		}

		// Calculate exact fraction bound using quadratic solution
		double fraction_bound_exact = CalculateFractionBoundExact(protein_total, ligand_total, kd);

		// Calculate actual concentrations at equilibrium
		double complex_conc = fraction_bound_exact * protein_total; // [PL]
		double free_protein = protein_total - complex_conc;			// [P]
		double free_ligand = ligand_total - complex_conc;			// [L]

		// Assess validity of approximation (ligand should be >>10x protein for approximation)
		double excess_ratio = ligand_total / protein_total;

		std::cout << "\n========== SINGLE-SITE BINDING ANALYSIS ==========" << std::endl;
		std::cout << "Initial Conditions:" << std::endl;
		std::cout << "  Total Protein [P0] = " << protein_total << " μM" << std::endl;
		std::cout << "  Total Ligand  [L0] = " << ligand_total << " μM" << std::endl;
		std::cout << "  Kd                  = " << kd << " μM" << std::endl;
		std::cout << "  Ligand excess ratio = " << excess_ratio << "x" << std::endl;

		std::cout << "\nEquilibrium Concentrations (exact solution):" << std::endl;
		std::cout << "  Complex [PL] = " << complex_conc << " μM" << std::endl;
		std::cout << "  Free Protein [P] = " << free_protein << " μM" << std::endl;
		std::cout << "  Free Ligand  [L] = " << free_ligand << " μM" << std::endl;

		std::cout << "\nBinding Results:" << std::endl;
		std::cout << "  Fraction Bound (exact)  = " << fraction_bound_exact
				  << " (" << (fraction_bound_exact * 100.0) << "%)" << std::endl;

		std::cout << "=================================================" << std::endl;
	}

	double TwoSiteIndependent(double protein_total,
							  double ligand_total,
							  double kd1,
							  double kd2)
	{
		if (kd1 <= 0 || kd2 <= 0)
			throw std::invalid_argument("Kd values must be positive");

		if (protein_total <= 0 || ligand_total <= 0)
			throw std::invalid_argument("Concentrations must be positive");

		// Independent binding site model
		double theta1 = ligand_total / (kd1 + ligand_total);
		double theta2 = ligand_total / (kd2 + ligand_total);

		// Fraction of binding sites occupied
		double fraction_sites_bound = 0.5 * (theta1 + theta2);

		std::cout << "\n========== TWO-SITE INDEPENDENT BINDING ANALYSIS ==========" << std::endl;

		std::cout << "Initial Conditions:" << std::endl;
		std::cout << "  Total Protein [P0] = " << protein_total << " μM" << std::endl;
		std::cout << "  Total Ligand  [L0] = " << ligand_total << " μM" << std::endl;
		std::cout << "  Kd1                = " << kd1 << " μM" << std::endl;
		std::cout << "  Kd2                = " << kd2 << " μM" << std::endl;

		std::cout << "\nPer-Site Occupancy:" << std::endl;
		std::cout << "  Site 1 occupancy (θ1) = " << theta1
				  << " (" << theta1 * 100.0 << "%)" << std::endl;
		std::cout << "  Site 2 occupancy (θ2) = " << theta2
				  << " (" << theta2 * 100.0 << "%)" << std::endl;

		std::cout << "\nBinding Results:" << std::endl;
		std::cout << "  Avg Fraction of Sites Bound = " << fraction_sites_bound
				  << " (" << fraction_sites_bound * 100.0 << "%)" << std::endl;

		std::cout << "============================================================" << std::endl;

		return fraction_sites_bound;
	}

	double ThreeSiteIndependent(double protein_total,
								double ligand_total,
								double kd1,
								double kd2,
								double kd3)
	{
		if (kd1 <= 0 || kd2 <= 0 || kd3 <= 0)
			throw std::invalid_argument("Kd values must be positive");

		if (protein_total <= 0 || ligand_total <= 0)
			throw std::invalid_argument("Concentrations must be positive");

		// Independent binding site model
		double theta1 = ligand_total / (kd1 + ligand_total);
		double theta2 = ligand_total / (kd2 + ligand_total);
		double theta3 = ligand_total / (kd3 + ligand_total);

		// Fraction of binding sites occupied (average over 3 sites)
		double fraction_sites_bound = (theta1 + theta2 + theta3) / 3.0;

		std::cout << "\n========== THREE-SITE INDEPENDENT BINDING ANALYSIS ==========" << std::endl;

		std::cout << "Initial Conditions:" << std::endl;
		std::cout << "  Total Protein [P0] = " << protein_total << " μM" << std::endl;
		std::cout << "  Total Ligand  [L0] = " << ligand_total << " μM" << std::endl;
		std::cout << "  Kd1                = " << kd1 << " μM" << std::endl;
		std::cout << "  Kd2                = " << kd2 << " μM" << std::endl;
		std::cout << "  Kd3                = " << kd3 << " μM" << std::endl;

		std::cout << "\nPer-Site Occupancy:" << std::endl;
		std::cout << "  Site 1 occupancy (θ1) = " << theta1
				  << " (" << theta1 * 100.0 << "%)" << std::endl;
		std::cout << "  Site 2 occupancy (θ2) = " << theta2
				  << " (" << theta2 * 100.0 << "%)" << std::endl;
		std::cout << "  Site 3 occupancy (θ3) = " << theta3
				  << " (" << theta3 * 100.0 << "%)" << std::endl;

		std::cout << "\nBinding Results:" << std::endl;
		std::cout << "  Avg Fraction of Sites Bound = " << fraction_sites_bound
				  << " (" << fraction_sites_bound * 100.0 << "%)" << std::endl;

		std::cout << "==============================================================" << std::endl;

		return fraction_sites_bound;
	}

	/**
	 * Calculate the intrinsic off-rate (koff) from a measured exchange rate (kex)
	 * at a specific anchor point where protein_total and ligand_total are known.
	 *
	 * This uses the exact quadratic solution for the free ligand concentration:
	 *   [PL] = (P₀ + L₀ + Kd - sqrt((P₀ + L₀ + Kd)² - 4P₀L₀)) / 2
	 *   L_free = L₀ - [PL]
	 *
	 * Then: koff = kex / (1 + L_free / Kd)
	 *
	 * @param protein_total Total protein concentration at the anchor (μM)
	 * @param ligand_total  Total ligand concentration at the anchor (μM)
	 * @param kd            Dissociation constant (μM)
	 * @param kex_meas      Measured exchange rate at the anchor (s⁻¹)
	 * @return              Intrinsic off-rate koff (s⁻¹)
	 * @throws std::invalid_argument if inputs are invalid (non‑positive or kex ≤ 0)
	 */

	double CalculateOffRateFromKex(double protein_total,
								   double ligand_total,
								   double kd,
								   double kex_meas)
	{
		// --- Input validation ---
		if (kd <= 0.0)
			throw std::invalid_argument("Kd must be positive");
		if (protein_total <= 0.0 || ligand_total <= 0.0)
			throw std::invalid_argument("Protein and ligand concentrations must be positive");
		if (kex_meas <= 0.0)
			throw std::invalid_argument("kex_meas must be positive");

		// --- Quadratic solution for [PL] (same as TitrationPoint::Compute) ---
		double term = protein_total + ligand_total + kd;
		double discriminant = term * term - 4.0 * protein_total * ligand_total;
		if (discriminant < 0.0)
			discriminant = 0.0; // clamp to avoid numerical issues

		double PL = 0.5 * (term - std::sqrt(discriminant));
		// Clamp PL to physical bounds
		if (PL < 0.0)
			PL = 0.0;
		if (PL > std::min(protein_total, ligand_total))
			PL = std::min(protein_total, ligand_total);

		// --- Free ligand concentration ---
		double L_free = ligand_total - PL;
		if (L_free < 0.0)
			L_free = 0.0;

		// --- Extract koff ---
		double koff = kex_meas / (1.0 + L_free / kd);

		// --- Print results ---
		std::cout << "\n========================================\n";
		std::cout << "  KINETIC RATE ESTIMATION (from kex)\n";
		std::cout << "========================================\n";
		std::cout << "Input parameters:\n";
		std::cout << "  Protein at anchor  = " << protein_total << " µM\n";
		std::cout << "  Ligand at anchor   = " << ligand_total << " µM\n";
		std::cout << "  Kd                 = " << kd << " µM\n";
		std::cout << "  kex (measured)     = " << kex_meas << " s⁻¹\n";
		std::cout << "\nComputed equilibrium:\n";
		std::cout << "  [PL] (complex)     = " << PL << " µM\n";
		std::cout << "  L_free             = " << L_free << " µM\n";
		std::cout << "\nResult:\n";
		std::cout << "  koff (intrinsic)   = " << koff << " s⁻¹\n";
		std::cout << "  kon                = " << koff / kd << " µM⁻¹s⁻¹\n";
		std::cout << "========================================\n\n";

		return koff;
	}

} // namespace BindingUtils
