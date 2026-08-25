/**
 * @file binding_utils.h
 * @brief Utilities for non-competitive protein–ligand binding models.
 *
 * This header defines a collection of self-contained functions for estimating
 * equilibrium protein–ligand binding in systems without ligand competition.
 *
 *  - All parameters in this file represent equilibrium concentrations (μM)
 *    or equilibrium dissociation constants (Kd, μM).
 *  - No kinetic rate constants (kon, koff) are used here.
 *    Time-dependent or kinetic simulations are implemented separately in
 *    `kinetics_utils.cpp`.
 *
 * The models implemented here assume:
 *  - no competing ligands,
 *  - no cooperativity between binding sites,
 *  - mass–action equilibrium,
 *  - known total concentrations and dissociation constants.
 */

#pragma once

#include <iostream>
#include <cmath>

namespace BindingUtils
{

	/**
	 * @brief Compute exact fraction of protein sites bound.
	 *
	 * Solves the quadratic equation for single-site binding.
	 *
	 * @param protein_total Total protein concentration.
	 * @param ligand_total Total ligand concentration.
	 * @param kd Dissociation constant.
	 * @return Fraction of protein bound.
	 */
	double CalculateFractionBoundExact(double protein_total, double ligand_total, double kd);

	/**
	 * @brief Perform single-site binding analysis with output.
	 *
	 * Computes bound fraction and prints detailed results.
	 *
	 * @param protein_total Total protein concentration.
	 * @param ligand_total Total ligand concentration.
	 * @param kd Dissociation constant.
	 */
	void SingleSiteBinding(double protein_total, double ligand_total, double kd);

	/**
	 * @brief Fraction of occupied sites for a protein with two independent ligand-binding sites.
	 *
	 * Uses the standard independent-site formula:
	 *   theta_i = L / (Kd_i + L)
	 * Assumes ligand_total >> protein_total (so L_free ≈ L_total).
	 *
	 * @param protein_total Total protein concentration (μM)
	 * @param ligand_total Total ligand concentration (μM)
	 * @param kd1 Dissociation constant for site 1 (μM)
	 * @param kd2 Dissociation constant for site 2 (μM)
	 * @return Fraction of total binding sites occupied (0.0 to 1.0)
	 */
	double TwoSiteIndependent(double protein_total, double ligand_total, double kd1, double kd2);

	/**
	 * @brief Fraction of occupied sites for a protein with three independent ligand-binding sites.
	 *
	 * Uses the standard independent-site formula:
	 *   theta_i = L / (Kd_i + L)
	 * Assumes ligand_total >> protein_total (so L_free ≈ L_total).
	 *
	 * @param protein_total Total protein concentration (μM)
	 * @param ligand_total Total ligand concentration (μM)
	 * @param kd1 Dissociation constant for site 1 (μM)
	 * @param kd2 Dissociation constant for site 2 (μM)
	 * @param kd3 Dissociation constant for site 3 (μM)
	 * @return Fraction of total binding sites occupied (0.0 to 1.0)
	 */
	double ThreeSiteIndependent(double protein_total, double ligand_total, double kd1, double kd2, double kd3);

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
								   double kex_meas);

} // namespace BindingUtils
