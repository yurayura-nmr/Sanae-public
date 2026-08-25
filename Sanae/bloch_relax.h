/*!
 * @file bloch_relax.h
 * @brief Demonstration of Bloch free precession and relaxation (no exchange).
 *
 * Implements free precession and transverse relaxation (R2 decay) for two
 * independent spins in the rotating frame. This simple function highlights how
 * differences in T2 and chemical shifts of 2 states A and B influence signal evolution.
 *
 * @note This is a relaxation-only function. Exchange is handled in the
 *       Bloch–McConnell modules.
 */

#pragma once

/*!
 * @brief Calculate Bloch relaxation for two independent spins
 *
 * @param R2_A       Transverse relaxation rate (s^-1) for state A.
 * @param R2_B       Transverse relaxation rate (s^-1) for state B.
 * @param CSppm_A    Chemical shift of state A in ppm.
 * @param CSppm_B    Chemical shift of state B in ppm.
 * @param M0_A       Initial transverse magnetization of state A.
 * @param M0_B       Initial transverse magnetization of state B.
 * @param nucleus    NMR nucleus ('H' = 1H, 'N' = 15N). Default is 'H'.
 */
void bloch_relax(
    double R2_A = 20.0,
    double R2_B = 10.0,
    double CSppm_A = -1.0,
    double CSppm_B = 1.0,
    double M0_A = 0.5,
    double M0_B = 0.5,
    char nucleus = 'H');
