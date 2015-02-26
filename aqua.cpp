// aqua.cpp
// This code is part of NH3H2O_PROPS.
// Copyright 2015 Nicholas W Fette.
// Please refer to doi:10.1063/1.556015
// Tillner-Roth, R., & Friend, D. G. (1998). A Helmholtz free energy
// formulation of the thermodynamic properties of the mixture {water+ ammonia}.
// Journal of Physical and Chemical Reference Data, 27(1), 63-96.

#include "aqua.h"
#include <math.h>

aqua::aqua()
{
    R_m = 8.314471;          /* J mol^-1 K^-1 */
    a_ideal[1] = -7.720435;
    a_ideal[2] = 8.649358;
    a_ideal[3] = 3.00632;
    a_ideal[4] = 0.012436;
    a_ideal[5] = 0.97315;
    a_ideal[6] = 1.27950;
    a_ideal[7] = 0.96956;
    a_ideal[8] = 0.24873;
    a_ideal[9] = -16.444285;
    a_ideal[10] = 4.036946;
    a_ideal[11] = -1.0;
    a_ideal[12] = 10.69955;
    a_ideal[13] = -1.775436;
    a_ideal[14] = 0.82374034;
    theta_ideal[4] = 1.666;
    theta_ideal[5] = 4.578;
    theta_ideal[6] = 10.018;
    theta_ideal[7] = 11.964;
    theta_ideal[8] = 35.600;
    t_ideal[12] = 1./3.;
    t_ideal[13] = -3./2.;
    t_ideal[14] = -7./4.;
    T_n_ideal = 500.;          /* K */
    V_n_ideal_inv = 15000.;    /* mol m^-3 */

    /* Fluid 1: Water */
    T_c_01 = 647.096;          /* K */
    rho_c_01 = 322.;           /* kg m^-3 */
    p_c_01 = 22.064E6;         /* Pa */
    M_1 = 0.018015268;         /* kg mol^-1 */

    /* Fluid 2: Water */
    T_c_02 = 405.40;    /* K */
    rho_c_02 = 225.;    /* kg m^-3 */
    p_c_02 = 11.36E6;   /* Pa */
    M_2 = 0.01703026;   /* kg mol^-1 */

    /* Reducing functions */
    k_V = 1.2395117;
    k_T = 0.9648407;
    alpha = 1.125455;
    beta = 0.8978069;

    /* Departure function */
    gamma = 0.5248379;
    a_depart[1] = -1.855822E-02;
    a_depart[2] = 5.258010E-02;
    a_depart[3] = 3.552874E-10;
    a_depart[4] = 5.451379E-06;
    a_depart[5] = -5.998546E-13;
    a_depart[6] = -3.687808E-06;
    a_depart[7] = 0.2586192;
    a_depart[8] = -1.368072E-08;
    a_depart[9] = 1.226146E-02;
    a_depart[10] = -7.181443E-02;
    a_depart[11] = 9.970849E-02;
    a_depart[12] = 1.0584086E-03;
    a_depart[13] = -0.1963687;
    a_depart[14] = -0.7777897;
    t_depart[1] = 3./2.;
    t_depart[2] = 1./2.;
    t_depart[3] = 13./2.;
    t_depart[4] = 7./4.;
    t_depart[5] = 15.;
    t_depart[6] = 6;
    t_depart[7] = -1;
    t_depart[8] = 4;
    t_depart[9] = 7./2.;
    t_depart[10] = 0;
    t_depart[11] = -1.;
    t_depart[12] = 8;
    t_depart[13] = 15./2.;
    t_depart[14] = 4.;
    d_depart[1] = 4.;
    d_depart[2] = 5.;
    d_depart[3] = 15.;
    d_depart[4] = 12.;
    d_depart[5] = 12.;
    d_depart[6] = 15.;
    d_depart[7] = -4.;
    d_depart[8] = 15.;
    d_depart[9] = 4.;
    d_depart[10] = 5.;
    d_depart[11] = 6.;
    d_depart[12] = 10.;
    d_depart[13] = 6.;
    d_depart[14] = 2.;
    e_depart[2] = 1;
    e_depart[3] = 1;
    e_depart[4] = 1;
    e_depart[5] = 1;
    e_depart[6] = 2;
    e_depart[7] = 1;
    e_depart[8] = 1;
    e_depart[9] = 1;
    e_depart[10] = 1;
    e_depart[11] = 2;
    e_depart[12] = 2;
    e_depart[13] = 2;
    e_depart[14] = 2;
}

double aqua::Phi_ideal_01(double tau_ideal, double delta_ideal)
{
    double result = a_ideal[1]
            + a_ideal[2] * tau_ideal
            + a_ideal[3] * log(tau_ideal)
            + a_ideal[4] * log(1 - exp(- theta_ideal[4] * tau_ideal))
            + a_ideal[5] * log(1 - exp(- theta_ideal[5] * tau_ideal))
            + a_ideal[6] * log(1 - exp(- theta_ideal[6] * tau_ideal))
            + a_ideal[7] * log(1 - exp(- theta_ideal[7] * tau_ideal))
            + a_ideal[8] * log(1 - exp(- theta_ideal[8] * tau_ideal));
    return result;
}

double aqua::Phi_ideal_02(double tau_ideal, double delta_ideal)
{
    double result = a_ideal[9]
            + a_ideal[10] * tau_ideal
            + a_ideal[11] * log(tau_ideal)
            + a_ideal[12] * pow(tau_ideal, t_ideal[12])
            + a_ideal[13] * pow(tau_ideal, t_ideal[13])
            + a_ideal[14] * pow(tau_ideal, t_ideal[14]);
    return result;
}

/** Returns reduced Helmholtz ideal part, units J mol^-1
 * TODO: the log(delta_ideal) term may be part of the individual fluid.
 * See the compile warnings about delta_ideal.
 */
double aqua::Phi_ideal(double tau_ideal, double delta_ideal, double x)
{
    double tau_01 = tau_ideal;
    double tau_02 = tau_ideal;
    double delta_01 = delta_ideal;
    double delta_02 = delta_ideal;
    double term1 = log(delta_ideal);
    double oneminusx = 1.0 - x;
    double term2_right = Phi_ideal_01(tau_01, delta_01);
    double term2 = oneminusx * term2_right;
    double term3_right = Phi_ideal_02(tau_02, delta_02);
    double term3 = x * term3_right;
    double term4 = oneminusx * log(oneminusx);
    double term5 = x * log(x);
    double result = term1 + term2 + term3 + term4 + term5;
    return result;
}

double aqua::DeltaPhi_residual(double tau, double delta, double x)
{
    double rhs, divisor;
    divisor = x * (1.0 - pow(x, gamma));
    double term1, term2, term3, term4;
    term1 = a_depart[1] * pow(tau, t_depart[1]) * pow(delta, d_depart[1]);
    term2 = a_depart[2] * exp(- pow(delta, e_depart[2])) * pow(tau, t_depart[2]) * pow(delta, d_depart[2])
            + a_depart[3] * exp(- pow(delta, e_depart[3])) * pow(tau, t_depart[3]) * pow(delta, d_depart[3])
            + a_depart[4] * exp(- pow(delta, e_depart[4])) * pow(tau, t_depart[4]) * pow(delta, d_depart[4])
            + a_depart[5] * exp(- pow(delta, e_depart[5])) * pow(tau, t_depart[5]) * pow(delta, d_depart[5])
            + a_depart[6] * exp(- pow(delta, e_depart[6])) * pow(tau, t_depart[6]) * pow(delta, d_depart[6]);
    term3 = x * (
            a_depart[7] * exp(- pow(delta, e_depart[7])) * pow(tau, t_depart[7]) * pow(delta, d_depart[7])
            + a_depart[8] * exp(- pow(delta, e_depart[8])) * pow(tau, t_depart[8]) * pow(delta, d_depart[8])
            + a_depart[9] * exp(- pow(delta, e_depart[9])) * pow(tau, t_depart[9]) * pow(delta, d_depart[9])
            + a_depart[10] * exp(- pow(delta, e_depart[10])) * pow(tau, t_depart[10]) * pow(delta, d_depart[10])
            + a_depart[11] * exp(- pow(delta, e_depart[11])) * pow(tau, t_depart[11]) * pow(delta, d_depart[11])
            + a_depart[12] * exp(- pow(delta, e_depart[12])) * pow(tau, t_depart[12]) * pow(delta, d_depart[12])
            + a_depart[13] * exp(- pow(delta, e_depart[13])) * pow(tau, t_depart[13]) * pow(delta, d_depart[13])
            );
    term4 = a_depart[14] * x * x * exp(- pow(delta, e_depart[14])) * pow(tau, t_depart[14]) * pow(delta, d_depart[14]);
    rhs = term1 + term2 + term3 + term4;
    double result = divisor * rhs;
    return result;
}

/** Returns reduced Helmholtz energy, residual part,
 * with units of J mol^-1
 */
double aqua::Phi_residual(double tau, double delta, double x)
{
    throw -1;
    return 0;
}
