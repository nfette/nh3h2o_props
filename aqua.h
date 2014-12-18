// aqua.h
// This code is part of NH3H2O_PROPS.
// Copyright 2014 Nicholas W Fette.
// Please refer to doi:10.1063/1.556015

#ifndef AQUA_H
#define AQUA_H

class aqua
{
public:
    aqua();
    double R_m;

    double a_ideal[15];
    double theta_ideal[15];
    double t_ideal[15];

    double T_n_ideal;
    double V_n_ideal_inv;

    double T_c_01;
    double rho_c_01;
    double p_c_01;
    double M_1;

    double T_c_02;
    double rho_c_02;
    double p_c_02;
    double M_2;

    double k_V;
    double k_T;
    double alpha;
    double beta;
    double gamma;

    double a_depart[15];
    double t_depart[15];
    double d_depart[15];
    double e_depart[15];

    double Phi_ideal_01(double tau_ideal, double delta_ideal);
    double Phi_ideal_02(double tau_ideal, double delta_ideal);
    double Phi_ideal(double tau_ideal, double delta_ideal, double x);
    double DeltaPhi_residual(double tau, double delta, double x);
    double Phi_residual(double tau, double delta, double x);

};

#endif // AQUA_H
