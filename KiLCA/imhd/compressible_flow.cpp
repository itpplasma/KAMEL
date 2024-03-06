/*! \file compressible_flow.cpp
    \brief The definitions of functions declared in compressible_flow.h.
*/

#include <stdlib.h>
#include <stdio.h>
#include <cmath>
#include <cstring>
#include <climits>
#include <ctime>

#include <gsl/gsl_math.h>
#include <gsl/gsl_deriv.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_odeiv.h>
#include <gsl/gsl_errno.h>

#include "constants.h"
#include "shared.h"
#include "eval_back.h"
#include "imhd_zone.h"
#include "compressible_flow.h"

double adiabat = 5.0/3.0;
//double adiabat = 1.0e6;

/************************** RWM *****************************************/
/******************from Bondenson et al. (1996)*************************************************/

// Doppler shifted frequency

inline complex<double> omega_D (double r, void * params)
{
imhd_zone const * zone = (const imhd_zone *) params;

int m = zone->wd->m;
double kz = (zone->wd->n)/(zone->sd->bs->rtor);
complex<double> omega = zone->wd->olab;

double Vt; eval_Vt (r, zone->bp, &Vt);

double Vz; eval_Vz (r, zone->bp, &Vz);

return omega - (m/r)*Vt - kz*Vz;
}

/*******************************************************************/

 // function T

inline complex<double> T_flow (double r, void * params)
{
imhd_zone const * zone = (const imhd_zone *) params;

double mdens; eval_mass_density (r, zone->bp, &mdens);

double Bt; eval_Bt (r, zone->bp, &Bt);

double F = Ffunc (r,params);

complex<double> w = omega_D (r, params);

double Vt; eval_Vt (r, zone->bp, &Vt);

return F*Bt/(4.0*pi) + mdens*w*Vt;
}

/*******************************************************************/

//function A

inline complex<double> A_flow (double r, void * params)
{
imhd_zone const * zone = (const imhd_zone *) params;

double mdens; eval_mass_density (r, zone->bp, &mdens);

double F = Ffunc (r,params);

complex<double> w = omega_D (r, params);

return mdens*w*w - F*F/(4.0*pi);
}

/*******************************************************************/

// function Q

inline complex<double> Q_flow (double r, void * params)
{
imhd_zone const * zone = (const imhd_zone *) params;

double mdens; eval_mass_density (r, zone->bp, &mdens);

double F = Ffunc (r,params);

double Vt; eval_Vt (r, zone->bp, &Vt);

double Bt; eval_Bt (r, zone->bp, &Bt);

complex<double>  T = T_flow (r, params);

complex<double> w = omega_D (r, params);

complex<double> A = A_flow (r,params);

return mdens*(w*w*(Bt*Bt/(4.0*pi) - mdens*Vt*Vt) + 1.0/(4.0*pi)*(Bt*w + F*Vt)*(Bt*w + F*Vt));
}

/*******************************************************************/

// function S

inline complex<double> S_flow (double r, void * params)
{
imhd_zone const * zone = (const imhd_zone *) params;

double mdens; eval_mass_density (r, zone->bp, &mdens);

double F = Ffunc (r,params);

double Bt; eval_Bt (r, zone->bp, &Bt);

double Bz; eval_Bz (r, zone->bp, &Bz);

complex<double> w = omega_D (r, params);

double R[2]; eval_p_dp (r, zone->bp, R); double press = R[0];

return  ((Bt*Bt + Bz*Bz)/(4.0*pi) + adiabat*press)*mdens*w*w - adiabat*press*F*F/(4.0*pi);
}

/*******************************************************************/

// function C11

inline complex<double> C11_flow (double r, void * params)
{
imhd_zone const * zone = (const imhd_zone *) params;

double mdens; eval_mass_density (r, zone->bp, &mdens);

int m = zone->wd->m;

complex<double> w = omega_D (r, params);

complex<double> Q = Q_flow (r, params);

complex<double> S = S_flow (r, params);

complex<double> T = T_flow (r, params);

return  mdens*w*w*Q/(r*r) - 2.0*m*S*T/(r*r*r);
}

/*******************************************************************/

// function C12

inline complex<double> C12_flow (double r, void * params)
{
imhd_zone const * zone = (const imhd_zone *) params;

double mdens; eval_mass_density (r, zone->bp, &mdens);

int m = zone->wd->m;
double kz = (zone->wd->n)/(zone->sd->bs->rtor);

complex<double> S = S_flow (r, params);

complex<double> w = omega_D (r, params);

return  mdens*mdens*w*w*w*w - (kz*kz + m*m/(r*r))*S;
}

/*******************************************************************/

// function C21

inline complex<double> C21_flow (double r, void * params)
{
imhd_zone const * zone = (const imhd_zone *) params;

complex<double> S = S_flow (r, params);

complex<double> A = A_flow (r, params);

complex<double> T = T_flow (r, params);

complex<double> Q = Q_flow (r, params);

complex<double> C4 = C4_flow (r, params);

return A*S*C4/r - 4.0*S*T*T/(r*r*r) + Q*Q/(r*r*r);
}

/*******************************************************************/

// function C22

inline complex<double> C22_flow (double r, void * params)
{
return r*C11_flow (r, params);
}

/*******************************************************************/

// functio to be derived inside C4

inline double func_deriv_in_C4 (double r, void * params)
{
imhd_zone const * zone = (const imhd_zone *) params;

double mdens; eval_mass_density (r, zone->bp, &mdens);

double Vt; eval_Vt (r, zone->bp, &Vt);

double Bt; eval_Bt (r, zone->bp, &Bt);

return ((Bt*Bt)/(4.0*pi) - mdens*Vt*Vt)/(r*r);
}

/*******************************************************************/

// function C4

inline complex<double> C4_flow (double r, void * params)
{
imhd_zone const * zone = (const imhd_zone *) params;

double der, h = 1.0e-3, abserr;

gsl_function F = {&func_deriv_in_C4, params};

gsl_deriv_central (&F, r, h, &der, &abserr);

complex<double> A = A_flow (r, params);

return A + r*der;
}

/*******************************************************************/

// 2nd order ODE for eigenmode solving

int rhs_flow (double r, const double y[], double dy[], void * params)
{
const imhd_zone * zone = (const imhd_zone *) params;

complex<double> C11 = C11_flow (r, params);
complex<double> C12 = C12_flow (r, params);
complex<double> C21 = C21_flow (r, params);
complex<double> C22 = C22_flow (r, params);

complex<double> A = A_flow (r, params);
complex<double> S = S_flow (r, params);

complex<double> f = y[0] + I*y[1]; //(r*dzeta)
complex<double> g = y[2] + I*y[3]; //(pressure)

complex<double> df = (C11*f - C12*g)/(A*S)*r;
complex<double> dg = (C21*f - C22*g)/(A*S);

dy[0] = real(df);
dy[1] = imag(df);
dy[2] = real(dg);
dy[3] = imag(dg);

return GSL_SUCCESS;
}

/*******************************************************************/

void imhd_zone::calculate_basis_flow_gsl (void)
{
//for tests only, use general cvode version!

if (!(bc1 == BOUNDARY_CENTER || bc2 == BOUNDARY_IDEALWALL))
{
    fprintf (stdout, "\nimhd::calculate_basis_incompressible_gsl: zone=%d: error!", index);
    exit (1);
}

size_t Neq = 4; //2 real eqs, 2 imag eqs for (r*zeta), pressure

const gsl_odeiv_step_type *T = gsl_odeiv_step_rk8pd; //the best I have found

gsl_odeiv_step *step = gsl_odeiv_step_alloc (T, Neq);

gsl_odeiv_control *control = gsl_odeiv_control_y_new (eps_abs, eps_rel);

gsl_odeiv_evolve *evolve = gsl_odeiv_evolve_alloc (Neq);

//gsl_odeiv_system sys = {&rhs_incompressible, NULL, Neq, fp};
gsl_odeiv_system sys = {&rhs_flow, &jac_flow, Neq, this};

double t, tfinal;

double y[4];

int ind;

if (bc1 == BOUNDARY_CENTER)
{
    t = r1;
    tfinal = r2;

    //y start for integrator: r*zeta ~ r^|m|
    y[0] = pow(t, abs(wd->m));
    y[1] = 0.0;

    complex<double> rzeta = pow (t, abs(wd->m));

    complex<double> drzeta = abs(wd->m)*pow (t, abs(wd->m)-1);

    complex<double> pstar = (C11_flow(t, this)*rzeta - S_flow(t, this)*A_flow(t, this)*drzeta/t)/
                            C12_flow(t, this);
    y[2] = real (pstar);
    y[3] = imag (pstar);

    ind = 0; //index of a fundamental solution which is regular at the boundary
}
else if (bc2 == BOUNDARY_IDEALWALL)
{
    t = r2;
    tfinal = r1;

    y[0] = 0.0;
    y[1] = 0.0;
    y[2] = 1.0;
    y[3] = 0.0;

    ind = 0; //index of a fundamental solution which is regular at the boundary
}
else
{
    fprintf (stdout, "\nimhd::calculate_basis_compressible_gsl: error!");
    exit (1);
}

double h = (1.0e-6)*signum(tfinal-t);

size_t iter, status;

//arrays:
double *grid = new double[max_dim];

double *syst = new double[max_dim*Nwaves*Ncomps*2];

//initial point:
iter = 0;
grid[iter] = t;

syst[ib(iter, ind, 6, 0)] = y[0]; //real r*zeta for iter 0
syst[ib(iter, ind, 6, 1)] = y[1]; //imag r*zeta for iter 0
syst[ib(iter, ind, 7, 0)] = y[2]; //real pressure for iter 0
syst[ib(iter, ind, 7, 1)] = y[3]; //imag pressure for iter 0

state_to_EB_compressible_flow (grid[iter], syst+ib(iter, ind, 6, 0), syst+ib(iter, ind, 0, 0));

while (t != tfinal)
{
    status = gsl_odeiv_evolve_apply (evolve, control, step, &sys, &t, tfinal, &h, y);
    
    if (status != GSL_SUCCESS)
    {
      fprintf (stderr, "\ncalculate_basis_flow_gsl: ODE solver failed at r = %le", t);
      exit (1);
    }

    //store values:
    iter++;
    grid[iter] = t;

    syst[ib(iter, ind, 6, 0)] = y[0];
    syst[ib(iter, ind, 6, 1)] = y[1];
    syst[ib(iter, ind, 7, 0)] = y[2];
    syst[ib(iter, ind, 7, 1)] = y[3];

    state_to_EB_compressible_flow (grid[iter], syst+ib(iter, ind, 6, 0), syst+ib(iter, ind, 0, 0));

    if (iter == max_dim-1)
    {
        fprintf (stderr, "\nMaximum allowed iteration number is reached: iter=%d r=%le", iter, t);
        exit (1);
    }
}

dim = iter + 1;

//zone class data allocation:
r = new double[dim];

int len = dim*Nwaves*Ncomps*2;

basis = new double[len];

for (int i=0; i<len; i++) basis[i] = 0.0e0;

//rearrangement of the points:
if (bc1 == BOUNDARY_CENTER)
{
    for (int j=0; j<dim; j++)
    {
        r[j] = grid[j];

        for (int comp=0; comp<Ncomps; comp++)
        {
            basis[ib(j, ind, comp, 0)] = syst[ib(j, ind, comp, 0)];
            basis[ib(j, ind, comp, 1)] = syst[ib(j, ind, comp, 1)];
        }
    }
    //normalize and scale basis functions:
    int j = dim-1;
    normalize_imhd_basis_ (&Ncomps, &Nwaves, &dim, basis, &ind, &j);
}
else if (bc2 == BOUNDARY_IDEALWALL)
{
    int k = 0;
    for (int j=dim-1; j>-1; j--)
    {
        r[k] = grid[j];

        for (int comp=0; comp<Ncomps; comp++)
        {
            basis[ib(k, ind, comp, 0)] = syst[ib(j, ind, comp, 0)];
            basis[ib(k, ind, comp, 1)] = syst[ib(j, ind, comp, 1)];
        }

        k++;
    }
    //normalize and scale basis functions:
    int j = 0;
    normalize_imhd_basis_ (&Ncomps, &Nwaves, &dim, basis, &ind, &j);
}
else
{
    fprintf (stdout, "\nimhd::calculate_basis_compressible_gsl: error!");
    exit (1);
}

gsl_odeiv_evolve_free (evolve);
gsl_odeiv_control_free (control);
gsl_odeiv_step_free (step);

delete [] grid;
delete [] syst;
}

/*******************************************************************/

void imhd_zone::state_to_EB_compressible_flow (double r, double * state, double * EB)
{
//state vector (r*zeta), pressure
//EB vector: Er, Et, Ez, Br, Bt, Bz, (r*zeta), (pressure)

complex<double> rzeta = state[0] + I*state[1]; //(r*zeta)
complex<double> press = state[2] + I*state[3]; //(pressure)

complex<double> zeta_r =  rzeta/r;

//wave:
double kt = (wd->m)/r;
double kz = (wd->n)/(sd->bs->rtor);
double k0k0 = kt*kt + kz*kz;

//background:
double R[4];

eval_Bt_dBt_Bz_dBz (r, bp, R);

double B0t  = R[0];
double dB0t = R[1];
double B0z  = R[2];
double dB0z = R[3];

double F = kt*B0t + kz*B0z; //double G = kt*B0z - kz*B0t;

double mdens; eval_mass_density (r, bp, &mdens);

double Vt; eval_Vt (r, bp, &Vt);

double Vz; eval_Vz (r, bp, &Vz);

// frequency omega:
complex<double> omega = (wd->olab);

complex<double> omega_Dop = omega - kt*Vt - kz*Vz;

double K[2]; eval_p_dp (r, bp, K); double press_0 = K[0]; double dpress_0 = K[1];

complex<double> drzeta = r*(C11_flow(r, this)*rzeta - C12_flow(r, this)*press)/
                           (S_flow(r, this)*A_flow(r, this));

//////////////////////////////////////////////////

complex<double> A = - kz*adiabat*press_0*kt + B0t*B0z*k0k0/4.0/pi;

complex<double> M = zeta_r*I*(F/4.0/pi*dB0z + kz*dpress_0 - kt/4.0/pi*dB0z*B0t -

                              kz/4.0/pi*B0t*(B0t/r - dB0t)) +

	            drzeta*I*(kz*B0t*B0t/4.0/pi/r - kt/4.0/pi*B0t*B0z/r + kz*adiabat*press_0/r);

complex<double> G = zeta_r*I*(F/4.0/pi*(B0t/r + dB0t) + 2.0*omega_Dop*mdens*Vt/r + dpress_0*kt +

                              kt/4.0/pi*B0z*dB0z + kz/4.0/pi*B0z*(B0t/r - dB0t)) +

  		    drzeta*I*(kt/4.0/pi*B0z*B0z/r - kz/4.0/pi/r*B0z*B0t + kt*adiabat*press_0/r);

complex<double> N = -mdens*omega_Dop*omega_Dop + kz*kz*adiabat*press_0 + B0t*B0t/4.0/pi*k0k0;

complex<double> H = -mdens*omega_Dop*omega_Dop + kt*kt*adiabat*press_0 + B0z*B0z/4.0/pi*k0k0;

///////////////////////////////////////////////////

complex<double> dzeta_r = (drzeta - zeta_r)/r;

complex<double> zeta_t = (M*A+G*N)/(N*H-A*A);

complex<double> zeta_z = (G*A+M*H)/(N*H-A*A);

// perturbated fields:
complex<double> alpha = zeta_t*B0z - zeta_z*B0t;

//Br:
complex<double> Br = I*F*zeta_r;

//Bt:
complex<double> Bt = - (dzeta_r*B0t + zeta_r*dB0t) + I*kz*alpha;

//Bz:
complex<double> Bz = - (drzeta*B0z + rzeta*dB0z)/r - I*kt*alpha;

//Er:
complex<double> Er = O;

//Et:
complex<double> Et = O;

//Ez:
complex<double> Ez = O;

//setup output values:
EB[0]  = real (Er);     EB[1]  = imag (Er);
EB[2]  = real (Et);     EB[3]  = imag (Et);
EB[4]  = real (Ez);     EB[5]  = imag (Ez);
EB[6]  = real (Br);     EB[7]  = imag (Br);
EB[8]  = real (Bt);     EB[9]  = imag (Bt);
EB[10] = real (Bz);     EB[11] = imag (Bz);
EB[12] = real (rzeta);  EB[13] = imag (rzeta);
EB[14] = real (drzeta); EB[15] = imag (drzeta);

//check:
/*complex<double> divzeta = drzeta/r + I*(kt*zeta_t + kz*zeta_z);
fprintf (stdout, "\nr = %le\tdrzeta/r = %le %+lei", r, real(drzeta/r), imag(drzeta/r));
fprintf (stdout, "\nzeta_t = %le %+lei", real(zeta_t), imag(zeta_t));
fprintf (stdout, "\nzeta_z = %le %+lei", real(zeta_z), imag(zeta_z));
fprintf (stdout, "\ndivzeta = %le %+lei", real(divzeta), imag(divzeta));*/
}

/*******************************************************************/

int jac_flow (double r, const double y[], double *dfdy, double dfdt[], void *params)
{
fprintf (stdout, "\nwarning: jac_flow is called!");

return GSL_SUCCESS;
}

/*******************************************************************/
