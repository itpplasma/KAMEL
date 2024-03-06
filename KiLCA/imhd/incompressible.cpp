/*! \file incompressible.cpp
    \brief The definitions of functions declared in incompressible.h.
*/

#include <stdlib.h>
#include <stdio.h>
#include <cmath>
#include <cstring>
#include <climits>

#include <gsl/gsl_math.h>
#include <gsl/gsl_deriv.h>

#include "constants.h"
#include "shared.h"
#include "eval_back.h"
#include "imhd_zone.h"
#include "incompressible.h"

/*******************************************************************/

int rhs_incompressible (double r, const double y[], double dy[], void *params)
{
const imhd_zone *zone = (const imhd_zone *) params;

complex<double> A = Afunc_hi_inc (r, params);
complex<double> B = Bfunc_hi_inc (r, params);

//take A derivatives:
diff_params ps = {zone, 0};

gsl_function F = {&Afunc_hi_inc_part, &ps};

double h = 1.0e-3, abserr;
double der_re, der_im;

ps.part = 0;
gsl_deriv_central (&F, r, h, &der_re, &abserr);

ps.part = 1;
gsl_deriv_central (&F, r, h, &der_im, &abserr);

complex<double> dA = der_re + I*der_im;

complex<double> f = y[0] + I*y[1]; //r*dzeta
complex<double> g = y[2] + I*y[3]; //(r*dzeta)'

complex<double> df = g;
complex<double> dg = -(dA*g + B*f)/A;

//result:
dy[0] = real(df);
dy[1] = imag(df);
dy[2] = real(dg);
dy[3] = imag(dg);

return GSL_SUCCESS;
}

/*******************************************************************/

inline complex<double> Afunc_hi_inc (double r, void *params)
{
const imhd_zone *zone = (const imhd_zone *) params;

int m = zone->wd->m;
double kz = (zone->wd->n)/(zone->sd->bs->rtor);

complex<double> omega2 = (zone->wd->olab)*(zone->wd->olab);

double mdens;

eval_mass_density (r, zone->bp, &mdens);

double F = Ffunc (r, params);

complex<double> omega2_a = F*F/(4.0*pi*mdens); //omega^2_a

double rk2 = r*(kz*kz + m*m/(r*r));

return mdens*(omega2 - omega2_a)/rk2;
}

/*******************************************************************/

inline complex<double> Bfunc_hi_inc (double r, void *params)
{
const imhd_zone *zone = (const imhd_zone *) params;

int m = zone->wd->m;
double kz = (zone->wd->n)/(zone->sd->bs->rtor);

complex<double> omega2 = (zone->wd->olab)*(zone->wd->olab);

double mdens;

eval_mass_density (r, zone->bp, &mdens);

double F = Ffunc (r, params);

double Bt;
eval_Bt (r, zone->bp, &Bt);

complex<double> omega2_a = F*F/(4.0*pi*mdens);

double r3k2 = r*r*r*(kz*kz + m*m/(r*r));

complex<double> t1 = - mdens*(omega2 - omega2_a)/r;

complex<double> t2 = 4.0*kz*kz*Bt*Bt/(4.0*pi*r3k2)*(omega2_a/(omega2 - omega2_a));

double t3, h = 1.0e-3, abserr;

gsl_function FD = {&func_der, params};

gsl_deriv_central (&FD, r, h, &t3, &abserr);

return t1 + t2 + t3;
}

/*******************************************************************/

inline double Afunc_hi_inc_part (double r, void *params)
{
diff_params *dp = (diff_params *) params;

complex<double> A = Afunc_hi_inc (r, (void *)dp->zone);

if (dp->part == 0) return real (A);
else               return imag (A);
}

/*******************************************************************/

inline double Bfunc_hi_inc_part (double r, void *params)
{
diff_params *dp = (diff_params *) params;

complex<double> B = Bfunc_hi_inc (r, (void *)dp->zone);

if (dp->part == 0) return real (B);
else               return imag (B);
}

/*******************************************************************/

inline double dBfunc_hi_inc_part (double r, void *params)
{
diff_params *dp = (diff_params *) params;

gsl_function F = {&Bfunc_hi_inc_part, dp};

double h = 1.0e-3, abserr, der;

gsl_deriv_central (&F, r, h, &der, &abserr);

return der;
}

/*******************************************************************/

inline double dAfunc_hi_inc_part (double r, void *params)
{
diff_params *dp = (diff_params *) params;

gsl_function F = {&Afunc_hi_inc_part, dp};

double h = 1.0e-3, abserr, der;

gsl_deriv_central (&F, r, h, &der, &abserr);

return der;
}

/*******************************************************************/

inline double ddAfunc_hi_inc_part (double r, void *params)
{
diff_params *dp = (diff_params *) params;

gsl_function F = {&dAfunc_hi_inc_part, dp};

double h = 1.0e-3, abserr, der;

gsl_deriv_central (&F, r, h, &der, &abserr);

return der;
}

/*******************************************************************/

inline double func_der (double r, void *params)
{
const imhd_zone *zone = (const imhd_zone *) params;

int m = zone->wd->m;
double kz = (zone->wd->n)/(zone->sd->bs->rtor);

double G = Gfunc(r, params);

double Bt;
eval_Bt (r, zone->bp, &Bt);

double r2k2 = r*r*(kz*kz + m*m/(r*r));

return Bt*Bt/(4.0*pi*r*r) + 2.0*kz*Bt*G/(4.0*pi*r2k2);
}

/*******************************************************************/

int jac_incompressible (double r, const double y[], double *dfdy, double dfdt[], void *params)
{
const imhd_zone *zone = (const imhd_zone *) params;

complex<double> A = Afunc_hi_inc (r, params);
complex<double> B = Bfunc_hi_inc (r, params);

diff_params ps = {zone, 0};

ps.part = 0;
double dA_re = dAfunc_hi_inc_part (r, &ps);

ps.part = 1;
double dA_im = dAfunc_hi_inc_part (r, &ps);

complex<double> dA = dA_re + I*dA_im;

complex<double> C = - dA/A, D = - B/A;

double Cre = real(C), Cim = imag(C);
double Dre = real(D), Dim = imag(D);

dfdy[0]  = 0.0;   dfdy[1]  = 0.0;   dfdy[2]  = 1.0;   dfdy[3]  = 0.0;

dfdy[4]  = 0.0;   dfdy[5]  = 0.0;   dfdy[6]  = 1.0;   dfdy[7]  = 1.0;

dfdy[8]  = Dre;   dfdy[9]  = -Dim;  dfdy[10] = Cre;   dfdy[11] = -Cim;

dfdy[12] = Dim;   dfdy[13] = Dre;   dfdy[14] = Cim;   dfdy[15] = Cre;

ps.part = 0;
double ddA_re = ddAfunc_hi_inc_part (r, &ps);

ps.part = 1;
double ddA_im = ddAfunc_hi_inc_part (r, &ps);

complex<double> ddA = ddA_re + I*ddA_im;

ps.part = 0;
double dB_re = dBfunc_hi_inc_part (r, &ps);

ps.part = 1;
double dB_im = dBfunc_hi_inc_part (r, &ps);

complex<double> dB = dB_re + I*dB_im;

complex<double> dC = - (ddA*A - dA*dA)/A/A;

complex<double> dD = - (dB*A - B*dA)/A/A;

double dCre = real(dC), dCim = imag(dC);
double dDre = real(dD), dDim = imag(dD);

dfdt[0] = 0.0;
dfdt[1] = 0.0;
dfdt[2] = dDre*y[0] - dDim*y[1] + dCre*y[2] - dCim*y[3];
dfdt[3] = dDim*y[0] + dDre*y[1] + dCim*y[2] + dCre*y[3];

return GSL_SUCCESS;
}

/*******************************************************************/

void imhd_zone::calculate_basis_incompressible_gsl (void)
{
//for tests only, use general cvode version!

if (!(bc1 == BOUNDARY_CENTER || bc2 == BOUNDARY_IDEALWALL))
{
    fprintf (stdout, "\nimhd::calculate_basis_incompressible_gsl: zone=%d: error!", index);
    exit (1);
}

size_t Neq = 4; //2 real eqs, 2 imag eqs for (r*zeta), (r*zeta)'

const gsl_odeiv_step_type *T = gsl_odeiv_step_rk8pd; //the best I have found

gsl_odeiv_step *step = gsl_odeiv_step_alloc (T, Neq);

gsl_odeiv_control *control = gsl_odeiv_control_y_new (eps_abs, eps_rel);

gsl_odeiv_evolve *evolve = gsl_odeiv_evolve_alloc (Neq);

//gsl_odeiv_system sys = {&rhs_incompressible, NULL, Neq, fp};
gsl_odeiv_system sys = {&rhs_incompressible, &jac_incompressible, Neq, this};

double t, tfinal;

double y[4];

int ind;

if (bc1 == BOUNDARY_CENTER)
{
    t = r1;
    tfinal = r2;

    //y start for integrator
    y[0] = pow (t, abs(wd->m));
    y[1] = 0.0;
    y[2] = abs(wd->m)*pow (t, abs(wd->m)-1);
    y[3] = 0.0;

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
    fprintf (stdout, "\nimhd::calculate_basis_incompressible_gsl: error!");
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

syst[ib(iter, ind, 6, 0)] = y[0];
syst[ib(iter, ind, 6, 1)] = y[1];
syst[ib(iter, ind, 7, 0)] = y[2];
syst[ib(iter, ind, 7, 1)] = y[3];

state_to_EB_incompressible (grid[iter], syst+ib(iter, ind, 6, 0), syst+ib(iter, ind, 0, 0));

while (t != tfinal)
{
    status = gsl_odeiv_evolve_apply (evolve, control, step, &sys, &t, tfinal, &h, y);
    
    if (status != GSL_SUCCESS)
    {
      fprintf (stderr, "\ncalculate_basis_incompressible_gsl: ODE solver failed at r = %le", t);
      exit (1);
    }

    //store values:
    iter++;
    grid[iter] = t;

    syst[ib(iter, ind, 6, 0)] = y[0];
    syst[ib(iter, ind, 6, 1)] = y[1];
    syst[ib(iter, ind, 7, 0)] = y[2];
    syst[ib(iter, ind, 7, 1)] = y[3];

    state_to_EB_incompressible (grid[iter], syst+ib(iter, ind, 6, 0), syst+ib(iter, ind, 0, 0));

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
    fprintf (stdout, "\nimhd::calculate_basis_incompressible_gsl: error!");
    exit (1);
}

gsl_odeiv_evolve_free (evolve);
gsl_odeiv_control_free (control);
gsl_odeiv_step_free (step);

delete [] grid;
delete [] syst;
}

/*******************************************************************/

void imhd_zone::state_to_EB_incompressible (double r, double *state, double *EB)
{
//state vector (r*zeta), (r*zeta)'
//EB vector: Er, Et, Ez, Br, Bt, Bz, (r*zeta), (r*zeta)'

complex<double> rzeta  = state[0] + I*state[1]; //(r*zeta)
complex<double> drzeta = state[2] + I*state[3]; //(r*zeta)'

complex<double> zeta_r = rzeta/r;
complex<double> dzeta_r = (drzeta - zeta_r)/r;

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

double rho; eval_mass_density (r, bp, &rho);

complex<double> tB = 2.0*I*kz*B0t/(4.0*pi*r);

complex<double> omega2 = (wd->olab)*(wd->olab);

complex<double> omegaa2 = F*F/(4.0*pi*rho);

complex<double> comfac1 = zeta_r*F*tB/(k0k0*rho*(omega2-omegaa2));

complex<double> comfac2 = drzeta*I/(r*k0k0);

complex<double> zeta_t = - comfac1*kz + comfac2*kt;
complex<double> zeta_z =   comfac1*kt + comfac2*kz;

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

//check for divergence:
double divzeta = abs(drzeta/r + I*(kt*zeta_t + kz*zeta_z));

double tst = abs(drzeta/r) + abs(kt*zeta_t) + abs(kz*zeta_z);

if (divzeta/tst > 1.0e-15)
{
    fprintf (stdout, "\nimhd_zone::state_to_EB_incompressible: warning: r = %le\t|divzeta| = %le\ttst = %le", r, divzeta, tst);
}
}

/*******************************************************************/
