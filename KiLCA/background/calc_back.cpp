/*! \file calc_back.cpp
    \brief The functions used to calculate equilibrium background profiles in cylindrical geometry.
*/

#include <stdlib.h>
#include <stdio.h>
#include <cmath>
#include <cstring>
#include <climits>
#include <inttypes.h>

#include <gsl/gsl_math.h>
#include <gsl/gsl_deriv.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_odeiv.h>
#include <gsl/gsl_errno.h>

#include "shared.h"
#include "settings.h"
#include "background.h"
#include "slatec.h"
#include "spline.h"
#include "inout.h"
#include "calc_back.h"
#include "eval_back.h"
#include "constants.h"
#include "conduct_parameters.h"

/*-----------------------------------------------------------------*/

int rhs_back (double r, const double y[], double dy[], void * params)
{
background * bp = (background *) params;

double rtor = bp->sd->bs->rtor;

spline_eval_ (bp->sid, 1, &r, 0, 1, bp->i_q, bp->i_T[1], bp->R); //R[n-Dmin+(j-Imin)*D1+k*D2], k=0 - new

double q = bp->R[0];

double n = bp->R[2];
double dn = bp->R[3];

double Ti = bp->R[4];
double dTi = bp->R[5];

double Te = bp->R[6];
double dTe = bp->R[7];

double dpress = boltz*((Ti+Te)*dn+(dTi+dTe)*n);

double g = 1.0 + r*r/rtor/rtor/q/q;

dy[0] = - 2.0*r*y[0]/(q*q*g*rtor*rtor) - 8.0*pi*dpress;

return GSL_SUCCESS;
}

/*-----------------------------------------------------------------*/

int jac_back (double r, const double y[], double *dfdy, double dfdt[], void *params)
{
fprintf (stderr, "\nwarning: impossible code flow: jac_back is called!");

return GSL_SUCCESS;
}

/*-----------------------------------------------------------------*/

int background::calculate_equilibrium (void)
{
/*! \fn calculate_equilibrium (void)
    \brief Calculates magnetic field and plasma current densities out of the input profiles.

    \param
*/

size_t Neq = 1;

const gsl_odeiv_step_type * T = gsl_odeiv_step_rk8pd; //the best I have found

gsl_odeiv_step * step = gsl_odeiv_step_alloc (T, Neq);

gsl_odeiv_control * control = gsl_odeiv_control_y_new (1.0e-16, 1.0e-16);

gsl_odeiv_evolve * evolve = gsl_odeiv_evolve_alloc (Neq);

gsl_odeiv_system sys = {&rhs_back, &jac_back, Neq, this};

double rtor = sd->bs->rtor;
double B0   = sd->bs->B0;

double rc = x[0], rf = x[1];

spline_eval_ (sid, 1, &rc, 0, 0, i_q, i_q, R);

double g = 1.0 + rc*rc/rtor/rtor/(R[0])/(R[0]); //starting value of u function

double * u = new double[dimx];

u[0] = B0*B0*g;

double uval[1] = {u[0]}; //current u value

double h = rf-rc; //initial step size

for (int i = 1; i < dimx; ++i)
{
    rf = x[i];

    size_t status, iter = 0;

    while (rc < rf)
    {
        status = gsl_odeiv_evolve_apply (evolve, control, step, &sys, &rc, rf, &h, uval);

        if (status != GSL_SUCCESS)
        {
            fprintf (stderr, "\ncalculate_equilibrium: ODE solver failed at r = %le", rc);
            exit (1);
        }

        if (iter == 100000)
        {
            fprintf (stderr, "\nMaximum allowed iteration number is reached: iter=%d r=%le", iter, rc);
            exit (1);
        }

        ++iter;
    }

    u[i] = uval[0];
}

gsl_odeiv_evolve_free (evolve);
gsl_odeiv_control_free (control);
gsl_odeiv_step_free (step);

if (sd->bs->flag_debug > 1) //save u if needed
{
    char *full_name = new char[1024];
    sprintf (full_name, "%s%s", path2background, "u.dat");
    save_real_array (dimx, x, u, full_name);
    delete [] full_name;
}

/*calc bth, bz, b0*/
//corresponding pointers:
double *q = y+(i_q)*dimx;
double *bth = y+(i_Bth)*dimx;
double *bz = y+(i_Bz)*dimx;
double *b0 = y+(i_B)*dimx;

for (int i = 0; i < dimx; ++i)
{
    bz[i] = signum(B0)*sqrt(u[i]/(1.0+(x[i])*(x[i])/rtor/rtor/q[i]/q[i]));

    bth[i] = bz[i]*(x[i])/q[i]/rtor;

    b0[i] = sqrt (bth[i]*bth[i]+bz[i]*bz[i]);
}

/*calculation hth hz:*/
double *hth = y+(i_hth)*dimx;
double *hz = y+(i_hz)*dimx;

for (int i = 0; i < dimx; ++i)
{
    hth[i] = bth[i]/b0[i];
    hz[i]  = bz[i]/b0[i];
}

int ierr;

/*splining bth, bz, b0, hth, hz: stored in a part of y array:*/
spline_calc_ (sid, bth, i_Bth, i_hz, NULL, &ierr);

/*calcs J0z, J0th*/
//corresponding pointers:
double *jth = y+(i_J0th[2])*dimx;
double *jz = y+(i_J0z[2])*dimx;

double dbz, dbth;

for (int i = 0; i < dimx; ++i)
{
    spline_eval_ (sid, 1, x+i, 1, 1, i_Bth, i_Bz, R);
    dbth = R[0];
    dbz = R[1];

    jth[i] = - c/4.0/pi*dbz;             /*J0th*/
    jz[i] = c/4.0/pi*(bth[i]/x[i]+dbth); /*J0z*/
}

/*splining jth, jz: stored in a part of y array:*/
spline_calc_ (sid, jth, i_J0th[2], i_J0z[2], NULL, &ierr);

delete [] u;

return 0;
}

/*-----------------------------------------------------------------*/

int background::check_and_spline_equilibrium (void)
{
/*splines input profiles for bth, bz, jth, jz and calcs the rest for equilibrium*/

/*calc bth, bz, b0*/
//corresponding pointers:
double *bth = y+(i_Bth)*dimx;
double *bz  = y+(i_Bz)*dimx;
double *b0  = y+(i_B)*dimx;
double *hth = y+(i_hth)*dimx;
double *hz  = y+(i_hz)*dimx;
double *jth = y+(i_J0th[2])*dimx;

int i;
for (i=0; i<dimx; i++)
{
    b0[i] = sqrt (bth[i]*bth[i]+bz[i]*bz[i]);
    hth[i] = bth[i]/b0[i];
    hz[i]  = bz[i]/b0[i];
}

int ierr;

/*splining bth, bz, b0, hth, hz, jth, jz: stored in a part of y array:*/
spline_calc_ (sid, bth, i_Bth, i_hz, NULL, &ierr);
spline_calc_ (sid, jth, i_J0th[2], i_J0z[2], NULL, &ierr);

//checks equilibrium:
double dev, Bt, dBt, Bz, dBz, p, dp;
for (i=0; i<dimx; i++)
{
    eval_Bt_dBt_Bz_dBz (x[i], this, R);
    Bt  = R[0];
    dBt = R[1];
    Bz  = R[2];
    dBz = R[3];

    eval_p_dp (x[i], this, R);
    p  = R[0];
    dp = R[1];

    dev = (4.0*pi*dp + Bt*dBt + Bz*dBz)/(Bt*Bt/x[i]) + 1.0;

    if (fabs(dev) > 1.0e-6)
    {
        fprintf (stderr, "\nwarning: check_and_spline_equilibrium: r = %le\tdev = %le", x[i], dev);
    }
}

return 0;
}

/*-----------------------------------------------------------------*/

int background::find_f0_parameters (void)
{
int ierr;

//renotation:
double *mass = sd->bs->mass;
double *charge = sd->bs->charge;
char *flag_back = sd->bs->flag_back;

double *n   = y + (i_n)*dimx;
double *Ti  = y + (i_T[0])*dimx;
double *Te  = y + (i_T[1])*dimx;
double *Vth = y + (i_Vth[2])*dimx;
double *Vz  = y + (i_Vz[2])*dimx;
double *Er  = y + (i_Er)*dimx;

double *b0  = y + (i_B)*dimx;
double *bth = y + (i_Bth)*dimx;
double *hth = y + (i_hth)*dimx;
double *hz  = y + (i_hz)*dimx;

double *jth = y + (i_J0th[2])*dimx;
double *jz  = y + (i_J0z[2])*dimx;

double *dPhi0 = y + (i_dPhi0)*dimx;

double *nui = y + (i_nu[0])*dimx;
double *nue = y + (i_nu[1])*dimx;

double *n_i_p = y + (i_n_p[0])*dimx;
double *n_e_p = y + (i_n_p[1])*dimx;

double *Vp_i_p = y + (i_Vp_p[0])*dimx;
double *Vp_e_p = y + (i_Vp_p[1])*dimx;

double *Vt_i_p = y + (i_Vt_p[0])*dimx;
double *Vt_e_p = y + (i_Vt_p[1])*dimx;

double Vs_tot, Vp_tot, Js_tot, Jp_tot;
double dpress, dn, dTi, dTe, Vsi_, Vse_;

/*estimations for the f0 parameters:*/
int i;

double nuee, nuei, nuie, nuii; //freqs
double Lee, Lei, Lie, Lii; //Coulomb logs
double vf = 1.0; //velocity factor

for (i=0; i<dimx; i++)
{
    /*collision frequences: NRL formulary: page 32 (in 2007), fast particles:*/
    //nui[i] = zion/sqrt(mass[0]/m_p)*pow(charge[0]/e, 4.0)*n[i]/pow(Ti[i], 1.5);
    //nue[i] = zele*pow(charge[0]/e, 2.0)*n[i]/pow(Te[i], 1.5);

    /*more correct collision frequences: NRL formulary: page 32 (in 2009)*/
    /*zion, zele is redefined so that true frequencies correspond to 1.0!*/
    Lee = 23.5 - log(sqrt(n[i])/pow(Te[i],1.25)) - (log(Te[i])-2.0)/4.0;
    Lei = 30.0 - log(sqrt(n[i])/pow(Ti[i],1.5)*pow(charge[0]/e, 2.0)/(mass[0]/m_p));
    Lie = Lei;
    Lii = 23.0 - log(pow(charge[0]/e,3.0)/pow(Ti[i],1.5)*sqrt(2.0*n[i]));

    //perpendicular coll freqs: (eps = T)
    nuee = (5.8e-6)*n[i]*Lee/sqrt(Te[i])/(vf*Te[i]);
    nuei = (7.7e-6)*n[i]*Lei*pow(charge[0]/e, 2.0)/pow(vf*Te[i], 1.5);
    nuie = (3.2e-9)*n[i]*Lie*pow(charge[0]/e, 2.0)/(mass[0]/m_p)/sqrt(Te[i])/(vf*Ti[i]);
    nuii = (1.4e-7)*n[i]*Lii*pow(charge[0]/e, 4.0)/sqrt(mass[0]/m_p)/sqrt(Ti[i])/(vf*Ti[i]);

    //sum up and add factors + limiting values:
    nui[i] = (sd->bs->zion)*(nuie + nuii + 10.0);
    nue[i] = (sd->bs->zele)*(nuee + nuei + 10.0);

    //density parameter:
    n_i_p[i] = n[i];
    n_e_p[i] = n[i];

    //thermal velocity:
    Vt_i_p[i] = sqrt (boltz*Ti[i]/mass[0]);
    Vt_e_p[i] = sqrt (boltz*Te[i]/mass[1]);

    //parallel velocity:
    Vp_tot = hth[i]*Vth[i] + hz[i]*Vz[i];
    Jp_tot = hth[i]*jth[i] + hz[i]*jz[i];

    Vp_i_p[i] = Vp_tot + Jp_tot*mass[1]/mass[0]/e/n[i]/(1.0+mass[1]/mass[0]);
    Vp_e_p[i] = Vp_tot - Jp_tot/e/n[i]/(1.0+mass[1]/mass[0]);

    if (flag_dPhi0_calc > 0) //recalculate dPhi0
    {
        //electric field: depends on a situation
        //Vs: velocity for either electrons or ions
        Vs_tot = hz[i]*Vth[i] - hth[i]*Vz[i];
        Js_tot = hz[i]*jth[i] - hth[i]*jz[i];

        Vsi_ = Vs_tot + Js_tot*mass[1]/mass[0]/e/n[i]/(1.0+mass[1]/mass[0]);
        Vse_ = Vs_tot - Js_tot/e/n[i]/(1.0+mass[1]/mass[0]);

        //one can only satisfy for the velocity of the one component or total:
        spline_eval_ (sid, 1, x+i, 1, 1, i_n, i_T[1], R);

        dn  = R[0];
        dTi = R[1];
        dTe = R[2];

        //from equation: Vsi(x[i]) = Vsi_;
        dpress = boltz*(Ti[i]*dn+dTi*n[i]);            //ions dpressure
        dPhi0[i] = Vsi_/c*b0[i]-dpress/charge[0]/n[i]; //el. field defined by ions
    }
    else //calculate dPhi0 from input Er profile
    {
        dPhi0[i] = - Er[i] + bth[i]*(sd->bs->V_gal_sys)/c;
    }
}

/*splining dPhi0, f0 parameters, nui, nue:*/
spline_calc_ (sid, dPhi0, i_dPhi0, i_nu[1], NULL, &ierr);

//density reestimation:
int spec;

for (i=0; i<dimx; i++)
{
    eval_and_set_background_parameters_spec_independent_ (x+i, flag_back, 1);

    for (spec=0; spec<2; spec++)
    {
        eval_and_set_background_parameters_spec_dependent_ (x+i, &spec, flag_back, 1);
        eval_and_set_f0_parameters_nu_and_derivs_ (x+i, &spec, flag_back, 1);
        dens_par_ (n+i, y+(i_n_p[spec])*dimx+i); //sets n_p out of n
    }
}

spline_calc_ (sid, n_i_p, i_n_p[0], i_n_p[0], NULL, &ierr);
spline_calc_ (sid, n_e_p, i_n_p[1], i_n_p[1], NULL, &ierr);

return ierr;
}

/*-----------------------------------------------------------------*/

int background::eval_and_save_f0_moments (void)
{
double *mass = sd->bs->mass;
double *charge = sd->bs->charge;
char *flag_back = sd->bs->flag_back;

const int num_moms = 16;
FILE **mom_files = new FILE *[num_moms];

char *full_name = new char[1024];

char moms_name[num_moms][10] = {{"ni_m"},  {"ne_m"},  {"Vth_m"}, {"Vz_m"},
				{"Vsi_m"}, {"Vse_m"}, {"Vpi_m"}, {"Vpe_m"},
				{"Vs0i_m"},{"Vs0e_m"},{"Ti_m"},  {"Te_m"},
				{"j0s_m"}, {"j0p_m"}, {"j0th_m"},{"j0z_m"}};

int i, k;

//fprintf (stdout, "\nReminder: all quantities are computed in the moving frame!!!");

/*opening all files*/
for (k=0; k<num_moms; k++)
{
	sprintf (full_name, "%s%s%s", path2background, moms_name[k], ".dat");
	if (!(mom_files[k] = fopen (full_name, "w")))
	{
		fprintf (stderr, "\nFailed to open file %s\a\n", full_name);
		continue;
	}
}

/*allocation: i,e,tot*/
double n_m[2], Vth_m[3], Vz_m[3], Vs_m[3], Vs_0_m[3], Vp_m[3], T_m[2];
double j0th_m[3], j0z_m[3], j0s_m[3], j0p_m[3];

int spec;

for (i=0; i<dimx; i++) /*over r points*/
{
    eval_and_set_background_parameters_spec_independent_ (&(x[i]), flag_back, 1);

    for (spec=0; spec<2; spec++) /*over specii*/
    {
        eval_and_set_background_parameters_spec_dependent_ (&(x[i]), &spec, flag_back, 1);

        eval_and_set_f0_parameters_nu_and_derivs_ (&(x[i]), &spec, flag_back, 1);

        //moments:
        dens_mom_(&n_m[spec]);

        vth_mom_(&Vth_m[spec]); /*it is contravariant component: r*Vth is physical*/

        Vth_m[spec] *= x[i]; /*this is physical*/

        vz_mom_(&Vz_m[spec]);

        vs_mom_(&Vs_m[spec]);

        vs_0 (x[i], spec, this, &Vs_0_m[spec]);

        vp_mom_(&Vp_m[spec]);

        eterm_mom_(&T_m[spec]);
        T_m[spec] = 2.0*T_m[spec]/3.0/n_m[spec]; /*up to const: kelvins*/
    }

    Vth_m[2] = (mass[0]*Vth_m[0]+mass[1]*Vth_m[1])/(mass[0]+mass[1]);
    Vz_m[2]  = (mass[0]*Vz_m[0] +mass[1]*Vz_m[1]) /(mass[0]+mass[1]);

    j0p_m[2] = e*(Vp_m[0]-Vp_m[1]);
    j0s_m[2] = e*(Vs_m[0]-Vs_m[1]);

    j0th_m[2]= e*(Vth_m[0]-Vth_m[1]);
    j0z_m[2] = e*(Vz_m[0]-Vz_m[1]);

    fprintf (mom_files[0], "%.15le\t%.15le\n", x[i], n_m[0]);
    fprintf (mom_files[1], "%.15le\t%.15le\n", x[i], n_m[1]);
    fprintf (mom_files[2], "%.15le\t%.15le\n", x[i], Vth_m[2]/n_m[0]);
    fprintf (mom_files[3], "%.15le\t%.15le\n", x[i], Vz_m[2]/n_m[0]);
    fprintf (mom_files[4], "%.15le\t%.15le\n", x[i], Vs_m[0]/n_m[0]);
    fprintf (mom_files[5], "%.15le\t%.15le\n", x[i], Vs_m[1]/n_m[1]);
    fprintf (mom_files[6], "%.15le\t%.15le\n", x[i], Vp_m[0]/n_m[0]);
    fprintf (mom_files[7], "%.15le\t%.15le\n", x[i], Vp_m[1]/n_m[1]);
    fprintf (mom_files[8], "%.15le\t%.15le\n", x[i], Vs_0_m[0]);
    fprintf (mom_files[9], "%.15le\t%.15le\n", x[i], Vs_0_m[1]);
    fprintf (mom_files[10], "%.15le\t%.15le\n", x[i], T_m[0]);
    fprintf (mom_files[11], "%.15le\t%.15le\n", x[i], T_m[1]);
    fprintf (mom_files[12], "%.15le\t%.15le\n", x[i], j0s_m[2]);
    fprintf (mom_files[13], "%.15le\t%.15le\n", x[i], j0p_m[2]);
    fprintf (mom_files[14], "%.15le\t%.15le\n", x[i], j0th_m[2]);
    fprintf (mom_files[15], "%.15le\t%.15le\n", x[i], j0z_m[2]);
}

for (k=0; k<num_moms; k++) fclose (mom_files[k]);

delete [] mom_files;

delete [] full_name;

return 0;
}

/*-----------------------------------------------------------------*/

int background::save_background (void)
{
/*saving data:*/
char *full_name = new char[1024];

sprintf (full_name, "%s%s", path2background, "b0th.dat");
save_real_array (dimx, x, y+(i_Bth)*dimx, full_name);

sprintf (full_name, "%s%s", path2background, "b0z.dat");
save_real_array (dimx, x, y+(i_Bz)*dimx, full_name);

sprintf (full_name, "%s%s", path2background, "b0.dat");
save_real_array (dimx, x, y+(i_B)*dimx, full_name);

sprintf (full_name, "%s%s", path2background, "hth.dat");
save_real_array (dimx, x, y+(i_hth)*dimx, full_name);

sprintf (full_name, "%s%s", path2background, "hz.dat");
save_real_array (dimx, x, y+(i_hz)*dimx, full_name);

sprintf (full_name, "%s%s", path2background, "j0th.dat");
save_real_array (dimx, x, y+(i_J0th[2])*dimx, full_name);

sprintf (full_name, "%s%s", path2background, "j0z.dat");
save_real_array (dimx, x, y+(i_J0z[2])*dimx, full_name);

sprintf (full_name, "%s%s", path2background, "nui.dat");
save_real_array (dimx, x, y+(i_nu[0])*dimx, full_name);

sprintf (full_name, "%s%s", path2background, "nue.dat");
save_real_array (dimx, x, y+(i_nu[1])*dimx, full_name);

sprintf (full_name, "%s%s", path2background, "dPhi0.dat");
save_real_array (dimx, x, y+(i_dPhi0)*dimx, full_name);

sprintf (full_name, "%s%s", path2background, "ni_p.dat");
save_real_array (dimx, x, y+(i_n_p[0])*dimx, full_name);

sprintf (full_name, "%s%s", path2background, "ne_p.dat");
save_real_array (dimx, x, y+(i_n_p[1])*dimx, full_name);

sprintf (full_name, "%s%s", path2background, "Vpi_p.dat");
save_real_array (dimx, x, y+(i_Vp_p[0])*dimx, full_name);

sprintf (full_name, "%s%s", path2background, "Vpe_p.dat");
save_real_array (dimx, x, y+(i_Vp_p[1])*dimx, full_name);

sprintf (full_name, "%s%s", path2background, "Vti_p.dat");
save_real_array (dimx, x, y+(i_Vt_p[0])*dimx, full_name);

sprintf (full_name, "%s%s", path2background, "Vte_p.dat");
save_real_array (dimx, x, y+(i_Vt_p[1])*dimx, full_name);

delete [] full_name;

return 0;
}

/*-----------------------------------------------------------------*/

inline void background::interp_basic_background_profiles (double r, double *q, double *n, double *Ti, double *Te, double *Vth, double *Vz, double *dPhi0)
{
spline_eval_ (sid, 1, &r, 0, 0, i_q, i_Vz[2], R);

*q   = R[0];
*n   = R[1];
*Ti  = R[2];
*Te  = R[3];
*Vth = R[4];
*Vz  = R[5];

spline_eval_ (sid, 1, &r, 0, 0, i_dPhi0, i_dPhi0, dPhi0);
}

/*-----------------------------------------------------------------*/

inline void background::transform_basic_background_profiles_to_lab_frame (double r, double *q, double *n, double *Ti, double *Te, double *Vth, double *Vz, double *dPhi0)
{
double Bth; spline_eval_ (sid, 1, &r, 0, 0, i_Bth, i_Bth, &Bth);

*Vz    += sd->bs->V_gal_sys;
*dPhi0 -= Bth*(sd->bs->V_gal_sys)/c;
}

/*-----------------------------------------------------------------*/

void background::interp_basic_background_profiles_in_lab_frame (int dim, double *r, double *q, double *n, double *Ti, double *Te, double *Vth, double *Vz, double *dPhi0)
{
for (int i=0; i<dim; i++)
{
    interp_basic_background_profiles (r[i], q+i, n+i, Ti+i, Te+i, Vth+i, Vz+i, dPhi0+i);

    transform_basic_background_profiles_to_lab_frame (r[i], q+i, n+i, Ti+i, Te+i, Vth+i, Vz+i, dPhi0+i);
}
}

/*-----------------------------------------------------------------*/
