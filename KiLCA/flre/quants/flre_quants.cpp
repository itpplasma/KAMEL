/*! \file flre_quants.cpp
    \brief The implementation of flre_quants class.
*/

#include <string.h>

#include "constants.h"
#include "flre_zone.h"
#include "flre_quants.h"
#include "calc_flre_quants.h"
#include "eval_cond.h"
#include "spline.h"
#include "shared.h"
#include "mode.h"

/*******************************************************************/

flre_quants::flre_quants (flre_zone *Z)
{
qloc = NULL;
qint = NULL;

C = NULL;
K = NULL;

Y = NULL;
S = NULL;
R = NULL;
sidY = 0;

jaE = NULL;
jaEi = NULL;

cdlab = NULL;

set_all_known_quants_parameters ();

zone = Z;

//begin external variables:

double rtor = zone->sd->bs->rtor;

flreo = zone->flre_order;

int num_quants = zone->sd->os->num_quants;

int NK = zone->cp->NK;
int dimK = zone->cp->dimK;

int NC = zone->cp->NC;
int dimC = zone->cp->dimC;

dimx = zone->dim;

x = zone->r;

//end external variables

//check:
if (num_quants > num_tot)
{
    fprintf (stderr, "\nerror: flre_quants: illegal value of num_quants: %d\n", num_quants);
    exit (1);
}

path2linear = new char[1024];

eval_path_to_linear_data (zone->sd->path2project, zone->wd->m, zone->wd->n, zone->wd->olab, path2linear);

//factor arising from integration over cylinder surface:
vol_fac = (2.0*pi)*(2.0*pi)*rtor;

//binomial coefficients used for computing some of the quantities:
bico = new double[(flreo+1)*(flreo+1)];
binomial_coefficients (flreo, bico);

//Set parameters in the computational list:
dni = new int[num_tot]; //global indices of the quants from computational list

numq = 0;
for (int i=0; i<num_quants; i++)
{
    if (zone->sd->os->flag_quants[i])
    {
        dni[numq] = i;
        ind[i] = numq++;
    }
}

if (numq == 0) //no quants to compute
{
    fprintf (stdout, "\nwarning: flre_quants: no quants to compute!\n");
    return;
}

//arrays allocation:
qloc = new double *[numq];
qint = new double *[numq];

for (int i=0; i<numq; i++)
{
    qloc[i] = new double[dimq[dni[i]]*dimx];
    qint[i] = new double[dimq[dni[i]]*dimx];
}

//conductivity matrices (+derivs) for current r point:
flagC = 0;
C = new double[(NC+1)*(dimC)]; //for C derivs: 0...NC

flagK = 0;
K = new double[(NK+1)*(dimK)]; //for K derivs: 0...NK (> flreo)!!!

N = NC;

//for jr splines:
nY = 2*1*1*3;

Y = new double[(dimx)*(nY)];
S = new double[(N+1)*(dimx)*(nY)];
R = new double[(N+1)*(nY)];

spline_alloc_ (N, 1, dimx, x, S, &sidY);

jaE  = new double[dimx];
jaEi = new double[dimx];

cdlab = new double[(dimq[CURRENT_DENS])*(dimx)];

#if DEBUG_FLAG
fprintf(stdout, "\nThere are %d additional quantities to compute:\n", numq);
for (int i = 0; i < numq; i++) {
  fprintf(stdout, "\n%s: dim=%d ind=%d", name[dni[i]], dimq[dni[i]],
          ind[dni[i]]);
}
fprintf(stdout, "\n");
#endif
}

/*******************************************************************/

void flre_quants::set_all_known_quants_parameters (void)
{
//global array of various quants: parameters are fixed
num_tot = 8;

dimq = new int[num_tot];
name = new char * [num_tot];

flag = new int[num_tot];
ind =  new int[num_tot];

calc_quant = new calc_func[num_tot];
save_quant = new save_func[num_tot];

int i;

//current density:
i = 0;
CURRENT_DENS = i;
dimq[i] = 2*3*2*3; //{{{re,im},comp={r,s,p},type={0,1},spec={i,e,t}}
calc_quant[i] = &calc_current_density;
save_quant[i] = &save_current_density;
name[i] = new char[16];
strcpy (name[i], "current");
i++;

//absorbed power density:
ABS_POWER_DENS = i;
dimq[i] = 1*2*3;  //{{re},type={0,1},spec={i,e,t}}
calc_quant[i] = &calc_absorbed_power_density;
save_quant[i] = &save_absorbed_power;
name[i] = new char[16];
strcpy (name[i], "abs_pow");
i++;

//dissipated power density:
DISS_POWER_DENS = i;
dimq[i] = 1*3*3;  //{{re},type={0,1,2=1-0},spec={i,e,t}}
calc_quant[i] = &calc_dissipated_power_density;
save_quant[i] = &save_dissipated_power;
name[i] = new char[16];
strcpy (name[i], "dis_pow");
i++;

//kinetic flux density:
KIN_FLUX = i;
dimq[i] = 1*3*3;  //{{re},type={0,1,2=1-0},spec={i,e,t}}
calc_quant[i] = &calc_kinetic_flux;
save_quant[i] = &save_kinetic_flux;
name[i] = new char[16];
strcpy (name[i], "kin_flux");
i++;

//Poyinting flux:
POY_FLUX = i;
dimq[i] = 1;  //{re}
calc_quant[i] = &calc_poynting_flux;
save_quant[i] = &save_poynting_flux;
name[i] = new char[16];
strcpy (name[i], "poy_flux");
i++;

//total flux:
TOT_FLUX = i;
dimq[i] = 1;  //{re}
calc_quant[i] = &calc_total_flux;
save_quant[i] = &save_total_flux;
name[i] = new char[16];
strcpy (name[i], "tot_flux");
i++;

//number density:
NUMBER_DENS = i;
dimq[i] = 2*3;  //{{re,im},spec={i,e,t}}
calc_quant[i] = &calc_number_density;
save_quant[i] = &save_number_density;
name[i] = new char[16];
strcpy (name[i], "density");
i++;

//Lorentz torque density:
LOR_TORQUE_DENS = i;
dimq[i] = 1*3*3;  //{{{re},comp={r,th,z}},spec={i,e,t}}
calc_quant[i] = &calc_lorentz_torque_density;
save_quant[i] = &save_lorentz_torque;
name[i] = new char[16];
strcpy (name[i], "torque");
i++;

//add more quants if needed:
//...

//check:
if (i != num_tot)
{
    fprintf (stderr, "\nerror: set_all_known_quants_parameters: illegal value of num_tot = %d\n", num_tot);
}

//more global setttings:
for (i=0; i<num_tot; i++)
{
    flag[i] =  0; //the quantity is not computed yet
    ind[i]  = -1; //the quantity has no index in the computational list yet
}
}

/*******************************************************************/

void flre_quants::calculate_local_profiles (void)
{
//The function omputes radial densities of the selected quantities

//begin external variables:
cond_profiles *cp = zone->cp;

//end external variables

if (numq == 0) return; //no quants to compute

Dmax = flreo; //actually flreo-1, must be <= NK!!!

if (Dmax > zone->cp->NK)
{
    fprintf (stdout, "\nerror: calculate_local_profiles: Dmax > NK!!!");
    exit (1);
}

for (int k=0; k<dimx; k++)
{
    //set the null state of node k for the the profiles struct:
    set_null_node_state (k);

    //calc conductivity:
    eval_all_K_matrices (cp, 0, Dmax, r, K);
    flagK = 1;

    eval_all_C_matrices (cp, 0, 0, r, C);
    flagC = 1;

    //calc local quants from the computational list: (depend on values at one r point)
    for (int i=CURRENT_DENS; i<=TOT_FLUX; i++)
    {
        if (zone->sd->os->flag_quants[i]) calc_quant[i](this);
    }
}

//compute the rest which depends on the previous quants computed on the grid: (derivs)
calc_splines_for_current_density (this);

calculate_field_profiles_poy_test (this);

for (int k=0; k<dimx; k++)
{
    //set the null state of node k for the the profiles struct:
    set_null_node_state (k);

    //calc quants from the computational list: (depend on values at several r points)
    for (int i=NUMBER_DENS; i<=LOR_TORQUE_DENS; i++)
    {
        if (zone->sd->os->flag_quants[i]) calc_quant[i](this);
    }
}
}

/*******************************************************************/

void flre_quants::set_null_node_state (int k)
{
node = k;

r = zone->r[k];

for (int i=0; i<num_tot; i++) flag[i] = 0;

flagC = 0;
flagK = 0;
}

/*******************************************************************/

void flre_quants::calculate_integrated_profiles (void)
{
//!The function integrates various densities over cylinder volume

if (zone->sd->os->flag_quants[ABS_POWER_DENS])
{
    calc_absorbed_power_in_cylinder (this);
}

if (zone->sd->os->flag_quants[DISS_POWER_DENS])
{
    calc_dissipated_power_in_cylinder (this);
}

if (zone->sd->os->flag_quants[LOR_TORQUE_DENS])
{
    calc_lorentz_torque_on_cylinder (this);
}
}

/*******************************************************************/

void flre_quants::save_profiles (void)
{
for (int i=0; i<num_tot; i++) if (zone->sd->os->flag_quants[i] == 2) save_quant[i](this);
}

/*******************************************************************/

inline int binary_search (double x, const double *xa, int ilo, int ihi)
{
//warning: returns dim-2 for x>=xa[dim-1] - Ok for splines!
while (ihi > ilo+1)
{
    int i = (ihi+ilo)/2;
    if (xa[i] > x) ihi = i;
    else           ilo = i;
}
return ilo;
}

/*******************************************************************/
