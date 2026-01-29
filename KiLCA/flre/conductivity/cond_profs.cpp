/*! \file cond_profs.cpp
    \brief The implementation of functions declared in cond_profs.h (version with simple interpolation used to generate adaptive radial grid).
*/

#include <stdlib.h>
#include <stdio.h>
#include <cmath>
#include <cstring>
#include <climits>

#include "cond_profs.h"
#include "calc_cond.h"
#include "core.h"
#include "calc_back.h"
#include "shared.h"
#include "flre_zone.h"

/*******************************************************************/

void alloc_conductivity_profiles_ (cond_profiles **cpptr)
{
cond_profiles *cp = new cond_profiles;
*cpptr = cp;
}

/*******************************************************************/

void cond_profiles::calc_and_spline_main_conductivity_profiles (const flre_zone *zone, int flag)
{
sd = zone->get_settings ();
bp = zone->get_background ();
wd = zone->get_wave_data ();

r1 = zone->r1;
r2 = zone->r2;

D = zone->D;
eps_out = zone->eps_out;
eps_res = zone->eps_res;

flag_back = new char[8];
strcpy (flag_back, sd->bs->flag_back);

path2linear = new char[1024];

eval_path_to_linear_data (zone->get_path_to_project(), wd->m, wd->n, wd->olab, path2linear);

flreo = zone->flre_order;

gal_corr = zone->gal_corr;

//binomial coefficients:
bico = new double[(flreo+1)*(flreo+1)];

binomial_coefficients (flreo, bico);

/*conductivity K profiles: spec=0:1, [type=0:dimt-1, p=0:flreo, q=0:flreo, i=0:2, j=0:2], (re, im)*/
NK = zone->N + flreo + 1; //splines degree for K matrices in rsp - odd!
dimt = 2; //2 types of K matrices
dimK = 2*(dimt)*(flreo+1)*(flreo+1)*3*3*2;

//arrays for temp adaptive grid:
int max_dim = zone->max_dim_c;

xt = new double[max_dim];
yt = new double[dimK*max_dim];

//arrays for adaptive grid:
double *xa = new double[max_dim];
double *ya = new double[max_dim];

//boundaries:
xa[0] = max (r1-1.0, bp->x[0]);
xa[2] = min (r2+1.0, bp->x[bp->dimx-1]);

if (xa[0] > r1 || xa[2] < r2)
{
    fprintf (stdout, "\nwarning: xa[0] or xa[2] is inside the zone.\n");
    fflush (stdout);
    exit(1);
}

if (wd->r_res != 0.0e0)
{
    xa[1] = 1.01*(wd->r_res);
}
else
{
    xa[1] = 0.5*(xa[0]+xa[2]);
}

int dimxa = 3;

ind = 0; //x grid index

for (int i=0; i<dimxa; i++) sample_cond_func (xa+i, ya+i, (void *)this);

double eps = eps_out;

calc_adaptive_1D_grid_4vector_ (sample_cond_func, (void *)this, &max_dim, &eps, &dimxa, xa, ya);

//check:
if (ind != dimxa)
{
    fprintf (stdout, "\nerror: calc_and_spline_main_conductivity_profiles: ind=%d\tdimxa=%d\n", ind, dimxa);
    fflush (stdout);
    exit (1);
}

delete [] xa;
delete [] ya;

dimx = dimxa;
x = new double[dimx];

K = new double[dimx*dimK];
CK = new double[dimx*dimK*(NK+1)];
RK = new double[(NK+1)*dimK];

if (flag == 0)
{
    calc_splines_for_K (this);
}
else
{
    //set_arrays_for_K (this);
    smooth_arrays_for_K (this);
}

delete [] xt;
delete [] yt;

if (zone->flag_debug > 1)
{
    save_K_matrices (this, -1, 0); //spec, type

    //save_K_matrices_fine (this, -1, 0, 10); //spec, type
}

if (flag) return;

NC = zone->N; //splines degree for C matrices

/*conductivity C profiles: spec=0:1, [type=0:dimt-1,s=0..2flreo, i=0:2, j=0:2], (re, im)*/
dimC = 2*dimt*(2*flreo+1)*3*3*2;
C = new double[dimx*dimC];
CC = new double[(NC+1)*dimx*dimC];
RC = new double[(NC+1)*dimC];

calc_splines_for_C (this);

if (zone->flag_debug > 1)
{
    save_C_matrices (this, -1, 0); //spec, type

    //save_C_matrices_fine (this, -1, 0, 10); //spec, type
}
}

/*******************************************************************/

void delete_conductivity_profiles_f_ (cond_profiles **cp_ptr)
{
delete *cp_ptr;
}

/*******************************************************************/

void calc_and_spline_conductivity_for_point_ (settings **sd_ptr, background **bp_ptr,
                          wave_data **wd_ptr, char *flag_back,
                          double *r, cond_profiles **cp_ptr)
{
//it is assumed that flre_sett and all other fortran modules are set properly!

alloc_conductivity_profiles_ (cp_ptr);

cond_profiles *cp = (cond_profiles *)(*cp_ptr);

cp->sd = (const settings *)  (*sd_ptr);
cp->bp = (const background *)(*bp_ptr);
cp->wd = (const wave_data *) (*wd_ptr);

cp->flag_back = new char[8];
cp->flag_back[0] = flag_back[0];

cp->path2linear = NULL; //not needed in this case

get_flre_order_ (&(cp->flreo));

get_gal_corr_ (&(cp->gal_corr));

//binomial coefficients:
cp->bico = new double[(cp->flreo+1)*(cp->flreo+1)];

binomial_coefficients (cp->flreo, cp->bico);

cp->NK = cp->sd->bs->N - (cp->flreo + 1); //flreo is assumed to be odd

cp->dimt = 2; //2 types of K matrices
cp->dimK = 2*(cp->dimt)*(cp->flreo+1)*(cp->flreo+1)*3*3*2;

int max_dim = max(3*(cp->NK + 1), 9); //minimum possible dimension for splines

cp->xt = new double[max_dim];
cp->yt = new double[max_dim*(cp->dimK)];

//arrays for adaptive grid:
double *xa = new double[max_dim];
double *ya = new double[max_dim];

//interval for splining around point *r:
double delta = 0.01*(cp->bp->x[cp->bp->dimx-1] - cp->bp->x[0]); //interval width
double rmin = fmax(cp->bp->x[0], *r-delta);                     //left boundary
double rmax = fmin(cp->bp->x[cp->bp->dimx-1], *r+delta);        //right boundary

//boundaries:
xa[0] = rmin;
xa[2] = rmax;
xa[1] = 0.5*(xa[0]+xa[2]);

int dimxa = 3;

cp->ind = 0; //x grid index

for (int i=0; i<dimxa; i++) sample_cond_func (xa+i, ya+i, (void *)cp);

double eps = 0.0e0;

calc_adaptive_1D_grid_4vector_ (sample_cond_func, (void *)cp, &max_dim, &eps, &dimxa, xa, ya);

//check:
if (cp->ind != dimxa)
{
    fprintf (stdout, "\nerror: calc_and_spline_conductivity_profiles_for_point: ind=%d\tdimxa=%d\n", cp->ind, dimxa);
    fflush (stdout);
    exit (1);
}

delete [] xa;
delete [] ya;

cp->dimx = dimxa;
cp->x = new double[cp->dimx];

cp->K = new double[(cp->dimx)*(cp->dimK)];
cp->CK = new double[(cp->dimx)*(cp->dimK)*(cp->NK+1)];
cp->RK = new double[(cp->NK+1)*(cp->dimK)];

calc_splines_for_K (cp);

delete [] cp->xt;
delete [] cp->yt;

/*conductivity C profiles: spec=0:1, [type=0:dimt-1,s=0..2flreo, i=0:2, j=0:2], (re, im)*/
cp->NC = cp->NK - (cp->flreo + 1); //flreo is assumed to be odd
cp->dimC = 2*(cp->dimt)*(2*(cp->flreo)+1)*3*3*2;
cp->C = new double[(cp->dimx)*(cp->dimC)];
cp->CC = new double[(cp->NC+1)*(cp->dimx)*(cp->dimC)];
cp->RC = new double[(cp->NC+1)*(cp->dimC)];

calc_splines_for_C (cp);
}

/*******************************************************************/
