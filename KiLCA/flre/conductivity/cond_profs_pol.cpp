/*! \file cond_profs_pol.cpp
    \brief The implementation of functions declared in cond_profs.h (version with polynomial interpolation used to generate adaptive radial grid).
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
#include "adaptive_grid_pol.h"
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

NK = zone->N + flreo + 1; //splines degree for K matrices in rsp - odd!

//binomial coefficients:
bico = new double[(flreo+1)*(flreo+1)];

binomial_coefficients (flreo, bico);

//conductivity K profiles: spec=0:1, [type=0:dimt-1, p=0:flreo, q=0:flreo, i=0:2, j=0:2], (re, im)
dimt = 2; //2 types of K matrices
dimK = 2*(dimt)*(flreo+1)*(flreo+1)*3*3*2;

//arrays for temp adaptive grid:
dimx = zone->max_dim_c;

xt = new double[dimx];
yt = new double[dimK*dimx];

double a = max (r1-1.0, bp->x[0]);
double b = min (r2+1.0, bp->x[bp->dimx-1]);

if (a > r1 || b < r2)
{
    fprintf (stdout, "\nwarning: a or b is inside the zone.\n");
    fflush (stdout);
    exit(1);
}

//indices for error control:
int l = 0, pmax = 0, qmax = 0;
int dim_err = 2*(pmax+1)*(qmax+1)*3*3*2;
int *ind_err = new int[dim_err];

for (int spec=0; spec<2; spec++)
{
    for (int type=0; type<1; type++)
    {
        for (int p=0; p<=pmax; p++)
       {
            for (int q=0; q<=qmax; q++)
            {
                for (int i=0; i<3; i++)
                {
                    for (int j=0; j<3; j++)
                    {
                        for (int part=0; part<2; part++)
                        {
                            ind_err[l++] = iKa (spec, type, p, q, i, j, part, 0);
                        }
                    }
                }
            }
        }
    }
}

if (l != dim_err)
{
    fprintf (stderr, "\nerror: wrong dimension of ind_err array: %d != %d", dim_err, l);
}

double epso = eps_out;
double epsi = eps_res;

//adaptive_grid_polynom (sample_cond_func_polynom, (void *)this, a, b, dimK, NK, &dimx, &epso, dim_err, ind_err, xt, yt);

adaptive_grid_polynom_res (sample_cond_func_polynom, (void *)this, a, b, dimK, NK, &dimx, &epso, wd->r_res, D, epsi, dim_err, ind_err, xt, yt);

delete [] ind_err;

//correct dimension dimx:
x = new double[dimx];
K = new double[(dimx)*(dimK)];
CK = new double[(dimx)*(dimK)*(NK+1)];
RK = new double[(NK+1)*(dimK)];

if (flag == 0)
{
    calc_splines_for_K_polynom (this);
}
else
{
    set_arrays_for_K_polynom (this);
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

    //save_C_matrices_fine (this, -1, 0, 5); //spec, type
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

cp->dimx = 3*(cp->NK + 1); //minimum possible dimension for splines

cp->xt = new double[(cp->dimx)];
cp->yt = new double[(cp->dimx)*(cp->dimK)];

//interval for splining around point *r:
double delta = 0.01*(cp->bp->x[cp->bp->dimx-1] - cp->bp->x[0]); //interval width
double a = fmax(cp->bp->x[0], *r-delta);                        //left boundary
double b = fmin(cp->bp->x[cp->bp->dimx-1], *r+delta);           //right boundary

//indices for error control:
int l = 0, pmax = 0, qmax = 0;
int dim_err = 2*(pmax+1)*(qmax+1)*3*3*2;
int *ind_err = new int[dim_err];

for (int spec=0; spec<2; spec++)
{
    for (int type=0; type<1; type++)
    {
        for (int p=0; p<=pmax; p++)
       {
            for (int q=0; q<=qmax; q++)
            {
                for (int i=0; i<3; i++)
                {
                    for (int j=0; j<3; j++)
                    {
                        for (int part=0; part<2; part++)
                        {
                            ind_err[l++] = cp->iKa (spec, type, p, q, i, j, part, 0);
                        }
                    }
                }
            }
        }
    }
}

if (l != dim_err)
{
    fprintf (stderr, "\nerror: wrong dimension of ind_err array: %d != %d", dim_err, l);
}

double eps = 0.0e0;

adaptive_grid_polynom (sample_cond_func_polynom, (void *)cp, a, b, cp->dimK, cp->NK, &(cp->dimx), &eps, dim_err, ind_err, cp->xt, cp->yt);

delete [] ind_err;

//correct dimension dimx:
cp->x = new double[cp->dimx];
cp->K = new double[(cp->dimx)*(cp->dimK)];
cp->CK = new double[(cp->dimx)*(cp->dimK)*(cp->NK+1)];
cp->RK = new double[(cp->NK+1)*(cp->dimK)];

calc_splines_for_K_polynom (cp);

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
