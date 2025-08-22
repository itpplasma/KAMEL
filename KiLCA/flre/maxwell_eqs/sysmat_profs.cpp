/*! \file sysmat_profs.cpp
    \brief The implementation of sysmat_profiles class.
*/

#include <stdlib.h>
#include <stdio.h>
#include <cmath>
#include <cstring>
#include <climits>
#include <cfloat>

#include <gsl/gsl_heapsort.h>

#include "sysmat_profs.h"
#include "eval_sysmat.h"
#include "shared.h"
#include "adaptive_grid.h"
#include "inout.h"
#include "lapack.h"
#include "flre_zone.h"

/*******************************************************************/

void alloc_sysmatrix_profiles_ (sysmat_profiles **spptr)
{
sysmat_profiles *sp = new sysmat_profiles;
*spptr = sp;
}

/*******************************************************************/

void sysmat_profiles::calc_and_spline_sysmatrix_profiles (const flre_zone *zone)
{
    //in constructor one should set all parameters depending on external world, and nowhere
    //else use a calls to external objects - this is especially apropriate if the class is small.

    //begin of settings section:

    flag_back = new char[8];
    strcpy (flag_back, zone->cp->flag_back); //background flag is the same

    path2linear = new char[1024];
    strcpy (path2linear, zone->cp->path2linear);

    Nwaves = zone->Nwaves;

    N = zone->cp->NC;

    int max_dim = zone->max_dim_c;

    double eps = zone->eps_out;

    int flag_debug = zone->flag_debug;

    double r1 = zone->r1;
    double r2 = zone->r2;
    double rm = zone->wd->r_res;

    //end of settings section

    dimM = 2*(Nwaves)*(Nwaves); //number of reals in the sys matrix

    R = new double[(N+1)*(dimM)]; //for evaluation

    ind = 0; //x grid index

    int dimxa = 3;
    double *xa = new double[dimxa];
    double *ya = new double[dimxa];

    //arrays for temp adaptive grid:
    xt = new double[max_dim];
    yt = new double[(dimM)*max_dim];

    //boundaries:
    xa[0] = r1;
    xa[2] = r2;



    if (rm != 0.0e0 && (rm > r1 && rm < r2)) xa[1] = 1.01*rm; else xa[1] = 0.5*(xa[0]+xa[2]);

    for (int i=0; i<dimxa; i++) sample_sysmat_func (xa+i, ya+i, (void *)this);

    calc_adaptive_1D_grid_4vector_ (sample_sysmat_func, (void *)this, &max_dim, &eps, &dimxa, xa, ya);

    //check:
    if (ind != dimxa)
    {
        fprintf (stdout, "\nerror: eval_and_spline_sysmatrix_profiles: ind=%d\tdimxa=%d\n", ind, dimxa);
    }

    delete [] xa;
    delete [] ya;

    dimx = dimxa;
    x = new double[dimx];
    M = new double[(dimx)*(dimM)];
    C = new double[(N+1)*(dimx)*(dimM)];

    /*sorting grid points*/
    size_t *perm = new size_t[dimx];

    gsl_heapsort_index (perm, xt, dimx, sizeof(double), compare_doubles);

    //rearanging arrays:
    for (int i=0; i<dimx; i++)
    {
        x[i] = xt[perm[i]];
        for (int j=0; j<dimM; j++) M[i+j*(dimx)] = yt[j+perm[i]*(dimM)];
    }

    delete [] perm;
    delete [] xt;
    delete [] yt;

    spline_alloc_ (N, 1, dimx, x, C, &sidM);

    if (DEBUG_FLAG) fprintf (stdout, "\neval_and_spline_sysmatrix_profiles: dimx = %d\tdimy = %d\tmax_err = %le\n", dimx, dimM, eps);

    int ierr;
    spline_calc_ (sidM, M, 0, dimM-1, NULL, &ierr);

    if (DEBUG_FLAG) fprintf (stdout, "\nM profiles are splined...\n");

    if (flag_debug > 1) save_M_matrix (10);
}

/*******************************************************************/

void sysmat_profiles::save_M_matrix (int dimf)
{
/*saving matrix on a fine grid:*/

char *filename = new char[1024];

int dimt = dimf*(dimx-1);

double *grid = new double[dimt];
double *vals = new double[dimM*dimt];

for (int i=0; i<dimx-1; i++)
{
    for (int k=0; k<dimf; k++)
    {
        int ind = k+dimf*i;

        double r = x[i] + k*(x[i+1] - x[i])/dimf;

        grid[ind] = r;

        //eval_diff_sys_matrix_ (sp, &r, vals+(dimM)*ind); //spline
        calc_diff_sys_matrix_ (&r, flag_back, vals+(dimM)*ind); //exact
    }
}

sprintf (filename, "%s%s", path2linear, "debug-data/amat");
save_cmplx_matrix (Nwaves, Nwaves, dimt, grid, vals, filename);

delete [] grid;
delete [] vals;
delete [] filename;
}

/*******************************************************************/

void sample_sysmat_func (double *r, double *f, void *p)
{
sysmat_profiles *sp = (sysmat_profiles *)(p);

sp->xt[sp->ind] = *r;

calc_diff_sys_matrix_ (r, sp->flag_back, sp->yt+(sp->dimM)*(sp->ind));

//target function:
*f = 0.0e0;

for (int j=0; j<sp->dimM; j++)
{
    (*f) += log(1.0+(sp->yt[j+(sp->dimM)*(sp->ind)])*(sp->yt[j+(sp->dimM)*(sp->ind)]));
}

(sp->ind)++;
}

/*******************************************************************/
