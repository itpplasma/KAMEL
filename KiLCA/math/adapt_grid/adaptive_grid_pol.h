/*! \file
    \brief The declarations of functions for generation of adaptive grids according to accuracy of polymomial interpolation.
*/

#include <stdlib.h>
#include <stdio.h>
#include <cmath>
#include <complex>
#include <cstring>
#include <climits>

#include "code_settings.h"

using namespace std;

#ifndef ADA_GRID_POLYNOM_INCLUDE

#define ADA_GRID_POLYNOM_INCLUDE

/*****************************************************************************/

int adaptive_grid_polynom (void (*func)(double *, double *, void *p), void *p,
                           double a, double b, int dimy, int deg, int *xdim, double *eps,
                           double *x1, double *y1);

int adaptive_grid_polynom (void (*func)(double *, double *, void *p), void *p,
                           double a, double b, int dimy, int deg, int *xdim, double *eps,
                           int dim_err, int *ind_err, double *x1, double *y1);

int adaptive_grid_polynom_res (void (*func)(double *, double *, void *p), void *p,
                               double a, double b, int dimy, int deg, int *xdim, double *eps,
                               double r_res, double D, double eps_res,
                               int dim_err, int *ind_err, double *x1, double *y1);

int sparse_grid_polynom (double *x, int dimx, double *y, int dimy, int deg, double *eps_a,
                         double *eps_r, double step, int *dim, double *x1, double *y1, int *ind1);

inline void eval_interp_polynom (int deg, double *xg, double *yg, int dimy, double x, double *y);

inline void find_index_for_interp (int deg, double x, int dimx, const double *xa, int *ind);

inline void search_array (double x, int dimx, const double *xa, int *ind);

inline int binary_search (double x, const double *xa, int ilo, int ihi);

inline int sign (double x);

void inline func_interp (int ind, double *yout, double *x, int dimx, double *y, int dimy);

inline void eval_neville_polynom (const double *xa, const double *ya,
                                  int yshift, int deg, double x, double *R);

/*****************************************************************************/

#endif

/*****************************************************************************/

inline void eval_interp_polynom (int deg, double *xg, double *yg, int dimy, double x, double *y)
{
for (int j=0; j<dimy; j++)
{
    y[j] = 0.0e0;

    for (int n=0; n<=deg; n++)
    {
        double fac = 1.0e0;
        for (int i=0; i<=deg; i++)
        {
            if (i != n) fac *= (x - xg[i])/(xg[n]-xg[i]);
        }
        y[j] += yg[n*dimy+j]*fac;
    }
}
}

/*****************************************************************************/

inline void find_index_for_interp (int deg, double x, int dimx, const double *xa, int *ind)
{
*ind = min (*ind+(deg-1)/2, dimx-2);

search_array (x, dimx, xa, ind);

*ind = min (max (0, *ind-(deg-1)/2), dimx-deg-1);
}

/*****************************************************************************/

inline void search_array (double x, int dimx, const double *xa, int *ind)
{
//0 <= ind <= dimx-2
//warning: returns dim-2 for x>=xa[dim-1] - Ok for splines and interpolation
//be careful to change something here:
if (x < xa[*ind])
{
    *ind = binary_search (x, xa, 0, *ind);
}
else if (x >= xa[*ind+1])
{
    *ind = binary_search (x, xa, *ind, dimx-1);
}
else
{}
}

/*****************************************************************************/

inline int binary_search (double x, const double *xa, int ilo, int ihi)
{
//warning: returns dim-2 for x>=xa[dim-1] - Ok for splines!
while (ihi > ilo+1)
{
    int i = (ihi+ilo)/2;
    if (xa[i] > x) ihi = i;
    else ilo = i;
}
return ilo;
}

/*****************************************************************************/

inline int sign (double x)
{
return (x<0.0) ? (-1):((x==0) ? 0:1);
}

/*****************************************************************************/

inline void eval_neville_polynom (const double *xa, const double *ya, int yshift, int deg, double x, double *R)
{
/*
if (x < xa[0] || x > xa[deg])
{
    fprintf (stdout, "\nwarning: eval_neville_polynom: x is outside the array [%le, %le]: x=%le", xa[0], xa[deg], x);
}
*/

double p[deg+1][deg+1];

for (int d=0; d<=deg; d++) p[d][d] = ya[d*yshift];

for (int d=1; d<=deg; d++) //over polymomial order
{
    for (int i=0; i<=deg-d; i++)
    {
        int j = i + d;

        p[i][j] = ((x-xa[j])/(xa[i]-xa[j]))*p[i][j-1]-((x-xa[i])/(xa[i]-xa[j]))*p[i+1][j];
    }
}

*R = p[0][deg];
}

/*****************************************************************************/
