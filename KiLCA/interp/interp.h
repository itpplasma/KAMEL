/*! \file
    \brief The declarations of functions used for polynomial interpolation and adaptive grids generation.
*/

#ifndef INTERP_COMMON_INCLUDE

#define INTERP_COMMON_INCLUDE

#include <stdlib.h>
#include <stdio.h>
#include <cmath>
#include <complex>

using namespace std;

/*****************************************************************************/

struct interpolant
{
    int N;     //polynomial degree

    int type;  //type of the interpolation

    int *ind;  //last search index for faster spline evaluation

    int dimx; /*dimension of a profile x grid (all the same)*/
    const double *x; /*x grid for profiles*/

    double *C; /*matrix of splines coefficients*/

    double *BC; /*array of binomial coefficients*/

    double *fac; /*factors used in the splines evaluation*/

    int calc_spline_boundaries (int dimy, const double *y, double *W);
    int calc_spline_coefficients (int Imin, int dimy, const double *y, double *W);
};

/*****************************************************************************/

int make_adaptive_grid (void (*func)(double *, double *, void *p), void *p, double a, double b, int deg, int *dim, double *eps, int order, int debug, double *x1, double *y1);

int sparse_grid_polynom (int dimx, double *x, double *y, int deg, double *eps_a, double *eps_r, double step, int *dim, double *x1, double *y1, int *ind1);

inline void eval_neville_polynom (const double *xa, const double *ya, int deg, double x, int Dmin, int Dmax, double *R);

inline void eval_neville_polynom (int dim, const double *xg, const double *yg, int deg, double x, int Dmin, int Dmax, int *ind, double *R);

inline void eval_lagrange_polynom (int dim, const double *xg, const double *yg, int deg, double x, int Dmin, int Dmax, int *ind, double *R);

void inline func_interp (int ind, double *y, double *yout);

inline void find_index_for_interp (int deg, double x, int dimx, const double *xa, int *ind);

inline void search_array (double x, int dimx, const double *xa, int *ind);

inline int binary_search (double x, const double *xa, int ilo, int ihi);

inline int sign (double x);

/*****************************************************************************/

extern "C"
{
void eval_neville_polynom_ (int * dim, const double * xg, const double * yg, int * deg,
                            double * x, int * Dmin, int * Dmax, int * ind, double * R);

void eval_neville_polynom_ready_ (const double * xa, const double * ya, int * deg,
                                  double * x, int * Dmin, int * Dmax, double * R);
}

/*****************************************************************************/

#endif

/*****************************************************************************/

inline void eval_lagrange_polynom (int dim, const double *xg, const double *yg, int deg, double x, int Dmin, int Dmax, int *ind, double *R)
{
if (Dmax > 0)
{
    fprintf (stdout, "\nwarning: derivatives are not implemented yet!"); fflush (stdout);
}

if (*ind < 0 || *ind > dim-1) *ind = dim/2;

find_index_for_interp (deg, x, dim, xg, ind);

const double *xa = xg + *ind;
const double *ya = yg + *ind;

//derivatives are not implemented!!!
R[0] = 0.0e0;

for (int n=0; n<=deg; n++)
{
    double fac = 1.0e0;
    for (int i=0; i<=deg; i++)
    {
        if (i != n) fac *= (x - xa[i])/(xa[n]-xa[i]);
    }
    R[0] += ya[n]*fac;
}
}

/*****************************************************************************/

inline void eval_neville_polynom (int dim, const double *xg, const double *yg, int deg, double x, int Dmin, int Dmax, int *ind, double *R)
{
constexpr double tol = 1e-9;
if (x < xg[0] - tol || x > xg[dim-1] + tol)
{
    fprintf (stdout, "\nwarning: eval_neville_polynom: x is outside the array: x=%le", x);
}

if (*ind < 0 || *ind > dim-1) *ind = dim/2;

find_index_for_interp (deg, x, dim, xg, ind);

eval_neville_polynom (xg + *ind, yg + *ind, deg, x, Dmin, Dmax, R);
}

/*****************************************************************************/

inline void eval_neville_polynom (const double *xa, const double *ya, int deg, double x, int Dmin, int Dmax, double *R)
{
constexpr double tol = 1e-9;
if (x < xa[0] - tol || x > xa[deg] + tol)
{
    fprintf (stdout, "\nwarning: eval_neville_polynom: x is outside the array [%.20le, %.20le]: x=%.20le", xa[0], xa[deg], x);
}

double p[Dmax+2][deg+1][deg+1]; //0th derivative corresponds to p[1][][]

for (int d=0; d<Dmax+2; d++)
{
    for (int i=0; i<=deg; i++)
    {
        for (int j=0; j<=deg; j++)
        {
            p[d][i][j] = 0.0e0;
        }
    }
}

for (int d=0; d<=deg; d++)
{
    p[1][d][d] = ya[d];
}

for (int n=0; n<=Dmax; n++) //over derivative order
{
    for (int d=1; d<=deg; d++) //over polymomial order
    {
        for (int i=0; i<=deg-d; i++)
        {
            int j = i + d;

            p[n+1][i][j] = n*(p[n][i][j-1] - p[n][i+1][j])/(xa[i]-xa[j]) +
                             ((x-xa[j])/(xa[i]-xa[j]))*p[n+1][i][j-1] -
                             ((x-xa[i])/(xa[i]-xa[j]))*p[n+1][i+1][j];
        }
    }
}

//store results:
for (int n=Dmin; n<=Dmax; n++) R[n-Dmin] = p[n+1][0][deg];
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
    else           ilo = i;
}
return ilo;
}

/*****************************************************************************/

inline int sign (double x)
{
return (x<0.0) ? (-1):((x==0) ? 0:1);
}

/*****************************************************************************/
