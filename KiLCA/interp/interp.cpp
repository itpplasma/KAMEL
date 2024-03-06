/*! \file
    \brief The definitions of functions declared in interp.h.
*/

#include <stdlib.h>
#include <stdio.h>
#include <cmath>
#include <complex>
#include <cstring>
#include <climits>

#include "interp.h"

/*****************************************************************************/

void eval_neville_polynom_ (int * dim, const double * xg, const double * yg, int * deg,
                            double * x, int * Dmin, int * Dmax, int * ind, double * R)
{
eval_neville_polynom(*dim, xg, yg, *deg, *x, *Dmin, *Dmax, ind, R);
}

/*****************************************************************************/

void eval_neville_polynom_ready_ (const double * xa, const double * ya, int * deg,
                                  double * x, int * Dmin, int * Dmax, double * R)
{
eval_neville_polynom(xa, ya, *deg, *x, *Dmin, *Dmax, R);
}

/*****************************************************************************/

int make_adaptive_grid (void (*func)(double *, double *, void *p), void *p, double a, double b, int deg, int *dim, double *eps, int order, int debug, double *x1, double *y1)
{
//makes an adaptive x grid for moving interpolating polynoma of degree deg
//better to use odd degrees

//additional space:
double *x2 = new double[(*dim)];
double *y2 = new double[(*dim)];

double *xold = x1, *yold = y1, *xnew = x2, *ynew = y2;

//initial grid: for 3 "moving" polynoma
int dimx = deg + 3;

xold[0] = a;
xold[dimx-1] = b;

double pi = 3.141592653589793238;

for (int k=1; k<=deg+1; k++) //Chebyshev nodes
{
    xold[k] = 0.5*(a+b) - 0.5*(b-a)*cos((2*k-1)*pi/(2*(deg+1)));
}

//initial y grid:
for (int k=0; k<dimx; k++) func (xold+k, yold+k, p);

//grid refinement is based on comparision of results of 3 interpolating polynoma at centers
int ind = dimx/2;

double errl, errr;
double xc, yl, ym, yr;

//for maximum error:
double xmerr, ymerr, merr = 0.0, errt, peps;

for (int iter=1; 1; iter++) //checks interpolation accuracy in the centers and adds points
{
    ind = 0; //index of a central interpolating polynom

    int node = 0; //points index for new grid

    //for maximum error estimation:
    peps = merr; //error from previous iteration
    merr = 0.0;
    xmerr = 0.0;
    ymerr = 0.0;

    for (int k=0; k<dimx-1; k++)
    {
        //new grid: copy from old grid
        xnew[node] = xold[k];
        ynew[node] = yold[k];
        node++;

        //interpolation in centers:
        xc = xold[k] + 0.5*(xold[k+1]-xold[k]);

        find_index_for_interp (deg, xc, dimx, xold, &ind);

        //middle polynom:
        eval_neville_polynom (xold+ind, yold+ind, deg, xc, order, order, &ym);

        //left polynom:
        int indl = max(0, ind-1);
        eval_neville_polynom (xold+indl, yold+indl, deg, xc, order, order, &yl);

        //right polynom:
        int indr = min(ind+1, dimx-deg-1);
        eval_neville_polynom (xold+indr, yold+indr, deg, xc, order, order, &yr);

        //check accuracy:
        int flag = 0;

        errl = abs(yl-ym);
        errr = abs(yr-ym);

        //determines maximum error:
        errt = max (errl, errr);
        if (abs(ym) > 1.0) errt = errt/abs(ym);

        if (errt > merr) //maximum error and its location
        {
           merr = errt;
           xmerr = xc;
           ymerr = ym;
        }

        if (errt > *eps) flag = 1; //add the point to the grid

        //check xnew for dimension:
        if (node >= *dim-2)
        {
            fprintf (stdout, "\n\nadaptive_grid_polynom: warning: maximum dimension is reached: node=%d: exit with the previous grid.", node);
            fflush (stdout);

            if (x1 != xold)
            {
                for (int m=0; m<dimx; m++)
                {
                    x1[m] = xold[m];
                    y1[m] = yold[m];
                }
            }

            *dim = dimx;
            *eps = peps;

            delete [] x2;
            delete [] y2;
            return 1;
        }

        if (!flag) continue; //accuracy is good: go to the next point in the old grid

        //adds new xc point to the new grid
        xnew[node] = xc;
        func (xnew+node, ynew+node, p);
        node++;
    }

    //add last grid point (b): k=dimx-1
    xnew[node] = xold[dimx-1];
    ynew[node] = yold[dimx-1];
    node++;

    //switch arrays:
    double *tmp;

    tmp = xold;
    xold = xnew;
    xnew = tmp;

    tmp = yold;
    yold = ynew;
    ynew = tmp;

    if (debug)
    {
        fprintf (stdout, "\n\nadaptive_grid_polynom:\nfor %d iteration grid dimension is %d points.", iter, node);

        fprintf (stdout, "\nmaximum error (%le) is found at (%le) for function value (%le).", merr, xmerr, ymerr);

        fflush (stdout);
    }

    if (node == dimx) break; //no new points have been added during iteration

    dimx = node; //new dimension of the "old" array
}

if (x1 != xold)
{
    for (int m=0; m<dimx; m++)
    {
        x1[m] = xold[m];
        y1[m] = yold[m];
    }
}

*dim = dimx;
*eps = merr;

delete [] x2;
delete [] y2;

return 0;
}

/*****************************************************************************/

int sparse_grid_polynom (int dimx, double *x, double *y, int deg, double *eps_a, double *eps_r, double step, int *dim, double *x1, double *y1, int *ind1)
{
//removes unnecessary points from initial grid (x1, y1) by error estimation for interpolating polynoma of degree deg

int xdim = dimx; //initial fine grid dimension

//additional space:
double *x2 = new double[xdim];
double *y2 = new double[xdim];
int *ind2 = new int[xdim];

double *xold = x1, *yold = y1, *xnew = x2, *ynew = y2;
int *iold = ind1, *inew = ind2;

//initial grid: for "moving" polynoma
dimx = deg + 1;

iold[0] = 0;
xold[0] = x[iold[0]];

iold[dimx-1] = xdim-1;
xold[dimx-1] = x[iold[dimx-1]];

double pi = 3.141592653589793238;

for (int k=1; k<deg; k++)
{
    iold[k] = (int) (0.5*xdim*(1.0-cos((2*k-1)*pi/(2*(deg-1)))));
    xold[k] = x[iold[k]];
}

//initial y grid:
for (int k=0; k<dimx; k++) func_interp (iold[k], y, &yold[k]);

//grid refinement based on comparision of results of interpolating polynoma at centers:
int ind, ic;
double xc, yc, errr, erra;
double ym;

//for maximum error:
int debug = 1; //flag for debugging
double xmerr, ymerr, merr = 0.0, errt, peps;

for (int iter=1; 1; iter++) //checks interpolation accuracy in the centers and adds points
{
    ind = 0; //index of a central interpolating polynom

    int node = 0; //points index for new grid

    //for maximum error estimation:
    peps = merr; //error from previous iteration
    merr = 0.0;
    xmerr = 0.0;
    ymerr = 0.0;

    for (int k=0; k<dimx-1; k++)
    {
        //new grid: copy from old grid
        xnew[node] = xold[k];
        inew[node] = iold[k];
        ynew[node] = yold[k];
        node++;

        //interpolation in index centers:
        ic = iold[k] + (int)(iold[k+1]-iold[k])/2;

        if (ic == iold[k]) continue;

        xc = x[ic];

        find_index_for_interp (deg, xc, dimx, xold, &ind);

        //middle polynom:
        eval_neville_polynom (xold+ind, yold+ind, deg, xc, 0, 0, &ym);

        //check accuracy for each component:
        int flag = 0;

        yc = y[ic];

        erra = abs(yc - ym);
        errr = (*eps_r) * abs(yc);

        //determines maximum error:
        errt = erra;
        if (abs(yc) > 1.0) errt = errt/abs(yc);

        if (errt > merr) //maximum error
        {
           merr = errt;
           xmerr = xc;
           ymerr = yc;
        }

        //checks values of errors:
        if (!(erra < (*eps_a) || erra < errr) || xc-xold[k] > step)
        {
            flag = 1;
        }

        if (!flag) continue; //accuracy is good: go to the next point in the old grid

        //adds new xc point to the new grid
        xnew[node] = xc;
        inew[node] = ic;
        func_interp (ic, y, &ynew[node]);
        node++;
    }

    //add last grid point (b): k = dimx-1
    xnew[node] = xold[dimx-1];
    inew[node] = iold[dimx-1];
    ynew[node] = yold[dimx-1];
    node++;

    //switch arrays:
    void *tmp;

    tmp = (void *)xold;
    xold = xnew;
    xnew = (double *)tmp;

    tmp = (void *)iold;
    iold = inew;
    inew = (int *)tmp;

    tmp = (void *)yold;
    yold = ynew;
    ynew = (double *)tmp;

    if (debug)
    {
        fprintf (stdout, "\n\nsparse_grid_polynom:\nfor %d iteration grid dimension is %d points.", iter, dimx);

        fprintf (stdout, "\nmaximum error (%le) is found at (%le) for function value (%le).", merr, xmerr, ymerr);

        fflush (stdout);
    }

    if (node == dimx) break; //no new points have been added during iteration

    dimx = node; //new dimension of the new "old" array
}

if (x1 != xold)
{
    for (int m=0; m<dimx; m++)
    {
        x1[m] = xold[m];
        y1[m] = yold[m];
        ind1[m] = iold[m];
    }
}

*dim = dimx;
*eps_a = merr;
*eps_r = merr;

if (dimx == xdim)
{
    fprintf (stdout, "\n\nsparse_grid_polynom: warning: no points have been removed from the grid: output dim = %d.", dimx);
}

delete [] x2;
delete [] ind2;
delete [] y2;

return 0;
}

/*****************************************************************************/

void inline func_interp (int ind, double *y, double *yout)
{
*yout = y[ind];
}

/*****************************************************************************/
