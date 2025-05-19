/*! \file
    \brief The implementation of functions declared in adaptive_grid_pol.h.
*/

#include <stdlib.h>
#include <stdio.h>
#include <cmath>
#include <complex>
#include <cstring>
#include <climits>

#include "adaptive_grid_pol.h"

/*****************************************************************************/

int adaptive_grid_polynom (void (*func)(double *, double *, void *p), void *p, 
                           double a, double b, int dimy, int deg, int *xdim, double *eps,
                           double *x1, double *y1)
{
//makes an adaptive x grid for moving interpolating polynoma of degree deg
//better to use odd degrees

//additional space:
double *x2 = new double[(*xdim)];
double *y2 = new double[(*xdim)*dimy];

double *xold = x1, *yold = y1, *xnew = x2, *ynew = y2;

//initial grid: for 3 "moving" polynoma
int dimx = deg + 3;

xold[0] = a;
xold[dimx-1] = b;

double pi = 3.141592653589793238;

for (int k=1; k<=deg+1; k++) //Chebishev nodes
{
    xold[k] = 0.5*(a+b) - 0.5*(b-a)*cos((2*k-1)*pi/(2*(deg+1)));
}

//initial y grid:
for (int k=0; k<dimx; k++)
{
    func (xold+k, yold+k*dimy, p);
}

//grid refinement based on comparision of results of 3 interpolating polynoma at centers:
int ind;
double xc;
double err, errl, errr;

//double yl[dimy], ym[dimy], yr[dimy];
double *yl = new double[dimy];
double *ym = new double[dimy];
double *yr = new double[dimy];

//for maximum error:
int jmerr, debug = 0; //flag for debugging
double xmerr, ymerr, merr = 0.0, errt, peps;

for (int iter=1; 1; iter++) //checks interpolation accuracy in the centers and adds points
{
    ind = 0; //index of a central interpolating polynom

    int node = 0; //points index for new grid

    //for maximum error estimation:
    peps = merr; //error from previous iteration
    merr = 0.0;
    jmerr = 0;
    xmerr = 0.0;
    ymerr = 0.0;

    for (int k=0; k<dimx-1; k++)
    {
        //new grid: copy from old grid
        xnew[node] = xold[k];
        for (int j=0; j<dimy; j++) ynew[node*dimy+j] = yold[k*dimy+j];
        node++;

        //interpolation in centers:
        xc = xold[k] + 0.5*(xold[k+1]-xold[k]);

        find_index_for_interp (deg, xc, dimx, xold, &ind);

        //middle polynom:
        eval_interp_polynom (deg, xold+ind, yold+ind*dimy, dimy, xc, ym);

        //left polynom:
        int indl = max(0, ind-1);
        eval_interp_polynom (deg, xold+indl, yold+indl*dimy, dimy, xc, yl);

        //right polynom:
        int indr = min(ind+1, dimx-deg-1);
        eval_interp_polynom (deg, xold+indr, yold+indr*dimy, dimy, xc, yr);

        //check accuracy for each component:
        int flag = 0;

        for (int j=0; j<dimy; j++)
        {
            errl = abs(yl[j]-ym[j]);
            errr = abs(yr[j]-ym[j]);
            err  = (*eps) * abs(ym[j]);

            //determines maximum error:
            errt = max (errl, errr);

            //if (abs(ym[j]) > 1.0e-16) errt = min (errt, errt/abs(ym[j]));
            if (abs(ym[j]) > 1.0) errt /= abs(ym[j]);

            if (errt > merr) //maximum error and its location
            {
               merr = errt;
               xmerr = xc;
               ymerr = ym[j];
               jmerr = j;
            }

            //checks values of errors:
            //if (!((errl < (*eps) || errl < err) && (errr < (*eps) || errr < err)))
            //{
            //    flag = 1;
            //}
            if (errt > *eps) flag = 1; //add the point to the grid
        }

        //check xnew for dimension:
        if (node >= *xdim-2)
        {
#if DEBUG_FLAG
          fprintf(stdout,
                  "\n\nadaptive_grid_polynom: warning: maximum dimension is "
                  "reached: node=%d: exit with the previous grid.",
                  node);
          fflush(stdout);
#endif

          if (x1 != xold) {
            for (int m = 0; m < dimx; m++) {
              x1[m] = xold[m];
              for (int j = 0; j < dimy; j++)
                y1[m * dimy + j] = yold[m * dimy + j];
            }
          }

          *xdim = dimx;
          *eps = peps;

          delete[] x2;
          delete[] y2;

          delete[] ym;
          delete[] yr;
          delete[] yl;

          return 1;
        }

        if (!flag) continue; //accuracy is good: go to the next point in the old grid

        //adds new xc point to the new grid
        xnew[node] = xc;
        func (xnew+node, ynew+node*dimy, p);
        node++;
    }

    //add last grid point (b): k=dimx-1
    xnew[node] = xold[dimx-1];
    for (int j=0; j<dimy; j++) ynew[node*dimy+j] = yold[(dimx-1)*dimy+j];
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

        fprintf (stdout, "\nmaximum error (%le) is found at (%le) for %d-th component (%le).", merr, xmerr, jmerr, ymerr);

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
        for (int j=0; j<dimy; j++) y1[m*dimy+j] = yold[m*dimy+j];
    }
}

*xdim = dimx;
*eps = merr;

delete [] x2;
delete [] y2;

delete [] ym;
delete [] yr;
delete [] yl;

return 0;
}

/*****************************************************************************/

int adaptive_grid_polynom (void (*func)(double *, double *, void *p), void *p,
                           double a, double b, int dimy, int deg, int *xdim, double *eps,
                           int dim_err, int *ind_err, double *x1, double *y1)
{
//makes an adaptive x grid for moving interpolating polynoma of degree deg
//better to use odd degrees

//additional space:
double *x2 = new double[(*xdim)];
double *y2 = new double[(*xdim)*dimy];

double *xold = x1, *yold = y1, *xnew = x2, *ynew = y2;

//initial grid: for 3 "moving" polynoma
int dimx = deg + 3;

xold[0] = a;
xold[dimx-1] = b;

double pi = 3.141592653589793238;

for (int k=1; k<=deg+1; k++) //Chebishev nodes
{
    xold[k] = 0.5*(a+b) - 0.5*(b-a)*cos((2*k-1)*pi/(2*(deg+1)));
}

//initial y grid:
for (int k=0; k<dimx; k++)
{
    func (xold+k, yold+k*dimy, p);
}

//grid refinement based on comparision of results of 3 interpolating polynoma at centers:
int ind;
double xc;
double errl, errr;

double *yl = new double[dim_err];
double *ym = new double[dim_err];
double *yr = new double[dim_err];

//for maximum error:
int jmerr, debug = 0;
double xmerr, ymerr, merr = 0.0, errt, peps;

double err;

for (int iter=1; 1; iter++) //checks interpolation accuracy in the centers and adds points
{
    ind = 0; //index of a central interpolating polynom

    int node = 0; //points index for new grid

    //for maximum error estimation:
    peps = merr; //error from previous iteration
    merr = 0.0;
    jmerr = 0;
    xmerr = 0.0;
    ymerr = 0.0;

    for (int k=0; k<dimx-1; k++)
    {
        //new grid: copy from old grid
        xnew[node] = xold[k];
        for (int j=0; j<dimy; j++) ynew[node*dimy+j] = yold[k*dimy+j];
        node++;

        //interpolation in centers:
        xc = xold[k] + 0.5*(xold[k+1]-xold[k]);

        find_index_for_interp (deg, xc, dimx, xold, &ind);

        int indl = max(0, ind-1);
        int indr = min(ind+1, dimx-deg-1);

        int flag = 0;

        for (int l=0; l<dim_err; l++)
        {
            int j = ind_err[l];

            //middle polynom:
            eval_neville_polynom (xold+ind,  yold+ind*dimy+j, dimy, deg, xc, &ym[l]);

            //left polynom:
            eval_neville_polynom (xold+indl, yold+indl*dimy+j, dimy, deg, xc, &yl[l]);

            //right polynom:
            eval_neville_polynom (xold+indr, yold+indr*dimy+j, dimy, deg, xc, &yr[l]);

            errl = abs(yl[l]-ym[l]);
            errr = abs(yr[l]-ym[l]);

            //determines maximum error:
            errt = max(errl, errr);

            if (abs(ym[l]) > 1.0) errt /= abs(ym[l]);

            if (errt > merr) //maximum error and its location
            {
               merr = errt;
               xmerr = xc;
               ymerr = ym[l];
               jmerr = j;
            }

            err = *eps;

            if (errt > err) flag = 1; //flag to add the point
        }

        //check xnew for dimension:
        if (node >= *xdim-2)
        {
#if DEBUG_FLAG
          fprintf(stdout,
                  "\n\nadaptive_grid_polynom: warning: maximum dimension is "
                  "reached: node=%d: exit with the previous grid.",
                  node);
          fflush(stdout);
#endif

          if (x1 != xold) {
            for (int m = 0; m < dimx; m++) {
              x1[m] = xold[m];
              for (int j = 0; j < dimy; j++)
                y1[m * dimy + j] = yold[m * dimy + j];
            }
          }

          *xdim = dimx;
          *eps = peps;

          delete[] x2;
          delete[] y2;

          delete[] ym;
          delete[] yr;
          delete[] yl;

          return 1;
        }

        if (!flag) continue; //accuracy is good: go to the next point in the old grid

        //adds new xc point to the new grid
        xnew[node] = xc;
        func (xnew+node, ynew+node*dimy, p);
        node++;
    }

    //add last grid point (b): k=dimx-1
    xnew[node] = xold[dimx-1];
    for (int j=0; j<dimy; j++) ynew[node*dimy+j] = yold[(dimx-1)*dimy+j];
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

        fprintf (stdout, "\nmaximum error (%le) is found at (%le) for %d-th component (%le).", merr, xmerr, jmerr, ymerr);

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
        for (int j=0; j<dimy; j++) y1[m*dimy+j] = yold[m*dimy+j];
    }
}

*xdim = dimx;
*eps = merr;

delete [] x2;
delete [] y2;

delete [] ym;
delete [] yr;
delete [] yl;

return 0;
}

/*****************************************************************************/

int adaptive_grid_polynom_res (void (*func)(double *, double *, void *p), void *p,
                               double a, double b, int dimy, int deg, int *xdim, double *eps,
                               double r_res, double D, double eps_res,
                               int dim_err, int *ind_err, double *x1, double *y1)
{
//makes an adaptive x grid for moving interpolating polynoma of degree deg
//better to use odd degrees

//additional space:
double *x2 = new double[(*xdim)];
double *y2 = new double[(*xdim)*dimy];

double *xold = x1, *yold = y1, *xnew = x2, *ynew = y2;

//initial grid: for 3 "moving" polynoma
int dimx = deg + 3;

xold[0] = a;
xold[dimx-1] = b;

double pi = 3.141592653589793238;

for (int k=1; k<=deg+1; k++) //Chebishev nodes
{
    xold[k] = 0.5*(a+b) - 0.5*(b-a)*cos((2*k-1)*pi/(2*(deg+1)));
}

//initial y grid:
for (int k=0; k<dimx; k++)
{
    func (xold+k, yold+k*dimy, p);
}

//grid refinement based on comparision of results of 3 interpolating polynoma at centers:
int ind;
double xc;
double errl, errr;

double *yl = new double[dim_err];
double *ym = new double[dim_err];
double *yr = new double[dim_err];

//for maximum error:
int jmerr, debug = 0;
double xmerr, ymerr, merr = 0.0, errt, peps;

double err;

for (int iter=1; 1; iter++) //checks interpolation accuracy in the centers and adds points
{
    ind = 0; //index of a central interpolating polynom

    int node = 0; //points index for new grid

    //for maximum error estimation:
    peps = merr; //error from previous iteration
    merr = 0.0;
    jmerr = 0;
    xmerr = 0.0;
    ymerr = 0.0;

    for (int k=0; k<dimx-1; k++)
    {
        //new grid: copy from old grid
        xnew[node] = xold[k];
        for (int j=0; j<dimy; j++) ynew[node*dimy+j] = yold[k*dimy+j];
        node++;

        //interpolation in centers:
        xc = xold[k] + 0.5*(xold[k+1]-xold[k]);

        find_index_for_interp (deg, xc, dimx, xold, &ind);

        int indl = max(0, ind-1);
        int indr = min(ind+1, dimx-deg-1);

        int flag = 0;

        for (int l=0; l<dim_err; l++)
        {
            int j = ind_err[l];

            //middle polynom:
            eval_neville_polynom (xold+ind,  yold+ind*dimy+j,  dimy, deg, xc, &ym[l]);

            //left polynom:
            eval_neville_polynom (xold+indl, yold+indl*dimy+j, dimy, deg, xc, &yl[l]);

            //right polynom:
            eval_neville_polynom (xold+indr, yold+indr*dimy+j, dimy, deg, xc, &yr[l]);

            errl = abs(yl[l]-ym[l]);
            errr = abs(yr[l]-ym[l]);

            //determines maximum error:
            errt = max(errl, errr);

            if (abs(ym[l]) > 1.0) errt /= abs(ym[l]);

            if (errt > merr) //maximum error and its location
            {
               merr = errt;
               xmerr = xc;
               ymerr = ym[l];
               jmerr = j;
            }

            //change the error in resonance zone:
            err = (*eps-eps_res)*(1.0 - exp(-(xc-r_res)*(xc-r_res)/D/D)) + eps_res;

            if (errt > err) flag = 1; //flag to add the point
        }

        //check xnew for dimension:
        if (node >= *xdim-2)
        {
#if DEBUG_FLAG
          fprintf(stdout,
                  "\n\nadaptive_grid_polynom: warning: maximum dimension is "
                  "reached: node=%d: exit with the previous grid.",
                  node);
          fflush(stdout);
#endif

          if (x1 != xold) {
            for (int m = 0; m < dimx; m++) {
              x1[m] = xold[m];
              for (int j = 0; j < dimy; j++)
                y1[m * dimy + j] = yold[m * dimy + j];
            }
          }

          *xdim = dimx;
          *eps = peps;

          delete[] x2;
          delete[] y2;

          delete[] ym;
          delete[] yr;
          delete[] yl;

          return 1;
        }

        if (!flag) continue; //accuracy is good: go to the next point in the old grid

        //adds new xc point to the new grid
        xnew[node] = xc;
        func (xnew+node, ynew+node*dimy, p);
        node++;
    }

    //add last grid point (b): k=dimx-1
    xnew[node] = xold[dimx-1];
    for (int j=0; j<dimy; j++) ynew[node*dimy+j] = yold[(dimx-1)*dimy+j];
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

        fprintf (stdout, "\nmaximum error (%le) is found at (%le) for %d-th component (%le).", merr, xmerr, jmerr, ymerr);

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
        for (int j=0; j<dimy; j++) y1[m*dimy+j] = yold[m*dimy+j];
    }
}

*xdim = dimx;
*eps = merr;

delete [] x2;
delete [] y2;

delete [] ym;
delete [] yr;
delete [] yl;

return 0;
}

/*****************************************************************************/

int sparse_grid_polynom (double *x, int dimx, double *y, int dimy, int deg,
                         double *eps_a, double *eps_r, double step,
                         int *dim, double *x1, double *y1, int *ind1)
{
//removes unnecessary points from initial grid (x1, y1) by error estimation for interpolating polynoma of degree deg

int xdim = dimx; //fine grid dimension

//additional space:
double *x2 = new double[xdim];
double *y2 = new double[xdim*dimy];
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
for (int k=0; k<dimx; k++)
{
    func_interp (iold[k], yold+k*dimy, x, xdim, y, dimy);
}

//grid refinement based on comparision of results of interpolating polynoma at centers:
int ind, ic;
double xc, yc, errr, erra;
//double ym[dimy];
double *ym = new double[dimy];

//for maximum error:
int jmerr, debug = 0; //flag for debugging
double xmerr, ymerr, merr = 0.0, errt, peps;

for (int iter=1; 1; iter++) //checks interpolation accuracy in the centers and adds points
{
    ind = 0; //index of a central interpolating polynom

    int node = 0; //points index for new grid

    //for maximum error estimation:
    peps = merr; //error from previous iteration
    merr = 0.0;
    jmerr = 0;
    xmerr = 0.0;
    ymerr = 0.0;

    for (int k=0; k<dimx-1; k++)
    {
        //new grid: copy from old grid
        xnew[node] = xold[k];
        inew[node] = iold[k];
        for (int j=0; j<dimy; j++) ynew[node*dimy+j] = yold[k*dimy+j];
        node++;

        //interpolation in index centers:
        ic = iold[k] + (int)(iold[k+1]-iold[k])/2;

        if (ic == iold[k]) continue;

        xc = x[ic];

        find_index_for_interp (deg, xc, dimx, xold, &ind);

        //middle polynom:
        eval_interp_polynom (deg, xold+ind, yold+ind*dimy, dimy, xc, ym);

        //check accuracy for each component:
        int flag = 0;

        for (int j=0; j<dimy; j++)
        {
            yc = y[ic*dimy+j];

            erra = abs(yc - ym[j]);
            errr = (*eps_r) * abs(yc);

            //determines maximum error:
            errt = erra;
            if (abs(yc) > 1.0e-16) errt = min (errt, errt/abs(yc));

            if (errt > merr) //maximum error and its location
            {
               merr = errt;
               xmerr = xc;
               ymerr = yc;
               jmerr = j;
            }

            //checks values of errors:
            if (!(erra < (*eps_a) || erra < errr) || xc-xold[k] > step)
            {
                flag = 1;
            }
        }

        if (!flag) continue; //accuracy is good: go to the next point in the old grid

        //adds new xc point to the new grid
        xnew[node] = xc;
        inew[node] = ic;
        func_interp (ic, ynew+node*dimy, x, xdim, y, dimy);
        node++;
    }

    //add last grid point (b): k=dimx-1
    xnew[node] = xold[dimx-1];
    inew[node] = iold[dimx-1];
    for (int j=0; j<dimy; j++) ynew[node*dimy+j] = yold[(dimx-1)*dimy+j];
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

        fprintf (stdout, "\nmaximum error (%le) is found at (%le) for %d-th component (%le).", merr, xmerr, jmerr, ymerr);

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
        ind1[m] = iold[m];
        for (int j=0; j<dimy; j++) y1[m*dimy+j] = yold[m*dimy+j];
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
delete [] ym;

return 0;
}

/*****************************************************************************/

void inline func_interp (int ind, double *yout, double *x, int dimx, double *y, int dimy)
{
for (int j=0; j<dimy; j++) yout[j] = y[ind*dimy+j];
}

/*****************************************************************************/
