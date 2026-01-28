/*! \file
    \brief The implementation of functions declared in adaptive_grid.h.
*/

#include <cmath>
#include <cfloat>

#include <gsl/gsl_heapsort.h>

#include "adaptive_grid.h"
#include "shared.h"

/******************************************************************************/

void set_interval (struct interval5 *I, double x1, double x2, double y1, double y2,
                   void (*f)(double *, double *, void *p), void *p)
{
I->xs = x1;
I->ys = y1;

I->xe = x2;
I->ye = y2;

I->xm = 0.5*(x1+x2);
f (&I->xm, &I->ym, p);

I->xl = (3*x1+x2)/4.0;
f (&I->xl, &I->yl, p);

I->xr = (x1+3*x2)/4.0;
f (&I->xr, &I->yr, p);
}

/******************************************************************************/

void add_new_interval_to_the_array (interval5 *Iarr, int max_ind, int Icount,
                                    void (*f)(double *, double *, void *p), void *p)
{
//adds new interval at the end of the array (which is right part of the max_ind interval)
Iarr[Icount].xs = Iarr[max_ind].xm;
Iarr[Icount].ys = Iarr[max_ind].ym;

Iarr[Icount].xe = Iarr[max_ind].xe;
Iarr[Icount].ye = Iarr[max_ind].ye;

Iarr[Icount].xm = Iarr[max_ind].xr;
Iarr[Icount].ym = Iarr[max_ind].yr;

Iarr[Icount].xl = (3*Iarr[Icount].xs+Iarr[Icount].xe)/4.0;
f (&Iarr[Icount].xl, &Iarr[Icount].yl, p);

Iarr[Icount].xr = (Iarr[Icount].xs+3*Iarr[Icount].xe)/4.0;
f (&Iarr[Icount].xr, &Iarr[Icount].yr, p);
}

/******************************************************************************/

void update_max_err_interval (interval5 *Iarr, int max_ind,
                              void (*f)(double *, double *, void *p), void *p)
{
//modifies max_ind interval (set it to the left part of the max_ind interval)
Iarr[max_ind].xe = Iarr[max_ind].xm;
Iarr[max_ind].ye = Iarr[max_ind].ym;

Iarr[max_ind].xm = Iarr[max_ind].xl;
Iarr[max_ind].ym = Iarr[max_ind].yl;

Iarr[max_ind].xl = (3*Iarr[max_ind].xs+Iarr[max_ind].xe)/4.0;
f (&Iarr[max_ind].xl, &Iarr[max_ind].yl, p);

Iarr[max_ind].xr = (Iarr[max_ind].xs+3*Iarr[max_ind].xe)/4.0;
f (&Iarr[max_ind].xr, &Iarr[max_ind].yr, p);
}

/******************************************************************************/

inline void eval_error (interval5 *I, double *err)
{
//computes absolute and relative error of the Simpson rule:
//double quad = (I->xe - I->xs)*fabs(I->ys + 4.0*I->ym + I->ye)/6.0;
//double rele = (I->xe - I->xs)*fabs(I->ys - 4.0*I->yl + 6.0*I->ym - 4.0*I->yr + I->ye);
//*err = fmin(rele/quad, rele); bad if quad == 0
//*err = fmin(rele/(quad+DBL_EPSILON), rele); also not very good

//computes absolute error of the Simpson rule: works good!
*err = fabs ((I->xe - I->xs)*fabs(I->ys - 4.0*I->yl + 6.0*I->ym - 4.0*I->yr + I->ye));
}

/******************************************************************************/

void calc_adaptive_1D_grid_ (void (*f)(double *, double *, void *p), void *p,
                             const int *max_dimx, double *eps, int *dimx, double *x, double *y)
{
/*
finds adaptive 1D grid for the function func on the interval (x1, x2)
maximal allowed grid dimension max_dim, error for stopping criterion eps.
on input: initial grid (x, y) has dimension dimx.
on output: generated grid stored in (x, y) and dimension is set in dimx.
*/

int max_dimI = (*max_dimx-1)/4; //maximal allowed number of intervals, each has 4 grid points (+1)

struct interval5 *Iarr = new struct interval5[max_dimI]; //array of intervals

double *err = new double[max_dimI]; //array of errors

int i;
for (i=0; i<*dimx-1; i++)
{
    set_interval (Iarr+i, x[i], x[i+1], y[i], y[i+1], f, p); //initial intervals
    eval_error (Iarr+i, err+i);
}

int dimI = *dimx-1; //number of elements in the array of intervals

double max_err = 0.0;
int max_ind;

while (dimI < max_dimI)
{
    //search for an interval with maximal error:
    max_err = err[0];
    max_ind = 0;

    for (i=1; i<dimI; i++)
    {
        if (err[i] > max_err)
        {
            max_err = err[i];
            max_ind = i;
        }
    }

    if (max_err < *eps) break;

    //split max_ind interval: add a new one to the end and modify the max_ind interval
    add_new_interval_to_the_array (Iarr, max_ind, dimI, f, p);
    eval_error (Iarr+dimI, err+dimI);
    dimI++;

    update_max_err_interval (Iarr, max_ind, f, p);
    eval_error (Iarr+max_ind, err+max_ind);
}

*eps = max_err; //the error reached dirung grid refinement

for(i=0; i<dimI; i++) err[i] = Iarr[i].xs; //err is not needed any more

/*sorting intervals: qsort is faster (uses pointers) but we need indices*/
size_t *perm = new size_t[dimI];

gsl_heapsort_index (perm, err, dimI, sizeof(double), compare_doubles);

//setting up output grid values:
for(i=0; i<dimI; i++)
{
    x[4*i] = Iarr[perm[i]].xs;
    y[4*i] = Iarr[perm[i]].ys;

    x[4*i+1] = Iarr[perm[i]].xl;
    y[4*i+1] = Iarr[perm[i]].yl;

    x[4*i+2] = Iarr[perm[i]].xm;
    y[4*i+2] = Iarr[perm[i]].ym;

    x[4*i+3] = Iarr[perm[i]].xr;
    y[4*i+3] = Iarr[perm[i]].yr;
}

//adds the last grid point:
x[4*dimI] = Iarr[perm[dimI-1]].xe;
y[4*dimI] = Iarr[perm[dimI-1]].ye;

*dimx = 4*dimI+1;

//check if there are condensations of points related to numerical noise and remove them from the grid:

//minimal allowed space between grid points dx: *% from average:
double dx = fmax((1.0e-6)*(x[*dimx-1]-x[0])/(*dimx), 10*DBL_EPSILON);

check_and_remove_grid_condensations_ (dx, dimx, x, y);

/*check output grid if needed:*/
/*
FILE *fout;

if (!(fout=fopen ("xgrid_clean", "w")))
{
    printf ("\nFailed to open file\a\a\a\n");
}

for (i=0; i<*dimx; i++)
{
    fprintf (fout, "%.16le\t%.16le\n", x[i], y[i]);
}
fclose (fout);
*/

delete [] Iarr;
delete [] err;
delete [] perm;
}

/******************************************************************************/

void calc_adaptive_1D_grid_4vector_ (void (*f)(double *, double *, void *p), void *p,
                                     const int *max_dimx, double *eps,
                                     int *dimx, const double *x, const double *y)
{
    /*
    finds adaptive 1D grid for the function func on the interval (x1, x2)
    maximal allowed grid dimension max_dim, error for stopping criterion eps.
    on input: initial grid (x, y) has dimension dimx.
    on output: grid dimension is set to dimx.
    during function evaluation a vector is computed and stored somewhere p points
    grid sorting should be made outside
    */

    int max_dimI = (*max_dimx-1)/4; //maximal allowed number of intervals, each has 4 grid points (+1)

    struct interval5 *Iarr = new struct interval5[max_dimI]; //array of intervals

    double *err = new double[max_dimI]; //array of errors

    int i;
    for (i=0; i<*dimx-1; i++)
    {
        set_interval (Iarr+i, x[i], x[i+1], y[i], y[i+1], f, p); //initial intervals
        eval_error (Iarr+i, err+i);
    }

    int dimI = *dimx-1; //number of elements in the array of intervals

    double max_err = 0.0;
    int max_ind;

    while (dimI < max_dimI)
    {
        //search for an interval with maximal error:
        max_err = err[0];
        max_ind = 0;

        for (i=1; i<dimI; i++)
        {
            if (err[i] > max_err)
            {
                max_err = err[i];
                max_ind = i;
            }
        }

        if (max_err < *eps) break;

        //split max_ind interval: add a new one to the end and modify the max_ind interval
        add_new_interval_to_the_array (Iarr, max_ind, dimI, f, p);
        eval_error (Iarr+dimI, err+dimI);
        dimI++;

        update_max_err_interval (Iarr, max_ind, f, p);
        eval_error (Iarr+max_ind, err+max_ind);
    }

    *eps = max_err; //the error reached dirung grid refinement

    *dimx = 4*dimI+1; //number of evaluated grid points

    delete [] Iarr;
    delete [] err;
}

/******************************************************************************/

void check_and_remove_grid_condensations_ (double dx, int *dimx, double *x, double *y)
{
//checks if there are condensations of points related with narrow peaks and remove them from the sorted grid
//(designed against the grid condensation in regions with a numerical noise in a data)

int i, j, step;

int k = 0; //cleaned array index

for (i=0; i<*dimx; i+=step)
{
    if (k-i) //if there were condensed points
    {
        x[k] = x[i];
        y[k] = y[i];
    }

    step = 0;
    for (j=i+1; j<*dimx; j++)
    {
        if (x[j]-x[i] > dx)
        {
            step = j-i;
            if (step > 1)
            {
                fprintf (stderr, "\nwarning: calc_adaptive_1D_grid: grid condensation detected:\n %d points are found in interval [%.16le, %.16le].", step+1, x[i], x[j]);
            }
            break;
        }
    }

    k++;

    if (step == 0)
    {
        if (i < (*dimx)-1) //for the last array point step = 0 is Ok.
        {
                fprintf (stderr, "\nwarning: calc_adaptive_1D_grid: grid condensation detected:\n %d points are found in interval [%.16le, %.16le].", j-i, x[i], x[j-1]);
        }
        break; //there are no more points within x[i]+dx interval
    }
}

//new grid dimension:
*dimx = k;
}

/******************************************************************************/
