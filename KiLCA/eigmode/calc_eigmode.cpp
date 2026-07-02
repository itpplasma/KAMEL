/*! \file calc_eigmode.cpp
    \brief The implementation of the functions declared in calc_eigmode.h
*/

#include "eigmode_sett.h"
#include "calc_eigmode.h"
#include "inout.h"
#include "mode.h"

#include "fortnum.h"

/**********************************************************************************/

// fortnum_vector_fn residual: x = (Re f, Im f), out = (Re det, Im det).
// Replaces the gsl_vector-based eval_det; the multiroot solver builds the
// Jacobian internally by central differences (the former eval_jac).
void eval_det (int n, const double *x, double *f, void *params)
{
(void) n;

int ind = ((det_params *) params)->ind;
int m = ((det_params *) params)->m;
int nn = ((det_params *) params)->n;

core_data *cd = ((det_params *) params)->cd;

const double fre = x[0];
const double fim = x[1];

complex<double> olab = 2.0*pi*(fre + I*fim);

cd->mda[ind] = new mode_data (m, nn, olab, (const settings *)cd->sd, (const background *)cd->bp);

cd->mda[ind]->calc_all_mode_data ();

complex<double> det = (cd->mda[ind])->wd->det;

//for debugging:
FILE *out;
if (!(out = fopen ("det.dat", "a")))
{
    fprintf (stderr, "\nFailed to open file %s\a\n", "det.dat");
}

fprintf (out, "\n%.20le %.20le\t%.20le %.20le", fre, fim, real(det), imag(det));

fclose (out);
//end debugging

f[0] = real(det);
f[1] = imag(det);

//clean up:
delete cd->mda[ind];
cd->mda[ind] = NULL;
clear_all_data_in_mode_data_module_ (); //clean up fortran module data
}

/**********************************************************************************/

int find_det_zeros (int ind, int m, int n, core_data *cd)
{
const eigmode_sett *es = cd->sd->es;

//output file:
char *full_name = new char[1024];
sprintf (full_name, "%s%s", cd->sd->path2project, es->fname);

FILE *out;
if (!(out = fopen (full_name, "w")))
{
    fprintf (stderr, "\nFailed to open file %s\a\n", full_name);
}

int status = 0;

const int N = 2;

det_params p = {ind, m, n, es->delta, cd};

double x_init[2];
double x_root[2];

complex<double> center, quad1, quad2;
double radius;

fprintf (out, "%%Re(f)\t\t\tIm(f)\t\t\tRe(det)\t\t\tIm(det)");
fclose (out);

int k;
for (k=es->kmin; k<=es->kmax; k++)
{
    if (!(out = fopen (full_name, "a")))
    {
        fprintf (stderr, "\nFailed to open file %s\a\n", full_name);
    }

    x_init[0] = real (es->fstart[k]);
    x_init[1] = imag (es->fstart[k]);

    // fortnum_multiroot_hybrid does a hybrid Newton solve building the Jacobian
    // by central differences, matching the former gsl_multiroot_fdfsolver
    // (hybridsj) usage with the finite-difference eval_jac. The residual and
    // step-size tolerances reuse the eigmode-options eps values.
    status = fortnum_multiroot_hybrid (&eval_det, N, x_init, es->eps_abs,
                                       es->eps_res, (int) 1e6, x_root, &p);

    double det_out[2];
    eval_det (N, x_root, det_out, &p);

    fprintf (out, "\n%.20le  %.20le\t%.20le  %.20le", x_root[0], x_root[1], det_out[0], det_out[1]);

    fprintf (out, "\n%%status = %d", status);

    //check if it is really root:
    if (es->test_roots == 1)
    {
        radius = 1.0e1;

        center = x_root[0] + I*x_root[1];

        quad1 = calc_circle_integral (center, radius, inv_det, &p);

        center += 5.0*radius*(1.0+I);

        quad2 = calc_circle_integral (center, radius, inv_det, &p);

        fprintf (stdout, "\nquad1 = %.20le %.20le\tquad2 = %.20le %.20le",
        real(quad1), imag(quad1), real(quad2), imag(quad2));

        fprintf (out, "\n%%quad1 = %.20le %.20le\tquad2 = %.20le %.20le\terr = %.20le",
        real(quad1), imag(quad1), real(quad2), imag(quad2), abs(quad1/quad2));
    }
    fflush (out);
    fclose (out);
}

delete [] full_name;

return status;
}

/**********************************************************************************/

complex<double> calc_circle_integral(complex<double> center, double radius, cmplx_func F, void *params)
{
int k, N = 100;
double dphi = 2.0*pi/N;
complex<double> expon, res;

res = O;
for (k=0; k<N; k++)
{
    expon = exp(k*dphi*I);
    res += F (center + radius*expon, params)*expon;
}
return dphi*I*radius*res;
}

/**********************************************************************************/

complex<double> inv_det(complex<double> z, void *params)
{
double x[2];
double f[2];

x[0] = real(z);
x[1] = imag(z);

eval_det (2, x, f, params);

complex<double> det = f[0] + I*f[1];

return E/det;
}

/**********************************************************************************/

complex<double> test_func(complex<double> z, void *params)
{
return exp(z)/(z-I)/(z-I)/(z-I);
//return (E+I)/z;
}

/**********************************************************************************/

int loop_over_frequences (int ind, int m, int n, core_data *cd)
{
const eigmode_sett *es = cd->sd->es;

//output file:
char *full_name = new char[1024];
sprintf (full_name, "%s%s", cd->sd->path2project, es->fname);

FILE *out;
if (!(out = fopen (full_name, "w")))
{
    fprintf (stderr, "\nFailed to open file %s\a\n", full_name);
}

delete [] full_name;

fprintf (out, "%%iter\tRe(f)\t\t\tIm(f)\t\t\tRe(det)\t\t\tIm(det)");

for (int i=0; i<es->rdim; i++)
{
    double fre = es->rfmin + i*(es->rfmax - es->rfmin)/max(es->rdim-1, 1);

    for (int k=0; k<es->idim; k++)
    {
        double fim = es->ifmin + k*(es->ifmax - es->ifmin)/max(es->idim-1, 1);

        complex<double> olab = 2.0*pi*(fre + I*fim);

        cd->mda[ind] = new mode_data (m, n, olab, (const settings *)cd->sd, (const background *)cd->bp);

        cd->mda[ind]->calc_all_mode_data ();

        fprintf (out, "\n%6u\t%.20le  %.20le\t%.20le  %.20le", i*(es->idim)+k, fre, fim,
                 real(cd->mda[ind]->wd->det), imag(cd->mda[ind]->wd->det));
        fflush (out);

        delete cd->mda[ind];
        cd->mda[ind] = NULL;
        clear_all_data_in_mode_data_module_ ();
    }
}
fclose (out);
return 0;
}

/**********************************************************************************/

void calc_det (double *freq, double *absdet, void *params)
{
double x[2];
double f[2];

x[0] = 0.0e0;
x[1] = *freq;

eval_det (2, x, f, params);

*absdet = f[1]; //imaginary part
}

/**********************************************************************************/
