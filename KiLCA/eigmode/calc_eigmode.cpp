/*! \file calc_eigmode.cpp
    \brief The implementation of the functions declared in calc_eigmode.h
*/

#include <math.h>

#include "eigmode_sett.h"
#include "calc_eigmode.h"
#include "inout.h"
#include "mode.h"

/**********************************************************************************/

int eval_det (const double x[2], void *params, double f[2])
{
int ind = ((det_params *) params)->ind;
int m = ((det_params *) params)->m;
int n = ((det_params *) params)->n;

core_data *cd = ((det_params *) params)->cd;

const double fre = x[0];
const double fim = x[1];

complex<double> olab = 2.0*pi*(fre + I*fim);

cd->mda[ind] = new mode_data (m, n, olab, (const settings *)cd->sd, (const background *)cd->bp);

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

return 0;
}

/**********************************************************************************/

int eval_jac (const double x[2], void *params, double J[4])
{
const double fre = x[0];
const double fim = x[1];

double freq[2];
double det_rm[2], det_rp[2], det_im[2], det_ip[2];

double delta = ((det_params *) params)->delta;

//set f-delta:
freq[0] = fre - delta;
freq[1] = fim;
eval_det (freq, params, det_rm);

//set f+delta:
freq[0] = fre + delta;
freq[1] = fim;
eval_det (freq, params, det_rp);

//set f-i*delta:
freq[0] = fre;
freq[1] = fim - delta;
eval_det (freq, params, det_im);

//set f+i*delta:
freq[0] = fre;
freq[1] = fim + delta;
eval_det (freq, params, det_ip);

J[0] = (det_rp[0] - det_rm[0])/2.0/delta;
J[1] = (det_ip[0] - det_im[0])/2.0/delta;
J[2] = (det_rp[1] - det_rm[1])/2.0/delta;
J[3] = (det_ip[1] - det_im[1])/2.0/delta;

return 0;
}

/**********************************************************************************/

int eval_det_jac (const double x[2], void *params, double f[2], double J[4])
{
eval_det (x, params, f);
eval_jac (x, params, J);

return 0;
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
size_t iter = 0;

det_params p = {ind, m, n, es->delta, cd};

double x[2], fv[2], J[4], dx[2];

complex<double> center, quad1, quad2;
double radius;

fprintf (out, "%%iter\tRe(f)\t\t\tIm(f)\t\t\tRe(det)\t\t\tIm(det)");
fclose (out);

int k;
for (k=es->kmin; k<=es->kmax; k++)
{
    if (!(out = fopen (full_name, "a")))
    {
        fprintf (stderr, "\nFailed to open file %s\a\n", full_name);
    }

    iter = 0;
    x[0] = real (es->fstart[k]);
    x[1] = imag (es->fstart[k]);

    eval_det_jac (x, &p, fv, J);

    if (DEBUG_FLAG)
    {
        print_search_state (iter, x, fv);
    }

    fprintf (out, "\n%5zu\t%.20le  %.20le\t%.20le  %.20le", iter, x[0], x[1], fv[0], fv[1]);

    const char *status_text = "did not converge";

    //Newton iteration on det(f) = 0 with the finite-difference Jacobian
    do
    {
        iter++;

        double detJ = J[0]*J[3] - J[1]*J[2];

        if (detJ == 0.0)
        {
            status = 1;
            status_text = "singular Jacobian, iteration is stuck";
            break;
        }

        dx[0] = -( J[3]*fv[0] - J[1]*fv[1])/detJ;
        dx[1] = -(-J[2]*fv[0] + J[0]*fv[1])/detJ;

        x[0] += dx[0];
        x[1] += dx[1];

        eval_det_jac (x, &p, fv, J);

        if (DEBUG_FLAG)
        {
            print_search_state (iter, x, fv);
        }

        fprintf (out, "\n%5zu\t%.20le  %.20le\t%.20le  %.20le", iter, x[0], x[1], fv[0], fv[1]);
        fflush (out);

        //stopping criteria:
        if (es->stop_flag == 0)
        {
            //residual test as in gsl_multiroot_test_residual
            status = (fabs (fv[0]) + fabs (fv[1]) < es->eps_res) ? 0 : -1;
        }
        else if (es->stop_flag == 1)
        {
            //step-size test as in gsl_multiroot_test_delta
            status = (fabs (dx[0]) < es->eps_abs + es->eps_rel*fabs (x[0])
                      && fabs (dx[1]) < es->eps_abs + es->eps_rel*fabs (x[1])) ? 0 : -1;
        }
        else
        {}

        if (status == 0) status_text = "converged";
    }
    while (status != 0 && iter < 1e6);

    if (DEBUG_FLAG)
    {
        printf ("\nstatus = %s", status_text);
    }
    fprintf (out, "\n%%status = %s", status_text);

    //check if it is really root:
    if (es->test_roots == 1)
    {
        radius = 1.0e1;

        center = x[0] + I*x[1];

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

void print_search_state (size_t iter, const double x[2], const double f[2])
{
printf ("\niter = %4zu f = (%le, %le) det = (%le, %le)\n", iter, x[0], x[1], f[0], f[1]);
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
double x[2], f[2];

x[0] = real(z);
x[1] = imag(z);

eval_det (x, params, f);

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
double x[2], f[2];

x[0] = 0.0e0;
x[1] = *freq;

eval_det (x, params, f);

*absdet = f[1]; //imaginary part
}

/**********************************************************************************/
