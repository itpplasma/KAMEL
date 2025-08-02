/*! \file
    \brief Confluent hypergeometric function 1F1(a,b,z) for a = 1 and complex b & z. Both kummer series and continued fractions are used.
*/

#include <stdlib.h>
#include <stdio.h>
#include <complex>
#include <cmath>

#include <gsl/gsl_integration.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_sum.h>
#include "constants.h"

using namespace std;

/*******************************************************************/

extern "C"
{
int hypergeometric1f1_quad_ (double *b_re, double *b_im, double *z_re, double *z_im, double *f_re, double *f_im);

int hypergeometric1f1_kummer_nmax_ (double *b_re, double *b_im, double *z_re, double *z_im, double *f_re, double *f_im);

int hypergeometric1f1_kummer_ada_ (double *b_re, double *b_im, double *z_re, double *z_im, double *f_re, double *f_im);

int hypergeometric1f1_kummer_modified_0_nmax_ (double *b_re, double *b_im, double *z_re, double *z_im, double *f_re, double *f_im);

int hypergeometric1f1_kummer_modified_0_ada_ (double *b_re, double *b_im, double *z_re, double *z_im, double *f_re, double *f_im);

int hypergeometric1f1_cont_fract_1_modified_0_ada_ (double *b_re, double *b_im, double *z_re, double *z_im, double *f_re, double *f_im);

int hypergeometric1f1_kummer_modified_0_accel_ (double *b_re, double *b_im, double *z_re, double *z_im, double *f_re, double *f_im);

int hypergeometric1f1_kummer_modified_1_ (double *b_re, double *b_im, double *z_re, double *z_im, double *f_re, double *f_im);

int hypergeometric1f1_cont_fract_1_dir_ (double *b_re, double *b_im, double *z_re, double *z_im, double *f_re, double *f_im);

int hypergeometric1f1_cont_fract_1_inv_nmax_ (double *b_re, double *b_im, double *z_re, double *z_im, double *f_re, double *f_im);

int hypergeometric1f1_cont_fract_1_inv_ada_ (double *b_re, double *b_im, double *z_re, double *z_im, double *f_re, double *f_im);
}

/*******************************************************************/

struct quad_params
{
    complex<double> b;
    complex<double> z;
    int part;
};

/*******************************************************************/

double exp1mt (double t, void *params)
{
quad_params *qp = (quad_params *)params;

complex<double> ans = ((qp->b)-1.0) * exp((qp->z)*t) * pow(1.0-t, (qp->b)-2.0);

if (qp->part == 0) return real(ans);
else               return imag(ans);
}

/*******************************************************************/

int hypergeometric1f1_quad_ (double *b_re, double *b_im, double *z_re, double *z_im, double *f_re, double *f_im)
{
//computes function 1F1(a,b,z) for a = 1 and complex b & z by quadrature
//must be optimized: avoid memory allocation!

gsl_set_error_handler_off ();

complex<double> b(*b_re, *b_im), z(*z_re, *z_im);

size_t limit = 100;
double epsabs = 1.0e-12, epsrel = 1.0e-12, err;

gsl_integration_workspace *w = gsl_integration_workspace_alloc (limit);

quad_params qp = {b,z};

gsl_function F;
F.function = &exp1mt;
F.params = &qp;

qp.part = 0;
gsl_integration_qag (&F, 0.0, 1.0, epsabs, epsrel, limit, GSL_INTEG_GAUSS21, w, f_re, &err);

qp.part = 1;
gsl_integration_qag (&F, 0.0, 1.0, epsabs, epsrel, limit, GSL_INTEG_GAUSS21, w, f_im, &err);

gsl_integration_workspace_free (w);

return 0;
}

/*******************************************************************/

int hypergeometric1f1_kummer_nmax_ (double *b_re, double *b_im, double *z_re, double *z_im, double *f_re, double *f_im)
{
//computes function 1F1(a,b,z) for a = 1 and complex b & z by kummer series

complex<double> b(*b_re, *b_im), z(*z_re, *z_im);

long int N = (int) ceil(-20.0/log10(abs(z/b))) + 5;

//fprintf (stdout, "\nhypergeometric1f1_kummer: N = %ld", N);

if (N < 1 || N > 1e6)
{
    fprintf (stdout, "\nwarning: hypergeometric1F1_kummer: N=%ld", N);
    fprintf (stdout, "\nb=%le %le\nz=%le %le\nabs(z/b)=%le", real(b), imag(b), real(z), imag(z),abs(z/b));
}

complex<double> term = z/(b+(double)N);

for (long int n=N-1; n>=0; n--)
{
    term = 1.0 + z/(b+(double)n)*term;
}

*f_re = real(term);
*f_im = imag(term);

return 0;
}

/*******************************************************************/

int hypergeometric1f1_kummer_ada_ (double *b_re, double *b_im, double *z_re, double *z_im, double *f_re, double *f_im)
{
//computes function 1F1(a,b,z) for a = 1 and complex b & z by kummer series

complex<double> b(*b_re, *b_im), z(*z_re, *z_im);

complex<double> term, S1, S2;
double err, eps = DBL_EPSILON;
long int n, Nmax, maxNmax = 1e8;

for (Nmax=4; Nmax<maxNmax; Nmax*=2)
{
    n = Nmax;
    term = z/(b+(double)n);

    for (n=Nmax-1; n>=0; n--)
    {
        term = 1.0 + z/(b+(double)n)*term;
    }

    S1 = term;

    n = Nmax+1;
    term = z/(b+(double)n);

    for (n=Nmax; n>=0; n--)
    {
        term = 1.0 + z/(b+(double)n)*term;
    }

    S2 = term;

    err = min (abs((S2-S1)/S2), abs(S2-S1));
    if (err < eps) break;
}

//fprintf (stdout, "\nhypergeometric1f1_kummer_ada: Nmax = %ld err = %le\nS2 = %.20le %.20le", Nmax, err, real(S2), imag(S2));

if (Nmax >= maxNmax)
{
    fprintf (stderr, "\nwarning: hypergeometric1f1_kummer_ada: Nmax = %ld err = %le\nS2 = %le %le", Nmax, err, real(S2), imag(S2));
    fprintf (stderr, "\nb = %le %le z = %le %le", real(b), imag(b), real(z), imag(z));
}

*f_re = real(S2);
*f_im = imag(S2);

return 0;
}

/*******************************************************************/

int hypergeometric1f1_kummer_modified_1_ (double *b_re, double *b_im, double *z_re, double *z_im, double *f_re, double *f_im)
{
//computes modified function 1F1m(a,b,z) for a = 1 and complex b & z by kummer series
//1F1 = 1 + z/b + z^2/b/(b+1)*1F1m

complex<double> b(*b_re, *b_im), z(*z_re, *z_im);

long int N = (int) ceil(-20.0/log10(abs(z/b))) + 5;

if (N < 1 || N > 1e6)
{
    fprintf (stdout, "\nwarning: hypergeometric1F1_kummer: N=%ld", N);
    fprintf (stdout, "\nb=%le %le\nz=%le %le\nabs(z/b)=%le", real(b), imag(b), real(z), imag(z),abs(z/b));
}

complex<double> term = z/(b+(double)N);

for (long int n=N-1; n>=2; n--)
{
    term = 1.0 + z/(b+(double)n)*term;
}

*f_re = real(term);
*f_im = imag(term);

return 0;
}

/*******************************************************************/

int hypergeometric1f1_kummer_modified_0_nmax_ (double *b_re, double *b_im, double *z_re, double *z_im, double *f_re, double *f_im)
{
//computes modified function 1F1m(a,b,z) for a = 1 and complex b & z by kummer series
//1F1 = 1 + z/b + z^2/b/(b+1)*(1 + 1F1m)

complex<double> b(*b_re, *b_im), z(*z_re, *z_im);

long int N = (long int) ceil(-20.0/log10(abs(z/b))) + 5;

if (N < 1 || N > 1e6)
{
    fprintf (stdout, "\nwarning: hypergeometric1F1_kummer: N=%ld", N);
    fprintf (stdout, "\nb=%le %le\nz=%le %le\nabs(z/b)=%le", real(b), imag(b), real(z), imag(z),abs(z/b));
}

complex<double> term = z/(b+(double)N);

for (long int n=N-1; n>2; n--)
{
    term = 1.0 + z/(b+(double)n)*term;
}

term *= z/(b+2.0);

*f_re = real(term);
*f_im = imag(term);

return 0;
}

/*******************************************************************/

int hypergeometric1f1_kummer_modified_0_ada_ (double *b_re, double *b_im, double *z_re, double *z_im, double *f_re, double *f_im)
{
//computes modified function 1F1m(a,b,z) for a = 1 and complex b & z by kummer series
//1F1 = 1 + z/b + z^2/b/(b+1)*(1 + 1F1m)

complex<double> b(*b_re, *b_im), z(*z_re, *z_im);

complex<double> term, S1, S2;
double err, eps = DBL_EPSILON;
long int n, Nmax, maxNmax = 1e8;

for (Nmax=4; Nmax<maxNmax; Nmax*=2)
{
    n = Nmax;
    term = z/(b+(double)n);

    for (n=Nmax-1; n>2; n--)
    {
        term = 1.0 + (z/(b+(double)n))*term;
    }

    S1 = term*(z/(b+2.0));

    n = Nmax+1;
    term = z/(b+(double)n);

    for (n=Nmax; n>2; n--)
    {
        term = 1.0 + (z/(b+(double)n))*term;
    }

    S2 = term*(z/(b+2.0));

    //err = min (abs((S2-S1)/S2), abs(S2-S1));
    err = abs((S2-S1)/S2);

    if (err < eps) break;
}

if (Nmax >= maxNmax)
{
    fprintf (stderr, "\nwarning: hypergeometric1f1_kummer_modified_0_ada: Nmax = %ld err = %le\nS2 = %le %le", Nmax, err, real(S2), imag(S2));
    fprintf (stderr, "\nb = %le %le z = %le %le", real(b), imag(b), real(z), imag(z));
}

*f_re = real(S2);
*f_im = imag(S2);

return 0;
}

/*******************************************************************/

int hypergeometric1f1_kummer_modified_0_accel_ (double *b_re, double *b_im, double *z_re, double *z_im, double *f_re, double *f_im)
{
//computes modified function 1F1m(a,b,z) for a = 1 and complex b & z by kummer series
//1F1 = 1 + z/b + z^2/b/(b+1)*(1 + 1F1m)
//do not use: gives unstable results with wrong error estimation!!!

complex<double> b(*b_re, *b_im), z(*z_re, *z_im);

int N = min(50, (int) ceil(-20.0/log10(abs(z/b))) + 5);

complex<double> term[N];
double t_re[N], t_im[N];

int n = 0;
term[n] = z/(b+3.0);
t_re[n] = real(term[n]);
t_im[n] = imag(term[n]);

for (n=1; n<N; n++)
{
    term[n] = term[n-1]*(z/(b+(double)(n+3)));
    t_re[n] = real(term[n]);
    t_im[n] = imag(term[n]);
}

gsl_sum_levin_u_workspace *w = gsl_sum_levin_u_alloc (N);

double err;

gsl_sum_levin_u_accel (t_re, N, w, f_re, &err);
if (err > 1.e-16)
{
    fprintf (stdout, "\nerr = %.16le sum_re = %.16le using %zu terms", err, *f_re, w->terms_used);
}

gsl_sum_levin_u_accel (t_im, N, w, f_im, &err);
if (err > 1.e-16)
{
    fprintf (stdout, "\nerr = %.16le sum_im = %.16le using %zu terms", err, *f_im, w->terms_used);
}

gsl_sum_levin_u_free (w);

return 0;
}

/*******************************************************************/

int hypergeometric1f1_cont_fract_1_modified_0_ada_ (double *b_re, double *b_im, double *z_re, double *z_im, double *f_re, double *f_im)
{
//computes modified function 1F1m(a,b,z) for a = 1 and complex b & z by continued fraction
//1F1 = 1 + z/b + z^2/b/(b+1)*(1 + 1F1m)

complex<double> b(*b_re, *b_im), z(*z_re, *z_im), F11m(0.0, 0.0);

if (abs(z/b) < 0.1e0)
{
    hypergeometric1f1_kummer_modified_0_ada_ (b_re, b_im, z_re, z_im, f_re, f_im);
}
else //big numbers substraction - better to implement direct continued fraction!
{
    hypergeometric1f1_cont_fract_1_inv_ada_ (b_re, b_im, z_re, z_im, f_re, f_im);
    F11m = ((*f_re) + (*f_im)*I - 1.0 - z/b)*(b/z)*((b + 1.0)/z) - 1.0;
    *f_re = real(F11m);
    *f_im = imag(F11m);
}

return 0;
}

/*******************************************************************/

int hypergeometric1f1_cont_fract_1_inv_ada_ (double *b_re, double *b_im, double *z_re, double *z_im, double *f_re, double *f_im)
{
//computes function 1F1(a,b,z) for a = 1 and complex b & z by inverse evaluation of continued fractions
complex<double> b(*b_re-1.0, *b_im), z(*z_re, *z_im);

complex<double> term, S1, S2;
double err, eps = DBL_EPSILON;
int n, Nmax, maxNmax = 1e6;

for (Nmax=4; Nmax<maxNmax; Nmax*=2)
{
    n = Nmax;
    term = ((double)n) * z / (b - z + (double)n);
    for (n=Nmax-1; n>0; n--)
    {
        term = ((double)n) * z / (b - z + (double)n + term);
    }

    S1 = b / (b - z + term);

    n = Nmax+1;
    term = ((double)n) * z / (b - z + (double)n);
    for (n=Nmax; n>0; n--)
    {
        term = ((double)n) * z / (b - z + (double)n + term);
    }

    S2 = b / (b - z + term);

    //err = min (abs((S2-S1)/S2), abs(S2-S1));
    err = abs((S2-S1)/S2);

    if (err < eps) break;
}

//fprintf (stdout, "\nhypergeometric1f1_cont_fract_1_inv_ada: Nmax = %d err = %le\nS2 = %.20le %.20le", Nmax, err, real(S2), imag(S2));

if (Nmax >= maxNmax)
{
  fprintf (stderr, "\nwarning: hypergeometric1f1_cont_fract_1_inv_ada: Nmax = %d err = %le\nS2 = %le %le", Nmax, err, real(S2), imag(S2));
    fprintf (stderr, "\nb = %le %le z = %le %le term = %le %le", real(b)+1.0, imag(b), real(z), imag(z), real(term), imag(term));
}

*f_re = real(S2);
*f_im = imag(S2);

return 0;
}

/*******************************************************************/

int hypergeometric1f1_cont_fract_1_inv_nmax_ (double *b_re, double *b_im, double *z_re, double *z_im, double *f_re, double *f_im)
{
//computes function 1F1(a,b,z) for a = 1 and complex b & z by inverse evaluation of continued fractions
const int Nmax = 1e3; //actually must be calculated from error estimation!

complex<double> b(*b_re-1.0, *b_im), z(*z_re, *z_im);

complex<double> term;

int n = Nmax;

term = ((double)n) * z / (b - z + (double)n);

for (n=Nmax-1; n>0; n--)
{
    term = ((double)n) * z / (b - z + (double)n + term);
}

term = b / (b - z + term);

*f_re = real(term);
*f_im = imag(term);

return 0;
}

/*******************************************************************/

int hypergeometric1f1_cont_fract_1_dir_ (double *b_re, double *b_im, double *z_re, double *z_im, double *f_re, double *f_im)
{
//computes function 1F1(a,b,z) for a = 1 and complex b & z by direct evaluation of continued fractions: warning: this method is UNSTABLE sometimes!
const long int Nmax = 1e6;
const double eps = DBL_EPSILON;
double err;

complex<double> b(*b_re-1.0, *b_im), z(*z_re, *z_im);

complex<double> An2, An1, An, Bn2, Bn1, Bn;
complex<double> an, bn;
complex<double> Sn1, Sn;

An2 = 1.0e0; An1 = 0.0e0;
Bn2 = 0.0e0; Bn1 = 1.0e0;
Sn1 = An1/Bn1;

long int n;

for (n=1; n<Nmax; n++)
{
    an = ((double)n) * z;
    bn = b - z + (double)n;

    An = bn*An1 + an*An2;
    Bn = bn*Bn1 + an*Bn2;
    Sn = An/Bn;

    err = min (abs((Sn-Sn1)/Sn), abs(Sn-Sn1));
    if (err < eps) break;

    An2 = An1; An1 = An;
    Bn2 = Bn1; Bn1 = Bn;
    Sn1 = Sn;
}

Sn = b/(b-z+Sn);

if (n == Nmax)
{
    fprintf (stderr, "\nwarning: hypergeometric1f1_cont_fract_1_dir: n = %ld err = %le\nSn = %le %le", n, err, real(Sn), imag(Sn));
    fprintf (stderr, "\nb = %le %le z = %le %le", real(b)+1.0, imag(b), real(z), imag(z));
    fprintf (stderr, "\nBn = %le %le An = %le %le", real(Bn), imag(Bn), real(An), imag(An));
}

*f_re = real(Sn);
*f_im = imag(Sn);

return 0;
}

/*******************************************************************/
