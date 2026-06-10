#include <stdlib.h>
#include <stdio.h>
#include <complex>
#include <cmath>
#include <cstring>
#include <climits>

#include "slatec.h"
#include "spline.h"
#include "constants.h"

#include "adaptive_quad.h"

/*-----------------------------------------------------------------*/

struct quad_func_params
{
    int k;
    int part;
    int sid;
    double k0;
    double x0;
};

/*-----------------------------------------------------------------*/

extern "C"
{
int calc_direct_fourier_ (int dimx, const double *x, const double *y, int M, double x0, double delta, double *f);

int calc_inverse_fourier_ (int dimx, const double *x, const double *f, int M, double x0, double delta, double *y);
}

/*-----------------------------------------------------------------*/

double func (double x, void *params)
{
quad_func_params *P = (quad_func_params *)params;

double y[2];
spline_eval_d_ (P->sid, 1, &x, 0, 0, 0, 1, (double *)y);

complex<double> arg = -(x - P->x0)*I*(P->k0)*double(P->k);
complex<double> Y(y[0], y[1]);

if (P->part == 0) return real(exp(arg)*Y);
else              return imag(exp(arg)*Y);
}

/*-----------------------------------------------------------------*/

int calc_direct_fourier_ (int dimx, const double *x, const double *y, int M, double x0, double delta, double *f)
{
int ierr;
int sid;
int N = 5; //spline order

double *C = new double[(N+1)*(dimx)*(2)];
double *yt = new double[(dimx)*(2)];

int k;

for (k=0; k<dimx; k++)
{
    yt[k]      = y[2*k];
    yt[dimx+k] = y[2*k+1];
}

spline_alloc_ (N, 1, dimx, x, C, &sid);

spline_calc_ (sid, yt, 0, 1, NULL, &ierr);

struct quad_func_params P = {0, 0, sid, pi/delta, x0};

double re, im, err_re, err_im;

for (k=-M; k<=M; k++)
{
    P.k = k;

    P.part = 0;
    quad_qag (&func, &P, x0-delta, x0+delta, 1.0e-8, 1.0e-8, 1000, 31, &re, &err_re);

    P.part = 1;
    quad_qag (&func, &P, x0-delta, x0+delta, 1.0e-8, 1.0e-8, 1000, 31, &im, &err_im);

    //fprintf (stdout, "\nk = %d:\tresult = %lg%+lgi\terror = %lg%+lgi", k, re, im, err_re, err_im);

    f[2*(M+k)+0] = re/2.0/delta;
    f[2*(M+k)+1] = im/2.0/delta;
}

delete [] C;
delete [] yt;

spline_free_ (sid);

return 0;
}

/*-----------------------------------------------------------------*/

int calc_inverse_fourier_ (int dimx, const double *x, const double *f, int M, double x0, double delta, double *y)
{
int k, i;

complex<double> arg;
complex<double> Y, tmp;

double k0 = pi/delta;

for (i=0; i<dimx; i++)
{
    tmp = O;
    for (k=-M; k<=M; k++)
    {
        arg = (x[i]-x0)*I*k0*double(k);
        Y = f[2*(M+k)+0] + f[2*(M+k)+1]*I;
        tmp += Y*exp(arg);
    }
    y[2*i+0] = real(tmp);
    y[2*i+1] = imag(tmp);
}
return 0;
}

/*-----------------------------------------------------------------*/
