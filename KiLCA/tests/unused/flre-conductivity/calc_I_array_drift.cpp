/*! \file
    \brief Special Imnl function.
*/

#include <complex>
#include <cmath>

#include <gsl/gsl_integration.h>
#include <gsl/gsl_errno.h>

#include "constants.h"

#include <iostream>
#include <iomanip>

/*******************************************************************/

extern "C"
{
int calc_imnl_quad_ (int * m_max, int * n_max, int * l_max,
                     double * x1_re, double * x1_im,
                     double * x2_re, double * x2_im,
                     double * x3_re, double * x3_im,
                     double * x4_re, double * x4_im,
                     double * Imnl);
}

/*******************************************************************/

struct quad_params
{
    complex<double> x1;
    complex<double> x2;
    complex<double> x3;
    complex<double> x4;

    int l;

    complex<double> gamma;

    double x;
    double y;

    int part;
};

/*******************************************************************/

double subint (double tau, void * params)
{
quad_params * p = static_cast<quad_params *>(params);

complex<double> gme = (p->gamma) - E;

complex<double> ExP = std::exp(-tau);

complex<double> gExP = (p->gamma)*ExP;

complex<double> gExP2 = gExP*gExP;

complex<double> den = gme*gme - gExP2;

complex<double> r2 = (p->x)*(p->x) + (p->y)*(p->y);

complex<double> ans = - sqrt2pi / pow(E + (p->x4)*tau*I, (p->l)+1) / sqrt(den) *

                        exp( (-(p->x1)*(p->x1) + (p->x2)*I)*tau +

                        (E - ExP) / (gme + gExP) * (p->x1) * ((2.0*(p->gamma) - 1.0)*(p->x1) + ((p->x) + (p->y))*I) +

                        (r2*(gExP*ExP - gme) + 2.0e0*(p->x)*(p->y)*ExP) / 2.0e0 / den);

//std::cout << "tau=" << tau << std::endl;

//std::cout << "pow=" << pow(E + (p->x4)*tau*I, (p->l)+1) << std::endl;

//std::cout << "arg=" << (-(p->x1)*(p->x1) + (p->x2)*I)*tau +
//
//                       (E - ExP) / (gme + gExP) * (p->x1) * ((2.0*(p->gamma) - 1.0)*(p->x1) + ((p->x) + (p->y))*I) +
//
//                        (r2*(gExP*ExP - gme) + 2.0e0*(p->x)*(p->y)*ExP) / 2.0e0 / den << std::endl;

//std::cout << ans << std::endl;

if (p->part == 0) return real(ans);
else              return imag(ans);
}

/*******************************************************************/

int calc_imnl_quad_ (int * m_max, int * n_max, int * l_max,
                     double * x1_re, double * x1_im,
                     double * x2_re, double * x2_im,
                     double * x3_re, double * x3_im,
                     double * x4_re, double * x4_im,
                     double * Imnl)
{
gsl_set_error_handler_off();

complex<double> x1(*x1_re, *x1_im);
complex<double> x2(*x2_re, *x2_im);
complex<double> x3(*x3_re, *x3_im);
complex<double> x4(*x4_re, *x4_im);

//bar x* values:
complex<double> alpha = sqrt(0.25e0*E - x3*I) - 0.5e0*E;

complex<double> ep2a = E + 2.0e0*alpha;

complex<double> gamma = alpha / ep2a;

x1 = x1 / pow(ep2a, 1.5e0);

x2 = x2 / ep2a + gamma * I;

x4 = x4 / ep2a;

quad_params qp = {x1, x2, x3, x4, 0, gamma, 0, 0, 0};

size_t limit = 1024;
double epsabs = 1.0e-14, epsrel = 1.0e-14, err;

gsl_integration_workspace * w = gsl_integration_workspace_alloc(limit);

gsl_function F;
F.function = &subint;
F.params = &qp;

double x = 0.0;
double y = 0.0;

int m = 0, n = 0;

complex<double> InT;

double tmax = 4.0e1;
double step = std::min(tmax, 1.0e2/std::abs(I*x2-x1*x1));
double t1 = 0.0e0, t2 = step;

std::cout << step;

for(int l = 0; l <= *l_max; ++l)
{
    qp.l = l;
    qp.x = x;
    qp.y = y;

    int ind = 2 * (m + (*m_max+1) * (n + (*n_max+1) * (l)));

    Imnl[ind+0] = 0.0e0;
    Imnl[ind+1] = 0.0e0;

    while(t2 <= tmax)
    {
        double I_re;
        qp.part = 0;
        //gsl_integration_qagiu(&F, 0.0e0, epsabs, epsrel, limit, w, &I_re, &err);
        gsl_integration_qag(&F, t1, t2, epsabs, epsrel, limit, GSL_INTEG_GAUSS21, w, &I_re, &err);

        double I_im;
        qp.part = 1;
        //gsl_integration_qagiu(&F, 0.0e0, epsabs, epsrel, limit, w, &I_im, &err);
        gsl_integration_qag(&F, t1, t2, epsabs, epsrel, limit, GSL_INTEG_GAUSS21, w, &I_im, &err);

        InT = (I_re + I_im * I) / pow(ep2a, 0.5e0*(3 + m + n));

        Imnl[ind+0] += real(InT);
        Imnl[ind+1] += imag(InT);

        t1 = t2;
        t2 += step;
    }
}

gsl_integration_workspace_free(w);

return 0;
}

/*******************************************************************/
