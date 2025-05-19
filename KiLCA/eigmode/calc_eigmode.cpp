/*! \file calc_eigmode.cpp
    \brief The implementation of the functions declared in calc_eigmode.h
*/

#include "eigmode_sett.h"
#include "calc_eigmode.h"
#include "inout.h"
#include "mode.h"

#include <gsl/gsl_errno.h>
#include <gsl/gsl_matrix.h>

/**********************************************************************************/

int eval_det (const gsl_vector *x, void *params, gsl_vector *f)
{
int ind = ((det_params *) params)->ind;
int m = ((det_params *) params)->m;
int n = ((det_params *) params)->n;

core_data *cd = ((det_params *) params)->cd;

const double fre = gsl_vector_get (x, 0);
const double fim = gsl_vector_get (x, 1);

complex<double> olab = 2.0*pi*(fre + I*fim);

cd->mda[ind] = new mode_data (m, n, olab, (const settings_t *)cd->sd, (const background *)cd->bp);

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

gsl_vector_set (f, 0, real(det));
gsl_vector_set (f, 1, imag(det));

//clean up:
delete cd->mda[ind];
cd->mda[ind] = NULL;
clear_all_data_in_mode_data_module_ (); //clean up fortran module data

return GSL_SUCCESS;
}

/**********************************************************************************/

int eval_jac (const gsl_vector *x, void *params, gsl_matrix *J)
{
const double fre = gsl_vector_get (x, 0);
const double fim = gsl_vector_get (x, 1);

gsl_vector *freq = gsl_vector_alloc (2);
gsl_vector *det_rm = gsl_vector_alloc (2);
gsl_vector *det_rp = gsl_vector_alloc (2);
gsl_vector *det_im = gsl_vector_alloc (2);
gsl_vector *det_ip = gsl_vector_alloc (2);

double delta = ((det_params *) params)->delta;

//set f-delta:
gsl_vector_set (freq, 0, fre-delta);
gsl_vector_set (freq, 1, fim);
eval_det (freq, params, det_rm);

//set f+delta:
gsl_vector_set (freq, 0, fre+delta);
gsl_vector_set (freq, 1, fim);
eval_det (freq, params, det_rp);

//set f-i*delta:
gsl_vector_set (freq, 0, fre);
gsl_vector_set (freq, 1, fim-delta);
eval_det (freq, params, det_im);

//set f+i*delta:
gsl_vector_set (freq, 0, fre);
gsl_vector_set (freq, 1, fim+delta);
eval_det (freq, params, det_ip);

gsl_matrix_set (J, 0, 0, (gsl_vector_get (det_rp, 0)-gsl_vector_get (det_rm, 0))/2.0/delta);
gsl_matrix_set (J, 0, 1, (gsl_vector_get (det_ip, 0)-gsl_vector_get (det_im, 0))/2.0/delta);
gsl_matrix_set (J, 1, 0, (gsl_vector_get (det_rp, 1)-gsl_vector_get (det_rm, 1))/2.0/delta);
gsl_matrix_set (J, 1, 1, (gsl_vector_get (det_ip, 1)-gsl_vector_get (det_im, 1))/2.0/delta);

gsl_vector_free (freq);
gsl_vector_free (det_rm);
gsl_vector_free (det_rp);
gsl_vector_free (det_im);
gsl_vector_free (det_ip);

return GSL_SUCCESS;
}

/**********************************************************************************/

int eval_det_jac (const gsl_vector *x, void *params, gsl_vector *f, gsl_matrix *J)
{
eval_det (x, params, f);
eval_jac (x, params, J);

return GSL_SUCCESS;
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

//solver initialization:
//const gsl_multiroot_fsolver_type *T;
//gsl_multiroot_fsolver *s;
const gsl_multiroot_fdfsolver_type *T;
gsl_multiroot_fdfsolver *s;

int status = 0;
size_t iter = 0;

const size_t N = 2;

det_params p = {ind, m, n, es->delta, cd};

//gsl_multiroot_function f = {&eval_det, N, &p};
gsl_multiroot_function_fdf f = {&eval_det, &eval_jac, &eval_det_jac, N, &p};

double x_init[2];
gsl_vector *x = gsl_vector_alloc (N);

//T = gsl_multiroot_fsolver_hybrids;
//s = gsl_multiroot_fsolver_alloc (T, 2);
T = gsl_multiroot_fdfsolver_hybridsj;
s = gsl_multiroot_fdfsolver_alloc (T, 2);

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
    x_init[0] = real (es->fstart[k]);
    x_init[1] = imag (es->fstart[k]);

    gsl_vector_set (x, 0, x_init[0]);
    gsl_vector_set (x, 1, x_init[1]);

    //gsl_multiroot_fsolver_set (s, &f, x);
    gsl_multiroot_fdfsolver_set (s, &f, x);

#if DEBUG_FLAG
    print_fdf_search_state (iter, s);
#endif

    fprintf (out, "\n%5u\t%.20le  %.20le\t%.20le  %.20le", iter, gsl_vector_get (s->x, 0), gsl_vector_get (s->x, 1), gsl_vector_get (s->f, 0), gsl_vector_get (s->f, 1));

    do
    {
        iter++;

        //status = gsl_multiroot_fsolver_iterate (s);
        status = gsl_multiroot_fdfsolver_iterate (s);

#if DEBUG_FLAG
        print_fdf_search_state (iter, s);
#endif

        fprintf (out, "\n%5u\t%.20le  %.20le\t%.20le  %.20le", iter, gsl_vector_get (s->x, 0), gsl_vector_get (s->x, 1), gsl_vector_get (s->f, 0), gsl_vector_get (s->f, 1));
        fflush (out);

        if (status) /* check if solver is stuck */
        {
            break;
        }

        //stopping criteria:
        if (es->stop_flag == 0)
        {
            status = gsl_multiroot_test_residual (s->f, es->eps_res);
        }
        else if (es->stop_flag == 1)
        {
            //status = gsl_multiroot_test_delta (gsl_multiroot_fsolver_dx(s),
            //                                   gsl_multiroot_fsolver_root(s),
            //                                   eps_abs, eps_rel);
            status = gsl_multiroot_test_delta (gsl_multiroot_fdfsolver_dx(s),
                                               gsl_multiroot_fdfsolver_root(s),
                                               es->eps_abs, es->eps_rel);
        }
        else
        {}
    }
    while (status == GSL_CONTINUE && iter < 1e6);

#if DEBUG_FLAG
        printf ("\nstatus = %s", gsl_strerror (status));
#endif
    fprintf (out, "\n%%status = %s", gsl_strerror (status));

    //check if it is really root:
    if (es->test_roots == 1)
    {
        gsl_set_error_handler_off ();

        radius = 1.0e1;

        center = gsl_vector_get (s->x, 0) + I*gsl_vector_get (s->x, 1);

        //quad1 = calc_circle_integral_gkq (center, radius, inv_det, &p);
        quad1 = calc_circle_integral (center, radius, inv_det, &p);

        center += 5.0*radius*(1.0+I);

        //quad2 = calc_circle_integral_gkq (center, radius, inv_det, &p);
        quad2 = calc_circle_integral (center, radius, inv_det, &p);

        fprintf (stdout, "\nquad1 = %.20le %.20le\tquad2 = %.20le %.20le",
        real(quad1), imag(quad1), real(quad2), imag(quad2));

        fprintf (out, "\n%%quad1 = %.20le %.20le\tquad2 = %.20le %.20le\terr = %.20le",
        real(quad1), imag(quad1), real(quad2), imag(quad2), abs(quad1/quad2));
    }
    fflush (out);
    fclose (out);
}

//gsl_multiroot_fsolver_free (s);
gsl_multiroot_fdfsolver_free (s);
gsl_vector_free (x);

delete [] full_name;

return status;
}

/**********************************************************************************/

void print_f_search_state (size_t iter, gsl_multiroot_fsolver *s)
{
printf ("\niter = %4u f = (%le, %le) det = (%le, %le)\n", iter, gsl_vector_get (s->x, 0), gsl_vector_get (s->x, 1), gsl_vector_get (s->f, 0), gsl_vector_get (s->f, 1));
}

/**********************************************************************************/

void print_fdf_search_state (size_t iter, gsl_multiroot_fdfsolver *s)
{
printf ("\niter = %4u f = (%le, %le) det = (%le, %le)\n", iter, gsl_vector_get (s->x, 0), gsl_vector_get (s->x, 1), gsl_vector_get (s->f, 0), gsl_vector_get (s->f, 1));
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

complex<double> calc_circle_integral_gkq (complex<double> center, double radius, cmplx_func f, void *params)
{
gsl_integration_workspace *w = gsl_integration_workspace_alloc (100);

double result, error;

size_t neval;

complex<double> ans;

quad_params p = {f, params, center, radius, 0};

gsl_function F;
F.function = &func_on_circle;
F.params = &p;

p.part = 0;

//gsl_integration_qags (&F, 0.0, 2.0*pi, 1.0e-3, 1.0e-3, 100, w, &result, &error);
gsl_integration_qng (&F, 0.0, 2.0*pi, 1.0e-3, 1.0e-3, &result, &error, &neval);

printf ("\nresult = %.16le\terror = %.16le", result, error);
ans = result;

p.part = 1;

//gsl_integration_qags (&F, 0.0, 2.0*pi, 1.0e-3, 1.0e-3, 100, w, &result, &error);
gsl_integration_qng (&F, 0.0, 2.0*pi, 1.0e-3, 1.0e-3, &result, &error, &neval);

printf ("\nresult = %.16le\terror = %.16le", result, error);
ans += result*I;

gsl_integration_workspace_free (w);

return ans*I*radius;
}

/**********************************************************************************/

double func_on_circle (double phi, void *params)
{
quad_params *p = (quad_params *)params;

complex<double> expon, res;

expon = exp(I*phi);
res = p->func(p->center + (p->radius)*expon, p->params)*expon;

if (p->part == 0) return real(res);
else              return imag(res);
}

/**********************************************************************************/

complex<double> inv_det(complex<double> z, void *params)
{
gsl_vector *x = gsl_vector_alloc (2);
gsl_vector *f = gsl_vector_alloc (2);

gsl_vector_set (x, 0, real(z));
gsl_vector_set (x, 1, imag(z));

eval_det (x, params, f);

complex<double> det = gsl_vector_get (f, 0) + I*gsl_vector_get (f, 1);

gsl_vector_free (x);
gsl_vector_free (f);

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

        cd->mda[ind] = new mode_data (m, n, olab, (const settings_t *)cd->sd, (const background *)cd->bp);

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
gsl_vector *x = gsl_vector_alloc (2);
gsl_vector *f = gsl_vector_alloc (2);

gsl_vector_set (x, 0, 0.0e0);
gsl_vector_set (x, 1, *freq);

eval_det (x, params, f);

*absdet = gsl_vector_get (f, 1); //imaginary part

gsl_vector_free (x);
gsl_vector_free (f);
}

/**********************************************************************************/
