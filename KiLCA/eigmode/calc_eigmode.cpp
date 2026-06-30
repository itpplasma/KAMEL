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

intptr_t cd = ((det_params *) params)->cd;

const double fre = x[0];
const double fim = x[1];

complex<double> olab = 2.0*pi*(fre + I*fim);

char cd_sd_path2project[1024];

settings_get_path2project_ (core_data_get_sd_(cd), cd_sd_path2project);

core_data_set_mda_element_ (cd, ind, mode_data_create_ (m, nn, real(olab), imag(olab), core_data_get_sd_(cd), core_data_get_bp_(cd), cd_sd_path2project));

mode_data_calc_all_mode_data_ (core_data_get_mda_element_(cd, ind), 0);

complex<double> det = complex<double>(wave_data_get_det_re_(mode_data_get_wd_(core_data_get_mda_element_(cd, ind))), wave_data_get_det_im_(mode_data_get_wd_(core_data_get_mda_element_(cd, ind))));

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
mode_data_destroy_ (core_data_get_mda_element_(cd, ind));
core_data_set_mda_element_ (cd, ind, 0);
clear_all_data_in_mode_data_module_ (); //clean up fortran module data
}

/**********************************************************************************/

int find_det_zeros (int ind, int m, int n, intptr_t cd)
{
//output file:
char *full_name = new char[1024];
char es_fname[1024];
get_eigmode_fname_ (es_fname);
char cd_sd_path2project[1024];
settings_get_path2project_ (core_data_get_sd_(cd), cd_sd_path2project);
sprintf (full_name, "%s%s", cd_sd_path2project, es_fname);

FILE *out;
if (!(out = fopen (full_name, "w")))
{
    fprintf (stderr, "\nFailed to open file %s\a\n", full_name);
}

int status = 0;

const int N = 2;

det_params p = {ind, m, n, get_eigmode_delta_ (), cd};

double x_init[2];
double x_root[2];

complex<double> center, quad1, quad2;
double radius;

fprintf (out, "%%Re(f)\t\t\tIm(f)\t\t\tRe(det)\t\t\tIm(det)");
fclose (out);

int k;
for (k=get_eigmode_kmin_ (); k<=get_eigmode_kmax_ (); k++)
{
    if (!(out = fopen (full_name, "a")))
    {
        fprintf (stderr, "\nFailed to open file %s\a\n", full_name);
    }

    x_init[0] = get_eigmode_fstart_re_ (k);
    x_init[1] = get_eigmode_fstart_im_ (k);

    // fortnum_multiroot_hybrid does a hybrid Newton solve building the Jacobian
    // by central differences, matching the former gsl_multiroot_fdfsolver
    // (hybridsj) usage with the finite-difference eval_jac. The residual and
    // step-size tolerances reuse the eigmode-options eps values.
    status = fortnum_multiroot_hybrid (&eval_det, N, x_init, get_eigmode_eps_abs_ (),
                                       get_eigmode_eps_res_ (), (int) 1e6, x_root, &p);

    double det_out[2];
    eval_det (N, x_root, det_out, &p);

    fprintf (out, "\n%.20le  %.20le\t%.20le  %.20le", x_root[0], x_root[1], det_out[0], det_out[1]);

    fprintf (out, "\n%%status = %d", status);

    //check if it is really root:
    if (get_eigmode_test_roots_ () == 1)
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

int loop_over_frequences (int ind, int m, int n, intptr_t cd)
{
//output file:
char *full_name = new char[1024];
char es_fname[1024];
get_eigmode_fname_ (es_fname);
char cd_sd_path2project[1024];
settings_get_path2project_ (core_data_get_sd_(cd), cd_sd_path2project);
sprintf (full_name, "%s%s", cd_sd_path2project, es_fname);

FILE *out;
if (!(out = fopen (full_name, "w")))
{
    fprintf (stderr, "\nFailed to open file %s\a\n", full_name);
}

delete [] full_name;

fprintf (out, "%%iter\tRe(f)\t\t\tIm(f)\t\t\tRe(det)\t\t\tIm(det)");

int es_rdim = get_eigmode_rdim_ ();
int es_idim = get_eigmode_idim_ ();
double es_rfmin = get_eigmode_rfmin_ ();
double es_rfmax = get_eigmode_rfmax_ ();
double es_ifmin = get_eigmode_ifmin_ ();
double es_ifmax = get_eigmode_ifmax_ ();

for (int i=0; i<es_rdim; i++)
{
    double fre = es_rfmin + i*(es_rfmax - es_rfmin)/max(es_rdim-1, 1);

    for (int k=0; k<es_idim; k++)
    {
        double fim = es_ifmin + k*(es_ifmax - es_ifmin)/max(es_idim-1, 1);

        complex<double> olab = 2.0*pi*(fre + I*fim);

        char cd_sd_path2project[1024];

        settings_get_path2project_ (core_data_get_sd_(cd), cd_sd_path2project);

        core_data_set_mda_element_ (cd, ind, mode_data_create_ (m, n, real(olab), imag(olab), core_data_get_sd_(cd), core_data_get_bp_(cd), cd_sd_path2project));

        mode_data_calc_all_mode_data_ (core_data_get_mda_element_(cd, ind), 0);

        fprintf (out, "\n%6u\t%.20le  %.20le\t%.20le  %.20le", i*es_idim+k, fre, fim,
                 wave_data_get_det_re_(mode_data_get_wd_(core_data_get_mda_element_(cd, ind))), wave_data_get_det_im_(mode_data_get_wd_(core_data_get_mda_element_(cd, ind))));
        fflush (out);

        mode_data_destroy_ (core_data_get_mda_element_(cd, ind));
        core_data_set_mda_element_ (cd, ind, 0);
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
