/*! \file imhd_zone.cpp
    \brief The implementation of imhd_zone class.
*/

#include "imhd_zone.h"
#include "mode.h"
#include "shared.h"
#include "inout.h"
#include "eval_back.h"

/*****************************************************************************/

void imhd_zone::read_settings (char * file)
{
read (file); //base class read function

//derived class specific:
FILE *in;

if ((in=fopen (file, "r"))==NULL)
{
    fprintf(stderr, "\nerror: hmedium_zone: read_settings: failed to open file %s\a\n", file);
    exit(0);
}

char *str_buf = new char[1024];

//skip lines (already readed and values are set):
read_line_2skip_it (in, &str_buf);
read_line_2skip_it (in, &str_buf);
read_line_2skip_it (in, &str_buf);
read_line_2skip_it (in, &str_buf);
read_line_2skip_it (in, &str_buf);
read_line_2skip_it (in, &str_buf);
read_line_2skip_it (in, &str_buf);
read_line_2skip_it (in, &str_buf);

//Solution settings:
read_line_2skip_it (in, &str_buf);
read_line_2get_int (in, &(max_dim));
read_line_2get_double (in, &(eps_rel));
read_line_2get_double (in, &(eps_abs));
read_line_2skip_it (in, &str_buf);

//Space out settings:
read_line_2skip_it (in, &str_buf);
read_line_2get_int (in, &(deg));
read_line_2get_double (in, &(reps));
read_line_2get_double (in, &(aeps));
read_line_2get_double (in, &(step));
read_line_2skip_it (in, &str_buf);

//Debugging settings:
read_line_2skip_it (in, &str_buf);
read_line_2get_int (in, &(flag_debug));
read_line_2skip_it (in, &str_buf);

fclose (in);

Nwaves = 2;
Ncomps = 8;

imhd_zone *pointer = this;

set_imhd_data_module_ (&pointer);

if (flag_debug) print_settings ();

delete [] str_buf;
}

/*****************************************************************************/

void imhd_zone::print_settings (void)
{
print (); //base class print function
}

/*****************************************************************************/

void imhd_zone::calc_basis_fields (int flag)
{
switch (version)
{
    case 0: //incompressible without background flows with gsl integrator
    calculate_basis_incompressible_gsl ();
    break;

    case 1: //incompressible with background flows with gsl integrator
    calculate_basis_flow_gsl ();
    break;

    default:
    fprintf (stderr, "\nerror: imhd_zone::calc_basis_fields: unknown code version.");
    exit (1);
}

//transform basis to the lab frame if needed!
}

/*****************************************************************************/

void imhd_zone::copy_E_and_B_fields (double *EB_p)
{
for (int node=0; node<dim; node++)
{
    for (int comp=0; comp<6; comp++) //Er, Et, Ez, Br, Bt, Bz
    {
        EB_p[iFFM(node, comp, 0)] = EB[iF(node, comp, 0)];
        EB_p[iFFM(node, comp, 1)] = EB[iF(node, comp, 1)];
    }
}
}

/*****************************************************************************/

double Ffunc (double r, void *params)
{
const imhd_zone *zone = (const imhd_zone *) params;

//wave staff:
double kt = (zone->wd->m)/r;
double kz = (zone->wd->n)/(zone->sd->bs->rtor);

double BtBz[2];

eval_Bt_Bz (r, zone->bp, BtBz);

return kt*BtBz[0] + kz*BtBz[1];
}

/*******************************************************************/

double Gfunc (double r, void *params)
{
const imhd_zone *zone = (const imhd_zone *) params;

//wave staff:
double kt = (zone->wd->m)/r;
double kz = (zone->wd->n)/(zone->sd->bs->rtor);

double BtBz[2];

eval_Bt_Bz (r, zone->bp, BtBz);

return kt*BtBz[1] - kz*BtBz[0];
}

/*******************************************************************/

complex<double> calc_k_vals (const imhd_zone *zone, double r, double *kt, double *kz,
                             double *ks, double *kp, double *k2, double *kB, complex<double> *kA)
{
int m = zone->wd->m;
int n = zone->wd->n;
complex<double> olab = zone->wd->olab;

//wave staff:
*kt = m/r;
*kz = n/(zone->sd->bs->rtor);
*k2 = (*kt)*(*kt) + (*kz)*(*kz);

//background stuff:
double R[5];

eval_B0_ht_hz_n0_Vz (r, zone->bp, R);

double B0 = R[0];
double ht = R[1];
double hz = R[2];
double n0 = R[3];
double Vz = R[4];

*ks = hz*(*kt) - ht*(*kz);
*kp = ht*(*kt) + hz*(*kz);

*kB = (*kp)*B0;

double VA = B0/sqrt(4.0*pi*(zone->sd->bs->mass[0])*n0);

*kA = (olab - (*kz)*Vz)/VA;

return E - ((*kA)*(*kA))/((*kp)*(*kp)); //kfac
}

/*******************************************************************/

void calc_k_vals_ (imhd_zone **zone_ptr, double *r, double *kt, double *kz,
                   double *ks, double *kp, double *k2, double *kB,
                   double *re_kA, double *im_kA, double *re_kfac, double *im_kfac)
{
const imhd_zone *zone = (const imhd_zone *)(*zone_ptr);

complex<double> kfac, kA;

kfac = calc_k_vals (zone, *r, kt, kz, ks, kp, k2, kB, &kA);

*re_kA = real(kA);
*im_kA = imag(kA);

*re_kfac = real(kfac);
*im_kfac = imag(kfac);
}

/*******************************************************************/
