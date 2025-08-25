/*! \file
    \brief The definitions of functions declared in wave_code_interface.h.
*/

#include <stdlib.h>
#include <stdio.h>
#include <cmath>
#include <cstring>
#include <climits>

#include "core.h"
#include "mode.h"
#include "background.h"
#include "transforms.h"
#include "spline.h"
#include "flre_zone.h"
#include "wave_code_interface.h"

/*******************************************************************/

void calc_wave_code_data_ (core_data ** cdptr, char const * run_path, int * pathlength)
{
    //!The function computes wave fields and other quantities which migth be obtained by subsequent calls of other get_* () interface functions
    //gets path to the project
    char * path = new char[1024];
    strcpy (path, run_path);
    path[*pathlength] = '\0'; //end of string symbol

    if (path[strlen(path)-1] != '/') strcat(path, "/");

    //!Allocates core data structure containing pointers to all important code data
    core_data * cd = new core_data (path);
    set_core_data_in_core_module_ (&cd);
    *cdptr = cd;

    cd->calc_and_set_mode_independent_core_data ();
    cd->calc_and_set_mode_dependent_core_data_antenna_interface ();

    delete [] path;
}

/*******************************************************************/

void calc_wave_code_data_for_mode_ (core_data ** cdptr, char const * run_path, int * pathlength, int * m, int * n)
{
    //!The function computes wave fields and other quantities which migth be obtained by subsequent calls of other get_* () interface functions
    //gets path to the project
    char * path = new char[1024];
    strcpy (path, run_path);
    path[*pathlength] = '\0'; //end of string symbol
    if (path[strlen(path)-1] != '/') strcat(path, "/");

    //!Allocates core data structure containing pointers to all important code data
    core_data * cd = new core_data (path);
    set_core_data_in_core_module_ (&cd);
    *cdptr = cd;

    cd->calc_and_set_mode_independent_core_data ();

    cd->calc_and_set_mode_dependent_core_data_antenna_interface (*m, *n);

    delete [] path;
}

/*******************************************************************/

void clear_wave_code_data_ (core_data ** cdptr)
{
//!The functions deletes all wave code data
delete *cdptr;
}

/*******************************************************************/

void get_basic_background_profiles_from_wave_code_ (core_data ** cdptr, int * dim_r, double * r,
                                                    double * q, double * n,
                                                    double * Ti, double * Te,
                                                    double * Vth, double * Vz, double * dPhi0)
{
core_data * cd = *cdptr;

cd->bp->interp_basic_background_profiles_in_lab_frame (*dim_r, r, q, n, Ti, Te, Vth, Vz, dPhi0);
}

/*******************************************************************/

void get_wave_vectors_from_wave_code_ (core_data ** cdptr, int * dim_r, double * r,
                                       int * m, int * n, double * ks, double * kp)
{
core_data *cd = *cdptr;

for (int i=0; i<*dim_r; i++)
{
    double kth = (*m)/r[i];
    double kz  = (*n)/(cd->sd->bs->rtor);

    spline_eval_ (cd->bp->sid, 1, r+i, 0, 0, cd->bp->i_hth, cd->bp->i_hz, cd->bp->R);

    double ht = cd->bp->R[0];
    double hz = cd->bp->R[1];

    kp[i] = ht*kth + hz*kz;
    ks[i] = hz*kth - ht*kz;
}
}

/*******************************************************************/

void get_background_magnetic_fields_from_wave_code_ (core_data ** cdptr, int * dim_r, double * r,
                                                     double * Bt, double * Bz, double * B0)
{
core_data * cd = *cdptr;

for (int i=0; i<*dim_r; i++)
{
    spline_eval_ (cd->bp->sid, 1, r+i, 0, 0, cd->bp->i_Bth, cd->bp->i_B, cd->bp->R);

    Bt[i] = cd->bp->R[0];
    Bz[i] = cd->bp->R[1];
    B0[i] = cd->bp->R[2];
}
}

/*******************************************************************/

void get_collision_frequences_from_wave_code_ (core_data ** cdptr, int * dim_r, double * r,
                                               double * nui, double * nue)
{
core_data * cd = *cdptr;

for (int i=0; i<*dim_r; i++)
{
    spline_eval_ (cd->bp->sid, 1, r+i, 0, 0, cd->bp->i_nu[0], cd->bp->i_nu[0], nui+i);
    spline_eval_ (cd->bp->sid, 1, r+i, 0, 0, cd->bp->i_nu[1], cd->bp->i_nu[1], nue+i);
}
}

/*******************************************************************/

void get_wave_fields_from_wave_code_ (core_data ** cdptr, int * dim_r, double * r,
                                      int * m, int * n,
                                      double * Er, double * Es, double * Ep,
                                      double * Et, double * Ez,
                                      double * Br, double * Bs, double * Bp,
                                      double * Bt, double * Bz)
{
core_data *cd = *cdptr;

//search for mode index
int num = -1;

for (int i=0; i<cd->dim; i++)
{
    if (cd->mda[i]->wd->m == *m && cd->mda[i]->wd->n == *n)
    {
        num = i;
        break;
    }
}

if (num == -1)
{
    fprintf (stdout, "\nwarning: get_wave_fields_from_wave_code: (%d, %d) mode is not in the spectrum.", *m, *n);
    return;
}

for (int i=0; i<*dim_r; i++)
{
    complex<double> EBcyl[6];
    complex<double> EBrsp[6];

    cd->mda[num]->eval_EB_fields (r[i], EBcyl);

    transform_EB_from_cyl_to_rsp (cd->bp, r[i], EBcyl, EBrsp);

    //copy values to fortran complex arrays passed as C double arrays:
    int ind = 2*i;

    Er[ind+0] = real (EBcyl[0]);
    Er[ind+1] = imag (EBcyl[0]);

    Es[ind+0] = real (EBrsp[1]);
    Es[ind+1] = imag (EBrsp[1]);

    Ep[ind+0] = real (EBrsp[2]);
    Ep[ind+1] = imag (EBrsp[2]);

    Et[ind+0] = real (EBcyl[1]);
    Et[ind+1] = imag (EBcyl[1]);

    Ez[ind+0] = real (EBcyl[2]);
    Ez[ind+1] = imag (EBcyl[2]);

    Br[ind+0] = real (EBcyl[3]);
    Br[ind+1] = imag (EBcyl[3]);

    Bs[ind+0] = real (EBrsp[4]);
    Bs[ind+1] = imag (EBrsp[4]);

    Bp[ind+0] = real (EBrsp[5]);
    Bp[ind+1] = imag (EBrsp[5]);

    Bt[ind+0] = real (EBcyl[4]);
    Bt[ind+1] = imag (EBcyl[4]);

    Bz[ind+0] = real (EBcyl[5]);
    Bz[ind+1] = imag (EBcyl[5]);
}
}

/*******************************************************************/

void get_diss_power_density_from_wave_code_ (core_data ** cdptr, int * dim_r, double * r,
                                             int * m, int * n,
                                             int * type, int * spec, double * pdis)
{
core_data *cd = *cdptr;

//search for mode index
int num = -1;

for (int i=0; i<cd->dim; i++)
{
    if (cd->mda[i]->wd->m == *m && cd->mda[i]->wd->n == *n)
    {
        num = i;
        break;
    }
}

if (num == -1)
{
    fprintf (stdout, "\nwarning: get_wave_fields_from_wave_code: (%d, %d) mode is not in the spectrum.", *m, *n);
    return;
}

for (int i=0; i<*dim_r; i++)
{
    cd->mda[num]->eval_diss_power_density (r[i], *type, *spec, pdis+i);
}
}

/*******************************************************************/

void get_antenna_spectrum_dim_ (core_data ** cdptr, int * dim_mn)
{
core_data *cd = *cdptr;

*dim_mn = cd->dim;
}

/*******************************************************************/

void get_antenna_spectrum_numbers_ (core_data ** cdptr, int * dim_mn, int * m_vals, int * n_vals)
{
core_data *cd = *cdptr;

if (*dim_mn != cd->dim)
{
   fprintf (stdout, "\nwarning: get_antenna_spectrum_numbers: spectrum dimensions must match");
   return;
}

for (int i=0; i<cd->dim; i++)
{
    m_vals[i] = cd->mda[i]->wd->m;
    n_vals[i] = cd->mda[i]->wd->n;
}
}

/*******************************************************************/

void get_current_densities_from_wave_code_ (core_data ** cdptr, int * dim_r, double * r,
                                            int * m, int * n,
                                            double * Jri, double * Jsi, double * Jpi,
                                            double * Jre, double * Jse, double * Jpe)
{
core_data *cd = *cdptr;

//search for mode index
int num = -1;

for (int i=0; i<cd->dim; i++)
{
    if (cd->mda[i]->wd->m == *m && cd->mda[i]->wd->n == *n)
    {
        num = i;
        break;
    }
}

if (num == -1)
{
    fprintf (stdout, "\nwarning: get_current_densities_from_wave_code: (%d, %d) mode is not in the spectrum.", *m, *n);
    return;
}

for (int i=0; i<*dim_r; i++)
{
    int ind = 2*i;
    cd->mda[num]->eval_current_density (r[i], 0, 0, 0, Jri+ind);
    cd->mda[num]->eval_current_density (r[i], 0, 0, 1, Jsi+ind);
    cd->mda[num]->eval_current_density (r[i], 0, 0, 2, Jpi+ind);
    cd->mda[num]->eval_current_density (r[i], 0, 1, 0, Jre+ind);
    cd->mda[num]->eval_current_density (r[i], 0, 1, 1, Jse+ind);
    cd->mda[num]->eval_current_density (r[i], 0, 1, 2, Jpe+ind);
}
}

/*******************************************************************/

void activate_kilca_modules_for_flre_zone_ (core_data ** cdptr)
{
//!The function activates the fortran modules for flre zone

flre_zone * fz = static_cast<flre_zone *>((*cdptr)->mda[0]->zones[0]);

activate_fortran_modules_for_zone_ (&fz);
}

/*******************************************************************/

void deactivate_kilca_modules_for_flre_zone_ (core_data ** cdptr)
{
//!The function activates the fortran modules for flre zone

flre_zone * fz = static_cast<flre_zone *>((*cdptr)->mda[0]->zones[0]);

deactivate_fortran_modules_for_zone_ (&fz);
}

/*******************************************************************/

void
get_kilca_conductivity_array_
(
core_data ** cdptr,
int * m,
int * n,
int * zone_ind,
int * spec,
int * flreo,
int * dim,
double ** r,
double ** cptr
)
{
core_data *cd = *cdptr;

//search for mode index
int num = -1;

for (int i=0; i<cd->dim; i++)
{
    if (cd->mda[i]->wd->m == *m && cd->mda[i]->wd->n == *n)
    {
        num = i;
        break;
    }
}

if (num == -1)
{
    fprintf (stdout, "\nwarning: get_kilca_conductivity_array: (%d, %d) mode is not in the spectrum.", *m, *n);
    return;
}

flre_zone * zone = static_cast<flre_zone *>(cd->mda[num]->zones[*zone_ind]);

*flreo = zone->flre_order;
*dim   = zone->cp->dimx;
*r     = zone->cp->x;
*cptr  = &(zone->cp->K[zone->cp->iKs(*spec, 0, 0, 0, 0, 0, 0, 0)]);
}

/*******************************************************************/

void
calc_conductivity_matrices_for_mode_
(core_data ** cdptr, char const * run_path, int * pathlength, int * m, int * n)
{
//!The function computes wave fields and other quantities which migth be obtained by subsequent calls of other get_* () interface functions

//gets path to the project
char * path = new char[1024];

strcpy (path, run_path);

path[*pathlength] = '\0'; //end of string symbol

if (path[strlen(path)-1] != '/') strcat(path, "/");

//!Allocates core data structure containing pointers to all important code data
core_data * cd = new core_data (path);
set_core_data_in_core_module_ (&cd);
*cdptr = cd;

cd->calc_and_set_mode_independent_core_data ();

cd->calc_and_set_mode_dependent_core_data_antenna_interface (*m, *n, 1);

delete [] path;
}

/*******************************************************************/

void get_mode_parameters_ (core_data ** cdptr, int * m, int * n, double * kz, double * omega_mov_re, double * omega_mov_im)
{
core_data *cd = *cdptr;

//search for mode index
int num = -1;

for (int i=0; i<cd->dim; i++)
{
    if (cd->mda[i]->wd->m == *m && cd->mda[i]->wd->n == *n)
    {
        num = i;
        break;
    }
}

if (num == -1)
{
    fprintf (stdout, "\nwarning: get_kilca_conductivity_array: (%d, %d) mode is not in the spectrum.", *m, *n);
    return;
}

*kz = cd->mda[num]->wd->n / cd->bp->sd->bs->rtor;

*omega_mov_re = real(cd->mda[num]->wd->omov);
*omega_mov_im = imag(cd->mda[num]->wd->omov);
}

/*******************************************************************/

void set_wave_parameters_ (core_data ** cdptr, int * m, int * n)
{
core_data *cd = *cdptr;

//search for mode index
int num = -1;

for (int i=0; i<cd->dim; i++)
{
    if (cd->mda[i]->wd->m == *m && cd->mda[i]->wd->n == *n)
    {
        num = i;
        break;
    }
}

if (num == -1)
{
    fprintf (stdout, "\nwarning: set_wave_parameters: (%d, %d) mode is not in the spectrum.", *m, *n);
    return;
}

//set_wave_parameters_in_mode_data_module_(&(cd->mda[num]->wd->m), &(cd->mda[num]->wd->n),
//                                         &(real(cd->mda[num]->wd->olab)), &(imag(cd->mda[num]->wd->olab)),
//                                         &(real(cd->mda[num]->wd->omov)), &(imag(cd->mda[num]->wd->omov)));

int mm = cd->mda[num]->wd->m, nn = cd->mda[num]->wd->n;

double olab_re = real(cd->mda[num]->wd->olab), olab_im = imag(cd->mda[num]->wd->olab);
double omov_re = real(cd->mda[num]->wd->omov), omov_im = imag(cd->mda[num]->wd->omov);

set_wave_parameters_in_mode_data_module_(&mm, &nn, &olab_re, &olab_im, &omov_re, &omov_im);
}

/*******************************************************************/

void unset_wave_parameters_ (void)
{
clear_all_data_in_mode_data_module_();
}

/*******************************************************************/
