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
#include "eval_back.h"
#include "transforms_dispatch.h"
#include "spline.h"
#include "flre_zone_dispatch.h"
#include "cond_profs.h"
#include "wave_code_interface.h"

/*******************************************************************/

void calc_wave_code_data_ (intptr_t * cdptr, char const * run_path, int * pathlength)
{
    //!The function computes wave fields and other quantities which migth be obtained by subsequent calls of other get_* () interface functions
    //gets path to the project
    char * path = new char[1024];
    strcpy (path, run_path);
    path[*pathlength] = '\0'; //end of string symbol

    if (path[strlen(path)-1] != '/') strcat(path, "/");

    //!Allocates core data structure containing pointers to all important code data
    intptr_t cd = core_data_create_ (path);
    set_core_data_in_core_module_ (&cd);
    *cdptr = cd;

    core_data_calc_and_set_mode_independent_ (cd);
    core_data_calc_and_set_mode_dependent_antenna_interface_ (cd);

    delete [] path;
}

/*******************************************************************/

void calc_wave_code_data_for_mode_ (intptr_t * cdptr, char const * run_path, int * pathlength, int * m, int * n)
{
    //!The function computes wave fields and other quantities which migth be obtained by subsequent calls of other get_* () interface functions
    //gets path to the project
    char * path = new char[1024];
    strcpy (path, run_path);
    path[*pathlength] = '\0'; //end of string symbol
    if (path[strlen(path)-1] != '/') strcat(path, "/");

    //!Allocates core data structure containing pointers to all important code data
    intptr_t cd = core_data_create_ (path);
    set_core_data_in_core_module_ (&cd);
    *cdptr = cd;

    core_data_calc_and_set_mode_independent_ (cd);

    core_data_calc_and_set_mode_dependent_antenna_interface_mn_ (cd, *m, *n, 0);

    delete [] path;
}

/*******************************************************************/

void clear_wave_code_data_ (intptr_t * cdptr)
{
//!The functions deletes all wave code data
core_data_destroy_ (*cdptr);
}

/*******************************************************************/

void get_basic_background_profiles_from_wave_code_ (intptr_t * cdptr, int * dim_r, double * r,
                                                    double * q, double * n,
                                                    double * Ti, double * Te,
                                                    double * Vth, double * Vz, double * dPhi0)
{
intptr_t cd = *cdptr;

interp_basic_background_profiles_in_lab_frame_ (*dim_r, r, q, n, Ti, Te, Vth, Vz, dPhi0);
}

/*******************************************************************/

void get_wave_vectors_from_wave_code_ (intptr_t * cdptr, int * dim_r, double * r,
                                       int * m, int * n, double * ks, double * kp)
{
intptr_t cd = *cdptr;

for (int i=0; i<*dim_r; i++)
{
    double kth = (*m)/r[i];
    double kz  = (*n)/(get_background_rtor_());

    double htz[2];
    eval_hthz (r[i], 0, 0, reinterpret_cast<const background *>(core_data_get_bp_(cd)), htz);

    double ht = htz[0];
    double hz = htz[1];

    kp[i] = ht*kth + hz*kz;
    ks[i] = hz*kth - ht*kz;
}
}

/*******************************************************************/

void get_background_magnetic_fields_from_wave_code_ (intptr_t * cdptr, int * dim_r, double * r,
                                                     double * Bt, double * Bz, double * B0)
{
intptr_t cd = *cdptr;

for (int i=0; i<*dim_r; i++)
{
    get_background_magnetic_fields_ (r[i], Bt+i, Bz+i, B0+i);
}
}

/*******************************************************************/

void get_collision_frequences_from_wave_code_ (intptr_t * cdptr, int * dim_r, double * r,
                                               double * nui, double * nue)
{
intptr_t cd = *cdptr;

for (int i=0; i<*dim_r; i++)
{
    get_background_collision_freqs_ (r[i], nui+i, nue+i);
}
}

/*******************************************************************/

void get_wave_fields_from_wave_code_ (intptr_t * cdptr, int * dim_r, double * r,
                                      int * m, int * n,
                                      double * Er, double * Es, double * Ep,
                                      double * Et, double * Ez,
                                      double * Br, double * Bs, double * Bp,
                                      double * Bt, double * Bz)
{
intptr_t cd = *cdptr;

//search for mode index
int num = -1;

for (int i=0; i<core_data_get_dim_(cd); i++)
{
    if (wave_data_get_m_(mode_data_get_wd_(core_data_get_mda_element_(cd, i))) == *m && wave_data_get_n_(mode_data_get_wd_(core_data_get_mda_element_(cd, i))) == *n)
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

    mode_data_eval_EB_fields_ (core_data_get_mda_element_(cd, num), r[i], reinterpret_cast<double*>(EBcyl));

    transform_EB_from_cyl_to_rsp (reinterpret_cast<const background *>(core_data_get_bp_(cd)), r[i], EBcyl, EBrsp);

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

void get_diss_power_density_from_wave_code_ (intptr_t * cdptr, int * dim_r, double * r,
                                             int * m, int * n,
                                             int * type, int * spec, double * pdis)
{
intptr_t cd = *cdptr;

//search for mode index
int num = -1;

for (int i=0; i<core_data_get_dim_(cd); i++)
{
    if (wave_data_get_m_(mode_data_get_wd_(core_data_get_mda_element_(cd, i))) == *m && wave_data_get_n_(mode_data_get_wd_(core_data_get_mda_element_(cd, i))) == *n)
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
    mode_data_eval_diss_power_density_ (core_data_get_mda_element_(cd, num), r[i], *type, *spec, pdis+i);
}
}

/*******************************************************************/

void get_antenna_spectrum_dim_ (intptr_t * cdptr, int * dim_mn)
{
intptr_t cd = *cdptr;

*dim_mn = core_data_get_dim_(cd);
}

/*******************************************************************/

void get_antenna_spectrum_numbers_ (intptr_t * cdptr, int * dim_mn, int * m_vals, int * n_vals)
{
intptr_t cd = *cdptr;

if (*dim_mn != core_data_get_dim_(cd))
{
   fprintf (stdout, "\nwarning: get_antenna_spectrum_numbers: spectrum dimensions must match");
   return;
}

for (int i=0; i<core_data_get_dim_(cd); i++)
{
    m_vals[i] = wave_data_get_m_(mode_data_get_wd_(core_data_get_mda_element_(cd, i)));
    n_vals[i] = wave_data_get_n_(mode_data_get_wd_(core_data_get_mda_element_(cd, i)));
}
}

/*******************************************************************/

void get_current_densities_from_wave_code_ (intptr_t * cdptr, int * dim_r, double * r,
                                            int * m, int * n,
                                            double * Jri, double * Jsi, double * Jpi,
                                            double * Jre, double * Jse, double * Jpe)
{
intptr_t cd = *cdptr;

//search for mode index
int num = -1;

for (int i=0; i<core_data_get_dim_(cd); i++)
{
    if (wave_data_get_m_(mode_data_get_wd_(core_data_get_mda_element_(cd, i))) == *m && wave_data_get_n_(mode_data_get_wd_(core_data_get_mda_element_(cd, i))) == *n)
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
    mode_data_eval_current_density_ (core_data_get_mda_element_(cd, num), r[i], 0, 0, 0, Jri+ind);
    mode_data_eval_current_density_ (core_data_get_mda_element_(cd, num), r[i], 0, 0, 1, Jsi+ind);
    mode_data_eval_current_density_ (core_data_get_mda_element_(cd, num), r[i], 0, 0, 2, Jpi+ind);
    mode_data_eval_current_density_ (core_data_get_mda_element_(cd, num), r[i], 0, 1, 0, Jre+ind);
    mode_data_eval_current_density_ (core_data_get_mda_element_(cd, num), r[i], 0, 1, 1, Jse+ind);
    mode_data_eval_current_density_ (core_data_get_mda_element_(cd, num), r[i], 0, 1, 2, Jpe+ind);
}
}

/*******************************************************************/

void activate_kilca_modules_for_flre_zone_ (intptr_t * cdptr)
{
//!The function activates the fortran modules for flre zone

intptr_t fz = mode_data_get_zone_handle_(core_data_get_mda_element_(*cdptr, 0), 0);

activate_fortran_modules_for_zone_ (&fz);
}

/*******************************************************************/

void deactivate_kilca_modules_for_flre_zone_ (intptr_t * cdptr)
{
//!The function activates the fortran modules for flre zone

intptr_t fz = mode_data_get_zone_handle_(core_data_get_mda_element_(*cdptr, 0), 0);

deactivate_fortran_modules_for_zone_ (&fz);
}

/*******************************************************************/

void
get_kilca_conductivity_array_
(
intptr_t * cdptr,
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
intptr_t cd = *cdptr;

//search for mode index
int num = -1;

for (int i=0; i<core_data_get_dim_(cd); i++)
{
    if (wave_data_get_m_(mode_data_get_wd_(core_data_get_mda_element_(cd, i))) == *m && wave_data_get_n_(mode_data_get_wd_(core_data_get_mda_element_(cd, i))) == *n)
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

intptr_t zone = mode_data_get_zone_handle_(core_data_get_mda_element_(cd, num), *zone_ind);
intptr_t zone_cp = flre_zone_get_cp_ (zone);

*flreo = flre_zone_get_flre_order_ (zone);
*dim   = get_cond_dimx_ (zone_cp);
*r     = get_cond_x_ptr_ (zone_cp);
*cptr  = get_cond_k_ptr_ (zone_cp) + get_cond_iks_ (zone_cp, *spec, 0, 0, 0, 0, 0, 0, 0);
}

/*******************************************************************/

void
calc_conductivity_matrices_for_mode_
(intptr_t * cdptr, char const * run_path, int * pathlength, int * m, int * n)
{
//!The function computes wave fields and other quantities which migth be obtained by subsequent calls of other get_* () interface functions

//gets path to the project
char * path = new char[1024];

strcpy (path, run_path);

path[*pathlength] = '\0'; //end of string symbol

if (path[strlen(path)-1] != '/') strcat(path, "/");

//!Allocates core data structure containing pointers to all important code data
intptr_t cd = core_data_create_ (path);
set_core_data_in_core_module_ (&cd);
*cdptr = cd;

core_data_calc_and_set_mode_independent_ (cd);

core_data_calc_and_set_mode_dependent_antenna_interface_mn_ (cd, *m, *n, 1);

delete [] path;
}

/*******************************************************************/

void get_mode_parameters_ (intptr_t * cdptr, int * m, int * n, double * kz, double * omega_mov_re, double * omega_mov_im)
{
intptr_t cd = *cdptr;

//search for mode index
int num = -1;

for (int i=0; i<core_data_get_dim_(cd); i++)
{
    if (wave_data_get_m_(mode_data_get_wd_(core_data_get_mda_element_(cd, i))) == *m && wave_data_get_n_(mode_data_get_wd_(core_data_get_mda_element_(cd, i))) == *n)
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

*kz = wave_data_get_n_(mode_data_get_wd_(core_data_get_mda_element_(cd, num))) / get_background_rtor_();

*omega_mov_re = get_wave_data_obj_omov_re_(mode_data_get_wd_(core_data_get_mda_element_(cd, num)));
*omega_mov_im = get_wave_data_obj_omov_im_(mode_data_get_wd_(core_data_get_mda_element_(cd, num)));
}

/*******************************************************************/

void set_wave_parameters_ (intptr_t * cdptr, int * m, int * n)
{
intptr_t cd = *cdptr;

//search for mode index
int num = -1;

for (int i=0; i<core_data_get_dim_(cd); i++)
{
    if (wave_data_get_m_(mode_data_get_wd_(core_data_get_mda_element_(cd, i))) == *m && wave_data_get_n_(mode_data_get_wd_(core_data_get_mda_element_(cd, i))) == *n)
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

//set_wave_parameters_in_mode_data_module_(&(mode_data_get_wd_(core_data_get_mda_element_(cd, num))->m), &(mode_data_get_wd_(core_data_get_mda_element_(cd, num))->n),
//                                         &(real(mode_data_get_wd_(core_data_get_mda_element_(cd, num))->olab)), &(imag(mode_data_get_wd_(core_data_get_mda_element_(cd, num))->olab)),
//                                         &(real(mode_data_get_wd_(core_data_get_mda_element_(cd, num))->omov)), &(imag(mode_data_get_wd_(core_data_get_mda_element_(cd, num))->omov)));

int mm = wave_data_get_m_(mode_data_get_wd_(core_data_get_mda_element_(cd, num))), nn = wave_data_get_n_(mode_data_get_wd_(core_data_get_mda_element_(cd, num)));

double olab_re = wave_data_get_olab_re_(mode_data_get_wd_(core_data_get_mda_element_(cd, num))), olab_im = wave_data_get_olab_im_(mode_data_get_wd_(core_data_get_mda_element_(cd, num)));
double omov_re = get_wave_data_obj_omov_re_(mode_data_get_wd_(core_data_get_mda_element_(cd, num))), omov_im = get_wave_data_obj_omov_im_(mode_data_get_wd_(core_data_get_mda_element_(cd, num)));

set_wave_parameters_in_mode_data_module_(&mm, &nn, &olab_re, &olab_im, &omov_re, &omov_im);
}

/*******************************************************************/

void unset_wave_parameters_ (void)
{
clear_all_data_in_mode_data_module_();
}

/*******************************************************************/
