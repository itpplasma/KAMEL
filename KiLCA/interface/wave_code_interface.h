/*! \file wave_code_interface.h
    \brief The declarations of functions implementing Fortran interface to KiLCA library functions (for balance code).
*/

#include <stdlib.h>
#include <stdio.h>
#include <cmath>
#include <cstring>
#include <climits>
#include <inttypes.h>

#include "core.h"

/*******************************************************************/

extern "C"
{
void calc_wave_code_data_ (core_data ** cdptr, char const * run_path, int * pathlength);

void calc_wave_code_data_for_mode_ (core_data ** cdptr, char const * run_path, int * pathlength, int * m, int * n);

void clear_wave_code_data_ (core_data ** cdptr);

void get_basic_background_profiles_from_wave_code_ (core_data ** cdptr, int * dim_r,
                                                    double * r, double * q, double * n,
                                                    double * Ti, double * Te, double * Vth,
                                                    double * Vz, double * dPhi0);

void get_wave_vectors_from_wave_code_ (core_data ** cdptr, int * dim_r, double * r,
                                       int * m, int * n, double * ks, double * kp);

void get_background_magnetic_fields_from_wave_code_ (core_data ** cdptr, int * dim_r,
                                                     double * r, double *Bt, double *Bz,
                                                     double *B0);

void get_wave_fields_from_wave_code_ (core_data ** cdptr, int * dim_r, double * r,
                                      int * m, int * n,
                                      double * Er, double * Es, double * Ep,
                                      double * Et, double * Ez,
                                      double * Br, double * Bs, double * Bp,
                                      double * Bt, double * Bz);

void get_diss_power_density_from_wave_code_ (core_data ** cdptr, int * dim_r, double * r,
                                             int * m, int * n,
                                             int * type, int * spec, double * pdis);

void get_antenna_spectrum_dim_ (core_data ** cdptr, int * dim_mn);

void get_antenna_spectrum_numbers_ (core_data ** cdptr, int * dim_mn, int * m_vals, int * nvals);

void get_collision_frequences_from_wave_code_ (core_data ** cdptr, int * dim_r, double * r,
                                               double * nui, double * nue);

void get_current_densities_from_wave_code_ (core_data ** cdptr, int * dim_r, double * r,
                                            int * m, int * n,
                                            double * Jri, double * Jsi, double * Jpi,
                                            double * Jre, double * Jse, double * Jpe);

//get profiles from balance code:
void get_background_dimension_from_balance_ (int *dimx);

void get_background_profiles_from_balance_ (int *dimx, double *x, double *q, double *n,
                                            double *Ti, double *Te,
                                            double *Vth, double *Vz, double *Er);

void activate_kilca_modules_for_flre_zone_ (core_data ** cdptr);

void deactivate_kilca_modules_for_flre_zone_ (core_data ** cdptr);

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
);

void
calc_conductivity_matrices_for_mode_
(core_data ** cdptr, char const * run_path, int * pathlength, int * m, int * n);

void get_mode_parameters_ (core_data ** cdptr, int * m, int * n, double * kz, double * omega_mov_re, double * omega_mov_im);

void set_wave_parameters_ (core_data ** cdptr, int * m, int * n);

void unset_wave_parameters_ (void);
}

/*******************************************************************/
