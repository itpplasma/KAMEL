/*! \file flre_quants.h
    \brief C entry points for the auxiliary FLRE quantities (current/power/
           flux/density/torque), now owned by the Fortran
           kilca_flre_quants_m module (each instance addressed by an opaque
           handle, since each flre_zone owns its own). The former C++
           flre_quants class has been translated away.
*/

#ifndef FLRE_QUANTS_INCLUDE

#define FLRE_QUANTS_INCLUDE

#include <cstdint>

extern "C"
{
intptr_t flre_quants_create_ (intptr_t zone_cp, intptr_t zone_me, void *zone_bp,
                              const char *path2linear_p, int flreo, int dimx,
                              double *x_ptr, int ncomps, double *eb_mov_ptr,
                              double wd_omov_re, double wd_omov_im,
                              int bc1, int bc2, int zone_index);

void flre_quants_destroy_ (intptr_t handle);

void flre_quants_calculate_local_profiles_ (intptr_t handle);

void flre_quants_calculate_integrated_profiles_ (intptr_t handle);

void flre_quants_save_profiles_ (intptr_t handle);

void flre_quants_calculate_jae_ (intptr_t handle);

void flre_quants_transform_quants_to_lab_cyl_frame_ (intptr_t handle);

void flre_quants_interp_diss_power_density_ (intptr_t handle, double x, int type,
                                             int spec, double *dpd);

void flre_quants_interp_current_density_ (intptr_t handle, double x, int type,
                                          int spec, int comp, double *J);
}

#endif
