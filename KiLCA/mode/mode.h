/*! \file
    \brief C entry points for the per-mode orchestrator, now owned by the
           Fortran kilca_mode_data_m module (KiLCA/mode/mode_data_m.f90).
           The former C++ mode_data class has been translated away; still-C++
           callers (core.cpp, eigmode/*.cpp, wave_code_interface.cpp) hold an
           opaque intptr_t handle instead of a `mode_data*` and dispatch
           through the functions below.
*/

#ifndef MODE_INCLUDE

#define MODE_INCLUDE

#include <complex>
#include <cstdint>

#include "settings.h"
#include "background.h"
#include "wave_data_dispatch.h"
#include "zone_dispatch.h"

/*****************************************************************************/

extern "C"
{
intptr_t mode_data_create_ (int m, int n, double olab_re, double olab_im,
                            intptr_t sd_ptr, intptr_t bp_ptr, const char *path2project);

void mode_data_destroy_ (intptr_t handle);

void mode_data_calc_all_mode_data_ (intptr_t handle, int flag);

intptr_t mode_data_get_wd_ (intptr_t handle);

intptr_t mode_data_get_zone_handle_ (intptr_t handle, int zone_ind);

void mode_data_eval_EB_fields_ (intptr_t handle, double x, double *EB_out);

void mode_data_eval_diss_power_density_ (intptr_t handle, double x, int type, int spec, double *dpd);

void mode_data_eval_current_density_ (intptr_t handle, double x, int type, int spec, int comp, double *J);

void set_settings_in_mode_data_module_ (const settings **);

void set_back_profiles_in_mode_data_module_ (const background **);

void set_wave_parameters_in_mode_data_module_ (int *, int *, double *, double *, double *, double *);

void set_current_density_in_antenna_module_ (void);

void clear_all_data_in_mode_data_module_ (void);
}

/*****************************************************************************/

#endif
