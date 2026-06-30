/*! \file core.h
    \brief C entry points for the top-level KiLCA orchestrator, now owned by
           the Fortran kilca_core_data_m module (KiLCA/core/core_data_m.f90).
           The former C++ core_data class has been translated away; still-C++
           callers (wave_code_interface.cpp, eigmode/*.cpp) hold an opaque
           intptr_t handle instead of a `core_data*`.
*/

#ifndef CORE_INCLUDE

#define CORE_INCLUDE

#include <cstring>
#include <cstdint>

#include "code_settings.h"

#include "settings.h"
#include "background.h"
#include "mode.h"

extern "C"
{
intptr_t core_data_create_ (const char *path);

void core_data_destroy_ (intptr_t handle);

void core_data_calc_and_set_mode_independent_ (intptr_t handle);

void core_data_calc_and_set_mode_dependent_antenna_ (intptr_t handle);

void core_data_calc_and_set_mode_dependent_eigmode_ (intptr_t handle);

void core_data_calc_and_set_mode_dependent_antenna_interface_ (intptr_t handle);

void core_data_calc_and_set_mode_dependent_antenna_interface_mn_ (intptr_t handle, int m, int n, int flag);

int core_data_get_dim_ (intptr_t handle);

intptr_t core_data_get_mda_element_ (intptr_t handle, int ind);

void core_data_set_mda_element_ (intptr_t handle, int ind, intptr_t val);

intptr_t core_data_get_bp_ (intptr_t handle);

intptr_t core_data_get_sd_ (intptr_t handle);

void core_data_get_path2project_ (intptr_t handle, char *buf);

void get_pointer_precision_ (int *);

void set_core_data_in_core_module_ (intptr_t *);

void set_background_in_core_module_ (background **);

void clear_all_data_in_mode_data_module_ (void);
}

/*******************************************************************/

#endif
