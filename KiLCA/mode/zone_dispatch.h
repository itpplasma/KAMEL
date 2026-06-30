/*! \file
    \brief Declarations for the Fortran zone_t dispatch shim (kilca_zone_m,
    KiLCA/mode/zone_m.f90) and its concrete-subtype factories. Replaces
    zone.h/hmedium_zone.h/imhd_zone.h/flre_zone.h's C++ class declarations
    now that the whole zone hierarchy is Fortran (F2003 type-extension
    polymorphism). Zone handles are opaque intptr_t values (1-based indices
    into a Fortran-side pool, not addresses); still-C++ callers (mode.cpp,
    calc_mode.cpp, wave_code_interface.cpp) treat them exactly as they
    treated the old `zone*`/`flre_zone*` pointers - construct via a
    *_zone_create_ factory, dispatch virtual-equivalent calls through the
    zone_<method>_ shims below, never dereference directly.

    This header is temporary: S6 removes it once mode_data itself becomes
    Fortran and can hold class(zone_t) natively, dropping the shim layer.
*/

#ifndef ZONE_DISPATCH_INCLUDE
#define ZONE_DISPATCH_INCLUDE

#include <cstdint>

#define PLASMA_MODEL_VACUUM 0
#define PLASMA_MODEL_MEDIUM 1
#define PLASMA_MODEL_IMHD 2
#define PLASMA_MODEL_RMHD 3
#define PLASMA_MODEL_FLRE 4

#define BOUNDARY_CENTER 0
#define BOUNDARY_INFINITY 1
#define BOUNDARY_IDEALWALL 2
#define BOUNDARY_INTERFACE 3
#define BOUNDARY_ANTENNA 4

extern "C"
{
intptr_t hmedium_zone_create_ (intptr_t sd_ptr, intptr_t bp_ptr, intptr_t wd_handle, const char *path, int index_p);
intptr_t imhd_zone_create_ (intptr_t sd_ptr, intptr_t bp_ptr, intptr_t wd_handle, const char *path, int index_p);
intptr_t flre_zone_create_ (intptr_t sd_ptr, intptr_t bp_ptr, intptr_t wd_handle, const char *path, int index_p);

void zone_read_settings_ (intptr_t handle, const char *file);
void zone_print_settings_ (intptr_t handle);
void zone_calc_basis_fields_ (intptr_t handle, int flag);
void zone_copy_E_and_B_fields_ (intptr_t handle, double *EB_out);
void zone_calc_final_fields_ (intptr_t handle);
void zone_calc_dispersion_ (intptr_t handle);
void zone_save_dispersion_ (intptr_t handle);
void zone_calc_all_quants_ (intptr_t handle);
void zone_save_all_quants_ (intptr_t handle);
void zone_eval_diss_power_density_ (intptr_t handle, double x, int ttype, int spec, double *dpd);
void zone_eval_current_density_ (intptr_t handle, double x, int ttype, int spec, int comp, double *J);

double zone_get_r1_ (intptr_t handle);
double zone_get_r2_ (intptr_t handle);
int zone_get_dim_of_basis_ (intptr_t handle);
int zone_get_dim_of_basis_vector_ (intptr_t handle);
int zone_get_radial_grid_dimension_ (intptr_t handle);
int zone_get_code_version_ (intptr_t handle);
int zone_get_medium_ (intptr_t handle);
int zone_get_bc1_ (intptr_t handle);
int zone_get_bc2_ (intptr_t handle);
double *zone_get_basis_at_left_boundary_ (intptr_t handle);
double *zone_get_basis_at_right_boundary_ (intptr_t handle);

void zone_save_basis_fields_ (intptr_t handle, const char *path2linear);
void zone_calc_superposition_of_basis_fields_ (intptr_t handle, double *S_p);
void zone_save_final_fields_ (intptr_t handle, const char *path2linear);
void zone_copy_radial_grid_ (intptr_t handle, double *r_p);
void zone_destroy_ (intptr_t handle);

void get_right_boundary_of_zone_ (intptr_t *ptr, double *r);
void get_left_boundary_of_zone_ (intptr_t *ptr, double *r);
}

#endif
