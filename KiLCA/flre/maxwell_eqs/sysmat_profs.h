/*! \file sysmat_profs.h
    \brief C entry points for the ODE system matrix profiles, now owned by the
           Fortran kilca_sysmat_profiles_m module (each instance addressed by
           an opaque handle, since each flre_zone owns its own). The former
           C++ sysmat_profiles class has been translated away.
*/

#ifndef SYSMAT_PROFS_INCLUDE

#define SYSMAT_PROFS_INCLUDE

#include <cstdint>

extern "C"
{
intptr_t sysmat_profiles_create_ (int Nwaves, const char *flag_back_p, const char *path2linear_p,
                                  int NC, int max_dim, double eps_out, int flag_debug,
                                  double r1, double r2, double rm);

void sysmat_profiles_destroy_ (intptr_t handle);

intptr_t get_sysmat_sidm_ (intptr_t handle);

int get_sysmat_dimm_ (intptr_t handle);

int get_sysmat_dimx_ (intptr_t handle);

double * get_sysmat_x_ptr_ (intptr_t handle);

char get_sysmat_flag_back_ (intptr_t handle);

void calc_diff_sys_matrix_ (double *r, char *flag_back, double *R, int flag_back_len);
}

#endif
