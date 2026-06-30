/*! \file disp_profs.h
    \brief C entry points for dispersion-relation profiles, now owned by the
           Fortran kilca_disp_profiles_m module (each instance addressed by an
           opaque handle, since each flre_zone owns its own). The former C++
           disp_profiles class has been translated away.
*/

#ifndef DISP_PROFS_INCLUDE

#define DISP_PROFS_INCLUDE

#include <cstdint>

extern "C"
{
intptr_t disp_profiles_create_ (int Nw, int dimx_p, double *x_p, const char *flag_back_p);

void disp_profiles_destroy_ (intptr_t handle);

void disp_profiles_calculate_ (intptr_t handle);

void disp_profiles_save_ (intptr_t handle, const char *filename);
}

#endif
