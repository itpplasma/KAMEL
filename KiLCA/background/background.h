/*! \file background.h
    \brief C entry points for the background equilibrium profiles, now owned
           by the Fortran kilca_background_data_m module. The former C++
           background class has been translated away; `background` is kept
           as a forward-declared, never-defined marker type purely so the
           still-C++ zone hierarchy (zone.h/mode.h/core.h/flre_zone.h/
           imhd_zone.h/hmedium_zone.h and the IMHD zone implementation
           files) can keep compiling unchanged - they only ever pass
           `const background *bp` around as an opaque token to the eval_*
           functions below (confirmed by an exhaustive grep: none of them
           dereference it), never allocate or inspect one directly. The
           Fortran side ignores the pointer value entirely, since
           background is a singleton (the only live `new background(...)`
           call site anywhere was core.cpp, called exactly once per run).
*/

#ifndef BACKGROUND_INCLUDE

#define BACKGROUND_INCLUDE

#include <inttypes.h>

class background;

extern "C"
{
intptr_t background_create_ (const char *path2project_p);

void background_set_profiles_from_files_ (void);

void background_set_profiles_from_interface_ (void);

void interp_basic_background_profiles_in_lab_frame_ (int dim, double *r,
                                                      double *q, double *n,
                                                      double *Ti, double *Te,
                                                      double *Vth, double *Vz, double *dPhi0);

double get_background_x0_ (void);

double get_background_xlast_ (void);

int get_background_dimx_ (void);

void get_background_magnetic_fields_ (double r, double *Bt, double *Bz, double *B0);

void get_background_collision_freqs_ (double r, double *nui, double *nue);
}

#endif
