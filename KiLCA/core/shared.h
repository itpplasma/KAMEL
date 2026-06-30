/*! \file shared.h
    \brief C entry points for the live subset of core/shared.{h,cpp}'s
           helpers, now owned by the Fortran kilca_shared_m module
           (KiLCA/core/shared_m.f90). Trimmed to the two functions still
           called from not-yet-translated C++ (progs/main_eig_param.cpp,
           math/adapt_grid/adaptive_grid.cpp) - both kept extern "C" so
           those call sites need no changes.
*/

#ifndef SHARED_INCLUDE

#define SHARED_INCLUDE

#include <stddef.h>

extern "C"
{
int signum (double x);

void sort_index_doubles (size_t *perm, const double *x, size_t n);
}

#endif
