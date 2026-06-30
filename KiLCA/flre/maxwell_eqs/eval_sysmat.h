#ifndef EVAL_SYSMAT_INCLUDE

#define EVAL_SYSMAT_INCLUDE

#include <cstdint>

#include "settings.h"
#include "constants.h"
#include "spline.h"
#include "sysmat_profs.h"

extern "C"
{
void eval_diff_sys_matrix_ (intptr_t sp, double *r, double *M);
}

#endif
