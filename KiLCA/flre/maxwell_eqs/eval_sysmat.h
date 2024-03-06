#ifndef EVAL_SYSMAT_INCLUDE

#define EVAL_SYSMAT_INCLUDE

#include "settings.h"
#include "constants.h"
#include "settings.h"
#include "spline.h"
#include "sysmat_profs.h"

extern "C"
{
void eval_diff_sys_matrix_ (const sysmat_profiles *sp, double *r, double *M);
void eval_diff_sys_matrix_f_ (sysmat_profiles **sp_ptr, double *r, double *M);
}

#endif
