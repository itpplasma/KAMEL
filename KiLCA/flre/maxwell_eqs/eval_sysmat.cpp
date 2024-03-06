#include <stdlib.h>
#include <stdio.h>
#include <cmath>
#include <cstring>

#include "eval_sysmat.h"

/*******************************************************************/

void eval_diff_sys_matrix_ (const sysmat_profiles *sp, double *r, double *M)
{
spline_eval_d_ (sp->sidM, 1, r, 0, 0, 0, sp->dimM-1, M);
}

/*******************************************************************/

void eval_diff_sys_matrix_f_ (sysmat_profiles **sp_ptr, double *r, double *M)
{
const sysmat_profiles *sp = *sp_ptr;
spline_eval_d_ (sp->sidM, 1, r, 0, 0, 0, sp->dimM-1, M);
}

/*******************************************************************/
