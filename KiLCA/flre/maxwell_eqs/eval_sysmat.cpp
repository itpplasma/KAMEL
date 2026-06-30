#include <stdlib.h>
#include <stdio.h>
#include <cmath>
#include <cstring>

#include "eval_sysmat.h"

/*******************************************************************/

void eval_diff_sys_matrix_ (intptr_t sp, double *r, double *M)
{
spline_eval_d_ ((uintptr_t) get_sysmat_sidm_ (sp), 1, r, 0, 0, 0, get_sysmat_dimm_ (sp)-1, M);
}

/*******************************************************************/
