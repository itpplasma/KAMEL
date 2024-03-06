/*! \file eval_cond.h
    \brief The declarations of functions used to evaluate conductivity matrices.
*/

#ifndef EVAL_COND_INCLUDE

#define EVAL_COND_INCLUDE

#include "cond_profs.h"

void eval_K_matrices (const cond_profiles *cp, int spec, int type, int Dmin, int Dmax, double r, double *K);

void eval_all_K_matrices (const cond_profiles *cp, int Dmin, int Dmax, double r, double *K);

void eval_C_matrices (const cond_profiles *cp, int spec, int type, int Dmin, int Dmax, double r, double *C);

void eval_all_C_matrices (const cond_profiles *cp, int Dmin, int Dmax, double r, double *C);

extern "C"
{
void eval_c_matrices_f_ (cond_profiles **cd_ptr, int *spec, int *type, int *Dmin, int *Dmax, double *r, double *ct);

void eval_k_matrices_f_ (cond_profiles **cp_ptr, int *spec, int *type, int *Dmin, int *Dmax, double *r, double *km);
}

#endif
