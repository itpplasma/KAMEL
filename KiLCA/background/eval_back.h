/*! \file eval_back.h
    \brief Functions used to evaluate the splined background profiles and parameters.
*/

#ifndef EVAL_BACK_INCLUDE

#define EVAL_BACK_INCLUDE

#include "background.h"

double q (double r, const background *bp);

void eval_hthz (double r, int omin, int omax, const background *bp, double *R);

void eval_Bt_Jt_Jz (double r, const background *bp, double *R);

void eval_Bt_dBz_dpress (double r, const background *bp, double *R);

void eval_B0_ht_hz_n0_Vz (double r, const background *bp, double *R);

void eval_Bt_dBt_Bz_dBz (double r, const background *bp, double *R);

void eval_p_dp (double r, const background *bp, double *R);

void eval_mass_density (double r, const background *bp, double *R);

void eval_dmass_density (double r, const background *bp, double *R);

void eval_Bt_Bz (double r, const background *bp, double *R);

void eval_Bt (double r, const background *bp, double *R);

void eval_Bz (double r, const background *bp, double *R);

void eval_Vt (double r, const background *bp, double *R);

void eval_Vz (double r, const background *bp, double *R);

extern "C"
{
void eval_background_spec_independent_ (double *r, background **bpro, double *R);

void eval_f0_parameters_nu_and_derivs_ (double *r, int *spec, background **bpro, double *R);

void vs_0 (double r, int spec, const background *bp, double *res);

void vs_0_f_ (double *r, int *spec, background **bp_ptr, double *vs);

void eval_hthz_ (double *r, int *omin, int *omax, background **bpro, double *R);
}

#endif
