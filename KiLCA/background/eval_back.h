/*! \file eval_back.h
    \brief C entry points for evaluating the splined background profiles and
           parameters, now owned by the Fortran kilca_background_data_m
           module. bp is accepted by every function below (matching the
           original signatures exactly, so callers need no changes) but
           ignored: background is a singleton (see background.h).
*/

#ifndef EVAL_BACK_INCLUDE

#define EVAL_BACK_INCLUDE

#include "background.h"

extern "C"
{
double q (double r, const background *bp);

void eval_hthz (double r, int omin, int omax, const background *bp, double *R);

void eval_Bt_dBt_Bz_dBz (double r, const background *bp, double *R);

void eval_p_dp (double r, const background *bp, double *R);

void eval_mass_density (double r, const background *bp, double *R);

void eval_Bt_Bz (double r, const background *bp, double *R);

void eval_Bt (double r, const background *bp, double *R);

void eval_Bz (double r, const background *bp, double *R);

void eval_Vt (double r, const background *bp, double *R);

void eval_Vz (double r, const background *bp, double *R);

void eval_B0_ht_hz_n0_Vz (double r, const background *bp, double *R);
}

#endif
