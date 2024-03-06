/*! \file eval_back.cpp
    \brief The implementation of the functions declared in eval_back.h
*/

#include <stdlib.h>
#include <stdio.h>
#include <cmath>
#include <cstring>
#include <climits>

#include "constants.h"
#include "eval_back.h"
#include "spline.h"

/*-----------------------------------------------------------------*/

void eval_background_spec_independent_ (double *r, background **bpro, double *R)
{
//evaluates B0, hth, hz, dPhi0 and 2 derivatives
background *bp = (background *)(*bpro);
spline_eval_ (bp->sid, 1, r, 0, 2, bp->i_B, bp->i_dPhi0, R);
}

/*-----------------------------------------------------------------*/

void eval_f0_parameters_nu_and_derivs_ (double *r, int *spec, background **bpro, double *R)
{
//evaluates n_p, Vp_p, Vt_p, nu + 2 derivatives
background *bp = (background *)(*bpro);
spline_eval_ (bp->sid, 1, r, 0, 2, bp->i_n_p[*spec], bp->i_nu[*spec], R);
}

/*-----------------------------------------------------------------*/

void vs_0 (double r, int spec, const background *bp, double *res)
{
double R[2];

spline_eval_ (bp->sid, 1, &r, 0, 1, bp->i_n, bp->i_n, R);

double n  = R[0];
double dn = R[1];

spline_eval_ (bp->sid, 1, &r, 0, 1, bp->i_T[spec], bp->i_T[spec], R);

double T  = R[0];
double dT = R[1];

double dpress = boltz*(dn*T+n*dT);

double B;
spline_eval_ (bp->sid, 1, &r, 0, 0, bp->i_B, bp->i_B, &B);

double dPhi0;
spline_eval_ (bp->sid, 1, &r, 0, 0, bp->i_dPhi0, bp->i_dPhi0, &dPhi0);

*res = c/B*(dPhi0+dpress/(bp->sd->bs->charge[spec])/n);
}

/*-----------------------------------------------------------------*/

void vs_0_f_ (double *r, int *spec, background **bp_ptr, double *vs)
{
vs_0 (*r, *spec, (const background *)(*bp_ptr), vs);
}

/*-----------------------------------------------------------------*/

double q (double r, const background *bp)
{
double ans;

spline_eval_ (bp->sid, 1, &r, 0, 0, bp->i_q, bp->i_q, &ans);

return ans;
}

/*-----------------------------------------------------------------*/

void eval_hthz_ (double *r, int *omin, int *omax, background **bpro, double *R)
{
//evaluates hth, hz derivs omin...omax

background *bp = (background *)(*bpro);

spline_eval_ (bp->sid, 1, r, *omin, *omax, bp->i_hth, bp->i_hz, R);
}

/*-----------------------------------------------------------------*/

void eval_hthz (double r, int omin, int omax, const background *bp, double *R)
{
//evaluates hth, hz derivs omin...omax
spline_eval_ (bp->sid, 1, &r, omin, omax, bp->i_hth, bp->i_hz, R);
}

/*-----------------------------------------------------------------*/

void eval_B0_ht_hz_n0_Vz (double r, const background *bp, double *R)
{
spline_eval_ (bp->sid, 1, &r, 0, 0, bp->i_B, bp->i_hz, R+0); //B, hth, hz
spline_eval_ (bp->sid, 1, &r, 0, 0, bp->i_n, bp->i_n, R+3); //n
spline_eval_ (bp->sid, 1, &r, 0, 0, bp->i_Vz[2], bp->i_Vz[2], R+4); //Vz - total
}

/*-----------------------------------------------------------------*/

void eval_Bt_dBz_dpress (double r, const background *bp, double *R)
{
spline_eval_ (bp->sid, 1, &r, 0, 0, bp->i_Bth, bp->i_Bth, R+0); //B0th
spline_eval_ (bp->sid, 1, &r, 1, 1, bp->i_Bz, bp->i_Bz, R+1);   //dBz

double S[4];

spline_eval_ (bp->sid, 1, &r, 0, 1, bp->i_n, bp->i_n, S);

double n  = S[0];
double dn = S[1];

spline_eval_ (bp->sid, 1, &r, 0, 1, bp->i_T[0], bp->i_T[1], S);

double Ti  = S[0], dTi = S[1];
double Te  = S[2], dTe = S[3];

R[2] = boltz*(dn*(Ti+Te)+n*(dTi+dTe)); //dpress
}

/*-----------------------------------------------------------------*/

void eval_Bt_Jt_Jz (double r, const background *bp, double *R)
{
spline_eval_ (bp->sid, 1, &r, 0, 0, bp->i_Bth, bp->i_Bth, R+0);  //B0th
spline_eval_ (bp->sid, 1, &r, 0, 0, bp->i_J0th[2], bp->i_J0z[2], R+1); //Jt and Jz
}

/*-----------------------------------------------------------------*/

void eval_Bt_dBt_Bz_dBz (double r, const background *bp, double *R)
{
spline_eval_ (bp->sid, 1, &r, 0, 1, bp->i_Bth, bp->i_Bz, R);
}

/*-----------------------------------------------------------------*/

void eval_p_dp (double r, const background *bp, double *R)
{
double S[4];

spline_eval_ (bp->sid, 1, &r, 0, 1, bp->i_n, bp->i_n, S);

double n  = S[0];
double dn = S[1];

spline_eval_ (bp->sid, 1, &r, 0, 1, bp->i_T[0], bp->i_T[1], S);

double Ti  = S[0], dTi = S[1];
double Te  = S[2], dTe = S[3];

R[0] = boltz*(n*(Ti+Te)); //pressure
R[1] = boltz*(dn*(Ti+Te)+n*(dTi+dTe)); //pressure'
}

/*-----------------------------------------------------------------*/

void eval_mass_density (double r, const background *bp, double *R)
{
spline_eval_ (bp->sid, 1, &r, 0, 0, bp->i_n, bp->i_n, R);

R[0] *= bp->sd->bs->mass[0];
}

/*-----------------------------------------------------------------*/

void eval_dmass_density (double r, const background *bp, double *R)
{
spline_eval_ (bp->sid, 1, &r, 1, 1, bp->i_n, bp->i_n, R);

R[0] *= bp->sd->bs->mass[0];
}

/*-----------------------------------------------------------------*/

void eval_Bt_Bz (double r, const background *bp, double *R)
{
spline_eval_ (bp->sid, 1, &r, 0, 0, bp->i_Bth, bp->i_Bz, R);
}

/*-----------------------------------------------------------------*/

void eval_Bt (double r, const background *bp, double *R)
{
spline_eval_ (bp->sid, 1, &r, 0, 0, bp->i_Bth, bp->i_Bth, R);
}

/*-----------------------------------------------------------------*/

void eval_Bz (double r, const background *bp, double *R)
{
spline_eval_ (bp->sid, 1, &r, 0, 0, bp->i_Bz, bp->i_Bz, R);
}

/*-----------------------------------------------------------------*/

void eval_Vt (double r, const background *bp, double *R)
{
spline_eval_ (bp->sid, 1, &r, 0, 0, bp->i_Vth[2], bp->i_Vth[2], R);
}

/*-----------------------------------------------------------------*/

void eval_Vz (double r, const background *bp, double *R)
{
spline_eval_ (bp->sid, 1, &r, 0, 0, bp->i_Vz[2], bp->i_Vz[2], R);
}

/*-----------------------------------------------------------------*/
