/*! \file eval_cond.cpp
    \brief The implementation of functions declared in eval_cond.h.
*/

#include <stdlib.h>
#include <stdio.h>
#include <cmath>
#include <cstring>

#include "eval_cond.h"
#include "cond_profs.h"
#include "spline.h"

/*******************************************************************/

void eval_K_matrices (const cond_profiles *cp, int spec, int type, int Dmin, int Dmax, double r, double *K)
{
int dimk = 2*9*(cp->flreo+1)*(cp->flreo+1); //number of K's to eval
int ind_ka = cp->iKa(spec, type, 0, 0, 0, 0, 0, 0);

spline_eval_d_ (cp->sidK, 1, &r, Dmin, Dmax, ind_ka, ind_ka+dimk-1, K);

if (cp->flag_back[0] != 'f') //r is multiplied on huge_factor
{
    for (int n=Dmin; n<=Dmax; n++)
    {
        double scale_fac = pow (cp->sd->bs->huge_factor, n);
        int ind = dimk*(n-Dmin);
        for (int j=0; j<dimk; j++) K[j+ind] /= scale_fac;
    }
}
}

/*******************************************************************/

void eval_all_K_matrices (const cond_profiles *cp, int Dmin, int Dmax, double r, double *K)
{
//evaluates all K matrices for different spec, types, etc...

int dimk = cp->dimK; //number of K's to eval

spline_eval_d_ (cp->sidK, 1, &r, Dmin, Dmax, 0, dimk-1, K);

int n, j, ind;
double scale_fac;

if (cp->flag_back[0] != 'f') //r is multiplied on a huge_factor
{
    for (n=Dmin; n<=Dmax; n++)
    {
        scale_fac = pow (cp->sd->bs->huge_factor, n);
        ind = dimk*(n-Dmin);
        for (j=0; j<dimk; j++) K[j+ind] /= scale_fac;
    }
}
}

/*******************************************************************/

void eval_C_matrices (const cond_profiles *cp, int spec, int type, int Dmin, int Dmax, double r, double *C)
{
int dimc = 18*(2*(cp->flreo)+1); //number of C's to eval
int ind_ca = dimc*(type + (cp->dimt)*spec);

spline_eval_d_ (cp->sidC, 1, &r, Dmin, Dmax, ind_ca, ind_ca+dimc-1, C);

int n, j, ind;
double scale_fac;

if (cp->flag_back[0] != 'f') //r is multiplied on huge_factor
{
    for (n=Dmin; n<=Dmax; n++)
    {
        scale_fac = pow (cp->sd->bs->huge_factor, n);
        ind = dimc*(n-Dmin);
        for (j=0; j<dimc; j++) C[j+ind] /= scale_fac;
    }
}
}

/*******************************************************************/

void eval_all_C_matrices (const cond_profiles *cp, int Dmin, int Dmax, double r, double *C)
{
int dimc = cp->dimC; //number of C's to eval

spline_eval_d_ (cp->sidC, 1, &r, Dmin, Dmax, 0, dimc-1, C);

int n, j, ind;
double scale_fac;

if (cp->flag_back[0] != 'f') //r is multiplied on huge_factor
{
    for (n=Dmin; n<=Dmax; n++)
    {
        scale_fac = pow (cp->sd->bs->huge_factor, n);
        ind = dimc*(n-Dmin);
        for (j=0; j<dimc; j++) C[j+ind] /= scale_fac;
    }
}
}

/*******************************************************************/

void eval_c_matrices_f_ (cond_profiles **cp_ptr, int *spec, int *type, int *Dmin, int *Dmax, double *r, double *ct)
{
const cond_profiles *cp = *cp_ptr;
eval_C_matrices (cp, *spec, *type, *Dmin, *Dmax, *r, ct);
}

/*******************************************************************/

void eval_k_matrices_f_ (cond_profiles **cp_ptr, int *spec, int *type, int *Dmin, int *Dmax, double *r, double *km)
{
const cond_profiles *cp = *cp_ptr;
eval_K_matrices (cp, *spec, *type, *Dmin, *Dmax, *r, km);
}

/*******************************************************************/
