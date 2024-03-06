/*! \file
    \brief The declarations of functions implementing ODE solver for stiff set of linear equations u' = A(t)*u(t).
*/

#include <stdio.h>
#include <limits.h>
#include <math.h>
#include <stdlib.h>

#include "shared.h"
#include "lapack.h"

#include <nvector/nvector_serial.h>       /* serial N_Vector types, fct. and macros */

struct solver_settings
{
    int Nort;          //max number of orthonormalization steps (ONS) for the solver
    double eps_rel;    //relative accuracy for the solver
    double eps_abs;    //absolute accuracy for the solver
    double norm_fac;   //controlling factor for ONS by QR: norm_max/norm_min > norm_fac
    int debug;         //debug flag
};

typedef void (*SysRHSFcn)(double, double *, double *, void *);

int func (realtype t, N_Vector y, N_Vector ydot, void *params);

int integrate_basis_vecs_ (SysRHSFcn f, int Nfs, int Nw, int dim, double *rvec, double *Smat, int *Nort, double *mem, void *params);

int integrate_basis_vecs (SysRHSFcn f, int Nfs, int Nw, int dim, double *rvec, double *Smat, solver_settings *ss, void *params);

int renorm_basis_vecs_ (int Nfs, int Nw, int dim, double *rvec, double *Smat, int Nort, double *rdata, double *ydata, double *taudata);

int superpose_basis_vecs_ (int Nfs, int Nw, int nsteps, double *Smat, double *Cvec, double *Svec);
