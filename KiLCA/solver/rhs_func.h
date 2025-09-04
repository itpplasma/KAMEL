/*! \file
    \brief The declaration of a right hand side function for ODE solver.
*/

#include "sysmat_profs.h"

#include <cvode/cvode.h>
#include <nvector/nvector_serial.h>  /* serial N_Vector types, fct. and macros */
#include <sundials/sundials_dense.h> /* use generic DENSE solver in preconditioning */
#include <sundials/sundials_types.h> /* definition of sunrealtype */

struct rhs_func_params {
    const int Nwaves;
    const int Nphys;
    const int Nfs;
    double* Dmat;
    const sysmat_profiles* sp;
};

/*-----------------------------------------------------------------*/

void rhs_func(double, double*, double*, void*);

void rhs_func_coeff(double, double*, double*, void*);

/*-----------------------------------------------------------------*/

int Jacobian(long int N, sunrealtype t, N_Vector y, N_Vector fy, SUNDlsMat Jac, void* user_data,
    N_Vector tmp1, N_Vector tmp2, N_Vector tmp3);

/*-----------------------------------------------------------------*/
