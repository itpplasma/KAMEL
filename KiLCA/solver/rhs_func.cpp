/*! \file
    \brief The definition of a right hand side function for ODE solver.
*/

#include "rhs_func.h"
#include "lapack.h"

/*-----------------------------------------------------------------*/

void rhs_func(double r, double* y, double* ydot, void* params) {
    rhs_func_params* fp = (rhs_func_params*)params;

    /* Dmat is a space of length = 2*Nw*Nw (for rhs matrix) + 2*Nw*Nfs (for yt)+ 2*Nw*Nfs (for dyt)
     */

#if USE_SPLINES_IN_RHS_EVALUATION == 1

    eval_diff_sys_matrix_(fp->sp, &r, fp->Dmat); // by spline

#else

    calc_diff_sys_matrix_(&r, fp->sp->flag_back, fp->Dmat, 1); // exact

#endif

    double alpha[2] = {1.0, 0.0}, beta[2] = {0.0, 0.0};
    char trans = 'N';

    int Nw = fp->Nwaves, Nfs = fp->Nfs;

    /* multiply Dmat on matrix y to get matrix ydot: ydot = matmul(Dmat, y) */
    zgemm_(&trans, &trans, &Nw, &Nfs, &Nw, alpha, fp->Dmat, &Nw, y, &Nw, beta, ydot, &Nw);
}

/*-----------------------------------------------------------------*/

int Jacobian(long int N, realtype t, N_Vector y, N_Vector fy, DlsMat J, void* user_data,
    N_Vector tmp1, N_Vector tmp2, N_Vector tmp3) {
    rhs_func_params* fp = (rhs_func_params*)user_data;

    int Nw = fp->Nwaves;

#if USE_SPLINES_IN_RHS_EVALUATION == 1

    eval_diff_sys_matrix_(fp->sp, &t, fp->Dmat); // by spline

#else

    calc_diff_sys_matrix_(&t, fp->sp->flag_back, fp->Dmat, 1); // exact

#endif

    double *col0, *col1;

    int k, j, i;

    // Jacobian of the real system obtained from the complex one: u' = Du + f
    // The (i, j)-th element of J is referenced by colj[i].

    for (k = 0; k < fp->Nfs; k++) // over fundamental solutions
    {
        for (j = 0; j < Nw; j++) // over columns of a complex matrix
        {
            col0 = DENSE_COL(J, 2 * Nw * k + 2 * j + 0);
            col1 = DENSE_COL(J, 2 * Nw * k + 2 * j + 1);

            for (i = 0; i < Nw; i++) // //over rows of a complex matrix
            {
                // Df[2i, 2j]    =  real(Dmat[i,j])
                col0[2 * Nw * k + 2 * i + 0] = fp->Dmat[2 * i + 0 + 2 * Nw * j];

                // Df[2i, 2j+1]  =- imag(Dmat[i,j])
                col1[2 * Nw * k + 2 * i + 0] = -fp->Dmat[2 * i + 1 + 2 * Nw * j];

                // Df[2i+1, 2j]  =  imag(Dmat[i,j])
                col0[2 * Nw * k + 2 * i + 1] = fp->Dmat[2 * i + 1 + 2 * Nw * j];

                // Df[2i+1, 2j+1]=  real(Dmat[i,j])
                col1[2 * Nw * k + 2 * i + 1] = fp->Dmat[2 * i + 0 + 2 * Nw * j];
            }
        }
    }

    return 0;
}

/*-----------------------------------------------------------------*/
