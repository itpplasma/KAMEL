/*! \file
    \brief The declarations of Fortran functions used to build and solve linear system of stitching equations (at zones boundaries).
*/

#ifndef STITCHING_INCLUDE

#define STITCHING_INCLUDE

/*****************************************************************************/

extern "C"
{
void update_system_matrix_and_rhs_vector_ (int *Nc, double *A, double *B, int *neq, int *ieq, int *nvar, int *ivar, double *M, double *J);

void calc_system_determinant_ (int *Nc, double *A, double *det);

void find_superposition_coeffs_ (int *Nc, double *A, double *B, double *S);

void eval_superposition_of_basis_functions_ (int *D, int *Nw, int *dim, double *basis, double *S, double *EB);
}

/*****************************************************************************/

#endif
