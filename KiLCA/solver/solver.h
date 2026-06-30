/*! \file
    \brief C entry points for the ODE solver implementing the FLRE
           basis-vector integration u' = A(r)*u(r), now owned by the Fortran
           kilca_solver_m module. The former solver.cpp implementation has
           been translated away.
*/

struct solver_settings
{
    int Nort;          //max number of orthonormalization steps (ONS) for the solver
    double eps_rel;    //relative accuracy for the solver
    double eps_abs;    //absolute accuracy for the solver
    double norm_fac;   //controlling factor for ONS by QR: norm_max/norm_min > norm_fac
    int debug;         //debug flag
};

typedef void (*SysRHSFcn)(double, double *, double *, void *);

extern "C"
{
int integrate_basis_vecs (SysRHSFcn f, int Nfs, int Nw, int dim, double *rvec, double *Smat, solver_settings *ss, void *params);
}
