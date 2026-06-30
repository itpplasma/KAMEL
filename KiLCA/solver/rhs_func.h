/*! \file
    \brief C entry point for the right hand side function of the ODE solver,
           now owned by the Fortran kilca_solver_m module. The former
           rhs_func.cpp implementation has been translated away.
*/

#include <cstdint>

struct rhs_func_params {
    const int Nwaves;
    const int Nphys;
    const int Nfs;
    double* Dmat;
    const intptr_t sp;
};

extern "C"
{
void rhs_func(double, double*, double*, void*);
}
