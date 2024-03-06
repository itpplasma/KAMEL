/*! \file compressible_flow.h
    \brief The declarations of functions implementing ideal compressible MHD equations with plasma flows.
*/

#ifndef IMHD_COMPRESSIBLE_FLOW_INCLUDE

#define IMHD_COMPRESSIBLE_FLOW_INCLUDE

#include "constants.h"
#include "imhd_zone.h"

/*******************************************************************/

struct diff_params
{
    const imhd_zone * zone;  //pointer to imhd object
    int part;               //real or imag
};

/*******************************************************************/

int rhs_flow (double r, const double y[], double dy[], void * params);

int jac_flow (double r, const double y[], double *dfdy, double dfdt[], void *params);

/*****************************RWM  Bondenson et al.(1987)**************************************/

inline complex<double> omega_D (double r, void * params);

inline complex<double> Q_flow (double r, void * params);

inline complex<double> A_flow (double r, void * params);

inline complex<double> T_flow (double r, void * params);

inline complex<double> S_flow (double r, void * params);

inline complex<double> C11_flow (double r, void * params);

inline complex<double> C12_flow (double r, void * params);

inline complex<double> C21_flow (double r, void * params);

inline complex<double> C22_flow (double r, void * params);

inline complex<double> C4_flow (double r, void * params);

inline double func_deriv_in_C4 (double r, void * params);

/*******************************************************************/

extern "C"
{
}

/*******************************************************************/

#endif
