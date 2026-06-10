/*! \file incompressible.h
    \brief The declarations of functions implementing incompressible flowless MHD equations.
*/

#ifndef IMHD_INCOMPRESSIBLE_INCLUDE

#define IMHD_INCOMPRESSIBLE_INCLUDE

#include <ctime>
#include "constants.h"
#include "imhd_zone.h"

/*******************************************************************/

struct diff_params
{
    const imhd_zone *zone;  //pointer to imhd object
    int part;               //real or imag
};

/*******************************************************************/

int rhs_incompressible (double r, const double y[], double f[], void *params);

inline complex<double> Afunc_hi_inc (double r, void *params);

inline double Afunc_hi_inc_part(double r, void *params);

inline complex<double> Bfunc_hi_inc (double r, void *params);

inline double func_der (double r, void *params);

/*******************************************************************/

extern "C"
{
}

/*******************************************************************/

#endif
