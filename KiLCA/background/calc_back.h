/*! \file calc_back.h
    \brief Some functions declarations related to the background data.
*/

#ifndef CALC_BACK_INCLUDE

#define CALC_BACK_INCLUDE

#include "background.h"

/*-----------------------------------------------------------------*/

extern "C"
{
/*from fortran*/
void dens_par_ (double *, double *);
void dens_mom_ (double *);
void vth_mom_ (double *);
void vz_mom_ (double *);
void vs_mom_ (double *);
void vp_mom_ (double *);
void eterm_mom_ (double *);

void eval_and_set_background_parameters_spec_independent_ (double *r, char *flag_back, int len);

void eval_and_set_background_parameters_spec_dependent_ (double *r, int *spec, char *flag_back, int len);

void eval_and_set_f0_parameters_nu_and_derivs_ (double *r, int *spec, char *flag_back, int len);
}

/*-----------------------------------------------------------------*/

#endif
