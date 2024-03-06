/*! \file calc_flre_quants.h
    \brief The functions used to calculate auxiliary quantities for FLRE zone.
*/

#ifndef CALC_FLRE_QUANTS_INCLUDE

#define CALC_FLRE_QUANTS_INCLUDE

#include "flre_quants.h"
#include "shared.h"

/*******************************************************************/

extern "C"
{
void current_density_ (double *);
void cyl2rsp_ (const double *, double *, double *, double *, double *);
void calc_current_density_r_s_p_ (double *, double *);

void eval_and_set_background_parameters_spec_independent_ (double *r, char *flag_back);
void eval_and_set_wave_parameters_ (double *r, char *flag_back);
void get_wave_parameters_ (double *kvals);
void get_magnetic_field_parameters_ (double *hvals);
}

void calc_splines_for_current_density (flre_quants *qp);

void calculate_JaE (flre_quants *qp);
void calculate_JaE_delta (flre_quants *qp);
void calculate_JaE_distributed (flre_quants *qp);

void save_total_absorbed_power (flre_quants *qp);

void calc_current_density (flre_quants *qp);
void save_current_density (const flre_quants *qp);

void calc_absorbed_power_density (flre_quants *qp);
void calc_absorbed_power_in_cylinder (flre_quants *qp);
void save_absorbed_power (const flre_quants *qp);

void calc_dissipated_power_density (flre_quants *qp);
void calc_dissipated_power_in_cylinder (flre_quants *qp);
void save_dissipated_power (const flre_quants *qp);

void calc_kinetic_flux (flre_quants *qp);
void save_kinetic_flux (const flre_quants *qp);

void calc_poynting_flux (flre_quants *qp);
void save_poynting_flux (const flre_quants *qp);

void calc_total_flux (flre_quants *qp);
void save_total_flux (const flre_quants *qp);

void calc_number_density (flre_quants *qp);
void save_number_density (const flre_quants *qp);

void calc_lorentz_torque_density (flre_quants *qp);
void calc_lorentz_torque_on_cylinder (flre_quants *qp);
void save_lorentz_torque (const flre_quants *qp);

void calculate_field_profiles_poy_test (flre_quants *qp);

inline int binary_search (double x, const double *xa, int ilo, int ihi);

void vec_product_3D (complex<double> *a, complex<double> *b, complex<double> *res);

void integrate_over_cylinder (int dim, double *x, double *q, double vol_fac, double *qi);

void transform_quants_to_lab_cyl_frame (const flre_quants *qp);

void eval_current_dens_in_lab_frame (const flre_quants *qp, double *cd);

void eval_current_dens_in_cyl_sys (const flre_quants *qp, double *cd);

void save_current_density (const flre_quants *qp, double *cd, const char *frame, const char *comp);

void interp_diss_power_density (const flre_quants *qp, double x, int type, int spec, double * dpd);

void interp_current_density (const flre_quants *qp, double x, int type, int spec, int comp, double * J);

/*******************************************************************/

#endif
