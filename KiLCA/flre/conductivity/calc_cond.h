/*! \file calc_cond.h
    \brief The declarations of functions used to calculate conductivity matrices and their splines.
*/

#ifndef CALC_COND_INCLUDE

#define CALC_COND_INCLUDE

#include "settings.h"
#include "adaptive_grid.h"
#include "constants.h"
#include "settings.h"
#include "background.h"
#include "cond_profs.h"
#include "spline.h"
#include "conduct_parameters.h"
#include "eval_cond.h"

void sample_cond_func (double *r, double *f, void *p);
void calc_splines_for_K (cond_profiles *cp);
void set_arrays_for_K (cond_profiles *cp);
void smooth_arrays_for_K (cond_profiles *cp);

void sample_cond_func_polynom (double *r, double *f, void *p);
void calc_splines_for_K_polynom (cond_profiles *cp);
void set_arrays_for_K_polynom (cond_profiles *cp);

void save_K_matrices (const cond_profiles *cp, int spec, int type);
void save_K_matrices_fine (const cond_profiles *cp, int spec, int type, int dimf);

void calc_splines_for_C (cond_profiles *cp);
void calc_C_matrices (const cond_profiles *cp, int spec, int type, double r, double *C);
void save_C_matrices (const cond_profiles *cp, int spec, int type);
void save_C_matrices_fine (const cond_profiles *cp, int spec, int type, int dimf);

extern "C"
{
void calc_and_add_galilelian_correction_ (double *r, int *spec, char *flag_back, double *ct, int len);

void eval_electric_drift_velocities_ (void);
void eval_fgi_arrays_ (void);
void eval_a_matrix_ (void);
void calc_dem_djmi_arrays_ (double *);

void calc_w1_array_ (void);
void calc_w2_array_ (int *);
void calc_d_array_ (void);
void calc_k_matrices_ (double *);
void calc_k1_matrices_ (double *);
void calc_dummy_matrices_ (double *);

void transform_c_matrices_to_rsp_ (double *, double *);

void print_conduct_params_ (void);
void print_conduct_arrays_ (void);

void smooth_array_gauss_ (int *, int *, double *);
}

#endif
