/*! \file conduct_parameters.h
    \brief The C declarations of Fortran subroutines used in evaluation of conductivity matrices.
*/

#ifndef CONDUCT_PARAMETERS_INCLUDE

#define CONDUCT_PARAMETERS_INCLUDE

extern "C"
{
void set_back_aliases_in_conductivity_parameters_ ();

void eval_and_set_conduct_parameters_spec_independent_ (double *r, char *flag_back, int len);
void eval_and_set_conduct_parameters_spec_dependent_ (double *r, int *spec, char *flag_back, int len);
void eval_and_set_all_conduct_parameters_ (double *r, int *spec, char *flag_back, int len);

void eval_and_set_background_parameters_spec_independent_ (double *r, char *flag_back, int len);
void eval_and_set_wave_parameters_ (double *r, char *flag_back, int len);

void eval_and_set_background_parameters_spec_dependent_ (double *r, int *spec, char *flag_back, int len);
void eval_and_set_f0_parameters_nu_and_derivs_ (double *r, int *spec, char *flag_back, int len);

void print_conduct_params_ (void);
void print_conduct_arrays_ (void);

void eval_z_and_oms_ (void);
void eval_w_array_ (void);

void get_wave_parameters_ (double *kvals);

//from conductivity_arrays:
void allocate_and_set_conductivity_arrays_ (void);
void eval_electric_drift_velocities_ (void);
void eval_FGI_arrays_ (void) ;
void eval_A_matrix_ (void);
}

#endif
