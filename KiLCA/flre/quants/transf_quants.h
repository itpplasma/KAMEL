/*! \file transf_quants.h
    \brief The declarations of functions used to perform Galilean transforms of quantities to and from a moving frame.
*/

#ifndef TRANSF_QUANTS_INCLUDE

#define TRANSF_QUANTS_INCLUDE

#include "quants_profs.h"
#include "field_profs.h"
#include "shared.h"

/*******************************************************************/

void transform_quants_to_lab_sys (const quants_profiles *qp);

void eval_current_dens_in_lab_frame (const quants_profiles *qp, double *cd);

void eval_current_dens_in_cyl_sys (const quants_profiles *qp, double *cd);

void save_current_density (const quants_profiles *qp, double *cd, const char *frame, const char *comp);

void eval_and_save_EB_fields_in_lab_frame (const quants_profiles *qp, double *field);

void save_EB_fields (const quants_profiles *qp, const double *fields, const char *frame, const char *comp);

void transform_EB_fields_to_mov_frame (const background *bp, int dim, double *x, double *field);

void transform_start_values_to_mov_frame (const field_profiles *fp, double r, double *Smat);

/*******************************************************************/

#endif
