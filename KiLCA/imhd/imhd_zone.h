/*! \file imhd_zone.h
    \brief The declaration of imhd_zone class.
*/

#ifndef IMHD_ZONE

#define IMHD_ZONE

#include <stdlib.h>
#include <stdio.h>

#include "zone.h"
#include "hmedium_zone.h"
#include "settings.h"
#include "background.h"
#include "wave_data.h"

/*****************************************************************************/

/*! \class imhd_zone
    \brief The class represents properties and data of a plasma zone described by ideal MHD equations.
*/
class imhd_zone : public zone
{
public:
    inline int iF(int comp)
    {
        //indexing function for Er, Et, Ez, Br, Bt, Bz fields in a state vector
        return comp;
    }

    inline int iF(int node, int comp, int part)
    {
        //indexing function for E, B fields in EB state array
        return part + 2*(iF(comp) + Ncomps*(node));
    }

public:
    imhd_zone (const settings *sd_p, const background *bp_p, const wave_data *wd_p,
               char *path_p, int index_p) : zone (sd_p, bp_p, wd_p, path_p, index_p) {}

    ~imhd_zone (void) {}

    void read_settings (char * file);

    void print_settings (void);

    void calc_basis_fields (int flag = 0);

    void copy_E_and_B_fields (double *);

    void state_to_EB_incompressible (double, double *, double *);

    void calculate_basis_incompressible_gsl (void);

    void calculate_basis_flow_gsl (void);

    void state_to_EB_compressible_flow (double, double *, double *);

    void calc_final_fields (void) {}

    void calc_dispersion (void) {}

    void save_dispersion (void) {}

    void calc_all_quants (void) {}

    void save_all_quants (void) {}

    void eval_diss_power_density (double x, int type, int spec, double * dpd)
    {
        fprintf (stdout, "\nerror: eval_diss_power_density() is not implemented for the imhd_zone");
    }

    void eval_current_density (double x, int type, int spec, int comp, double * J)
    {
        fprintf (stdout, "\nerror: eval_current_density() is not implemented for the imhd_zone");
    }
};

/*****************************************************************************/

double Ffunc (double r, void *params);
double Gfunc (double r, void *params);

complex<double> calc_k_vals (const imhd_zone *zone, double r, double *kt, double *kz, double *ks, double *kp, double *k2, double *kB, complex<double> *kA);

/*****************************************************************************/

extern "C"
{
void set_imhd_data_module_ (imhd_zone **);

void calc_k_vals_ (imhd_zone **zone_ptr, double *r, double *kt, double *kz,
                   double *ks, double *kp, double *k2, double *kB,
                   double *re_kA, double *im_kA, double *re_kfac, double *im_kfac);

void center_equations_imhd_ (int *Nw, int *len, imhd_zone **code, double *EB, int *neq, int *nvar, double *M, double *J);

void infinity_equations_imhd_ (int *Nw, int *len, imhd_zone **code, double *EB, int *neq, int *nvar, double *M, double *J);

void ideal_wall_equations_imhd_ (int *Nw, int *len, imhd_zone **code, double *EB, int *neq, int *nvar, double *M, double *J);

void stitching_equations_imhd_hommed_ (int *Nw1, int *len1, imhd_zone **code1, double *EB1,
                                       int *Nw2, int *len2, hmedium_zone **code2, double *EB2,
                                       int *flg_ant, int *neq, int *nvar, double *M, double *J);

void stitching_equations_imhd_imhd_ (int *Nw1, int *len1, imhd_zone **code1, double *EB1,
                                     int *Nw2, int *len2, imhd_zone **code2, double *EB2,
                                     int *flg_ant, int *neq, int *nvar, double *M, double *J);

void normalize_imhd_basis_ (int * Ncomp, int * Nwaves, int * dim, double * basis, int * ind, int * node);
}

/*****************************************************************************/

#endif
