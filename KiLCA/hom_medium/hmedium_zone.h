/*! \file hmedium_zone.h
    \brief The declaration of hmedium_zone class.
*/

#ifndef HOM_MEDIUM_ZONE

#define HOM_MEDIUM_ZONE

#include <cstring>
#include <complex>
#include <stdlib.h>
#include <stdio.h>

using namespace :: std;

#include "zone.h"

/*****************************************************************************/

/*! \class hmedium_zone
    \brief The class represents properties and data of a plasma zone described by ideal MHD equations.
*/
class hmedium_zone : public zone
{
public:
    complex<double> sigma; //!<homogenius medium conductivity

    inline int iF(int comp)
    {
        //indexing function for E, B fields in a state vector
        return comp;
    }

    inline int iF(int node, int comp, int part)
    {
        //indexing function for E, B fields in EB state array (real part)
        return part + 2*(iF(comp) + Ncomps*(node));
    }

public:
    hmedium_zone (const settings *sd_p, const background *bp_p, const wave_data *wd_p,
                  char *path_p, int index_p) : zone (sd_p, bp_p, wd_p, path_p, index_p) {}

    ~hmedium_zone (void) {}

    void read_settings (char * file);

    void print_settings (void);

    void calc_basis_fields (int flag = 0);

    void copy_E_and_B_fields (double *);

    void calc_final_fields (void) {}

    void calc_dispersion (void) {}

    void save_dispersion (void) {}

    void calc_all_quants (void) {}

    void save_all_quants (void) {}

    void eval_diss_power_density (double x, int type, int spec, double * dpd)
    {
        fprintf (stdout, "\nerror: eval_diss_power_density() is not implemented for the hmedium_zone");
    }

    void eval_current_density (double x, int type, int spec, int comp, double * J)
    {
        J[0] = 0.0;
        J[1] = 0.0;
    }
};

/*****************************************************************************/

extern "C"
{
void eval_basis_in_hom_media_ (int *, double *, double *, double *, double *, double *);

void center_equations_hommed_ (int *Nw, int *len, hmedium_zone **code, double *EB,
			       int *neq, int *nvar, double *M, double *J);

void infinity_equations_hommed_ (int *Nw, int *len, hmedium_zone **code, double *EB,
				 int *neq, int *nvar, double *M, double *J);

void ideal_wall_equations_hommed_ (int *Nw, int *len, hmedium_zone **code, double *EB,
				   int *neq, int *nvar, double *M, double *J);

void stitching_equations_hommed_hommed_ (int *Nw1, int *len1, hmedium_zone **code1, double *EB1,
					 int *Nw2, int *len2, hmedium_zone **code2, double *EB2,
					 int *flg_ant, int *neq, int *nvar, double *M, double *J);

void calc_coeffs_vacuum_antenna_vacuum_ (int *, double *, double *, double *, double *);
}

/*****************************************************************************/

#endif
