/*! \file flre_zone.h
    \brief The declaration of flre_zone class.
*/

#ifndef FLRE_ZONE

#define FLRE_ZONE

#include "code_settings.h"

#include "zone.h"
#include "hmedium_zone.h"
#include "maxwell_eqs_data.h"
#include "cond_profs.h"
#include "sysmat_profs.h"
#include "disp_profs.h"
#include "flre_quants.h"

/*****************************************************************************/

/*! \class flre_zone
    \brief The class represents properties and data of a plasma zone described by FLRE.
*/
class flre_zone : public zone
{
public:
    maxwell_eqs_data *me;  //!<pointer to maxwell equations data
    cond_profiles *cp;     //!<pointer to conductivity structure
    sysmat_profiles *sp;   //!<pointer to system matrix profiles
    disp_profiles *dp;     //!<dispersion profiles
    flre_quants *qp;       //!<pointer to flre quantities

    //Conductivity settings:
    int flre_order; //!<order of FLR expansion
    int Nmax;       //!<highest cyclotron harmonic
    int gal_corr;   //!<flag if use correction term in conductivity: for flre_order=1 only!
    int N;          //!<splines degree: >= 2N+1, where N - order of flr expansion
    int max_dim_c;  //!<maximum dimension of the radial grid for conductivity martices
    double D;       //!<width of resonant layer (RL)
    double eps_out; //!<error parameter used in the adaptive radial grid generation outside RL
    double eps_res; //!<error parameter used in the adaptive radial grid generation in RL

    //ME output grid settings:
    double dr_out; //!<grid step outside the resonance region for the solutions
    double dr_res; //!<grid step in the resonance region for the solutions
    double del;    //!<width of the resonant layer

    //additional Maxwell equations settings with flre current density:
    int hom_sys;    //!<flag if system of equations should be used in homogenious limit

    //additional settings for ODE solver:
    int Nort;        //!<max number of orthonormalization steps (ONS) for the solver
    double norm_fac; //!<controlling factor for ONS by QR: norm_max/norm_min > norm_fac

    int Nfs;         //!<number of fundamental solutions for integration

    double *EB_mov;  //!<system vector in a moving frame

    int collmod[2];  //!<collisions model flag for ions and electrons

    int rsp;         //!<use cylindrical components for conductivity if rsp = 0 and rsp components otherwise

    //indexing function for state vectors (to be used after integration)
    inline int is(int node, int sol, int comp, int part)
    {
        //indexing function for state vectors array
        //must correspond to the following fortran definition:
        //complex(8), dimension(Nwaves,Nfs,dim) :: state;

        //node = 0:dim-1, sol  = 0:Nfs-1, comp = 0:Nwaves-1, part = 0:1 (real, imag)
        return part + 2*(comp + Nwaves*(sol + Nfs*(node)));
    }

    //warning: EB array contains not a system but state vectors up to calc_final_fields ()
    inline int iF(int comp)
    {
        //indexing function for Er, Es, Ep, Br, Bs, Bp fields in a system EB vector
        if (comp < 3) return me->iErsp_sys[comp];
        else          return me->iBrsp_sys[comp-3];
    }

    inline int iF(int node, int comp, int part)
    {
        //indexing function for E, B fields in EB system array
        return part + 2*(iF(comp) + Ncomps*(node));
    }

public:
    flre_zone (const settings *sd_p, const background *bp_p, const wave_data *wd_p,
               char *path_p, int index_p) : zone (sd_p, bp_p, wd_p, path_p, index_p)
    {
        me = NULL;
        cp = NULL;
        sp = NULL;
        dp = NULL;
        qp = NULL;

        EB_mov = NULL;
    }

    ~flre_zone (void)
    {
        if (me) delete me;
        if (cp) delete cp;
        if (sp) delete sp;
        if (dp) delete dp;
        if (qp) delete qp;

        delete [] EB_mov;
    }

    void read_settings (char * file);

    void print_settings (void);

    void calc_basis_fields (int flag = 0);

    void copy_E_and_B_fields (double *);

    void calc_final_fields (void);

    void calc_dispersion (void);

    void save_dispersion (void);

    void calculate_field_profiles_orth (void);

    void save_system_vector_in_mov_frame (void);

    void calc_all_quants (void);

    void save_all_quants (void);

    void eval_diss_power_density (double x, int type, int spec, double * dpd);

    void eval_current_density (double x, int type, int spec, int comp, double * J);
};

/*****************************************************************************/

inline void galilean_transform_of_flre_state_vector (const flre_zone *zone, double V, complex<double> omega, double r, const double *E1, double *E2);

inline void galilean_transform_of_flre_system_vector (const flre_zone *zone, double V, complex<double> omega, double r, const double *EB1, double *EB2);

inline void system_to_state_copy (const flre_zone *zone, double *system, double *state);

inline void state_to_system_copy (const flre_zone *zone, double *state, double *system);

inline void transform_of_flre_system_vector_to_cyl_coordinates (const flre_zone *zone, double r, const double *EB1, double *EB2);

/*****************************************************************************/

extern "C"
{
void allocate_and_set_conductivity_arrays_ (void);

void deallocate_conductivity_arrays_ (void);

void set_cond_profiles_in_mode_data_module_ (cond_profiles **);

void set_sysmat_profiles_in_mode_data_module_ (sysmat_profiles **);

void set_conductivity_settings_c_ (flre_zone **ptr, int *flre_order, int *Nmax, int *gal_corr, int *rsp);

void set_equations_settings_c_ (flre_zone **ptr, int *hom_sys, int *Nwaves, int *Nfs, int *Nphys, int *flag_debug);

void set_collisions_settings_c_ (flre_zone **ptr, int *collmod);

void setup_flre_data_module_ (flre_zone **);

void clean_flre_data_module_ ();

void calc_and_set_maxwell_system_parameters_module_ (void);

void clean_maxwell_system_parameters_module_ (void);

void center_equations_flre_ (int *Nw, int *len, flre_zone **code, double *EB,
			     int *neq, int *nvar, double *M, double *J);

void infinity_equations_flre_ (int *Nw, int *len, flre_zone **code, double *EB,
			       int *neq, int *nvar, double *M, double *J);

void ideal_wall_equations_flre_ (int *Nw, int *len, flre_zone **code, double *EB,
				 int *neq, int *nvar, double *M, double *J);

void stitching_equations_flre_hommed_ (int *Nw1, int *len1, flre_zone **code1, double *EB1,
				       int *Nw2, int *len2, hmedium_zone **code2, double *EB2,
				       int *flg_ant, int *neq, int *nvar, double *M, double *J);

void stitching_equations_flre_flre_ (int *Nw1, int *len1, flre_zone **code1, double *EB1,
				     int *Nw2, int *len2, flre_zone **code2, double *EB2,
				     int *flg_ant, int *neq, int *nvar, double *M, double *J);

void calc_start_values_anywhere_low_derivs_ (double *, double *);

void calc_start_values_center_with_correct_asymptotic_ (double *, double *);

void calc_start_values_ideal_wall_low_derivs_ (double *, double *);

void state2sys_ (double *, char *, double *, double *, double *, int len);

void get_sys_ind_array_ (flre_zone **ptr, int *sys_ind);

void activate_fortran_modules_for_zone_ (flre_zone **ptr);

void deactivate_fortran_modules_for_zone_ (flre_zone **ptr);

void calc_flre_basis_in_lab_cyl_frame_with_full_system_vectors_ (flre_zone **ptr, double *r, double *EB1, double *EB2);

void get_iersp_sys_array_ (flre_zone **ptr, int *iErsp_sys);

void get_ibrsp_sys_array_ (flre_zone **ptr, int *iBrsp_sys);

void get_flre_order_ (int *flre_order);

void get_gal_corr_ (int *gal_corr);

void normalize_flre_basis_ (int * Ncomp, int * Nwaves, int * dim, double * basis,
			    int * iErsp_sys, int * ind1, int * ind2, int * node);
}

/*****************************************************************************/

#endif
