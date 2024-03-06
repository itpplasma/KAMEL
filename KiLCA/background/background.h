/*! \file background.h
    \brief The declaration of background class representing background profiles and other data.
*/

#ifndef BACKGROUND_INCLUDE

#define BACKGROUND_INCLUDE

#include <inttypes.h>

#include "settings.h"

/*! \class background
    \brief The class for background profiles and data.
*/
class background
{
public:
    const settings *sd; //!<pointer to settings data

    ///below are local data calculated and stored in this structure:
    ///profiles:
    int Nprofiles;          //!<number of basic (input) background profiles
    char **profile_names;   //!<names of basic background profiles
    char *path2background;  //!<path to output background dir

    ///data:
    int dimx;   //!<dimension of a profile x grid (all the same)
    double *x;  //!<x grid for profiles

    int dimy;   //!<number of all background profiles
    double *y;  //!<y grids: all values are stored in one-dimensional array

    ///splines:
    int ind;    //!<search index
    double *C;  //!<matrix of all splines coefficients
    double *R;  //!<array of all spline values (+ derivs) at some point

    uintptr_t sid; //!<spline_data id

    ///Indices: used to get acces to a particular value in (y, C, R) arrays
    ///basic profiles:
    int i_q;
    int i_n;
    int *i_T;   /*i,e*/
    int *i_Vth; /*i,e,t*/
    int *i_Vz;  /*i,e,t*/
    int i_Er;

    /*from equlibrium:*/
    int i_hth;
    int i_hz;
    int i_Bth;
    int i_Bz;
    int i_B;

    int *i_J0th; /*i,e,t*/
    int *i_J0z;  /*i,e,t*/

    /*from f0 parameters search*/
    int i_dPhi0;
    int *i_n_p;  /*parameter of f0 for i,e*/
    int *i_Vp_p; /*parameter of f0 for i,e*/
    int *i_Vt_p; /*parameter of f0 for i,e*/

    /*additional:*/
    int *i_nu; /*i,e*/

    int flag_dPhi0_calc; //flag to recalculate dPhi0 or not

public:
    background (const settings *s);

    ~background (void);

    void clean_background_arrays (void);

    void set_profiles_indices (void);

    void set_background_profiles_from_files (void);

    void set_background_profiles_from_interface (void);

    int load_input_background_profiles (void);

    int spline_input_profiles (void);

    int find_f0_parameters (void);

    int calculate_equilibrium (void);

    int check_and_spline_equilibrium (void);

    int eval_and_save_f0_moments (void);

    int save_background (void);

    inline void interp_basic_background_profiles (double r, double *q, double *n,
                                                  double *Ti, double *Te,
                                                  double *Vth, double *Vz, double *dPhi0);

    inline void transform_basic_background_profiles_to_lab_frame (double r, double *q, double *n,
                                                                  double *Ti, double *Te,
                                                                  double *Vth, double *Vz, double *dPhi0);

    void interp_basic_background_profiles_in_lab_frame (int dim, double *r,
                                                        double *q, double *n,
                                                        double *Ti, double *Te,
                                                        double *Vth, double *Vz, double *dPhi0);
};

#endif
