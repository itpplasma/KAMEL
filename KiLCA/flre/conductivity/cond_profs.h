/*! \file cond_profs.h
    \brief The declaration of cond_profiles class.
*/

#ifndef COND_PROFS_INCLUDE

#define COND_PROFS_INCLUDE

#include <ctime>
#include <inttypes.h>

#include "constants.h"
#include "settings.h"
#include "background.h"
#include "wave_data.h"
#include "spline.h"

/*! \class cond_profiles
    \brief Class containes conductivity profiles.
*/
class cond_profiles
{
public:
    ///external data which are used for calculations:
    const settings *sd;   //!<pointer to settings data
    const background *bp; //!<pointer to background
    const wave_data *wd;  //!<wave data

    char *flag_back;   //!<background flag for this conductivity profiles
    char *path2linear; //!<path to linear-data

    int NC; //!<spline degree for C
    int NK; //!<spline degree for K

    int flreo;    //!<order of finite Larmour radius expansion
    int gal_corr; //!<flag whether to use Galilean correction term

    int ind; //!<index of a radial point for spline evaluation

    int dimt; //!<number of K types

    double r1, r2; //!<flre zone boundaries

    double D;       //!<width of resonant layer (RL)
    double eps_out; //!<error parameter used in the adaptive radial grid generation outside RL
    double eps_res; //!<error parameter used in the adaptive radial grid generation in RL

    ///local data to compute:
    int dimx;  //!<dimension of a profile x grid (all the same)
    double *x; //!<x grid for profiles

    int dimK;       //!<number of K splines
    double *K;      //!<K grid
    double *CK;     //!<matrix of splines coefficients
    double *RK;     //!<array of spline values (+derivs) at some point
    uintptr_t sidK; //!<spline_data id

    int dimC;       //!<number of C splines
    double *C;      //!<C grid
    double *CC;     //!<matrix of splines coefficients
    double *RC;     //!<array of spline values (+derivs) at some point
    uintptr_t sidC; //!<spline_data id

    double *bico;   //!<binomial coefficients

    double *xt;     //!<temp space
    double *yt;     //!<temp space

    cond_profiles (void) : sd(0), bp(0), wd(0),
                           flag_back(0), path2linear(0),
                           x(0),
                           K(0), CK(0), RK(0), sidK(0),
                           C(0), CC(0), RC(0), sidC(0),
                           bico(0),
                           xt(0), yt(0)
    {
    }

    ~cond_profiles (void)
    {
        delete [] x;

        delete [] K;
        delete [] CK;
        delete [] RK;

        delete [] C;
        delete [] CC;
        delete [] RC;

        delete [] bico;

        delete [] flag_back;
        delete [] path2linear;

        spline_free_ (sidK);
        spline_free_ (sidC);
    }

    void calc_and_spline_main_conductivity_profiles (const class flre_zone *, int flag = 0);

    /*
    K profiles memory aligment:
    {spec=0:1, {[type=0:1, {p=0:flreo, {q=0:flreo, {i=0:2, {j=0:2], {(re, im), {i=0:dimx-1}}}}}}}}
    */

    /*
    conductivity C profiles memory aligment:
    {spec=0:1, {type=0:1, {s=0..2flreo, {i=0:2, {j=0:2, {(re, im), {i=0:dimx-1}}}}}}}
    */

    int iKs (int spec, int type, int p, int q, int i, int j, int part, int node) const
    {
    ///index of the element in K array for splines:
    return node + (this->dimx)*(part + 2*(j + 3*(i + 3*(q + (this->flreo+1)*(p +
           (this->flreo+1)*(type + (this->dimt)*(spec)))))));
    }

    int iKa (int spec, int type, int p, int q, int i, int j, int part, int node) const
    {
    ///index of the element in initial K array of matrices:
    return part + 2*(j + 3*(i + 3*(q + (this->flreo+1)*(p + (this->flreo+1)*
           (type + (this->dimt)*(spec + 2*node))))));
    }

    int iCs (int spec, int type, int s, int i, int j, int part, int node) const
    {
    ///index of the element in C array for splines:
    return node + (this->dimx)*(part + 2*(j + 3*(i + 3*(s + (2*(this->flreo)+1)*
           (type + (this->dimt)*(spec))))));
    }

    int iCa (int s, int i, int j, int part) const
    {
    ///index of the element in C matrix:
    return part + 2*(j + 3*(i + 3*s));
    }
};

extern "C"
{
void alloc_conductivity_profiles_ (cond_profiles **cpptr);

void calc_and_spline_conductivity_for_point_ (settings **sd_ptr, background **bp_ptr,
                          wave_data **wd_ptr, char *flag_back, double *r,
                          cond_profiles **cp_ptr);

void delete_conductivity_profiles_f_ (cond_profiles **cp_ptr);
}

#endif
