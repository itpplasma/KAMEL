/*! \file sysmat_profs.h
    \brief The declaration of sysmat_profiles class.
*/

#ifndef SYSMAT_PROFS_INCLUDE

#define SYSMAT_PROFS_INCLUDE

#include <ctime>
#include <inttypes.h>

#include "settings.h"
#include "background.h"
#include "wave_data.h"
#include "spline.h"

/*******************************************************************/

/*! \class sysmat_profiles
    \brief Class represents ODE system (u' = A*u) matrix evaluated on a grid of radial nodes.
*/
class sysmat_profiles
{
public:
    char *flag_back;   //!<background flag for this profiles
    char *path2linear; //!<path to linear-data

    int N;   //!<spline degree
    int ind; //!<index

    int dimx;  //!<dimension of a profile x grid (all the same)
    double *x; //!<x grid for profiles

    int Nwaves;

    ///system matrix:
    int dimM;  //!<number of splines
    double *M; //!<M grid
    double *C; //!<matrix of splines coefficients
    double *R; //!<array for spline values (+derivs) at some point
    uintptr_t sidM; //!<spline id

    double *xt; //!<temp space
    double *yt; //!<temp space

    sysmat_profiles (void) {}

    ~sysmat_profiles (void)
    {
        delete [] flag_back;
        delete [] path2linear;

        delete [] x;

        delete [] M;
        delete [] C;
        delete [] R;

        spline_free_ (sidM);
    }

    void calc_and_spline_sysmatrix_profiles (const class flre_zone *zone);

    void save_M_matrix (int dimf);
};

/*******************************************************************/

void sample_sysmat_func (double *r, double *f, void *p);

/*******************************************************************/

extern "C"
{
void alloc_sysmatrix_profiles_ (sysmat_profiles **spptr);
void calc_diff_sys_matrix_ (double *r, char *flag_back, double *R, int flag_back_len);
}

#endif
