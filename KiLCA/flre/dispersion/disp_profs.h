/*! \file disp_profs.h
    \brief The declaration of disp_profiles class.
*/

#ifndef DISP_PROFS_INCLUDE

#define DISP_PROFS_INCLUDE

#include <inttypes.h>

#include "spline.h"

/*! \class disp_profiles
    \brief Class stores various data related to dispersion profiles.
*/
class disp_profiles
{
public:
    int Nwaves; //!<number of waves

    char *flag_back; //!<background flag

    int dimx;   //!<dimension of the x grid
    double *x;  //!<x grid for profiles

    int dimk;  //!<number of k dispersion profiles
    double *k; //!<k dispersion array

    int dimp;  //!<number of polarization vector profiles
    double *p; //!<p dispersion array

    ///splines:
    int N;          //!<spline degree
    uintptr_t sid;  //!<spline id
    double *S;      //!<matrix of splines coefficients
    double *R;      //!<array for all spline values (+derivs) at some point

    disp_profiles (int, int, double *, const char *);

    ~disp_profiles (void)
    {
        if (flag_back) delete [] flag_back;

        if (k) delete [] k;
        if (p) delete [] p;

        if (S) delete [] S;
        if (R) delete [] R;

        if (sid) spline_free_ (sid);
    }

    void calculate_dispersion_profiles (void);

    void sort_dispersion_profiles (void);

    void save_dispersion_profiles (char *);
};

extern "C"
{
void calc_dispersion_ (double *, char *, int *, double *, double *, int);

void dpsort_ (double *, int *, int *, int *, int *);
}

#endif
