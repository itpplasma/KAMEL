/*! \file
    \brief The declaration of spline_data class.
*/

#ifndef SPLINE_INCLUDE

#define SPLINE_INCLUDE

#include <inttypes.h>

/*! \class spline_data
    \brief Class implements splines of arbitrary degree for arrays of data.
*/
struct spline_data
{
    int N;              //!<spline degree
    int ind;            //!<last search index for the faster spline evaluation

    int type;           //!<type of spline boundary

    int dimx;           //!<dimension of a profile x grid (all the same)
    const double *x;    //!<x grid for profiles

    double *C;          //!<matrix of splines coefficients

    double *BC;         //!<array of binomial coefficients

    double *fac;        //!<factors used in the splines evaluation

    int calc_spline_boundaries (int dimy, const double *y, double *W);
    int calc_spline_coefficients (int Imin, int dimy, const double *y, double *W);
};

//interface subs:
extern "C"
{
void spline_alloc_ (int N, int type, int dimx, const double *x, double *C, uintptr_t *sid);

void spline_calc_ (uintptr_t sid, const double *y, int Imin, int Imax, double *W, int *ierr);

void spline_eval_ (uintptr_t sid, int dimz, double *z, int Dmin, int Dmax, int Imin, int Imax, double *R);

void spline_eval_d_ (uintptr_t sid, int dimz, double *z, int Dmin, int Dmax, int Imin, int Imax, double *R);

void spline_free_ (uintptr_t sid);
}

//external subs:
extern "C"
{
void dgesv_ (int *, int *, double *, int *, int *, double *, int *, int *);
void dgbsv_ (int *, int *, int *, int *, double *, int *, int *, double *, int *, int *);
}

//internal subs:
inline void search_array (double x, int dimx, const double *xa, int *ind);

inline int binary_search (double x, const double *xa, int ilo, int ihi);

inline int sign (double x);

inline void set_bc_array (int N, double *BC);

inline void set_fac_array (int N, double *fac);

#endif
