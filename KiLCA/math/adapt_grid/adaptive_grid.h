/*! \file
    \brief The declarations of functions for generation of adaptive grids according to simple estimate of accuracy.
*/

#ifndef ADA_GRID_INCLUDE

#define ADA_GRID_INCLUDE

struct interval5
{
    double xs, ys;
    double xl, yl;
    double xm, ym;
    double xr, yr;
    double xe, ye;
    double err;
};

//external:
extern "C"
{
void calc_adaptive_1D_grid_ (void (*f)(double *, double *, void *p), void *p,
                             const int *max_dimx, double *eps, int *dimx, double *x, double *y);

void calc_adaptive_1D_grid_4vector_ (void (*f)(double *, double *, void *p), void *p,
                                     const int *max_dimx, double *eps,
                                     int *dimx, const double *x, const double *y);
}

//internal:
void set_interval (interval5 *I, double x1, double x2, double y1, double y2,
                   void (*f)(double *, double *, void *p), void *p);

void add_new_interval_to_the_array (interval5 *Iarr, int max_ind, int Icount,
                                    void  (*f)(double *, double *, void *p), void *p);

void update_max_err_interval (interval5 *Iarr, int max_ind,
                              void (*f)(double *, double *, void *p), void *p);

inline void eval_error (interval5 *I, double *err);

void check_and_remove_grid_condensations_ (double eps, int *dimx, double *x, double *y);

#endif
