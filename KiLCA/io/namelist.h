#pragma once

#include <complex>

extern "C" {

void read_namelist(double *rtor, double *rp, double *B0, char *path2profiles, int *calc_back,
                   char *flag_back, int *N, double *V_gal_sys, double *V_scale, double *m_i,
                   double *zele, double *zion,
                   char *fname_buffer, bool *search_flag, int *rdim, double *rfmin, double *rfmax,
                   int *idim, double *ifmin, double *ifmax, bool *stop_flag, double *eps_res,
                   double *eps_abs, double *eps_rel, double *delta, bool *test_roots, int *Nguess,
                   int *kmin, int *kmax, std::complex<double> *fstart,
                   double *ra, double *wa, double *I0, std::complex<double> *flab, int *dma,
                   int *modes, bool *flag_eigmode,
                   bool *flag_debug);

void read_namelist_unit_test(int *a, double *b, char *c, char *d, int *e, std::complex<double> *f);

}
