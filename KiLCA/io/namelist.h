#pragma once

#include <complex>

extern "C" {

void read_namelist(double *rtor, double *rp, double *B0, char *path2profiles,
                   int *calc_back, char *flag_back, int *N, double *V_gal_sys,
                   double *V_scale, double *m_i, double *zele, double *zion,
                   bool *flag_debug);
void read_namelist_unit_test(int *a, double *b, char *c, char *d, int *e, std::complex<double> *f);

}
