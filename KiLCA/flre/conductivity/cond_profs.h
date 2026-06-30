/*! \file cond_profs.h
    \brief C entry points for the conductivity (K/C matrix) profiles, now owned
           by the Fortran kilca_cond_profiles_m module (each instance addressed
           by an opaque handle, since each flre_zone owns its own). The former
           C++ cond_profiles class has been translated away.
*/

#ifndef COND_PROFS_INCLUDE

#define COND_PROFS_INCLUDE

#include <cstdint>

extern "C"
{
intptr_t cond_profiles_create_ (const char *path2linear_p, int flreo, int gal_corr,
                                int N, int max_dim_c, double r1, double r2, double D,
                                double eps_out, double eps_res, double a, double b,
                                double r_res, double omov_re, double omov_im,
                                int flag_debug, int flag);

void cond_profiles_destroy_ (intptr_t handle);

int get_cond_nk_ (intptr_t handle);

int get_cond_dimk_ (intptr_t handle);

int get_cond_nc_ (intptr_t handle);

int get_cond_dimc_ (intptr_t handle);

int get_cond_dimx_ (intptr_t handle);

double * get_cond_x_ptr_ (intptr_t handle);

double * get_cond_k_ptr_ (intptr_t handle);

int get_cond_iks_ (intptr_t handle, int spec, int type, int p, int q, int i, int j, int part, int node);

char get_cond_flag_back_ (intptr_t handle);

void get_cond_path2linear_ (intptr_t handle, char *buf);

void eval_all_k_matrices_ (intptr_t handle, int Dmin, int Dmax, double r, double *K);

void eval_all_c_matrices_ (intptr_t handle, int Dmin, int Dmax, double r, double *C);
}

#endif
