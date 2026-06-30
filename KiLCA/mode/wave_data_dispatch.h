/*! \file
    \brief Declarations for the Fortran wave_data_t accessors (kilca_wave_data_m,
    KiLCA/mode/wave_data_m.f90), replacing wave_data.h's C++ class now that
    wave_data is a Fortran per-instance handle. Still-C++ callers (mode.cpp,
    calc_mode.cpp, wave_code_interface.cpp) hold the handle as an opaque
    intptr_t and read/write fields through these getters/setters instead of
    `wd->field`.
*/

#ifndef WAVE_DATA_DISPATCH_INCLUDE
#define WAVE_DATA_DISPATCH_INCLUDE

#include <cstdint>

extern "C"
{
intptr_t wave_data_create_ (int m, int n, double olab_re, double olab_im, double omov_re, double omov_im);
void wave_data_destroy_ (intptr_t handle);

int wave_data_get_m_ (intptr_t handle);
int wave_data_get_n_ (intptr_t handle);
double wave_data_get_olab_re_ (intptr_t handle);
double wave_data_get_olab_im_ (intptr_t handle);
double get_wave_data_obj_omov_re_ (intptr_t handle);
double get_wave_data_obj_omov_im_ (intptr_t handle);
double wave_data_get_r_res_ (intptr_t handle);
void wave_data_set_r_res_ (intptr_t handle, double val);
double wave_data_get_det_re_ (intptr_t handle);
double wave_data_get_det_im_ (intptr_t handle);

void set_det_in_wd_struct_ (intptr_t *wd_ptr, double *re_det, double *im_det);
}

#endif
