/*! \file maxwell_eqs_data.h
    \brief C entry points for the per-zone Maxwell-equations layout snapshot,
           now owned by the Fortran kilca_maxwell_eqs_data_m module (each
           instance addressed by an opaque handle, since multiple zones can be
           alive at once). The former C++ maxwell_eqs_data class has been
           translated away.
*/

#ifndef MAXWELL_EQS_INCLUDE

#define MAXWELL_EQS_INCLUDE

#include <cstdint>
#include <cstdio>

extern "C"
{
intptr_t maxwell_eqs_data_create_ (int Nwaves);

void maxwell_eqs_data_destroy_ (intptr_t handle);

int get_me_nwaves_ (intptr_t handle);

int get_me_num_vars_ (intptr_t handle);

int get_me_num_eqs_ (intptr_t handle);

int get_me_der_order_ (intptr_t handle, int i, int j);

int get_me_dim_ersp_state_ (intptr_t handle, int k);

int get_me_iersp_state_ (intptr_t handle, int k);

int get_me_dim_ersp_sys_ (intptr_t handle, int k);

int get_me_iersp_sys_ (intptr_t handle, int k);

int get_me_dim_brsp_sys_ (intptr_t handle, int k);

int get_me_ibrsp_sys_ (intptr_t handle, int k);

int get_me_sys_ind_ (intptr_t handle, int k);
}

/*! \brief Mirrors the former maxwell_eqs_data::print(), reading via the getters. */
inline void print_maxwell_eqs_data (intptr_t me)
{
fprintf (stdout, "\ncheck Maxwell system parameters:");
fprintf (stdout, "\nnum_vars=%d", get_me_num_vars_ (me));
fprintf (stdout, "\nnum_eqs=%d", get_me_num_eqs_ (me));

fprintf (stdout, "\ndim_Ersp_state: %d %d %d", get_me_dim_ersp_state_ (me, 0), get_me_dim_ersp_state_ (me, 1), get_me_dim_ersp_state_ (me, 2));

fprintf (stdout, "\niErsp_state: %d %d %d", get_me_iersp_state_ (me, 0), get_me_iersp_state_ (me, 1), get_me_iersp_state_ (me, 2));

fprintf (stdout, "\ndim_Ersp_sys: %d %d %d", get_me_dim_ersp_sys_ (me, 0), get_me_dim_ersp_sys_ (me, 1), get_me_dim_ersp_sys_ (me, 2));

fprintf (stdout, "\niErsp_sys: %d %d %d", get_me_iersp_sys_ (me, 0), get_me_iersp_sys_ (me, 1), get_me_iersp_sys_ (me, 2));

fprintf (stdout, "\ndim_Brsp_sys: %d %d %d", get_me_dim_brsp_sys_ (me, 0), get_me_dim_brsp_sys_ (me, 1), get_me_dim_brsp_sys_ (me, 2));

fprintf (stdout, "\niBrsp_sys: %d %d %d", get_me_ibrsp_sys_ (me, 0), get_me_ibrsp_sys_ (me, 1), get_me_ibrsp_sys_ (me, 2));

fprintf (stdout, "\njr der_orders: %d %d %d", get_me_der_order_ (me, 0, 0), get_me_der_order_ (me, 0, 1), get_me_der_order_ (me, 0, 2));

fprintf (stdout, "\njs der_orders: %d %d %d", get_me_der_order_ (me, 1, 0), get_me_der_order_ (me, 1, 1), get_me_der_order_ (me, 1, 2));

fprintf (stdout, "\njp der_orders: %d %d %d", get_me_der_order_ (me, 2, 0), get_me_der_order_ (me, 2, 1), get_me_der_order_ (me, 2, 2));

int Nwaves = get_me_nwaves_ (me);
for (int k=0; k<Nwaves; k++)
{
    fprintf (stdout, "\nsys_ind[%d] = %d", k, get_me_sys_ind_ (me, k));
}
}

#endif
