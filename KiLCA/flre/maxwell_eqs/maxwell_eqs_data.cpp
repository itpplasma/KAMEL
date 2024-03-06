/*! \file maxwell_eqs_data.cpp
    \brief The implementation of maxwell_eqs_data class.
*/

#include <stdlib.h>
#include <stdio.h>
#include <cmath>
#include <cstring>
#include <climits>

#include "maxwell_eqs_data.h"

/*******************************************************************/

void copy_module_data_to_maxwell_eqs_data_struct (maxwell_eqs_data *me)
{
copy_module_data_to_maxwell_eqs_data_struct_f_ (&(me->num_vars), &(me->num_eqs), me->dim_Ersp_sys, me->iErsp_sys, me->dim_Brsp_sys, me->iBrsp_sys, (int *)me->der_order);

get_ersp_state_indices_and_dims_f_ ((int *)me->dim_Ersp_state, (int *)me->iErsp_state);

get_sys_ind_array_f_ (me->sys_ind);

//in C index starts from 0:
for (int k=0; k<3; k++)
{
    me->iErsp_state[k]--;
    me->iErsp_sys[k]--;
    me->iBrsp_sys[k]--;
}

for (int k=0; k<me->Nwaves; k++) me->sys_ind[k]--;
}

/*******************************************************************/

void maxwell_eqs_data::print (void)
{
fprintf (stdout, "\ncheck Maxwell system parameters:");
fprintf (stdout, "\nnum_vars=%d", num_vars);
fprintf (stdout, "\nnum_eqs=%d", num_eqs);

fprintf (stdout, "\ndim_Ersp_state: %d %d %d", dim_Ersp_state[0], dim_Ersp_state[1], dim_Ersp_state[2]);

fprintf (stdout, "\niErsp_state: %d %d %d", iErsp_state[0], iErsp_state[1], iErsp_state[2]);

fprintf (stdout, "\ndim_Ersp_sys: %d %d %d", dim_Ersp_sys[0], dim_Ersp_sys[1], dim_Ersp_sys[2]);

fprintf (stdout, "\niErsp_sys: %d %d %d", iErsp_sys[0], iErsp_sys[1], iErsp_sys[2]);

fprintf (stdout, "\ndim_Brsp_sys: %d %d %d", dim_Brsp_sys[0], dim_Brsp_sys[1], dim_Brsp_sys[2]);

fprintf (stdout, "\niBrsp_sys: %d %d %d", iBrsp_sys[0], iBrsp_sys[1], iBrsp_sys[2]);

fprintf (stdout, "\njr der_orders: %d %d %d", der_order[0][0], der_order[0][1], der_order[0][2]);

fprintf (stdout, "\njs der_orders: %d %d %d", der_order[1][0], der_order[1][1], der_order[1][2]);

fprintf (stdout, "\njp der_orders: %d %d %d", der_order[2][0], der_order[2][1], der_order[2][2]);

for (int k=0; k<Nwaves; k++)
{
    fprintf (stdout, "\nsys_ind[%d] = %d", k, sys_ind[k]);
}
}

/*******************************************************************/
