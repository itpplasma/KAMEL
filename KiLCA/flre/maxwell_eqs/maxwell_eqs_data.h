/*! \file maxwell_eqs_data.h
    \brief The declaration of maxwell_eqs_data class.
*/

#ifndef MAXWELL_EQS_INCLUDE

#define MAXWELL_EQS_INCLUDE

/*! \class maxwell_eqs_data
    \brief Class containes data needed to build the system of Maxwell equations.
*/
class maxwell_eqs_data
{
public:
    int Nwaves;                 //!<number of waves
    int num_vars;               //!<number of all unknowns
    int num_eqs;                //!<number of equations
    int max_der_order_Ersp[3];  //!<maximum order of derivs
    int der_order[3][3];        //!<order of derivs of E in a current (rsp)

    ///names of the variables:
    char **names_sys;
    char **names_state;
    char names_comp[3];

    //indices:
    int iEr, iEs, iEp, iBr, idBr, iBs, idBs, iBp, idBp, iddBp;

    ///number of an appropriate E components in a state vector:
    int dim_Ersp_state[3];
    int iErsp_state[3];

    ///number of an appropriate E and B components in a system vector:
    int dim_Ersp_sys[3];
    int iErsp_sys[3];
    int dim_Brsp_sys[3];
    int iBrsp_sys[3];

    int *sys_ind; //!<indices of a state components in a system array

    maxwell_eqs_data (int Nwaves_p)
    {
        Nwaves = Nwaves_p;
        sys_ind = new int[Nwaves];
    }

    ~maxwell_eqs_data (void)
    {
        delete [] sys_ind;
    }

    void print (void);
};

void copy_module_data_to_maxwell_eqs_data_struct (maxwell_eqs_data *me);

extern "C"
{
void copy_module_data_to_maxwell_eqs_data_struct_f_ (int *num_vars, int *num_eqs,
int *dim_Ersp_sys, int *iErsp_sys, int *dim_Brsp_sys, int *iBrsp_sys, int *der_order);

void get_ersp_state_indices_and_dims_f_ (int *dim_Ersp_state_p, int *iErsp_state_p);

void get_sys_ind_array_f_ (int *sys_ind_p);
}

#endif
