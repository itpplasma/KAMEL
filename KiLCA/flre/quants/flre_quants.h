/*! \file flre_quants.h
    \brief The declaration of flre_quants class.
*/

#ifndef FLRE_QUANTS_INCLUDE

#define FLRE_QUANTS_INCLUDE

#include "constants.h"
#include "settings.h"
#include "cond_profs.h"

/*******************************************************************/

typedef void (*calc_func)(class flre_quants *);
typedef void (*save_func)(const class flre_quants *);

/*******************************************************************/

/*! \class flre_quants
    \brief The class represents auxiliary quantities (like current densities, dissipated power, etc) in a FLRE zone.
*/
class flre_quants
{
public:
    class flre_zone *zone;      //!<pointer to zone

    char *path2linear;          //!<path to linear data

    double vol_fac;             //!<factor arising in integration over cylinder surface

    double *bico;               //!<binomial coefficients

    int num_tot;                //!<fixed total number of quantities

    int *dimq;                  //!<array of dimensions for each quantity

    char **name;                //!<names of the quantities

    int *ind;                   //!<index of the quantity in the computational list

    calc_func *calc_quant;      //!<array of functions to calculate quants

    save_func *save_quant;      //!<array of functions to save quants

    //quants indices in the typical order of calculation
    int CURRENT_DENS;
    int ABS_POWER_DENS;
    int DISS_POWER_DENS;
    int KIN_FLUX;
    int POY_FLUX;
    int TOT_FLUX;

    //depend on previous quants on the grid:
    int NUMBER_DENS;
    int LOR_TORQUE_DENS;

    int numq; //!<actual number of quantities to compute

    int *dni; //!<index of the needed quantity in the global array

    //state variables:
    double r;           //!<current r value

    int node;           //!<current grid node index

    int *flag;          //!<flag showing whether the quantity is already computed for given r and node

    double **qloc;      //!<all local quants array

    double **qint;      //!<all integrated quants array

    int flagC;          //!<flag if C is computed for the node

    double *C;          //!<C matrices storage

    int flagK;          //!<flag if K is computed for the node

    double *K;          //!<K matrices storage

    int Dmax;           //!<maximum order of K derivatives used for evaluation

    int flreo;          //!<FLRE order

    int dimx;           //!<radial grid dimension

    double *x;          //!<radial grid

    double JaE;         //!<work of the perturbation field over antenna current density

    int N;              //!<spline degree

    //jr splines
    uintptr_t sidY;   //!<spline id
    int nY;           //!<number of splines
    double *Y;        //!<y grid of values
    double *S;        //!<matrix of spline coefficients
    double *R;        //!<array of spline values (+derivs) at some point
    int flagS;        //!<flag if spline is computed

    double *jaE;
    double *jaEi;

    double * cdlab;   //!current densities in lab frame used for interpolation

    flre_quants (flre_zone *Z);

    ~flre_quants (void)
    {
        delete [] path2linear;

        delete [] bico;

        for (int i=0; i<num_tot; i++) delete [] name[i];

        delete [] name;
        delete [] dimq;
        delete [] flag;

        delete [] calc_quant;
        delete [] save_quant;

        delete [] ind;
        delete [] dni;

        for (int i=0; i<numq; i++)
        {
            delete [] qloc[i];
            delete [] qint[i];
        }

        if (qloc) delete [] qloc;
        if (qint) delete [] qint;

        if (C) delete [] C;
        if (K) delete [] K;

        //splines:
        if (Y) delete [] Y;
        if (S) delete [] S;
        if (R) delete [] R;

        if (sidY) spline_free_ (sidY);

        if (jaE) delete [] jaE;
        if (jaEi) delete [] jaEi;

        if (cdlab) delete [] cdlab;
    }

    //member functions:
    void set_all_known_quants_parameters (void);

    void set_null_node_state (int k);

    void calculate_local_profiles (void);

    void calculate_integrated_profiles (void);

    void save_profiles (void);
};

/*******************************************************************/

inline int binary_search (double x, const double *xa, int ilo, int ihi);

/*******************************************************************/

extern "C"
{

}

/*******************************************************************/

#endif
