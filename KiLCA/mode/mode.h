/*! \file
    \brief The declaration of mode_data class.
*/

#ifndef MODE_INCLUDE

#define MODE_INCLUDE

#include <complex>

#include "settings.h"
#include "background.h"
#include "wave_data.h"
#include "zone.h"

/*****************************************************************************/

inline int iFFM(int node, int comp, int part)
{
    //indexing function for E, B fields in final EB array in mode class
    return part + 2*(comp + 6*node);
}

/*****************************************************************************/

/*! \class mode_data
    \brief Class represents data related to a particular mode of a perturbation.
*/
class mode_data
{
public:
    //pointers to common structures:
    const settings_t *sd;         //!<pointer to settings
    const background *bp;       //!<pointer to background

    //pointers to structures for the given mode:
    wave_data *wd;              //!<pointer to wave_data

    //Directories for a given mode:
    char *path2linear;          //!<path to linear-data
    char *path2dispersion;      //!<path to dispersion-data
    char *path2poincare;        //!<path to poincare-data

    int  Nzones;                //!<number of zones
    zone **zones;               //!<zones array (pointers)

    int dim;         //!<final radial grid dimension
    double *r;       //!<final radial grid
    double *EB;      //!<final EB fields on the r grid in lab frame in cyl coordinates
    double *EB_int;  //!<EB_lab fields ready for the interpolation
    int *index;      //!<zone first point index in the final grids

    int Nc;      //!<number of superposition coefficients
    double *A;   //!<system matrix
    double *B;   //!<system rhs vector
    double *S;   //!<superposition coefficients

public:
    mode_data (int m, int n, complex<double> olab, const settings_t *, const background *);

    ~mode_data (void)
    {
        delete [] path2linear;
        delete [] path2dispersion;
        delete [] path2poincare;

        delete wd;

        delete [] r;
        delete [] EB;
        delete [] EB_int;
        delete [] index;

        delete [] A;
        delete [] B;
        delete [] S;

        //loop over zones:
        for (int iz=0; iz<Nzones; iz++)
        {
            if (zones[iz]) delete zones[iz];
        }
        delete [] zones;
    }

    inline int iFint(int node, int comp, int part)
    {
        //indexing function for E, B fields in EB_int array used for interpolation
        return node + dim*(part + 2*(comp));
    }

    void calc_all_mode_data (int flag = 0);

    void allocate_and_setup_zones (void);

    void set_and_make_mode_data_directories (void);

    int determine_zone_index_for_point (double);

    void divEB (double x, complex<double> *divEB);

    void calc_and_save_divEB (void);

    void eval_EB_fields (double, complex<double> *);

    void eval_diss_power_density (double, int type, int spec, double *pdd);

    void eval_current_density (double x, int type, int spec, int comp, double * J);

private:
    int find_resonance_location (void);

    void check_zones_parameters (void);

    void calc_basis_fields_in_zones (int flag = 0);

    void calc_dispersion_in_zones (void);

    void save_dispersion_in_zones (void);

    void calc_stitching_equations (void);

    void calc_stitching_equations_determinant (void);

    void save_mode_det_data (void);

    void solve_stitching_equations (void);

    void calc_superposition_of_basis_fields_in_zones (void);

    void space_out_fields_in_zones (void);

    void combine_final_wave_fields (void);

    void setup_wave_fields_for_interpolation (void);

    void save_final_wave_fields (void);

    void calc_quants_in_zones (void);

    void save_quants_in_zones (void);

    void combine_final_quants (void);

    void save_final_quants (void);
};

/*****************************************************************************/

struct func_params
{
    const background *bp;
    double q_res;
};

/*****************************************************************************/

double qminusq0 (double x, void *p);

/*****************************************************************************/

void eval_path_to_linear_data (const char *path2project, int m, int n, complex<double> olab, char *path2linear);

void eval_path_to_dispersion_data (const char *path2project, int m, int n, complex<double> olab, char *path2dispersion);

void eval_path_to_poincare_data (const char *path2project, int m, int n, complex<double> olab, char *path2poincare);

/*****************************************************************************/

extern "C"
{
void set_settings_in_mode_data_module_ (const settings_t **);

void set_back_profiles_in_mode_data_module_ (const background **);

void set_wave_data_in_mode_data_module_ (wave_data **);

void set_wave_parameters_in_mode_data_module_ (int *, int *, double *, double *, double *, double *);

void set_resonance_location_in_mode_data_module_ (double *r_res);

void set_current_density_in_antenna_module_ (void);

void copy_mode_paths_to_mode_data_module_ (mode_data **);

void copy_mode_paths_from_mode_data_struct_ (mode_data **md, char *path2linear, char *path2dispersion, char *path2poincare);

void clear_all_data_in_mode_data_module_ (void);

void set_det_in_wd_struct_ (wave_data **wd_ptr, double *re_det, double *im_det);
}

/*****************************************************************************/

#endif
