/*! \file
    \brief The declaration of zone class - the basic class for all types of zones.
*/

#ifndef ZONE_INCLUDE

#define ZONE_INCLUDE

/*******************************************************************/

#include "settings.h"
#include "background.h"
#include "wave_data.h"

/*******************************************************************/

#define PLASMA_MODEL_VACUUM 0
#define PLASMA_MODEL_MEDIUM 1
#define PLASMA_MODEL_IMHD 2
#define PLASMA_MODEL_RMHD 3
#define PLASMA_MODEL_FLRE 4

#define BOUNDARY_CENTER 0
#define BOUNDARY_INFINITY 1
#define BOUNDARY_IDEALWALL 2
#define BOUNDARY_INTERFACE 3
#define BOUNDARY_ANTENNA 4

/*******************************************************************/

const int Nbc = 5;
const char bc_str[Nbc][64] = {{"center"}, {"infinity"}, {"idealwall"}, {"interface"}, {"antenna"}};

const int Nmed = 5;
const char med_str[Nmed][64] = {{"vacuum"}, {"medium"}, {"imhd"}, {"rmhd"}, {"flre"}};

/*******************************************************************/

/*! \class zone
    \brief This basic class for all types of zones describes the zone - an interval over radius where a particular plasma model is used.
*/
class zone
{
public:
    double r1;      //!<left  zone boundary
    double r2;      //!<right zone boundary

    int bc1;        //!<type of left boundary
    int bc2;        //!<type of right boundary

    int medium;     //!<type of the plasma model
    int version;    //!<version of the code for the model

    int index;      //!<zone index

    const settings   *sd;       //!<pointer to settings structure
    const background *bp;       //!<pointer to background structure
    const wave_data  *wd;       //!<pointer to wave data

    char *path;                 //!<path2project

    int dim;                    //!<radial grid dimension
    double *r;                  //!<radial grid
    double *basis;              //!<independent basis solutions
    double *EB;                 //!<superposition field (system vector in the lab frame)
    double *S;                  //!<superposition coefficients

    int Nwaves;                 //!<number of waves
    int Ncomps;                 //!<number of (E, B) components in a sytem vector

    //Settings for the ME solution:
    int max_dim;    //!<max dimension of the radial grid for the solution
    double eps_rel; //!<relative accuracy for the solution
    double eps_abs; //!<absolute accuracy for the solution

    //ME solution space out settings:
    int deg;        //!<degree of the polynomial used to space out the solution
    double reps;    //!<relative accuracy of the sparse solution
    double aeps;    //!<absolute accuracy of the sparse solution
    double step;    //!<max grid step in the solution

    int flag_debug; //!<flag for debugging mode (additional checks are performed in the code)

    //indexing functions:
    inline int ib(int node, int sol, int comp, int part) const
    {
        //indexing function for basis array
        //must correspond to the following fortran definition:
        //complex(8), dimension(Ncomps,Nwaves,dim), intent(in) :: basis;

        //node = 0:dim-1, sol  = 0:Nwaves-1, comp = 0:Ncomps-1, part = 0:1 (real, imag)
        return part + 2*(comp + Ncomps*(sol + Nwaves*(node)));
    }

    inline int iEB(int node, int comp, int part) const
    {
        //indexing function for EB field state array
        //node = 0:dim-1, comp = 0:Ncomps-1, part = 0:1 (real, imag)
        return part + 2*(comp + Ncomps*(node));
    }

public:

    zone (const settings *, const background *, const wave_data *, char *, int);

    virtual ~zone (void); //very important to have it virtual: otherwise...

    void read (char *file);

    void print (void) const;

    double get_r1 (void) const { return r1; }

    double get_r2 (void) const { return r2; }

    int get_dim_of_basis (void) const { return Nwaves; }

    int get_dim_of_basis_vector (void) const { return Ncomps; }

    int get_code_version (void) const { return version; }

    int get_radial_grid_dimension (void) const { return dim; }

    double * get_basis_at_left_boundary (void) const { return &basis[ib(0, 0, 0, 0)]; }

    double * get_basis_at_right_boundary (void) const { return &basis[ib(dim-1, 0, 0, 0)]; }

    const settings * get_settings (void) const { return sd; }

    const background * get_background (void) const { return bp; }

    const wave_data * get_wave_data (void) const { return wd; }

    char * get_path_to_project (void) const { return path; }

    void save_basis_fields (char *);

    void calc_superposition_of_basis_fields (double *S);

    void save_final_fields (char *);

    void copy_radial_grid (double *);

    //for derived classes:
    virtual void read_settings (char * file) = 0;

    virtual void print_settings (void) = 0;

    virtual void calc_basis_fields (int flag = 0) = 0;

    virtual void copy_E_and_B_fields (double *) = 0;

    virtual void calc_final_fields (void) = 0;

    virtual void calc_dispersion (void) = 0;

    virtual void save_dispersion (void) = 0;

    virtual void calc_all_quants (void) = 0;

    virtual void save_all_quants (void) = 0;

    virtual void eval_diss_power_density (double, int, int, double *) = 0;

    virtual void eval_current_density (double, int, int, int, double *) = 0;
};

/*******************************************************************/

int selector (const struct dirent *ent);
int determine_number_of_zones (char *path2project);
int determine_zone_type (char *file);
char * get_zone_file_name (char * path2project, int zone_index);

/*******************************************************************/

extern "C"
{
void get_right_boundary_of_zone_ (zone **ptr, double *r);
void get_left_boundary_of_zone_ (zone **ptr, double *r);
}

/*****************************************************************************/

#endif
