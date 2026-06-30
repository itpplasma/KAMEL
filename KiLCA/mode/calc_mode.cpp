/*! \file
    \brief The implementation of mode_data class and other functions declared in mode.h.
*/

#include "fortnum.h"

#include <stdlib.h>
#include <stdio.h>
#include <cmath>
#include <cstring>
#include <dirent.h>
#include <fnmatch.h>

#include "constants.h"
#include "mode.h"
#include "eval_back.h"
#include "inout.h"
#include "stitching.h"
#include "interp.h"

/*****************************************************************************/

//zone settings-file discovery, formerly in zone.cpp: scans path2project for
//zone_N.in files and determines each one's medium type. allocate_and_setup_zones
//(below) is the sole caller.

static int selector(const struct dirent *ent)
{
    if (!fnmatch("zone_*.in", ent->d_name, 0)) return 1;
    else                                       return 0;
}

static int determine_number_of_zones(char *path2project)
{
    int Nzones = 0;

    DIR *dp;
    struct dirent *ep;

    if (dp = opendir(path2project))
    {
        while (ep = readdir(dp))
        {
            if (!fnmatch("zone_*.in", ep->d_name, 0)) Nzones++;
        }
        closedir(dp);
    }
    else
    {
        fprintf(stderr, "\ndetermine_number_of_zones: faled to open the project directory %s", path2project);
        exit(1);
    }

    if (DEBUG_FLAG)
    {
        fprintf(stdout, "\nNzones = %d", Nzones); fflush(stdout);
    }

    return Nzones;
}

static char *get_zone_file_name(char *path2project, int zone_index)
{
    DIR *dp;
    struct dirent *ep;

    char *file_name = new char[1024];

    strcpy(file_name, path2project);

    int count = 0;

    if (dp = opendir(path2project))
    {
        char *file_pattern = new char[1024];

        sprintf(file_pattern, "*zone_%d*.in", zone_index + 1);

        while (ep = readdir(dp))
        {
            if (!fnmatch(file_pattern, ep->d_name, 0))
            {
                ++count;
                strcat(file_name, ep->d_name);
                break;
            }
        }
        closedir(dp);
        delete[] file_pattern;
    }
    else
    {
        fprintf(stderr, "\nget_zone_file_name: faled to open the project directory %s", path2project);
        exit(1);
    }

    if (DEBUG_FLAG)
    {
        fprintf(stdout, "\nget_zone_file_name: file name for zone %d is: %s", zone_index, file_name);
        fflush(stdout);
    }

    if (!count)
    {
        fprintf(stderr, "\nget_zone_file_name: failed to find the file name for zone %d!", zone_index);
        delete[] file_name;
        exit(1);
    }

    return file_name;
}

static int determine_zone_type(char *file)
{
    FILE *in;

    if ((in = fopen(file, "r")) == NULL)
    {
        fprintf(stderr, "\nerror: determine_zone_type: failed to open file %s\a\n", file);
        exit(0);
    }

    char *str_buf = new char[1024];
    char *mstr = new char[64];

    read_line_2skip_it(in, &str_buf);
    read_line_2skip_it(in, &str_buf);
    read_line_2skip_it(in, &str_buf);
    read_line_2get_string(in, &mstr);

    fclose(in);

    static const int Nmed = 5;
    static const char med_str[Nmed][64] = {{"vacuum"}, {"medium"}, {"imhd"}, {"rmhd"}, {"flre"}};

    int medium = -1;

    for (int k = 0; k < Nmed; k++)
    {
        if (!strcmp(med_str[k], mstr))
        {
            medium = k;
            break;
        }
    }

    if (medium == -1)
    {
        fprintf(stderr, "\nerror: determine_zone_type: medium type is unknown: %s", mstr);
        exit(1);
    }

    delete[] str_buf;
    delete[] mstr;

    return medium;
}

/*****************************************************************************/

extern "C"
{
void center_equations_hommed_ (int *Nw, int *len, intptr_t *code, double *EB, int *neq, int *nvar, double *M, double *J);
void infinity_equations_hommed_ (int *Nw, int *len, intptr_t *code, double *EB, int *neq, int *nvar, double *M, double *J);
void ideal_wall_equations_hommed_ (int *Nw, int *len, intptr_t *code, double *EB, int *neq, int *nvar, double *M, double *J);
void stitching_equations_hommed_hommed_ (int *Nw1, int *len1, intptr_t *code1, double *EB1,
                     int *Nw2, int *len2, intptr_t *code2, double *EB2,
                     int *flg_ant, int *neq, int *nvar, double *M, double *J);

void center_equations_imhd_ (int *Nw, int *len, intptr_t *code, double *EB, int *neq, int *nvar, double *M, double *J);
void infinity_equations_imhd_ (int *Nw, int *len, intptr_t *code, double *EB, int *neq, int *nvar, double *M, double *J);
void ideal_wall_equations_imhd_ (int *Nw, int *len, intptr_t *code, double *EB, int *neq, int *nvar, double *M, double *J);
void stitching_equations_imhd_hommed_ (int *Nw1, int *len1, intptr_t *code1, double *EB1,
                                       int *Nw2, int *len2, intptr_t *code2, double *EB2,
                                       int *flg_ant, int *neq, int *nvar, double *M, double *J);
void stitching_equations_imhd_imhd_ (int *Nw1, int *len1, intptr_t *code1, double *EB1,
                                     int *Nw2, int *len2, intptr_t *code2, double *EB2,
                                     int *flg_ant, int *neq, int *nvar, double *M, double *J);

void center_equations_flre_ (int *Nw, int *len, intptr_t *code, double *EB, int *neq, int *nvar, double *M, double *J);
void infinity_equations_flre_ (int *Nw, int *len, intptr_t *code, double *EB, int *neq, int *nvar, double *M, double *J);
void ideal_wall_equations_flre_ (int *Nw, int *len, intptr_t *code, double *EB, int *neq, int *nvar, double *M, double *J);
void stitching_equations_flre_hommed_ (int *Nw1, int *len1, intptr_t *code1, double *EB1,
                                       int *Nw2, int *len2, intptr_t *code2, double *EB2,
                                       int *flg_ant, int *neq, int *nvar, double *M, double *J);
void stitching_equations_flre_flre_ (int *Nw1, int *len1, intptr_t *code1, double *EB1,
                                     int *Nw2, int *len2, intptr_t *code2, double *EB2,
                                     int *flg_ant, int *neq, int *nvar, double *M, double *J);

void update_system_matrix_and_rhs_vector_ (int *Nc, double *A, double *B, int *neq, int *ieq, int *nvar, int *ivar, double *M, double *J);
void calc_system_determinant_ (int *NoC, double *A, double *det);
void find_superposition_coeffs_ (int *NoC, double *A, double *B, double *S);
}

/*****************************************************************************/

double qminusq0(double x, void *p)
{
    func_params *P = (func_params *)p;

    return q(x, P->bp) - P->q_res;
}

/*****************************************************************************/

int mode_data::find_resonance_location(void)
{
    double q_res = - ((double)(wave_data_get_m_(wd))) / ((double)(wave_data_get_n_(wd)));

    double r1 = get_background_x0_ (), r2 = get_background_xlast_ ();

    if ((q(r1, bp) - q_res) * (q(r2, bp) - q_res) > 0)
    {
        if (DEBUG_FLAG)
        {
            fprintf(stdout, "\nfind_resonance_location: resonant surface for the mode m=%d n=%d is absent", wave_data_get_m_(wd), wave_data_get_n_(wd));
        }
        wave_data_set_r_res_ (wd, 0.0e0);
        return 0;
    }

    int max_iter = 100;

    struct func_params params = {bp, q_res};

    double root = 0.0e0;
    int status = fortnum_root_brent(&qminusq0, r1, r2, 0.0e0, 0.0e0, max_iter,
                                    &root, &params);

    if (status != FORTNUM_OK)
    {
        fprintf(stdout, "\nfind_resonance_location: failed to find resonant surface for the mode m=%d n=%d", wave_data_get_m_(wd), wave_data_get_n_(wd));
        wave_data_set_r_res_ (wd, 0.0e0);
        return 0;
    }

    wave_data_set_r_res_ (wd, root);

    if (DEBUG_FLAG)
    {
        fprintf(stdout, "\nresonant surface for the mode m=%d n=%d is found at:\nr=%.16le,  q(r)=%.16le\n",
                wave_data_get_m_(wd), wave_data_get_n_(wd), wave_data_get_r_res_(wd), q(wave_data_get_r_res_(wd), bp));
    }
    return 1;
}

/*****************************************************************************/

void mode_data::check_zones_parameters(void)
{
    for (int k = 0; k < Nzones - 1; k++)
    {
        if (zone_get_r2_(zones[k]) != zone_get_r1_(zones[k + 1]) || zone_get_bc2_(zones[k]) != zone_get_bc1_(zones[k + 1]))
        {
            fprintf(stdout, "\ncheck_zones: boundaries are different: k = %d", k);
            exit(1);
        }

        if (!(zone_get_bc2_(zones[k]) == BOUNDARY_INTERFACE || zone_get_bc2_(zones[k]) == BOUNDARY_ANTENNA))
        {
            fprintf(stdout, "\ncheck_zones: improper type of the right boundary for zone %d: %d", k, zone_get_bc2_(zones[k]));
            exit(1);
        }
    }

    int iz;

    // check that the first boundary is either center or ideal wall:
    iz = 0;
    if (!(zone_get_bc1_(zones[iz]) == BOUNDARY_CENTER || zone_get_bc1_(zones[iz]) == BOUNDARY_IDEALWALL))
    {
        fprintf(stdout, "\ncheck_zones: improper type of the first boundary: %d", zone_get_bc1_(zones[iz]));
        exit(1);
    }

    // check that the last boundary is either infinity or ideal wall:
    iz = Nzones - 1;
    if (!(zone_get_bc2_(zones[iz]) == BOUNDARY_INFINITY || zone_get_bc2_(zones[iz]) == BOUNDARY_IDEALWALL))
    {
        fprintf(stdout, "\ncheck_zones: improper type of the last boundary: %d", zone_get_bc2_(zones[iz]));
        exit(1);
    }

    if (DEBUG_FLAG)
    {
        fprintf(stdout, "\nzones consistency check passed...\n");
    }
}

/*******************************************************************/

void mode_data::allocate_and_setup_zones(void)
{
    Nzones = determine_number_of_zones(sd->path2project);

    if (Nzones < 2)
    {
        fprintf(stderr, "\nallocate_and_setup_zones: Nzones must be >= 2: %d", Nzones);
        exit(1);
    }

    zones = new intptr_t[Nzones]; // array of handles to Fortran zone_t instances

    for (int k = 0; k < Nzones; k++)
    {
        char *filename = get_zone_file_name(sd->path2project, k);

        int type = determine_zone_type(filename);

        if (type > 1 && get_output_flag_background_() == 0)
        {
            fprintf(stderr, "\nerror: allocate_and_setup_zones: background is not set!");
            exit(1);
        }

        switch (type)
        {
        case PLASMA_MODEL_VACUUM:

            zones[k] = hmedium_zone_create_((intptr_t)sd, (intptr_t)bp, wd, sd->path2project, k);
            zone_read_settings_(zones[k], filename);
            break;

        case PLASMA_MODEL_MEDIUM:

            zones[k] = hmedium_zone_create_((intptr_t)sd, (intptr_t)bp, wd, sd->path2project, k);
            zone_read_settings_(zones[k], filename);
            break;

        case PLASMA_MODEL_IMHD:

            zones[k] = imhd_zone_create_((intptr_t)sd, (intptr_t)bp, wd, sd->path2project, k);
            zone_read_settings_(zones[k], filename);
            break;

        case PLASMA_MODEL_RMHD:

            fprintf(stdout, "\nThe plasma model for zone %d is not implemented.", k);
            exit(1);

        case PLASMA_MODEL_FLRE:

            zones[k] = flre_zone_create_((intptr_t)sd, (intptr_t)bp, wd, sd->path2project, k);
            zone_read_settings_(zones[k], filename);
            break;

        default:
            fprintf(stdout, "\nThe plasma model for zone %d is unknown.", k);
            exit(1);
        }

        delete[] filename;
    }

    check_zones_parameters();
}

/*******************************************************************/

void mode_data::calc_all_mode_data(int flag)
{
    calc_basis_fields_in_zones(flag);

    if (flag)
        return;

    calc_stitching_equations();

    calc_stitching_equations_determinant();

    if (get_output_flag_emfield_() > 1)
        save_mode_det_data();

    solve_stitching_equations();

    calc_superposition_of_basis_fields_in_zones();

    space_out_fields_in_zones();

    combine_final_wave_fields();

    setup_wave_fields_for_interpolation();

    if (get_output_flag_emfield_() > 1)
        save_final_wave_fields();

    if (DEBUG_FLAG)
    {
        calc_and_save_divEB();
    }

    if (get_output_flag_additional_() > 0)
    {
        calc_quants_in_zones();

        if (get_output_flag_additional_() > 1)
            save_quants_in_zones();

        combine_final_quants();

        if (get_output_flag_additional_() > 1)
            save_final_quants();
    }
}

/*******************************************************************/

void mode_data::calc_basis_fields_in_zones(int flag)
{
    for (int iz = 0; iz < Nzones; iz++)
    {
        zone_calc_basis_fields_(zones[iz], flag);
    }
}

/*****************************************************************************/

void mode_data::calc_dispersion_in_zones(void)
{
    for (int iz = 0; iz < Nzones; iz++)
    {
        zone_calc_dispersion_(zones[iz]);
    }
}

/*****************************************************************************/

void mode_data::save_dispersion_in_zones(void)
{
    for (int iz = 0; iz < Nzones; iz++)
    {
        zone_save_dispersion_(zones[iz]);
    }
}

/*****************************************************************************/

void mode_data::calc_stitching_equations(void)
{
    Nc = 0;
    for (int iz = 0; iz < Nzones; iz++)
    {
        Nc += zone_get_dim_of_basis_(zones[iz]);
    }

    if (DEBUG_FLAG)
    {
        fprintf(stdout, "\ncalc_stitching_equations: Nc = %d\n", Nc);
    }

    A = new double[2 * Nc * 2 * Nc]; // complex system matrix
    B = new double[2 * Nc];          // complex system rhs vector

    double *M = new double[2 * Nc * 2 * Nc]; // complex matrix (for each boundary)
    double *J = new double[2 * Nc];          // complex rhs vector (for each boundary)

    int neq = 1;  // number of equations
    int ieq = 0;  // equation index
    int nvar = 1; // number of variables
    int ivar = 0; // variable index

    // set A, B by zeros:
    update_system_matrix_and_rhs_vector_(&Nc, A, B, &neq, &ieq, &nvar, &ivar, M, J);

    // initial indices:
    ieq = 1;
    ivar = 1;

    // first boundary (usually, but not necessary, it is a center):
    int iz = 0;

    // dimensions of basis and a basis vector:
    int Nw = zone_get_dim_of_basis_(zones[iz]);
    int len = zone_get_dim_of_basis_vector_(zones[iz]);
    intptr_t code = zones[iz];

    // pointers to basis solutions at the boundary:
    double *EB = zone_get_basis_at_left_boundary_(zones[iz]);

    if (zone_get_medium_(zones[iz]) == PLASMA_MODEL_VACUUM || zone_get_medium_(zones[iz]) == PLASMA_MODEL_MEDIUM)
    {
        switch (zone_get_bc1_(zones[iz]))
        {
        case BOUNDARY_CENTER:

            center_equations_hommed_(&Nw, &len, &code, EB, &neq, &nvar, M, J);
            break;

        case BOUNDARY_IDEALWALL:

            ideal_wall_equations_hommed_(&Nw, &len, &code, EB, &neq, &nvar, M, J);
            break;

        default:

            fprintf(stderr, "\nerror: calc_stitching_equations: first boundary is unknown.");
            exit(1);
        }
    }
    else if (zone_get_medium_(zones[iz]) == PLASMA_MODEL_IMHD)
    {
        switch (zone_get_bc1_(zones[iz]))
        {
        case BOUNDARY_CENTER:

            center_equations_imhd_(&Nw, &len, &code, EB, &neq, &nvar, M, J);
            break;

        case BOUNDARY_IDEALWALL:

            ideal_wall_equations_imhd_(&Nw, &len, &code, EB, &neq, &nvar, M, J);
            break;

        default:

            fprintf(stderr, "\nerror: calc_stitching_equations: first boundary is unknown.");
            exit(1);
        }
    }
    else if (zone_get_medium_(zones[iz]) == PLASMA_MODEL_RMHD)
    {
        switch (zone_get_bc1_(zones[iz]))
        {
        case BOUNDARY_CENTER:

            // center_equations_rmhd_ (&Nw, &len, &code, EB, &neq, &nvar, M, J);
            // break;
            fprintf(stderr, "\nerror: calc_stitching_equations: not implemented.");
            exit(1);

        case BOUNDARY_IDEALWALL:

            // ideal_wall_equations_rmhd_ (&Nw, &len, &code, EB, &neq, &nvar, M, J);
            // break;
            fprintf(stderr, "\nerror: calc_stitching_equations: not implemented.");
            exit(1);

        default:

            fprintf(stderr, "\nerror: calc_stitching_equations: first boundary is unknown.");
            exit(1);
        }
    }
    else if (zone_get_medium_(zones[iz]) == PLASMA_MODEL_FLRE)
    {
        switch (zone_get_bc1_(zones[iz]))
        {
        case BOUNDARY_CENTER:

            center_equations_flre_(&Nw, &len, &code, EB, &neq, &nvar, M, J);
            break;

        case BOUNDARY_IDEALWALL:

            ideal_wall_equations_flre_(&Nw, &len, &code, EB, &neq, &nvar, M, J);
            break;

        default:

            fprintf(stderr, "\nerror: calc_stitching_equations: first boundary is unknown.");
            exit(1);
        }
    }
    else
    {
        fprintf(stderr, "\nerror: calc_stitching_equations: unknown plasma model.");
        exit(1);
    }

    // updates equations matrix and rhs vector:
    update_system_matrix_and_rhs_vector_(&Nc, A, B, &neq, &ieq, &nvar, &ivar, M, J);

    // updates indices:
    ieq += neq;
    ivar += 0; // the same zone and variables, so that + 0.

    // fprintf (stdout, "\niz = %d: ieq = %d ivar = %d", iz, ieq, ivar);

    // internal boundaries:
    for (iz = 0; iz < Nzones - 1; iz++)
    {
        int Nw1 = zone_get_dim_of_basis_(zones[iz + 0]);
        int Nw2 = zone_get_dim_of_basis_(zones[iz + 1]);

        int len1 = zone_get_dim_of_basis_vector_(zones[iz + 0]);
        int len2 = zone_get_dim_of_basis_vector_(zones[iz + 1]);

        intptr_t code1 = zones[iz + 0];
        intptr_t code2 = zones[iz + 1];

        double *EB1 = zone_get_basis_at_right_boundary_(zones[iz + 0]);
        double *EB2 = zone_get_basis_at_left_boundary_(zones[iz + 1]);

        int flg_ant;
        if (zone_get_bc2_(zones[iz]) == BOUNDARY_ANTENNA)
            flg_ant = 1;
        else
            flg_ant = 0;

        if ((zone_get_medium_(zones[iz + 0]) == PLASMA_MODEL_VACUUM ||
             zone_get_medium_(zones[iz + 0]) == PLASMA_MODEL_MEDIUM) &&
            (zone_get_medium_(zones[iz + 1]) == PLASMA_MODEL_VACUUM ||
             zone_get_medium_(zones[iz + 1]) == PLASMA_MODEL_MEDIUM))
        {
            stitching_equations_hommed_hommed_(&Nw1, &len1, &code1, EB1,
                                               &Nw2, &len2, &code2, EB2,
                                               &flg_ant, &neq, &nvar, M, J);
        }
        else if ((zone_get_medium_(zones[iz + 0]) == PLASMA_MODEL_IMHD) &&
                 (zone_get_medium_(zones[iz + 1]) == PLASMA_MODEL_VACUUM ||
                  zone_get_medium_(zones[iz + 1]) == PLASMA_MODEL_MEDIUM))
        {
            stitching_equations_imhd_hommed_(&Nw1, &len1, &code1, EB1,
                                             &Nw2, &len2, &code2, EB2,
                                             &flg_ant, &neq, &nvar, M, J);
        }
        else if (zone_get_medium_(zones[iz + 0]) == PLASMA_MODEL_IMHD &&
                 zone_get_medium_(zones[iz + 1]) == PLASMA_MODEL_IMHD)
        {
            stitching_equations_imhd_imhd_(&Nw1, &len1, &code1, EB1,
                                           &Nw2, &len2, &code2, EB2,
                                           &flg_ant, &neq, &nvar, M, J);
        }
        else if ((zone_get_medium_(zones[iz + 0]) == PLASMA_MODEL_RMHD) &&
                 (zone_get_medium_(zones[iz + 1]) == PLASMA_MODEL_VACUUM ||
                  zone_get_medium_(zones[iz + 1]) == PLASMA_MODEL_MEDIUM))
        {
            // stitching_equations_rmhd_hommed_ (&Nw1, &len1, &code1, EB1, &Nw2, &len2, &code2,
            //                                   EB2, &flg_ant, &neq, &nvar, M, J);
            fprintf(stderr, "\nerror: calc_stitching_equations: not implemented.");
            exit(1);
        }
        else if (zone_get_medium_(zones[iz + 0]) == PLASMA_MODEL_RMHD &&
                 zone_get_medium_(zones[iz + 1]) == PLASMA_MODEL_RMHD)
        {
            // stitching_equations_rmhd_rmhd_ (&Nw1, &len1, &code1, EB1, &Nw2, &len2, &code2,
            //                                 EB2, &flg_ant, &neq, &nvar, M, J);
            fprintf(stderr, "\nerror: calc_stitching_equations: not implemented.");
            exit(1);
        }
        else if ((zone_get_medium_(zones[iz + 0]) == PLASMA_MODEL_FLRE) &&
                 (zone_get_medium_(zones[iz + 1]) == PLASMA_MODEL_VACUUM ||
                  zone_get_medium_(zones[iz + 1]) == PLASMA_MODEL_MEDIUM))
        {
            stitching_equations_flre_hommed_(&Nw1, &len1, &code1, EB1,
                                             &Nw2, &len2, &code2, EB2,
                                             &flg_ant, &neq, &nvar, M, J);
        }
        else if (zone_get_medium_(zones[iz + 0]) == PLASMA_MODEL_FLRE &&
                 zone_get_medium_(zones[iz + 1]) == PLASMA_MODEL_FLRE)
        {
            stitching_equations_flre_flre_(&Nw1, &len1, &code1, EB1,
                                           &Nw2, &len2, &code2, EB2,
                                           &flg_ant, &neq, &nvar, M, J);
        }
        else
        {
            fprintf(stderr, "\nerror: calc_stitching_equations: unknown type of boundary.");
            exit(1);
        }

        // updates equations matrix and rhs vector:
        update_system_matrix_and_rhs_vector_(&Nc, A, B, &neq, &ieq, &nvar, &ivar, M, J);

        // updates indices:
        ieq += neq;
        ivar += Nw1;

        // fprintf (stdout, "\niz = %d: ieq = %d ivar = %d", iz, ieq, ivar);
    }

    // last boundary (usually, but not necessary it is "infinity"):
    iz = Nzones - 1;

    // dimensions of basis and a basis vector:
    Nw = zone_get_dim_of_basis_(zones[iz]);
    len = zone_get_dim_of_basis_vector_(zones[iz]);
    code = zones[iz];

    // pointers to basis solutions at the boundary:
    EB = zone_get_basis_at_right_boundary_(zones[iz]);

    if (zone_get_medium_(zones[iz]) == PLASMA_MODEL_VACUUM || zone_get_medium_(zones[iz]) == PLASMA_MODEL_MEDIUM)
    {
        switch (zone_get_bc2_(zones[iz]))
        {
        case BOUNDARY_INFINITY:

            infinity_equations_hommed_(&Nw, &len, &code, EB, &neq, &nvar, M, J);
            break;

        case BOUNDARY_IDEALWALL:

            ideal_wall_equations_hommed_(&Nw, &len, &code, EB, &neq, &nvar, M, J);
            break;

        default:

            fprintf(stderr, "\nerror: calc_stitching_equations: last boundary is unknown.");
            exit(1);
        }
    }
    else if (zone_get_medium_(zones[iz]) == PLASMA_MODEL_IMHD)
    {
        switch (zone_get_bc2_(zones[iz]))
        {
        case BOUNDARY_INFINITY:

            infinity_equations_imhd_(&Nw, &len, &code, EB, &neq, &nvar, M, J);
            break;

        case BOUNDARY_IDEALWALL:

            ideal_wall_equations_imhd_(&Nw, &len, &code, EB, &neq, &nvar, M, J);
            break;

        default:

            fprintf(stderr, "\nerror: calc_stitching_equations: last boundary is unknown.");
            exit(1);
        }
    }
    else if (zone_get_medium_(zones[iz]) == PLASMA_MODEL_RMHD)
    {
        switch (zone_get_bc2_(zones[iz]))
        {
        case BOUNDARY_INFINITY:

            // infinity_equations_rmhd_ (&Nw, &len, &code, EB, &neq, &nvar, M, J);
            // break;
            fprintf(stderr, "\nerror: calc_stitching_equations: not implemented.");
            exit(1);

        case BOUNDARY_IDEALWALL:

            // ideal_wall_equations_rmhd_ (&Nw, &len, &code, EB, &neq, &nvar, M, J);
            // break;
            fprintf(stderr, "\nerror: calc_stitching_equations: not implemented.");
            exit(1);

        default:

            fprintf(stderr, "\nerror: calc_stitching_equations: last boundary is unknown.");
            exit(1);
        }
    }
    else if (zone_get_medium_(zones[iz]) == PLASMA_MODEL_FLRE)
    {
        switch (zone_get_bc2_(zones[iz]))
        {
        case BOUNDARY_INFINITY:

            infinity_equations_flre_(&Nw, &len, &code, EB, &neq, &nvar, M, J);
            break;

        case BOUNDARY_IDEALWALL:

            ideal_wall_equations_flre_(&Nw, &len, &code, EB, &neq, &nvar, M, J);
            break;

        default:

            fprintf(stderr, "\nerror: calc_stitching_equations: last boundary is unknown.");
            exit(1);
        }
    }
    else
    {
        fprintf(stderr, "\nerror: calc_stitching_equations: unknown plasma model.");
        exit(1);
    }

    // updates equations matrix and rhs vector:
    update_system_matrix_and_rhs_vector_(&Nc, A, B, &neq, &ieq, &nvar, &ivar, M, J);

    // updates indices:
    ieq += neq;
    ivar += Nw;

    // fprintf (stdout, "\niz = %d: ieq = %d ivar = %d", iz, ieq, ivar);

    // check for indices:
    if (ieq - 1 != Nc || ivar - 1 != Nc)
    {
        fprintf(stderr, "\nerror: calc_stitching_equations: wrong final indices: ieq = %d, ivar =%d", ieq, ivar);
        exit(1);
    }

    // clean up:
    delete[] M;
    delete[] J;
}

/*****************************************************************************/

void mode_data::calc_stitching_equations_determinant(void)
{
    double det[2]; // fortran complex number

    calc_system_determinant_(&Nc, A, det);

    set_det_in_wd_struct_ (&wd, &det[0], &det[1]);
}

/*****************************************************************************/

void mode_data::solve_stitching_equations(void)
{
    S = new double[2 * Nc]; // complex coefficients (solution)

    find_superposition_coeffs_(&Nc, A, B, S);

    if (DEBUG_FLAG)
    {
        fprintf(stdout, "\ndeterminat = %.20le %+.20lei\n", wave_data_get_det_re_(wd), wave_data_get_det_im_(wd));
        fprintf(stdout, "\ncheck for superposition coefficients below:");
        for (int i = 0; i < Nc; i++)
        {
            fprintf(stdout, "\ni = %d: S = %.20le %+.20lei", i, S[2 * i + 0], S[2 * i + 1]);
        }
    }
}

/*****************************************************************************/

void mode_data::calc_superposition_of_basis_fields_in_zones(void)
{
    int ind = 0;

    for (int iz = 0; iz < Nzones; iz++)
    {
        zone_calc_superposition_of_basis_fields_(zones[iz], &S[ind]);

        ind += 2 * zone_get_dim_of_basis_(zones[iz]);
    }
}

/*****************************************************************************/

void mode_data::space_out_fields_in_zones(void)
{
    for (int iz = 0; iz < Nzones; iz++)
    {
        zone_calc_final_fields_(zones[iz]);

        if (DEBUG_FLAG)
            zone_save_final_fields_(zones[iz], path2linear);
    }
}

/*****************************************************************************/

void mode_data::combine_final_wave_fields(void)
{
    dim = 0;
    for (int iz = 0; iz < Nzones; iz++)
    {
        dim += zone_get_radial_grid_dimension_(zones[iz]);
    }

    r = new double[dim];

    EB = new double[dim * 2 * 6];

    index = new int[Nzones];

    int ind = 0;

    for (int iz = 0; iz < Nzones; iz++)
    {
        index[iz] = ind;

        zone_copy_radial_grid_(zones[iz], r + ind);

        zone_copy_E_and_B_fields_(zones[iz], EB + ind * 12);

        ind += zone_get_radial_grid_dimension_(zones[iz]);
    }
}

/*****************************************************************************/

void mode_data::setup_wave_fields_for_interpolation(void)
{
    EB_int = new double[dim * 2 * 6];

    for (int node = 0; node < dim; node++)
    {
        for (int comp = 0; comp < 6; comp++)
        {
            for (int part = 0; part < 2; part++)
            {
                EB_int[iFint(node, comp, part)] = EB[iFFM(node, comp, part)];
            }
        }
    }
}

/*****************************************************************************/

void mode_data::save_final_wave_fields(void)
{
    char *fname = new char[1024];
    FILE *out = NULL;

    sprintf(fname, "%sEB.dat", path2linear);

    if (!(out = fopen(fname, "w")))
    {
        fprintf(stderr, "\nFailed to open file %s\a\n", fname);
    }

    for (int i = 0; i < dim; i++)
    {
        fprintf(out, "%.16le", r[i]);

        for (int comp = 0; comp < 6; comp++)
        {
            fprintf(out, "\t%.16le\t%.16le", EB[iFFM(i, comp, 0)], EB[iFFM(i, comp, 1)]);
        }
        fprintf(out, "\n");
    }

    fclose(out);

    delete[] fname;
}

/*****************************************************************************/

void mode_data::calc_quants_in_zones(void)
{
    for (int iz = 0; iz < Nzones; iz++)
    {
        zone_calc_all_quants_(zones[iz]);
    }
}

/*****************************************************************************/

void mode_data::save_quants_in_zones(void)
{
    for (int iz = 0; iz < Nzones; iz++)
    {
        zone_save_all_quants_(zones[iz]);
    }
}

/*****************************************************************************/

void mode_data::combine_final_quants(void)
{
}

/*****************************************************************************/

void mode_data::save_final_quants(void)
{
}

/*****************************************************************************/

inline int mode_data::determine_zone_index_for_point(double x)
{
    // determines zone index for the x value:
    int ind = -1;

    for (int iz = 0; iz < Nzones; iz++)
    {
        if (x >= zone_get_r1_(zones[iz]) && x <= zone_get_r2_(zones[iz]))
        {
            ind = iz;
            break;
        }
    }

    if (ind == -1)
    {
        fprintf(stderr, "\nwarning: determine_zone_index_for_point: x is outside the grid: x = %le", x);
        return 0;
    }

    return ind;
}

/*****************************************************************************/

void mode_data::divEB(double x, complex<double> *div)
{
    int ind = determine_zone_index_for_point(x); // determines zone index for the x value

    int D = zone_get_radial_grid_dimension_(zones[ind]);

    int deg = 5; // degree of the interpolating polynom

    int i = 0;

    complex<double> EB[6], dEB[6];

    double *xg = r + index[ind]; // radial grid for the zone

    for (int comp = 0; comp < 6; comp++)
    {
        double *yg;

        double re[2], im[2];

        // real part:
        yg = &EB_int[iFint(index[ind], comp, 0)]; // EB grid for the zone

        eval_neville_polynom(D, xg, yg, deg, x, 0, 1, &i, re);

        // imag part:
        yg = &EB_int[iFint(index[ind], comp, 1)]; // EB grid for the zone

        eval_neville_polynom(D, xg, yg, deg, x, 0, 1, &i, im);

        // set value:
        EB[comp] = re[0] + I * im[0];
        dEB[comp] = re[1] + I * im[1];
    }

    double kt = (wave_data_get_m_(wd)) / x, kz = (wave_data_get_n_(wd)) / (get_background_rtor_());

    div[0] = EB[0] / x + dEB[0] + I * kt * EB[1] + I * kz * EB[2];
    div[1] = EB[3] / x + dEB[3] + I * kt * EB[4] + I * kz * EB[5];
}

/*******************************************************************/

void mode_data::calc_and_save_divEB(void)
{
    char *fname = new char[1024];
    FILE *out = NULL;

    sprintf(fname, "%sdivEB.dat", path2linear);

    if (!(out = fopen(fname, "w")))
    {
        fprintf(stderr, "\nFailed to open file %s\a\n", fname);
    }

    delete[] fname;

    complex<double> div[2];

    for (int i = 0; i < dim; i++)
    {
        divEB(r[i], div);

        fprintf(out, "%.16le", r[i]);

        fprintf(out, "\t%.16le\t%.16le\t%.16le\t%.16le", real(div[0]), imag(div[0]),
                real(div[1]), imag(div[1]));

        fprintf(out, "\n");
    }

    fclose(out);
}

/*******************************************************************/

void mode_data::eval_EB_fields(double x, complex<double> *EB)
{
    int ind = determine_zone_index_for_point(x); // determines zone index for the x value

    int D = zone_get_radial_grid_dimension_(zones[ind]);

    int deg = 5; // degree of the interpolating polynom

    int i = 0;

    double *xg = r + index[ind]; // radial grid for the zone

    for (int comp = 0; comp < 6; comp++)
    {
        double *yg;

        double re, im;

        // real part:
        yg = &EB_int[iFint(index[ind], comp, 0)]; // EB grid for the zone

        eval_neville_polynom(D, xg, yg, deg, x, 0, 0, &i, &re);

        // imag part:
        yg = &EB_int[iFint(index[ind], comp, 1)]; // EB grid for the zone

        eval_neville_polynom(D, xg, yg, deg, x, 0, 0, &i, &im);

        // set value:
        EB[comp] = re + I * im;
    }
}

/*******************************************************************/

void mode_data::eval_diss_power_density(double x, int type, int spec, double *dpd)
{
    int ind = determine_zone_index_for_point(x); // determines zone index for the x value

    zone_eval_diss_power_density_(zones[ind], x, type, spec, dpd);
}

/*******************************************************************/

void mode_data::eval_current_density(double x, int type, int spec, int comp, double *J)
{
    int ind = determine_zone_index_for_point(x); // determines zone index for the x value

    zone_eval_current_density_(zones[ind], x, type, spec, comp, J);
}

/*******************************************************************/
