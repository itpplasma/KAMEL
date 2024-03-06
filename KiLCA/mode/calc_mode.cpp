/*! \file
    \brief The implementation of mode_data class and other functions declared in mode.h.
*/

#include <gsl/gsl_errno.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_roots.h>

#include <stdlib.h>
#include <stdio.h>
#include <cmath>
#include <cstring>

#include "constants.h"
#include "mode.h"
#include "eval_back.h"
#include "zone.h"
#include "stitching.h"
#include "interp.h"
#include "transforms.h"

#include "hmedium_zone.h"
#include "imhd_zone.h"
#include "flre_zone.h"

/*****************************************************************************/

double qminusq0 (double x, void *p)
{
func_params *P = (func_params *) p;

return q(x, P->bp) - P->q_res;
}

/*****************************************************************************/

int mode_data::find_resonance_location (void)
{
double q_res = - ((double)(wd->m))/((double)(wd->n));

double r1 = bp->x[0], r2 = bp->x[bp->dimx-1];

if ((q(r1, bp)-q_res)*(q(r2, bp)-q_res) > 0)
{
    if (DEBUG_FLAG)
    {
        fprintf (stdout, "\nfind_resonance_location: resonant surface for the mode m=%d n=%d is absent", wd->m, wd->n);
    }
    wd->r_res = 0.0e0;
    return 0;
}

int status;
int iter = 0, max_iter = 100;
const gsl_root_fsolver_type *T;
gsl_root_fsolver *s;

struct func_params params = {bp, q_res};

gsl_function F;
F.function = &qminusq0;
F.params = &params;

T = gsl_root_fsolver_brent;
s = gsl_root_fsolver_alloc (T);

gsl_root_fsolver_set (s, &F, r1, r2);

do
{
    iter++;

    status = gsl_root_fsolver_iterate (s);

    r1 = gsl_root_fsolver_x_lower (s);
    r2 = gsl_root_fsolver_x_upper (s);

    status = gsl_root_test_interval (r1, r2, 1.0e-8, 1.0e-8);

} while (status == GSL_CONTINUE && iter < max_iter);

if (status != GSL_SUCCESS)
{
    fprintf (stdout, "\nfind_resonance_location: failed to find resonant surface for the mode m=%d n=%d", wd->m, wd->n);
    wd->r_res = 0.0e0;
    return 0;
}

wd->r_res = gsl_root_fsolver_root (s);

gsl_root_fsolver_free (s);

if (DEBUG_FLAG)
{
    fprintf(stdout, "\nresonant surface for the mode m=%d n=%d is found at:\nr=%.16le,  q(r)=%.16le\n", wd->m, wd->n, wd->r_res, q(wd->r_res, bp));
}

return 1;
}

/*****************************************************************************/

void mode_data::check_zones_parameters (void)
{
for (int k=0; k<Nzones-1; k++)
{
    if (zones[k]->r2 != zones[k+1]->r1 || zones[k]->bc2 != zones[k+1]->bc1)
    {
        fprintf (stdout, "\ncheck_zones: boundaries are different: k = %d", k);
        exit (1);
    }

    if (!(zones[k]->bc2 == BOUNDARY_INTERFACE || zones[k]->bc2 == BOUNDARY_ANTENNA))
    {
        fprintf (stdout, "\ncheck_zones: improper type of the right boundary for zone %d: %d", k, zones[k]->bc2);
        exit (1);
    }
}

int iz;

//check that the first boundary is either center or ideal wall:
iz = 0;
if (!(zones[iz]->bc1 == BOUNDARY_CENTER || zones[iz]->bc1 == BOUNDARY_IDEALWALL))
{
    fprintf (stdout, "\ncheck_zones: improper type of the first boundary: %d", zones[iz]->bc1);
    exit (1);
}

//check that the last boundary is either infinity or ideal wall:
iz = Nzones-1;
if (!(zones[iz]->bc2 == BOUNDARY_INFINITY || zones[iz]->bc2 == BOUNDARY_IDEALWALL))
{
    fprintf (stdout, "\ncheck_zones: improper type of the last boundary: %d", zones[iz]->bc2);
    exit (1);
}

if (DEBUG_FLAG)
{
    fprintf (stdout, "\nzones consistency check passed...\n");
}
}

/*******************************************************************/

void mode_data::allocate_and_setup_zones (void)
{
Nzones = determine_mumber_of_zones (sd->path2project);

if (Nzones < 2)
{
    fprintf (stderr, "\nallocate_and_setup_zones: Nzones must be >= 2: %d", Nzones);
    exit (1);
}

zones = new zone * [Nzones]; //array of pointers

for (int k=0; k<Nzones; k++)
{
    char * filename = get_zone_file_name (sd->path2project, k);

    int type = determine_zone_type (filename);

    if (type > 1 && sd->os->flag_background == 0)
    {
        fprintf (stderr, "\nerror: allocate_and_setup_zones: background is not set!");
        exit (1);
    }

    switch (type)
    {
        case PLASMA_MODEL_VACUUM:

        zones[k] = new hmedium_zone (sd, bp, (const wave_data *)wd, sd->path2project, k);
        zones[k]->read_settings (filename);
        break;

        case PLASMA_MODEL_MEDIUM:

        zones[k] = new hmedium_zone (sd, bp, (const wave_data *)wd, sd->path2project, k);
        zones[k]->read_settings (filename);
        break;

        case PLASMA_MODEL_IMHD:

        zones[k] = new imhd_zone (sd, bp, (const wave_data *)wd, sd->path2project, k);
        zones[k]->read_settings (filename);
        break;

        case PLASMA_MODEL_RMHD:

        fprintf (stdout, "\nThe plasma model for zone %d is not implemented.", k);
        exit (1);

        case PLASMA_MODEL_FLRE:

        zones[k] = new flre_zone (sd, bp, (const wave_data *)wd, sd->path2project, k);
        zones[k]->read_settings (filename);
        break;

        default:
        fprintf (stdout, "\nThe plasma model for zone %d is unknown.", k);
        exit (1);
    }

    delete [] filename;
}

check_zones_parameters ();
}

/*******************************************************************/

void mode_data::calc_all_mode_data (int flag)
{
calc_basis_fields_in_zones (flag);

if (flag) return;

calc_stitching_equations ();

calc_stitching_equations_determinant ();

if (sd->os->flag_emfield > 1) save_mode_det_data ();

solve_stitching_equations ();

calc_superposition_of_basis_fields_in_zones ();

space_out_fields_in_zones ();

combine_final_wave_fields ();

setup_wave_fields_for_interpolation ();

if (sd->os->flag_emfield > 1) save_final_wave_fields ();

if (DEBUG_FLAG)
{
    calc_and_save_divEB ();
}

if (sd->os->flag_additional > 0)
{
    calc_quants_in_zones ();

    if (sd->os->flag_additional > 1) save_quants_in_zones ();

    combine_final_quants ();

    if (sd->os->flag_additional > 1) save_final_quants ();
}
}

/*******************************************************************/

void mode_data::calc_basis_fields_in_zones (int flag)
{
for (int iz=0; iz<Nzones; iz++)
{
    zones[iz]->calc_basis_fields (flag);
}
}

/*****************************************************************************/

void mode_data::calc_dispersion_in_zones (void)
{
for (int iz=0; iz<Nzones; iz++)
{
    zones[iz]->calc_dispersion ();
}
}

/*****************************************************************************/

void mode_data::save_dispersion_in_zones (void)
{
for (int iz=0; iz<Nzones; iz++)
{
    zones[iz]->save_dispersion ();
}
}

/*****************************************************************************/

void mode_data::calc_stitching_equations (void)
{
Nc = 0;
for (int iz=0; iz<Nzones; iz++)
{
    Nc += zones[iz]->get_dim_of_basis ();
}

if (DEBUG_FLAG)
{
    fprintf (stdout, "\ncalc_stitching_equations: Nc = %d\n", Nc);
}

A = new double[2*Nc*2*Nc];  //complex system matrix
B = new double[2*Nc];       //complex system rhs vector

double *M = new double[2*Nc*2*Nc];  //complex matrix (for each boundary)
double *J = new double[2*Nc];       //complex rhs vector (for each boundary)

int neq = 1;  //number of equations
int ieq = 0;  //equation index
int nvar = 1; //number of variables
int ivar = 0; //variable index

//set A, B by zeros:
update_system_matrix_and_rhs_vector_ (&Nc, A, B, &neq, &ieq, &nvar, &ivar, M, J);

//initial indices:
ieq = 1;
ivar = 1;

//first boundary (usually, but not necessary, it is a center):
int iz = 0;

//dimensions of basis and a basis vector:
int Nw = zones[iz]->get_dim_of_basis ();
int len = zones[iz]->get_dim_of_basis_vector ();
zone * code = zones[iz];

//pointers to basis solutions at the boundary:
double *EB = zones[iz]->get_basis_at_left_boundary ();

if (zones[iz]->medium == PLASMA_MODEL_VACUUM || zones[iz]->medium == PLASMA_MODEL_MEDIUM)
{
    switch (zones[iz]->bc1)
    {
        case BOUNDARY_CENTER:

        center_equations_hommed_ (&Nw, &len, (hmedium_zone **)(&code), EB, &neq, &nvar, M, J);
        break;

        case BOUNDARY_IDEALWALL:

        ideal_wall_equations_hommed_ (&Nw, &len, (hmedium_zone **)(&code), EB, &neq, &nvar, M, J);
        break;

        default:

        fprintf (stderr, "\nerror: calc_stitching_equations: first boundary is unknown.");
        exit (1);
    }
}
else if (zones[iz]->medium == PLASMA_MODEL_IMHD)
{
    switch (zones[iz]->bc1)
    {
        case BOUNDARY_CENTER:

        center_equations_imhd_ (&Nw, &len, (imhd_zone **)(&code), EB, &neq, &nvar, M, J);
        break;

        case BOUNDARY_IDEALWALL:

        ideal_wall_equations_imhd_ (&Nw, &len, (imhd_zone **)(&code), EB, &neq, &nvar, M, J);
        break;

        default:

        fprintf (stderr, "\nerror: calc_stitching_equations: first boundary is unknown.");
        exit (1);
    }
}
else if (zones[iz]->medium == PLASMA_MODEL_RMHD)
{
    switch (zones[iz]->bc1)
    {
        case BOUNDARY_CENTER:

        //center_equations_rmhd_ (&Nw, &len, &code, EB, &neq, &nvar, M, J);
        //break;
        fprintf (stderr, "\nerror: calc_stitching_equations: not implemented.");
        exit (1);

        case BOUNDARY_IDEALWALL:

        //ideal_wall_equations_rmhd_ (&Nw, &len, &code, EB, &neq, &nvar, M, J);
        //break;
        fprintf (stderr, "\nerror: calc_stitching_equations: not implemented.");
        exit (1);

        default:

        fprintf (stderr, "\nerror: calc_stitching_equations: first boundary is unknown.");
        exit (1);
    }
}
else if (zones[iz]->medium == PLASMA_MODEL_FLRE)
{
    switch (zones[iz]->bc1)
    {
        case BOUNDARY_CENTER:

        center_equations_flre_ (&Nw, &len, (flre_zone **)(&code), EB, &neq, &nvar, M, J);
        break;

        case BOUNDARY_IDEALWALL:

        ideal_wall_equations_flre_ (&Nw, &len, (flre_zone **)(&code), EB, &neq, &nvar, M, J);
        break;

        default:

        fprintf (stderr, "\nerror: calc_stitching_equations: first boundary is unknown.");
        exit (1);
    }
}
else
{
    fprintf (stderr, "\nerror: calc_stitching_equations: unknown plasma model.");
    exit (1);
}

//updates equations matrix and rhs vector:
update_system_matrix_and_rhs_vector_ (&Nc, A, B, &neq, &ieq, &nvar, &ivar, M, J);

//updates indices:
ieq += neq;
ivar += 0; //the same zone and variables, so that + 0.

//fprintf (stdout, "\niz = %d: ieq = %d ivar = %d", iz, ieq, ivar);

//internal boundaries:
for (iz=0; iz<Nzones-1; iz++)
{
    int Nw1 = zones[iz+0]->get_dim_of_basis ();
    int Nw2 = zones[iz+1]->get_dim_of_basis ();

    int len1 = zones[iz+0]->get_dim_of_basis_vector ();
    int len2 = zones[iz+1]->get_dim_of_basis_vector ();

    zone * code1 = zones[iz+0];
    zone * code2 = zones[iz+1];

    double *EB1 = zones[iz+0]->get_basis_at_right_boundary ();
    double *EB2 = zones[iz+1]->get_basis_at_left_boundary ();

    int flg_ant;
    if (zones[iz]->bc2 == BOUNDARY_ANTENNA) flg_ant = 1;
    else                                    flg_ant = 0;

    if ((zones[iz+0]->medium == PLASMA_MODEL_VACUUM ||
         zones[iz+0]->medium == PLASMA_MODEL_MEDIUM) &&
        (zones[iz+1]->medium == PLASMA_MODEL_VACUUM ||
         zones[iz+1]->medium == PLASMA_MODEL_MEDIUM))
    {
        stitching_equations_hommed_hommed_ (&Nw1, &len1, (hmedium_zone **)(&code1), EB1,
					    &Nw2, &len2, (hmedium_zone **)(&code2), EB2,
					    &flg_ant, &neq, &nvar, M, J);
    }
    else if ((zones[iz+0]->medium == PLASMA_MODEL_IMHD) &&
             (zones[iz+1]->medium == PLASMA_MODEL_VACUUM ||
              zones[iz+1]->medium == PLASMA_MODEL_MEDIUM))
    {
        stitching_equations_imhd_hommed_ (&Nw1, &len1, (imhd_zone **)(&code1), EB1,
					  &Nw2, &len2, (hmedium_zone **)(&code2), EB2,
					  &flg_ant, &neq, &nvar, M, J);
    }
    else if (zones[iz+0]->medium == PLASMA_MODEL_IMHD &&
             zones[iz+1]->medium == PLASMA_MODEL_IMHD)
    {
        stitching_equations_imhd_imhd_ (&Nw1, &len1, (imhd_zone **)(&code1), EB1,
					&Nw2, &len2, (imhd_zone **)(&code2), EB2,
					&flg_ant, &neq, &nvar, M, J);
    }
    else if ((zones[iz+0]->medium == PLASMA_MODEL_RMHD) &&
             (zones[iz+1]->medium == PLASMA_MODEL_VACUUM ||
              zones[iz+1]->medium == PLASMA_MODEL_MEDIUM))
    {
        //stitching_equations_rmhd_hommed_ (&Nw1, &len1, &code1, EB1, &Nw2, &len2, &code2,
        //                                  EB2, &flg_ant, &neq, &nvar, M, J);
        fprintf (stderr, "\nerror: calc_stitching_equations: not implemented.");
        exit (1);
    }
    else if (zones[iz+0]->medium == PLASMA_MODEL_RMHD &&
             zones[iz+1]->medium == PLASMA_MODEL_RMHD)
    {
        //stitching_equations_rmhd_rmhd_ (&Nw1, &len1, &code1, EB1, &Nw2, &len2, &code2,
        //                                EB2, &flg_ant, &neq, &nvar, M, J);
        fprintf (stderr, "\nerror: calc_stitching_equations: not implemented.");
        exit (1);
    }
    else if ((zones[iz+0]->medium == PLASMA_MODEL_FLRE) &&
             (zones[iz+1]->medium == PLASMA_MODEL_VACUUM ||
              zones[iz+1]->medium == PLASMA_MODEL_MEDIUM))
    {
        stitching_equations_flre_hommed_ (&Nw1, &len1, (flre_zone **)(&code1), EB1,
					  &Nw2, &len2, (hmedium_zone **)(&code2), EB2,
					  &flg_ant, &neq, &nvar, M, J);
    }
    else if (zones[iz+0]->medium == PLASMA_MODEL_FLRE &&
             zones[iz+1]->medium == PLASMA_MODEL_FLRE)
    {
        stitching_equations_flre_flre_ (&Nw1, &len1, (flre_zone **)(&code1), EB1,
					&Nw2, &len2, (flre_zone **)(&code2), EB2,
					&flg_ant, &neq, &nvar, M, J);
    }
    else
    {
        fprintf (stderr, "\nerror: calc_stitching_equations: unknown type of boundary.");
        exit (1);
    }

    //updates equations matrix and rhs vector:
    update_system_matrix_and_rhs_vector_ (&Nc, A, B, &neq, &ieq, &nvar, &ivar, M, J);

    //updates indices:
    ieq += neq;
    ivar += Nw1;

    //fprintf (stdout, "\niz = %d: ieq = %d ivar = %d", iz, ieq, ivar);
}

//last boundary (usually, but not necessary it is "infinity"):
iz = Nzones-1;

//dimensions of basis and a basis vector:
Nw = zones[iz]->get_dim_of_basis ();
len = zones[iz]->get_dim_of_basis_vector ();
code = zones[iz];

//pointers to basis solutions at the boundary:
EB = zones[iz]->get_basis_at_right_boundary ();

if (zones[iz]->medium == PLASMA_MODEL_VACUUM || zones[iz]->medium == PLASMA_MODEL_MEDIUM)
{
    switch (zones[iz]->bc2)
    {
        case BOUNDARY_INFINITY:

        infinity_equations_hommed_ (&Nw, &len, (hmedium_zone **)(&code), EB, &neq, &nvar, M, J);
        break;

        case BOUNDARY_IDEALWALL:

        ideal_wall_equations_hommed_ (&Nw, &len, (hmedium_zone **)(&code), EB, &neq, &nvar, M, J);
        break;

        default:

        fprintf (stderr, "\nerror: calc_stitching_equations: last boundary is unknown.");
        exit (1);
    }
}
else if (zones[iz]->medium == PLASMA_MODEL_IMHD)
{
    switch (zones[iz]->bc2)
    {
        case BOUNDARY_INFINITY:

        infinity_equations_imhd_ (&Nw, &len, (imhd_zone **)(&code), EB, &neq, &nvar, M, J);
        break;

        case BOUNDARY_IDEALWALL:

        ideal_wall_equations_imhd_ (&Nw, &len, (imhd_zone **)(&code), EB, &neq, &nvar, M, J);
        break;

        default:

        fprintf (stderr, "\nerror: calc_stitching_equations: last boundary is unknown.");
        exit (1);
    }
}
else if (zones[iz]->medium == PLASMA_MODEL_RMHD)
{
    switch (zones[iz]->bc2)
    {
        case BOUNDARY_INFINITY:

        //infinity_equations_rmhd_ (&Nw, &len, &code, EB, &neq, &nvar, M, J);
        //break;
        fprintf (stderr, "\nerror: calc_stitching_equations: not implemented.");
        exit (1);

        case BOUNDARY_IDEALWALL:

        //ideal_wall_equations_rmhd_ (&Nw, &len, &code, EB, &neq, &nvar, M, J);
        //break;
        fprintf (stderr, "\nerror: calc_stitching_equations: not implemented.");
        exit (1);

        default:

        fprintf (stderr, "\nerror: calc_stitching_equations: last boundary is unknown.");
        exit (1);
    }
}
else if (zones[iz]->medium == PLASMA_MODEL_FLRE)
{
    switch (zones[iz]->bc2)
    {
        case BOUNDARY_INFINITY:

        infinity_equations_flre_ (&Nw, &len, (flre_zone **)(&code), EB, &neq, &nvar, M, J);
        break;

        case BOUNDARY_IDEALWALL:

        ideal_wall_equations_flre_ (&Nw, &len, (flre_zone **)(&code), EB, &neq, &nvar, M, J);
        break;

        default:

        fprintf (stderr, "\nerror: calc_stitching_equations: last boundary is unknown.");
        exit (1);
    }
}
else
{
    fprintf (stderr, "\nerror: calc_stitching_equations: unknown plasma model.");
    exit (1);
}

//updates equations matrix and rhs vector:
update_system_matrix_and_rhs_vector_ (&Nc, A, B, &neq, &ieq, &nvar, &ivar, M, J);

//updates indices:
ieq += neq;
ivar += Nw;

//fprintf (stdout, "\niz = %d: ieq = %d ivar = %d", iz, ieq, ivar);

//check for indices:
if (ieq-1 != Nc || ivar-1 != Nc)
{
    fprintf (stderr, "\nerror: calc_stitching_equations: wrong final indices: ieq = %d, ivar =%d", ieq, ivar);
    exit (1);
}

//clean up:
delete [] M;
delete [] J;
}

/*****************************************************************************/

void mode_data::calc_stitching_equations_determinant (void)
{
double det[2]; //fortran complex number

calc_system_determinant_ (&Nc, A, det);

wd->det = det[0] + I*det[1];
}

/*****************************************************************************/

void mode_data::solve_stitching_equations (void)
{
S = new double[2*Nc]; //complex coefficients (solution)

find_superposition_coeffs_ (&Nc, A, B, S);

if (DEBUG_FLAG)
{
    fprintf (stdout, "\ndeterminat = %.20le %+.20lei\n", real(wd->det), imag(wd->det));
    fprintf (stdout, "\ncheck for superposition coefficients below:");
    for (int i=0; i<Nc; i++)
    {
        fprintf (stdout, "\ni = %d: S = %.20le %+.20lei", i, S[2*i+0], S[2*i+1]);
    }
}
}

/*****************************************************************************/

void mode_data::calc_superposition_of_basis_fields_in_zones (void)
{
int ind = 0;

for (int iz=0; iz<Nzones; iz++)
{
    zones[iz]->calc_superposition_of_basis_fields (&S[ind]);

    ind += 2*zones[iz]->get_dim_of_basis ();
}
}

/*****************************************************************************/

void mode_data::space_out_fields_in_zones (void)
{
for (int iz=0; iz<Nzones; iz++)
{
    zones[iz]->calc_final_fields ();

    if (DEBUG_FLAG) zones[iz]->save_final_fields (path2linear);
}
}

/*****************************************************************************/

void mode_data::combine_final_wave_fields (void)
{
dim = 0;
for (int iz=0; iz<Nzones; iz++)
{
    dim += zones[iz]->get_radial_grid_dimension ();
}

r = new double[dim];

EB = new double[dim*2*6];

index = new int [Nzones];

int ind = 0;

for (int iz=0; iz<Nzones; iz++)
{
    index[iz] = ind;

    zones[iz]->copy_radial_grid (r+ind);

    zones[iz]->copy_E_and_B_fields (EB + ind*12);

    ind += zones[iz]->get_radial_grid_dimension ();
}
}

/*****************************************************************************/

void mode_data::setup_wave_fields_for_interpolation (void)
{
EB_int = new double[dim*2*6];

for (int node=0; node<dim; node++)
{
    for (int comp=0; comp<6; comp++)
    {
        for (int part=0; part<2; part++)
        {
            EB_int[iFint(node, comp, part)] = EB[iFFM(node, comp, part)];
        }
    }
}
}

/*****************************************************************************/

void mode_data::save_final_wave_fields (void)
{
char *fname = new char[1024];
FILE *out = NULL;

sprintf (fname, "%sEB.dat", path2linear);

if (!(out = fopen (fname, "w")))
{
    fprintf (stderr, "\nFailed to open file %s\a\n", fname);
}

for (int i=0; i<dim; i++)
{
    fprintf (out, "%.16le", r[i]);

    for (int comp=0; comp<6; comp++)
    {
       fprintf (out, "\t%.16le\t%.16le", EB[iFFM(i, comp, 0)], EB[iFFM(i, comp, 1)]);
    }
     fprintf (out, "\n");
}

fclose (out);

delete [] fname;
}

/*****************************************************************************/

void mode_data::calc_quants_in_zones (void)
{
for (int iz=0; iz<Nzones; iz++)
{
    zones[iz]->calc_all_quants ();
}
}

/*****************************************************************************/

void mode_data::save_quants_in_zones (void)
{
for (int iz=0; iz<Nzones; iz++)
{
    zones[iz]->save_all_quants ();
}
}

/*****************************************************************************/

void mode_data::combine_final_quants (void)
{

}

/*****************************************************************************/

void mode_data::save_final_quants (void)
{

}

/*****************************************************************************/

inline int mode_data::determine_zone_index_for_point (double x)
{
//determines zone index for the x value:
int ind = -1;

for (int iz=0; iz<Nzones; iz++)
{
    if (x >= zones[iz]->get_r1 () && x <= zones[iz]->get_r2 ())
    {
        ind = iz;
        break;
    }
}

if (ind == -1)
{
    fprintf (stderr, "\nwarning: determine_zone_index_for_point: x is outside the grid: x = %le", x);
    return 0;
}

return ind;
}

/*****************************************************************************/

void mode_data::divEB (double x, complex<double> *div)
{
int ind = determine_zone_index_for_point (x); //determines zone index for the x value

int D = zones[ind]->get_radial_grid_dimension ();

int deg = 5; //degree of the interpolating polynom

int i = 0;

complex<double> EB[6], dEB[6];

double *xg = r + index[ind]; //radial grid for the zone

for (int comp=0; comp<6; comp++)
{
    double *yg;

    double re[2], im[2];

    //real part:
    yg = &EB_int[iFint(index[ind], comp, 0)]; //EB grid for the zone

    eval_neville_polynom (D, xg, yg, deg, x, 0, 1, &i, re);

    //imag part:
    yg = &EB_int[iFint(index[ind], comp, 1)]; //EB grid for the zone

    eval_neville_polynom (D, xg, yg, deg, x, 0, 1, &i, im);

    //set value:
    EB[comp]  = re[0] + I*im[0];
    dEB[comp] = re[1] + I*im[1];
}

double kt = (wd->m)/x, kz = (wd->n)/(sd->bs->rtor);

div[0] = EB[0]/x + dEB[0] + I*kt*EB[1] + I*kz*EB[2];
div[1] = EB[3]/x + dEB[3] + I*kt*EB[4] + I*kz*EB[5];
}

/*******************************************************************/

void mode_data::calc_and_save_divEB (void)
{
char *fname = new char[1024];
FILE *out = NULL;

sprintf (fname, "%sdivEB.dat", path2linear);

if (!(out = fopen (fname, "w")))
{
    fprintf (stderr, "\nFailed to open file %s\a\n", fname);
}

delete [] fname;

complex<double> div[2];

for (int i=0; i<dim; i++)
{
    divEB (r[i], div);

    fprintf (out, "%.16le", r[i]);

    fprintf (out, "\t%.16le\t%.16le\t%.16le\t%.16le", real(div[0]), imag(div[0]),
                                                      real(div[1]), imag(div[1]));

    fprintf (out, "\n");
}

fclose (out);
}

/*******************************************************************/

void mode_data::eval_EB_fields (double x, complex<double> *EB)
{
int ind = determine_zone_index_for_point (x); //determines zone index for the x value

int D = zones[ind]->get_radial_grid_dimension ();

int deg = 5; //degree of the interpolating polynom

int i = 0;

double *xg = r + index[ind]; //radial grid for the zone

for (int comp=0; comp<6; comp++)
{
    double *yg;

    double re, im;

    //real part:
    yg = &EB_int[iFint(index[ind], comp, 0)]; //EB grid for the zone

    eval_neville_polynom (D, xg, yg, deg, x, 0, 0, &i, &re);

    //imag part:
    yg = &EB_int[iFint(index[ind], comp, 1)]; //EB grid for the zone

    eval_neville_polynom (D, xg, yg, deg, x, 0, 0, &i, &im);

    //set value:
    EB[comp] = re + I*im;
}
}

/*******************************************************************/

void mode_data::eval_diss_power_density (double x, int type, int spec, double *dpd)
{
int ind = determine_zone_index_for_point (x); //determines zone index for the x value

zones[ind]->eval_diss_power_density (x, type, spec, dpd);
}

/*******************************************************************/

void mode_data::eval_current_density (double x, int type, int spec, int comp, double * J)
{
int ind = determine_zone_index_for_point (x); //determines zone index for the x value

zones[ind]->eval_current_density (x, type, spec, comp, J);
}

/*******************************************************************/
