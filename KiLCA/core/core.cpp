/*! \file core.cpp
    \brief The implementation of core_data class.
*/

#include <stdlib.h>
#include <stdio.h>
#include <cmath>
#include <cstring>
#include <iostream>

#include "core.h"
#include "settings.h"
#include "calc_eigmode.h"

/*******************************************************************/

core_data::core_data (char *path)
{
/*! \fn core_data (char *path)
    \brief Constructor.

    \param path path to top level run dir of a project.
*/

//checks if pointer size and Fortran integer size are compatibile

int pp; //fortran integer size to store the C pointers

get_pointer_precision_ (&pp); //gets size of integers used to store C pointers in Fortran

//technical check: important for passing pointers to fortran as integers
if (sizeof(uintptr_t) != sizeof(core_data *) || pp != sizeof(uintptr_t))
{
  fprintf (stdout, "\nwarning: core_data: sizeof(uintptr_t)=%ld != sizeof(ptr)=%ld != pp=%d\n", sizeof(uintptr_t), sizeof(core_data *), pp);
  fprintf (stdout, "\nset apropriate value for integer, parameter :: pp in the file constants_m.f90");
  exit(1);
}

path2project = new char[1024];
strcpy (path2project, path);
}

/*******************************************************************/

core_data::~core_data (void)
{
/*! \fn ~core_data (void)
    \brief Destructor. Frees all allocated memory hierarchically.
*/

// Do NOT delete sd intentionally, as the static object, this member points to, will be reused
// in the next iteration of the time evolution.

// background is now a Fortran singleton (kilca_background_data_m) with no
// heap allocation behind bp's sentinel value - nothing to delete.

delete [] path2project;

if (mda == 0) return;

for (int ind=0; ind<dim; ind++) if (mda[ind]) mode_data_destroy_ (mda[ind]);

delete [] mda;

mda = 0;
}

/*******************************************************************/

void core_data::delete_modes_array (void)
{
if (mda == 0) return;

for (int ind=0; ind<dim; ind++) if (mda[ind]) mode_data_destroy_ (mda[ind]);

delete [] mda;

mda = 0;
}

/*******************************************************************/

void core_data::calc_and_set_mode_independent_core_data (void)
{
    static intptr_t static_settings_vacuum = 0;
    static intptr_t static_settings_flre = 0;

    auto const p2pstr = std::string{path2project};

    if (p2pstr.find("vacuum") != std::string::npos)
    {
        if (static_settings_vacuum == 0)
        {
            static_settings_vacuum = settings_create_ (path2project);
            settings_read_settings_ (static_settings_vacuum);
        }
        sd = static_settings_vacuum;
    }
    else if (p2pstr.find("flre") != std::string::npos)
    {
        if (static_settings_flre == 0)
        {
            static_settings_flre = settings_create_ (path2project);
            settings_read_settings_ (static_settings_flre);
        }
        sd = static_settings_flre;
    }
    else
    {
        std::cerr << "\nError: calc_and_set_mode_independent_core_data: unknown project type in path: " << path2project << '\n';
        exit(1);
    }

    set_settings_in_core_module_(&sd);

char sd_path2project[1024];
settings_get_path2project_ (sd, sd_path2project);

bp = (background *) background_create_ (sd_path2project);

set_background_in_core_module_ (&bp);

//loads and splines initial background profiles, computes equilibrium magnetic field,
//currents, f0 parameters and other stuff:
if (get_background_calc_back_() > 0)
{
    background_set_profiles_from_files_ ();
}
else if(get_background_calc_back_() < 0)
{
    background_set_profiles_from_interface_ ();
}
else
{
    fprintf (stderr, "\nwarning: calc_and_set_mode_independent_core_data: unknown flag in background.in!\n");
    exit(1);
}
}

/*******************************************************************/

void core_data::calc_and_set_mode_dependent_core_data_antenna (void)
{
//allocates modes array in core struct:
dim = get_antenna_dma_ ();
mda = new intptr_t [dim];

double flab_re, flab_im;
get_antenna_flab_ (&flab_re, &flab_im);
complex<double> olab = (2.0*pi)*complex<double>(flab_re, flab_im);

//loop over modes array:
for (int ind=0; ind<dim; ind++)
{
    int m, n;
    get_antenna_mode_ (ind, &m, &n);

    char sd_path2project[1024];
    settings_get_path2project_ (sd, sd_path2project);
    mda[ind] = mode_data_create_ (m, n, real(olab), imag(olab), sd, (intptr_t)bp, sd_path2project);

    mode_data_calc_all_mode_data_ (mda[ind], 0);

    mode_data_destroy_ (mda[ind]); //delete mode_data object (if only needed)
    mda[ind] = 0;

    clear_all_data_in_mode_data_module_ (); //clean up fortran module data
}
}

/*******************************************************************/

void core_data::calc_and_set_mode_dependent_core_data_eigmode (void)
{
//allocates modes array in core struct:
dim = get_antenna_dma_ ();
mda = new intptr_t [dim];

//loop over modes array:
for (int ind=0; ind<dim; ind++)
{
    int m, n;
    get_antenna_mode_ (ind, &m, &n);

    if (get_eigmode_search_flag_ () == 1)
    {
        //loop over frequences:
        loop_over_frequences (ind, m, n, this);
    }
    else if (get_eigmode_search_flag_ () == 0)
    {
        //complex zero search by gsl real solver:
        find_det_zeros (ind, m, n, this);
    }
    else if (get_eigmode_search_flag_ () == -1)
    {
        //all complex zeros search by complex zero finder:
        find_eigmodes (ind, m, n, this);
    }
    else
    {
        fprintf (stdout, "\nError: unknown search_flag in eigmode options file: %d.", get_eigmode_search_flag_ ());
        exit(1);
    }
}
}

/*******************************************************************/

void core_data::calc_and_set_mode_dependent_core_data_antenna_interface (void)
{
    //allocates modes array in core struct:
    dim = get_antenna_dma_ ();
    mda = new intptr_t [dim];

    double flab_re, flab_im;
get_antenna_flab_ (&flab_re, &flab_im);
complex<double> olab = (2.0*pi)*complex<double>(flab_re, flab_im);

    //loop over modes array:
    for (int ind=0; ind<dim; ind++)
    {
        int m, n;
        get_antenna_mode_ (ind, &m, &n);

        char sd_path2project[1024];
    settings_get_path2project_ (sd, sd_path2project);
    mda[ind] = mode_data_create_ (m, n, real(olab), imag(olab), sd, (intptr_t)bp, sd_path2project);

        mode_data_calc_all_mode_data_ (mda[ind], 0);

        clear_all_data_in_mode_data_module_ (); //clean up fortran module data
    }
}

/*******************************************************************/

void core_data::calc_and_set_mode_dependent_core_data_antenna_interface (int m, int n, int flag)
{
//allocates modes array in core struct:
dim = 1;
mda = new intptr_t [dim];

double flab_re, flab_im;
get_antenna_flab_ (&flab_re, &flab_im);
complex<double> olab = (2.0*pi)*complex<double>(flab_re, flab_im);

//loop over modes array:
for (int ind=0; ind<dim; ind++)
{
    char sd_path2project[1024];
    settings_get_path2project_ (sd, sd_path2project);
    mda[ind] = mode_data_create_ (m, n, real(olab), imag(olab), sd, (intptr_t)bp, sd_path2project);

    mode_data_calc_all_mode_data_ (mda[ind], flag);

    clear_all_data_in_mode_data_module_ (); //clean up fortran module data
}
}

/*******************************************************************/
