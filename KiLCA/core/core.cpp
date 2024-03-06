/*! \file core.cpp
    \brief The implementation of core_data class.
*/

#include <stdlib.h>
#include <stdio.h>
#include <cmath>
#include <cstring>

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

delete sd;

delete bp;

delete [] path2project;

if (mda == 0) return;

for (int ind=0; ind<dim; ind++) if (mda[ind]) delete mda[ind];

delete [] mda;

mda = 0;
}

/*******************************************************************/

void core_data::delete_modes_array (void)
{
if (mda == 0) return;

for (int ind=0; ind<dim; ind++) if (mda[ind]) delete mda[ind];

delete [] mda;

mda = 0;
}

/*******************************************************************/

void core_data::calc_and_set_mode_independent_core_data (void)
{
sd = new settings (path2project);

set_settings_in_core_module_ (&sd);

sd->read_settings ();

bp = new background (sd);

set_background_in_core_module_ (&bp);

//loads and splines initial background profiles, computes equilibrium magnetic field,
//currents, f0 parameters and other stuff:
if (sd->bs->calc_back > 0)
{
    bp->set_background_profiles_from_files ();
}
else if(sd->bs->calc_back < 0)
{
    bp->set_background_profiles_from_interface ();
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
dim = sd->as->dma;
mda = new mode_data * [dim];

complex<double> olab = (2.0*pi)*(sd->as->flab);

//loop over modes array:
for (int ind=0; ind<dim; ind++)
{
    int m = sd->as->modes[2*ind+0];
    int n = sd->as->modes[2*ind+1];

    mda[ind] = new mode_data (m, n, olab, (const settings *)sd, (const background *)bp);

    mda[ind]->calc_all_mode_data ();

    delete mda[ind]; //delete mode_data object (if only needed)
    mda[ind] = 0;

    clear_all_data_in_mode_data_module_ (); //clean up fortran module data
}
}

/*******************************************************************/

void core_data::calc_and_set_mode_dependent_core_data_eigmode (void)
{
//allocates modes array in core struct:
dim = sd->as->dma;
mda = new mode_data * [dim];

//loop over modes array:
for (int ind=0; ind<dim; ind++)
{
    int m = sd->as->modes[2*ind+0];
    int n = sd->as->modes[2*ind+1];

    if (sd->es->search_flag == 1)
    {
        //loop over frequences:
        loop_over_frequences (ind, m, n, this);
    }
    else if (sd->es->search_flag == 0)
    {
        //complex zero search by gsl real solver:
        find_det_zeros (ind, m, n, this);
    }
    else if (sd->es->search_flag == -1)
    {
        //all complex zeros search by complex zero finder:
        find_eigmodes (ind, m, n, this);
    }
    else
    {
        fprintf (stdout, "\nError: unknown search_flag in eigmode options file: %d.", sd->es->search_flag);
        exit(1);
    }
}
}

/*******************************************************************/

void core_data::calc_and_set_mode_dependent_core_data_antenna_interface (void)
{
//allocates modes array in core struct:
dim = sd->as->dma;
mda = new mode_data * [dim];

complex<double> olab = (2.0*pi)*(sd->as->flab);

//loop over modes array:
for (int ind=0; ind<dim; ind++)
{
    int m = sd->as->modes[2*ind+0];
    int n = sd->as->modes[2*ind+1];

    mda[ind] = new mode_data (m, n, olab, (const settings *)sd, (const background *)bp);

    mda[ind]->calc_all_mode_data ();

    clear_all_data_in_mode_data_module_ (); //clean up fortran module data
}
}

/*******************************************************************/

void core_data::calc_and_set_mode_dependent_core_data_antenna_interface (int m, int n, int flag)
{
//allocates modes array in core struct:
dim = 1;
mda = new mode_data * [dim];

complex<double> olab = (2.0*pi)*(sd->as->flab);

//loop over modes array:
for (int ind=0; ind<dim; ind++)
{
    mda[ind] = new mode_data (m, n, olab, (const settings *)sd, (const background *)bp);

    mda[ind]->calc_all_mode_data (flag);

    clear_all_data_in_mode_data_module_ (); //clean up fortran module data
}
}

/*******************************************************************/
