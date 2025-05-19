/*! \file background.cpp
    \brief The implementation of background class.
*/

#include <stdlib.h>
#include <stdio.h>
#include <cmath>
#include <cstring>
#include <climits>
#include <ctime>

#include "background.h"
#include "back_sett.h"
#include "calc_back.h"
#include "inout.h"
#include "spline.h"

#include "shared.h"
#include "constants.h"
#include "conduct_parameters.h"
#include "wave_code_interface.h"

/*-----------------------------------------------------------------*/

background::background (const settings_t *s)
{
sd = s;

set_profiles_indices ();

//for moments computation:
set_back_aliases_in_conductivity_parameters_ ();

path2background = new char[1024];

sprintf (path2background, "%s%s", sd->path2project, "background-data/");

if (sd->os->flag_background > 1)
{
    //make dir for background data:
    char *sys_command = new char[1024];

    sprintf (sys_command, "%s%s", "mkdir -p ", path2background);

    if (system (sys_command)==-1)
    {
        fprintf (stderr, "\nerror: set_background_profiles: system()!");
    }

    delete [] sys_command;
}
}

/*-----------------------------------------------------------------*/

background::~background (void)
{
for (int k=0; k<Nprofiles; k++) delete [] profile_names[k];

delete [] profile_names;
delete [] path2background;

delete [] i_T;
delete [] i_Vth;
delete [] i_Vz;

delete [] i_n_p;
delete [] i_Vp_p;
delete [] i_Vt_p;
delete [] i_nu;

delete [] i_J0th;
delete [] i_J0z;

if (x) delete [] x;
if (y) delete [] y;

if (C) delete [] C;
if (R) delete [] R;

if (sid) spline_free_ (sid);
}

/*-----------------------------------------------------------------*/

void background::clean_background_arrays (void)
{
if (x) delete [] x; x = NULL;
if (y) delete [] y; y = NULL;

if (C) delete [] C; C = NULL;
if (R) delete [] R; R = NULL;

spline_free_ (sid); sid = 0;
}

/*-----------------------------------------------------------------*/

void background::set_profiles_indices (void)
{
if      (abs(sd->bs->calc_back) == 1 || abs(sd->bs->calc_back) == 3) Nprofiles = 7;  //q, n, Ti, Te, Vth, Vz, Er
else if (abs(sd->bs->calc_back) == 2)                                Nprofiles = 11; //q, n, Ti, Te, Vth, Vz, Er, Bth, Bz, Jth, Jz
else
{
    fprintf (stderr, "\nwarning: set_profiles_indices: unknown flag!\n");
    exit(1);
}

profile_names = new char *[Nprofiles];

int k;

for (k=0; k<Nprofiles; k++) profile_names[k] = new char[8];

if (abs(sd->bs->calc_back) == 1 || abs(sd->bs->calc_back) == 3)
{
    strcpy (profile_names[0], "q");
    strcpy (profile_names[1], "n");
    strcpy (profile_names[2], "Ti");
    strcpy (profile_names[3], "Te");
    strcpy (profile_names[4], "Vth");
    strcpy (profile_names[5], "Vz");
    strcpy (profile_names[6], "Er");
}
else if (abs(sd->bs->calc_back) == 2)
{
    strcpy (profile_names[0], "q");
    strcpy (profile_names[1], "n");
    strcpy (profile_names[2], "Ti");
    strcpy (profile_names[3], "Te");
    strcpy (profile_names[4], "Vth");
    strcpy (profile_names[5], "Vz");
    strcpy (profile_names[6], "Er");
    strcpy (profile_names[7], "Bth");
    strcpy (profile_names[8], "Bz");
    strcpy (profile_names[9], "Jth");
    strcpy (profile_names[10], "Jz");
}

/*background splines indices: the order of indices below is VERY important and offers the best convenience: DO NOT CHANGE, ONLY ADD NEW IF NEEDED!!!*/
k = 0;

i_q = k++;

i_n = k++;

i_T = new int[2]; /*i,e*/
i_T[0] = k++;
i_T[1] = k++;

i_Vth = new int[3]; /*i,e,t*/
i_Vz  = new int[3]; /*i,e,t*/

i_Vth[2] = k++;
i_Vz[2]  = k++;

i_Er = k++;

/*from equlibrium:*/
i_Bth = k++;

i_Bz = k++;

i_B = k++;

i_hth = k++;

i_hz = k++;

i_dPhi0 = k++;

/*from f0 parameters search:*/
i_n_p  = new int[2]; /*i,e*/
i_Vp_p = new int[2]; /*i,e*/
i_Vt_p = new int[2]; /*i,e*/
i_nu   = new int[2]; /*i,e*/

i_n_p[0]  = k++;
i_Vp_p[0] = k++;
i_Vt_p[0] = k++;
i_nu[0]   = k++;

i_n_p[1]  = k++;
i_Vp_p[1] = k++;
i_Vt_p[1] = k++;
i_nu[1]   = k++;

/*from equlibrium:*/
i_J0th = new int[3]; /*i,e,t*/
i_J0z  = new int[3]; /*i,e,t*/

i_J0th[2] = k++;
i_J0z[2]  = k++;

i_J0th[1] = k++;
i_J0z[1]  = k++;

i_J0th[0] = k++;
i_J0z[0]  = k++;

i_Vth[1] = k++;
i_Vz[1]  = k++;

i_Vth[0] = k++;
i_Vz[0]  = k++;

dimy = k; /*total number of background profiles*/
}

/*-----------------------------------------------------------------*/

void background::set_background_profiles_from_files (void)
{
char *file_name = new char[1024];
sprintf (file_name, "%sn.dat", sd->bs->path2profiles);

/*arrays allocation for background profiles:*/
dimx = count_lines_in_file (file_name, 0);

delete [] file_name;

int N = sd->bs->N;

ind = (int)0.5*(dimx);

x = new double[dimx];

y = new double[(dimx)*(dimy)];

C = new double[(N+1)*(dimx)*(dimy)];

R = new double[(N+1)*(dimy)];

spline_alloc_ (N, 1, dimx, x, C, &(sid));

//load and spline input profiles:
int ierr;
ierr = load_input_background_profiles ();
//fprintf (stdout, "\nbackground profiles are loaded...\n");

ierr = spline_input_profiles (); //first 7 profs, see later
//fprintf (stdout, "\nbackground profiles are splined...\n");

if (abs(sd->bs->calc_back) == 1)
{
    ierr = calculate_equilibrium ();
    //fprintf (stdout, "\nequilibrium parameters are computed...\n");
}
else if (abs(sd->bs->calc_back) == 2)
{
    ierr = check_and_spline_equilibrium ();
    //fprintf (stdout, "\ninput equilibrium is checked...\n");
}
else
{
    fprintf (stderr, "\nset_background_profiles_from_files: this feature is not implemented yet!\n");
    exit(1);
}

flag_dPhi0_calc = 1; //to recalculate dPhi0

//to get background electric field never use Er array - it is in the lab frame. Always use dPhi0!

ierr = find_f0_parameters ();
//fprintf (stdout, "\nf0 parameters are computed...\n");

if (sd->os->flag_background > 1)
{
    save_background ();
    //fprintf (stdout, "\nbackground quantities are stored...\n");

    //independent computation of f0 moments for checking:
    eval_and_save_f0_moments ();
    //fprintf (stdout, "\nf0 moments are computed and stored...\n");
}
}

/*-----------------------------------------------------------------*/

void background::set_background_profiles_from_interface (void)
{
get_background_dimension_from_balance_ (&dimx);

ind = (int)0.5*(dimx);

int N = sd->bs->N;

x = new double[dimx];

y = new double[(dimx)*(dimy)];

C = new double[(N+1)*(dimx)*(dimy)];

R = new double[(N+1)*(dimy)];

spline_alloc_ (N, 1, dimx, x, C, &(sid));

double *q   = y + (i_q)*dimx;
double *n   = y + (i_n)*dimx;
double *Ti  = y + (i_T[0])*dimx;
double *Te  = y + (i_T[1])*dimx;
double *Vth = y + (i_Vth[2])*dimx;
double *Vz  = y + (i_Vz[2])*dimx;
double *Er  = y + (i_Er)*dimx;

get_background_profiles_from_balance_ (&dimx, x, q, n, Ti, Te, Vth, Vz, Er);

//transform background profiles to a moving frame
//Er is not transformed, but dPhi0 is transformed later in find_f0_parameters
for (int i=0; i<dimx; i++) Vz[i] -= sd->bs->V_gal_sys;

int ierr;

ierr = spline_input_profiles ();

//to skip calculation of dPhi0 in find_f0_parameters and use dPhi0 = -Er transformed to a moving frame
flag_dPhi0_calc = 0;

if (abs(sd->bs->calc_back) == 1)
{
    ierr = calculate_equilibrium ();
    //fprintf (stdout, "\nequilibrium parameters are computed...\n");
}
else if (abs(sd->bs->calc_back) == 2)
{
    //ierr = check_and_spline_equilibrium ();
    fprintf (stderr, "\nset_background_profiles_from_interface: this feature is not implemented yet!\n");
    exit(1);
}
else if (abs(sd->bs->calc_back) == 3)
{
    flag_dPhi0_calc = 1;
    ierr = calculate_equilibrium ();
    //fprintf (stdout, "\nequilibrium parameters are computed...\n");
}
else
{
    fprintf (stderr, "\nset_background_profiles_from_files: this feature is not implemented yet!\n");
    exit(1);
}

ierr = find_f0_parameters ();

if (sd->os->flag_background > 1)
{
    save_background ();

    eval_and_save_f0_moments (); //independent computation of f0 moments for checking
}
}

/*-----------------------------------------------------------------*/

int background::load_input_background_profiles (void)
{
//all profiles are assumed to be in the lab frame, Vz profile is shifted to a moving frame, but Er not!

int i, j;
double tmp1 = 0, tmp2 = 0;

char *file_name = new char[1024];
char *sys_command = new char[1024];

int ind[11] = {i_q, i_n, i_T[0], i_T[1], i_Vth[2], i_Vz[2], i_Er, i_Bth, i_Bz, i_J0th[2], i_J0z[2]};

/*loading profiles*/
for (i=0; i<Nprofiles; i++)
{
    sprintf (file_name, "%s%s.dat", sd->bs->path2profiles, profile_names[i]);
    load_data_file (file_name, dimx, 1, x, y+ind[i]*dimx);

    if (i==5) //shift (negative) and scale Vz profile:
    {
        for (j=0; j<dimx; j++)
        {
            *(y+ind[i]*dimx+j) *= sd->bs->V_scale;   //scale Vz profile
            *(y+ind[i]*dimx+j) -= sd->bs->V_gal_sys; //shift Vz profile
        }
    }

    /*check for r grid: must be the same for all profiles:*/
    tmp2 = 0.0e0;
    for (j=0; j<dimx; j++) tmp2 += x[j];
    if (tmp2!=tmp1 && i!=0)
    {
        fprintf (stdout, "\nwarning: load_input_profiles: r grids look different: tmp1=%.25le tmp2=%.25le i=%d!\n", tmp1, tmp2, i);
    }
    tmp1 = tmp2;

    if (sd->os->flag_background > 1)
    {
        /*copy input profile to background-data dir*/
        sprintf (sys_command, "%s%s%s%s%s%s", "cp ", file_name, " ", path2background, profile_names[i], "_i.dat");
        system (sys_command);
    }
}

if (x[0]==0.0)
{
    fprintf (stderr, "\nwarning: load_input_profiles: first point of profiles grid must be not zero!\n");
}

delete [] file_name;
delete [] sys_command;



return 0;
}

/*-----------------------------------------------------------------*/

int background::spline_input_profiles (void)
{
/*splines first 7 profiles: q, n, Ti, Te, Vth, Vz (in mov frame), Er (in lab frame)*/
int ierr;

spline_calc_ (sid, y + i_q, i_q, i_Er, NULL, &ierr);

return ierr;
}

/*-----------------------------------------------------------------*/
