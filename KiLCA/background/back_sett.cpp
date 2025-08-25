/*! \file back_sett.cpp
    \brief The implementation of back_sett class.
*/

#include <stdlib.h>
#include <stdio.h>
#include <cmath>
#include <cstring>
#include <climits>
#include <ctime>

#include "constants.h"
#include "back_sett.h"
#include "inout.h"
#include "shared.h"

/*-----------------------------------------------------------------*/

void back_sett::read_settings (char *path)
{
char *file_set = new char[1024];

sprintf (file_set, "%s/background.in", path);

FILE *in;

if ((in=fopen (file_set, "r"))==NULL)
{
    fprintf(stderr, "\nerror: read_background_settings: failed to open file %s\a\n", file_set);
    exit(0);
}

//reading: check for consistence with file!
char *str_buf = new char[1024]; //str buffer

//Machine settings:
read_line_2skip_it (in, &str_buf);
read_line_2get_double (in, &(rtor));
read_line_2get_double (in, &(rp));
read_line_2get_double (in, &(B0));
read_line_2skip_it (in, &str_buf);

//Backround field and plasma settings:
read_line_2skip_it (in, &str_buf);
path2profiles = new char[1024];
read_line_2get_string (in, &(path2profiles));
read_line_2get_int (in, &(calc_back));
flag_back = new char[8];
read_line_2get_string (in, &(flag_back));
read_line_2get_int (in, &(N));
read_line_2get_double (in, &(V_gal_sys));
read_line_2get_double (in, &(V_scale));
read_line_2get_double (in, &(m_i));
read_line_2get_double (in, &(zele));
read_line_2get_double (in, &(zion));
read_line_2skip_it (in, &str_buf);

//Checkings setting:
read_line_2skip_it (in, &str_buf);
read_line_2get_int (in, &(flag_debug));
read_line_2skip_it (in, &str_buf);

fclose (in);

delete [] str_buf;
delete [] file_set;

//Particles settings:
mass = new double[2];
charge = new double[2];

mass[0] = (m_i)*m_p;    /*ions mass*/
mass[1] = m_e;              /*electrons mass*/

charge[0] = e;  //ions charge
charge[1] = -e; //electrons charge

huge_factor = 1.0e20;

if (flag_debug) print_settings ();
}

/*-----------------------------------------------------------------*/

void back_sett::print_settings (void)
{
fprintf(stdout, "\nCheck for background parameters below:\n");

//Machine settings:
fprintf (stdout, "\ntorus big radius: %lg cm", rtor);
fprintf (stdout, "\nplasma radius: %lg cm", rp);
fprintf (stdout, "\ntoroidal magnetic field at the center: %lg G", B0);

//Backround field and plasma settings:
fprintf (stdout, "\npath to background profiles: %s", path2profiles);
fprintf (stdout, "\nflag if recalculate background: %d", calc_back);
fprintf (stdout, "\nflag for background: %s", flag_back);
fprintf (stdout, "\nsplines degree: %d", N);
fprintf (stdout, "\nvelocity of the moving frame: %lg cm/s", V_gal_sys);
fprintf (stdout, "\nscale factor for the Vz velocity profile: %lg", V_scale);
fprintf (stdout, "\nions mass in units of proton mass: %lg", m_i);
fprintf (stdout, "\ncollision coefficient for electrons: %lg", zele);
fprintf (stdout, "\ncollision coefficient for ions: %lg", zion);
fprintf (stdout, "\nflag for debugging mode: %d", flag_debug);

fprintf(stdout, "\n");
}

/*-----------------------------------------------------------------*/

void set_background_settings_c_ (back_sett **bsett, double *rtor, double *rp, double *B0, char *flag_back, double *V_gal_sys, double *V_scale, double *zele, double *zion, int *flag_debug)
{
back_sett *bs = (back_sett *)(*bsett);

*rtor = bs->rtor;
*rp   = bs->rp;
*B0   = bs->B0;

*flag_back = bs->flag_back[0];

*V_gal_sys = bs->V_gal_sys;
*V_scale = bs->V_scale;

*zele = bs->zele;
*zion = bs->zion;

*flag_debug = bs->flag_debug;
}

/*-----------------------------------------------------------------*/

void set_particles_settings_c_ (back_sett **bsett, double *mass, double *charge)
{
back_sett *bs = (back_sett *)(*bsett);

mass[0] = bs->mass[0];
mass[1] = bs->mass[1];

charge[0] = bs->charge[0];
charge[1] = bs->charge[1];
}

/*-----------------------------------------------------------------*/

void set_huge_factor_c_ (back_sett **bsett, double *fac)
{
back_sett *bs = (back_sett *)(*bsett);

*fac = bs->huge_factor;
}

/*-----------------------------------------------------------------*/
