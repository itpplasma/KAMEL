#include "back_sett.h"

#include <climits>
#include <cmath>
#include <cstring>
#include <ctime>
#include <stdio.h>
#include <stdlib.h>

void back_sett::print_settings()
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

strcpy (flag_back, bs->flag_back);

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
