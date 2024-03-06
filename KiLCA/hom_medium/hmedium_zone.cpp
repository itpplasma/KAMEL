/*! \file hmedium_zone.cpp
    \brief The implementation of hmedium_zone class.
*/

#include "hmedium_zone.h"
#include "mode.h"
#include "shared.h"
#include "inout.h"

/*****************************************************************************/

void hmedium_zone::read_settings (char * file)
{
read (file); //base class read function

//derived class specific:
FILE *in;

if ((in=fopen (file, "r"))==NULL)
{
    fprintf(stderr, "\nerror: hmedium_zone: read_settings: failed to open file %s\a\n", file);
    exit(0);
}

char *str_buf = new char[1024];

//skip lines (already readed and values are set):
read_line_2skip_it (in, &str_buf);
read_line_2skip_it (in, &str_buf);
read_line_2skip_it (in, &str_buf);
read_line_2skip_it (in, &str_buf);
read_line_2skip_it (in, &str_buf);
read_line_2skip_it (in, &str_buf);
read_line_2skip_it (in, &str_buf);
read_line_2skip_it (in, &str_buf);

//Medium settings:
read_line_2skip_it (in, &str_buf);
read_line_2get_complex (in, &(sigma));
read_line_2skip_it (in, &str_buf);

//Solution settings:
read_line_2skip_it (in, &str_buf);
read_line_2get_int (in, &(max_dim));
read_line_2get_double (in, &(eps_rel));
read_line_2get_double (in, &(eps_abs));
read_line_2skip_it (in, &str_buf);

//Space out settings:
read_line_2skip_it (in, &str_buf);
read_line_2get_int (in, &(deg));
read_line_2get_double (in, &(reps));
read_line_2get_double (in, &(aeps));
read_line_2get_double (in, &(step));
read_line_2skip_it (in, &str_buf);

//Debugging settings:
read_line_2skip_it (in, &str_buf);
read_line_2get_int (in, &(flag_debug));
read_line_2skip_it (in, &str_buf);

fclose (in);

Nwaves = 4;
Ncomps = 6;

if (flag_debug) print_settings ();

delete [] str_buf;
}

/*****************************************************************************/

void hmedium_zone::print_settings (void)
{
print (); //base class print function

fprintf (stdout, "\nmedium conductivity: (%le, %le)", real(sigma), imag(sigma));

fprintf(stdout, "\n");
}

/*****************************************************************************/

void hmedium_zone::calc_basis_fields (int flag)
{
dim = max_dim;

r = new double[dim];

basis = new double[dim*Nwaves*Ncomps*2];

int m = wd->m;
double kz = (wd->n)/(sd->bs->rtor);

double comega[2] = {real(wd->olab),  imag(wd->olab)}; //lab frame

double csigma[2] = {real(sigma), imag(sigma)}; //lab frame

//makes homogenious radial grid:
double dr = (r2-r1)/(dim-1);

for (int i=0; i<dim; i++)
{
    r[i] = r1 + dr*i;

    eval_basis_in_hom_media_ (&m, &kz, comega, csigma, &r[i], &basis[ib(i, 0, 0, 0)]);
}

//last point (to have exactly r2):
int i = dim-1;

r[i] = r2;

eval_basis_in_hom_media_ (&m, &kz, comega, csigma, &r[i], &basis[ib(i, 0, 0, 0)]);

//transform basis to the lab frame if needed!

//normalize and scale basis functions if needed!
}

/*****************************************************************************/

void hmedium_zone::copy_E_and_B_fields (double *EB_p)
{
for (int node=0; node<dim; node++)
{
    for (int comp=0; comp<6; comp++) //Er, Et, Ez, Br, Bt, Bz
    {
        EB_p[iFFM(node, comp, 0)] = EB[iF(node, comp, 0)];
        EB_p[iFFM(node, comp, 1)] = EB[iF(node, comp, 1)];
    }
}
}

/*****************************************************************************/
