/*! \file
    \brief The implementation of zone class.
*/

#include <stdlib.h>
#include <stdio.h>
#include <cmath>
#include <cstring>
#include <dirent.h>
#include <fnmatch.h>

#include "zone.h"
#include "shared.h"
#include "inout.h"
#include "constants.h"
#include "mode.h"
#include "stitching.h"

/*****************************************************************************/

zone::zone (const settings *sd_p, const background *bp_p, const wave_data *wd_p, char *path_p, int index_p)
{
sd = sd_p;
bp = bp_p;
wd = wd_p;

path = new char[1024];
strcpy (path, path_p);

index = index_p;

r = 0;
basis = 0;
EB = 0;
S = 0;
}

/*****************************************************************************/

zone::~zone ()
{
delete [] path;

delete [] r;
delete [] basis;
delete [] EB;
delete [] S;
}

/*****************************************************************************/

void zone::read (char *file)
{
FILE *in;

if ((in=fopen (file, "r"))==NULL)
{
    fprintf(stderr, "\nerror: read_settings: failed to open file %s\a\n", file);
    exit(0);
}

//string buffers:
char *str_buf = new char[1024];
char *bstr1 = new char[64];
char *bstr2 = new char[64];
char *mstr = new char[64];

read_line_2skip_it (in, &str_buf);
read_line_2get_double (in, &r1);
read_line_2get_string (in, &bstr1);
read_line_2get_string (in, &mstr);
read_line_2get_int (in, &version);
read_line_2get_string (in, &bstr2);
read_line_2get_double (in, &r2);
read_line_2skip_it (in, &str_buf);

fclose (in);

//Determination of BC type index:
bc1 = -1;
bc2 = -1;
for (int k=0; k<Nbc; k++)
{
    if (!strcmp (bc_str[k], bstr1)) bc1 = k;
    if (!strcmp (bc_str[k], bstr2)) bc2 = k;
}

if (bc1 == -1 || bc2 == -1)
{
    fprintf (stderr, "\nzone::read: BC type is not known: %s %s", bstr1, bstr2);
    exit (1);
}

//Determination of medium type index:
medium = -1;
for (int k=0; k<Nmed; k++)
{
    if (!strcmp (med_str[k], mstr))
    {
        medium = k;
        break;
    }
}

if (medium == -1)
{
    fprintf (stderr, "\nzone::read: medium type is not known: %s", mstr);
    exit (1);
}

delete [] str_buf;
delete [] mstr;
delete [] bstr1;
delete [] bstr2;
}

/*****************************************************************************/

void zone::print (void) const
{
fprintf (stdout, "\nzone index = %d", index);

fprintf (stdout, "\nr1 = %le", r1);
fprintf (stdout, "\nr2 = %le", r2);

fprintf (stdout, "\nbc1 = %d", bc1);
fprintf (stdout, "\nbc2 = %d", bc2);

fprintf (stdout, "\nmedium = %d", medium);
fprintf (stdout, "\nversion = %d", version);

fprintf (stdout, "\nnumber of waves: %d", Nwaves);
fprintf (stdout, "\nnumber of field components: %d", Ncomps);

fprintf (stdout, "\nmax dimension of the radial grid for the solution: %d", max_dim);
fprintf (stdout, "\nrelative accuracy for the solver: %lg", eps_rel);
fprintf (stdout, "\nabsolute accuracy for the solver: %lg", eps_abs);

fprintf (stdout, "\ndegree of the polynomial used to space out the solution: %d", deg);
fprintf (stdout, "\nrelative accuracy of of the sparse solution: %lg", reps);
fprintf (stdout, "\nabsolute accuracy of of the sparse solution: %lg", aeps);
fprintf (stdout, "\nmax grid step in the solution: %lg", step);

fprintf (stdout, "\nflag for debugging mode: %d", flag_debug);
}

/*****************************************************************************/

int determine_zone_type (char *file)
{
FILE *in;

if ((in=fopen (file, "r"))==NULL)
{
    fprintf(stderr, "\nerror: determine_zone_type: failed to open file %s\a\n", file);
    exit(0);
}

//string buffers:
char *str_buf = new char[1024];
char *mstr = new char[64];

read_line_2skip_it (in, &str_buf);
read_line_2skip_it (in, &str_buf);
read_line_2skip_it (in, &str_buf);
read_line_2get_string (in, &mstr);

fclose (in);

//Determination of medium type index:
int medium = -1;

for (int k=0; k<Nmed; k++)
{
    if (!strcmp (med_str[k], mstr))
    {
        medium = k;
        break;
    }
}

if (medium == -1)
{
    fprintf (stderr, "\nerror: determine_zone_type: medium type is unknown: %s", mstr);
    exit (1);
}

delete [] str_buf;
delete [] mstr;

return medium;
}

/*****************************************************************************/

int selector (const struct dirent *ent)
{
if (!fnmatch ("zone_*.in", ent->d_name, 0)) return 1;
else                                        return 0;
}

/*******************************************************************/

int determine_number_of_zones (char *path2project)
{
//determines number of zones: 2 versions: scandir() and readdir()

//struct dirent **all;

//int Nzones = scandir (path2project, &all, selector, alphasort);

int Nzones = 0;

DIR *dp;
struct dirent *ep;

if (dp = opendir (path2project))
{
    while (ep = readdir (dp))
    {
        if (!fnmatch ("zone_*.in", ep->d_name, 0)) Nzones++;
    }
    closedir (dp);
}
else
{
    fprintf(stderr, "\ndetermine_number_of_zones: faled to open the project directory %s", path2project);
    exit(1);
}

if (DEBUG_FLAG)
{
    fprintf (stdout, "\nNzones = %d", Nzones); fflush (stdout);
}

return Nzones;
}

/*******************************************************************/

char * get_zone_file_name (char * path2project, int zone_index)
{
DIR * dp;
struct dirent * ep;

char * file_name = new char[1024];

strcpy (file_name, path2project);

int count = 0;

if (dp = opendir (path2project))
{
  char * file_pattern = new char[1024];

  sprintf (file_pattern, "*zone_%d*.in", zone_index+1);

  while (ep = readdir (dp))
  {
    if (!fnmatch (file_pattern, ep->d_name, 0))
    {
      ++count;
      strcat (file_name, ep->d_name);
      break;
    }
  }
  closedir (dp);
  delete [] file_pattern;
}
else
{
  fprintf (stderr, "\nget_zone_file_name: faled to open the project directory %s", path2project);
  exit (1);
}

if (DEBUG_FLAG)
{
  fprintf (stdout, "\nget_zone_file_name: file name for zone %d is: %s", zone_index, file_name);
  fflush (stdout);
}

if (!count)
{
  fprintf (stderr, "\nget_zone_file_name: failed to find the file name for zone %d!", zone_index);
  delete [] file_name;
  exit (1);
}

return file_name;
}

/*******************************************************************/

void zone::calc_superposition_of_basis_fields (double *S_p)
{
EB = new double[dim*Ncomps*2];

S = new double[Nwaves*2];

for (int i=0; i<Nwaves*2; i++) S[i] = S_p[i];

eval_superposition_of_basis_functions_ (&Ncomps, &Nwaves, &dim, basis, S, EB);
}

/*****************************************************************************/

void zone::copy_radial_grid (double *r_p)
{
for (int i=0; i<dim; i++) r_p[i] = r[i];
}

/*****************************************************************************/

void zone::save_basis_fields (char *path2linear)
{
char *fname = new char[1024];

sprintf (fname, "%sdebug-data/zone_%d_basis.dat", path2linear, index);

FILE *out;
if (!(out = fopen (fname, "w")))
{
    fprintf (stderr, "\nFailed to open file %s\a\n", fname);
}

for (int i=0; i<dim; i++)
{
    fprintf (out, "%.16le", r[i]);

    for (int sol=0; sol<Nwaves; sol++)
    {
        for (int comp=0; comp<Ncomps; comp++)
        {
            fprintf (out, "\t%.16le\t%.16le", basis[ib(i, sol, comp, 0)], basis[ib(i, sol, comp, 1)]);
        }
    }
     fprintf (out, "\n");
}

fclose (out);

delete [] fname;
}

/*****************************************************************************/

void zone::save_final_fields (char *path2linear)
{
char *fname = new char[1024];

sprintf (fname, "%szone_%d_EB.dat", path2linear, index);

FILE *out;
if (!(out = fopen (fname, "w")))
{
    fprintf (stderr, "\nFailed to open file %s\a\n", fname);
}

for (int i=0; i<dim; i++)
{
    fprintf (out, "%.16le", r[i]);

    for (int comp=0; comp<Ncomps; comp++)
    {
       fprintf (out, "\t%.16le\t%.16le", EB[iEB(i, comp, 0)], EB[iEB(i, comp, 1)]);
    }
     fprintf (out, "\n");
}

fclose (out);

delete [] fname;
}

/*****************************************************************************/

void get_right_boundary_of_zone_ (zone **ptr, double *r)
{
zone *Z = (zone *)(*ptr);

*r = Z->get_r2 ();
}

/*****************************************************************************/

void get_left_boundary_of_zone_ (zone **ptr, double *r)
{
zone *Z = (zone *)(*ptr);

*r = Z->get_r1 ();
}

/*****************************************************************************/
