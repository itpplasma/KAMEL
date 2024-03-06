/*! \file output_sett.h
    \brief The implementation of output_sett class.
*/

#include <stdlib.h>
#include <stdio.h>
#include <cmath>
#include <cstring>

#include "output_sett.h"
#include "shared.h"
#include "inout.h"
#include "constants.h"

/*****************************************************************************/

void output_sett::read_settings (char *path)
{
char *file_set = new char[1024];

sprintf (file_set, "%s/output.in", path);

FILE *in;

if ((in=fopen (file_set, "r"))==NULL)
{
    fprintf(stderr, "\nerror: read_settings: failed to open file %s\a\n", file_set);
    exit(0);
}

//reading: check for consistence with file!
char *str_buf = new char[1024]; //str buffer

//Run settings:
read_line_2skip_it (in, &str_buf);
read_line_2get_int (in, &(flag_background));
read_line_2get_int (in, &(flag_emfield));
read_line_2get_int (in, &(flag_additional));
read_line_2get_int (in, &(flag_dispersion));
read_line_2skip_it (in, &str_buf);

read_line_2skip_it (in, &str_buf);
read_line_2get_int (in, &(flag_debug));
read_line_2skip_it (in, &str_buf);

num_quants = 8;

flag_quants = new int[num_quants];

read_line_2skip_it (in, &str_buf);

int i;
for (i=0; i<num_quants; i++)
{
    if (feof(in) || ferror(in)) break;
    read_line_2get_int (in, &(flag_quants[i]));
}

if(i != num_quants)
{
    fprintf(stderr, "\nerror: read_settings: read error or false flag_quants array dimension!\a\n");
    exit(0);
}

fclose (in);

delete [] str_buf;
delete [] file_set;

if (flag_debug) print_settings ();
}

/*****************************************************************************/

void output_sett::print_settings ()
{
fprintf(stdout, "\nCheck for output parameters below:\n");

//Run settings:
fprintf (stdout, "\nif compute background data: %d", flag_background);
fprintf (stdout, "\nif compute linear data: %d", flag_emfield);
fprintf (stdout, "\nif compute additional quants: %d", flag_additional);
fprintf (stdout, "\nif compute dispersion: %d", flag_dispersion);
fprintf (stdout, "\nflag for debugging: %d", flag_debug);

fprintf (stdout, "\ndimension of flags array for quantities: %d", num_quants);

fprintf (stdout, "\narray of flags for additional quantities: ");

int i;

for (i=0; i<num_quants; i++)
{
    fprintf (stdout, "%d ", flag_quants[i]);
}

fprintf(stdout, "\n");
}

/*****************************************************************************/
