/*! \file antenna.cpp
    \brief The implementation of antenna class.
*/

#include <stdlib.h>
#include <stdio.h>
#include <cmath>
#include <cstring>

#include "antenna.h"
#include "inout.h"

/*-----------------------------------------------------------------*/

void antenna::read_settings (char *path)
{
/*! \fn read_settings (char *path)
    \brief

    \param path path to top level run dir of a project.
*/

char *file_set = new char[1024];

sprintf (file_set, "%s/antenna.in", path);

FILE *in;

if ((in=fopen (file_set, "r"))==NULL)
{
    fprintf(stderr, "\nerror: read_antenna_settings: failed to open file %s\a\n", file_set);
    exit(0);
}

//reading: check for consistence with the input file!
char *str_buf = new char[1024]; //str buffer

//Antenna settings:
read_line_2skip_it (in, &str_buf);
read_line_2get_double (in, &(ra));
read_line_2get_double (in, &(wa));
read_line_2get_double (in, &(I0));
read_line_2get_complex (in, &(flab));
read_line_2get_int (in, &(dma));
read_line_2get_int (in, &(flag_debug));
read_line_2get_int (in, &(flag_eigmode));
read_line_2skip_it (in, &str_buf);

fclose (in);

int i;

//Reading modes.in file:
sprintf (file_set, "%s/modes.in", path);

if ((in=fopen (file_set, "r"))==NULL)
{
    fprintf(stderr, "\nerror: read_antenna_settings: failed to open file %s\a\n", file_set);
    exit(0);
}

modes = new int[2*dma]; //(m,n)
for (i=0; i<dma; i++)
{
    if (feof(in) || ferror(in)) break;
    fscanf (in, "(%d, %d)\n", &(modes[2*i]), &(modes[2*i+1]));
}

if(i != dma)
{
    fprintf(stderr, "\nerror: read_antenna_settings: read error or false modes array dimension!\a\n");
    exit(0);
}

fclose (in);

delete [] file_set;
delete [] str_buf;

if (flag_debug) print_settings ();
}

/*-----------------------------------------------------------------*/

void antenna::print_settings (void)
{
fprintf(stdout, "\nCheck for antenna parameters below:\n");

//Antenna settings:
fprintf (stdout, "\nantenna radius: %lg cm", ra);
fprintf (stdout, "\nantenna current layer width: %lg cm", wa);
fprintf (stdout, "\nantenna coils current: %lg statamps", I0);
fprintf (stdout, "\nantenna lab frequency: (%lg, %lg) 1/s", real(flab), imag(flab));
fprintf (stdout, "\ndimension of modes array: %d", dma);
fprintf (stdout, "\nflag for debugging mode: %d", flag_debug);
fprintf (stdout, "\nflag for eigmode search: %d", flag_eigmode);

fprintf (stdout, "\narray of mode numbers (m, n): ");

int i;
for (i=0; i<dma; i++)
{
    fprintf (stdout, "(%d, %d) ", modes[2*i], modes[2*i+1]);
}

fprintf(stdout, "\n");
}

/*-----------------------------------------------------------------*/

void set_antenna_settings_c_ (antenna **ptr, double *ra, double *wa, double *I0, double *flab_re, double *flab_im, int *dma, int *flag_debug)
{
/*! \fn set_antenna_settings_c_ (antenna **ptr, double *ra, double *wa, double *I0, double *flab_re, double *flab_im, int *dma, int *flag_debug)
    \brief The function copies class variables to the parameters which are in fact members of the antenna module. The parameter ptr is an address of an address of the antenna object passed.

    \param
*/

antenna *ap = *ptr;

*ra = ap->ra;
*wa = ap->wa;
*I0 = ap->I0;
*flab_re = real(ap->flab);
*flab_im = imag(ap->flab);
*dma = ap->dma;
*flag_debug = ap->flag_debug;
}

/*-----------------------------------------------------------------*/
