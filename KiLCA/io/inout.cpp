/*! \file
    \brief The definitions of functions declared in inout.h.
*/

#include <stdlib.h>
#include <stdio.h>
#include <cmath>
#include <cstring>

#include "inout.h"
#include "shared.h"
#include "constants.h"

/*******************************************************************/

int save_cmplx_matrix (int Nrows, int Ncols, int Npoints, const double *xgrid, const double *arr, const char *path_name)
{
FILE *outfile;

size_t nchar = 1024;

char *full_name = new char[nchar];

int i, j, k;

for (i=0; i<Ncols; i++)
{
    sprintf (full_name, "%s_%d.dat", path_name, i);
    /*printf ("\nsave_cmplx_matrix: file=%s", full_name);*/

    if (!(outfile = fopen (full_name, "w")))
    {
        fprintf (stderr, "\nFailed to open file %s\a\a\a\n", full_name);
        continue;
    }

    for (k=0; k<Npoints; k++)
    {
        fprintf (outfile, "%.16le", xgrid[k]);
        for (j=0; j<Nrows; j++) fprintf (outfile, "\t%.16le\t%.16le", arr[2*Nrows*Ncols*k+2*Nrows*i+2*j], arr[2*Nrows*Ncols*k+2*Nrows*i+2*j+1]);
        fprintf (outfile, "\n");
    }
    fclose (outfile);
}

delete [] full_name;
return 0;
}

/*******************************************************************/

int save_cmplx_matrix_to_one_file (int Nrows, int Ncols, int Npoints, const double *xgrid, const double *arr, const char *full_name)
{
//stores a complex matrix packed in 1D array (fortran order) to a file in the C order:
//x1, a[1,1], a[1,2],...,a[2,1],...
//x2, a[1,1], a[1,2],...,a[2,1],...
//.................................

FILE *outfile;

int i, j, k;

if (!(outfile = fopen (full_name, "w")))
{
    fprintf (stderr, "\nFailed to open file %s\a\a\a\n", full_name);
    return 1;
}

for (k=0; k<Npoints; k++)
{
    fprintf (outfile, "%.16le", xgrid[k]);
    //columns first:
    for (i=0; i<Nrows; i++)
    {
        for (j=0; j<Ncols; j++)
        {
            fprintf (outfile, "\t%.16le\t%.16le", arr[2*Nrows*Ncols*k+2*Nrows*j+2*i], arr[2*Nrows*Ncols*k+2*Nrows*j+2*i+1]);
        }
    }
    fprintf (outfile, "\n");
}
fclose (outfile);
return 0;
}

/*******************************************************************/

int save_real_matrix_to_one_file (int order, int Nrows, int Ncols, int Npoints, const double *xgrid, const double *arr, const char *full_name)
{
//stores a real matrix packed in 1D array (as defined by order) to a file as:
//x1, a[1,1], a[1,2],...,a[2,1],...
//x2, a[1,1], a[1,2],...,a[2,1],...
//.................................
//order = 0 means fortran order
//order = 1 means C order

FILE *outfile;

int i, j, k;

if (!(outfile = fopen (full_name, "w")))
{
    fprintf (stderr, "\nFailed to open file %s\a\a\a\n", full_name);
    return 1;
}

for (k=0; k<Npoints; k++)
{
    fprintf (outfile, "%.16le", xgrid[k]);
    //saves rows first (C order):
    for (i=0; i<Nrows; i++)
    {
        for (j=0; j<Ncols; j++)
        {
            if (order == 0) fprintf (outfile, "\t%.16le", arr[Nrows*Ncols*k+Nrows*j+i]);
            if (order == 1) fprintf (outfile, "\t%.16le", arr[Nrows*Ncols*k+Ncols*i+j]);
        }
    }
    fprintf (outfile, "\n");
}

fclose (outfile);
return 0;
}

/*******************************************************************/

int save_real_array (int dim, const double *xgrid, const double *arr, const char *full_name)
{
FILE *outfile;

if (!(outfile = fopen (full_name, "w")))
{
    fprintf (stderr, "\nFailed to open file %s\a\a\a\n", full_name);
    return 1;
}

int k;

for (k=0; k<dim; k++) fprintf (outfile, "%.16le\t%.16le\n", xgrid[k], arr[k]);

fclose (outfile);
return 0;
}

/*******************************************************************/

int save_complex_array (int dim, const double *xgrid, const double *arr, const char *full_name)
{
//re, im parts follow each other like in Fortran complex data
FILE *outfile;

if (!(outfile = fopen (full_name, "w")))
{
    fprintf (stderr, "\nFailed to open file %s\a\a\a\n", full_name);
    return 1;
}

int k;

for (k=0; k<dim; k++) fprintf (outfile, "%.16le\t%.16le\t%.16le\n", xgrid[k], arr[2*k], arr[2*k+1]);

fclose (outfile);
return 0;
}

/*******************************************************************/

char * trim (char *str)
{
const char delimiters[4] = {'\t',' ','\n','\0'};
return (char *) strtok (str, delimiters);
}

/*-----------------------------------------------------------------*/

void read_line_2get_double (FILE *in, double *value)
{
/*reads current line and stores value before # as a double value*/
size_t nchar = 1024;
char *char_buf = new char[nchar];
const char delimiters[4] = {'#','\n','\0'};
char *info;

if (getline (&char_buf, &nchar, in) < 2) /*reads line*/
{
   fprintf (stderr, "\nread_line_2get_double: file reading error: %s\n", char_buf);
   fflush (stderr);
}

*value = strtod((char *) strtok (char_buf, delimiters), NULL); /*takes what is before #*/
info = (char *) strtok (NULL, delimiters); /*takes what is after #*/

delete [] char_buf;
}

/*-----------------------------------------------------------------*/

void read_line_2get_complex (FILE *in, complex<double> *value)
{
/*reads current line and stores value before # as a complex<double> value*/
size_t nchar = 1024;
char *char_buf = new char[nchar];
const char delimiters[4] = {'#','\n','\0'};
const char del[3] = {'(',',',')'};
char *info;
double re, im;

if (getline (&char_buf, &nchar, in) < 2) /*reads line*/
{
    fprintf (stderr, "\nread_line_2get_complex: file reading error!\a\n");
    fflush (stderr);
}

info = (char *) strtok (char_buf, delimiters); /*takes what is before #*/
re = strtod((char *) strtok (info, del), NULL);/*takes what is before ,*/
im = strtod((char *) strtok (NULL, del), NULL);/*takes what is before )*/

*value = re + I*im;

/*printf ("\nvalue = (%le, %le)\n", re, im);*/

delete [] char_buf;
}

/*-----------------------------------------------------------------*/

void read_line_2get_int (FILE *in, int *value)
{
/*reads current line and stores value before # as an integer value*/
size_t nchar = 1024;
char *char_buf = new char[nchar];
const char delimiters[4] = {'#','\n','\0'};
char *info;

if (getline (&char_buf, &nchar, in) < 2) /*reads line*/
{
    fprintf (stderr, "\nread_line_2get_int: file reading error!\a\n");
    fflush (stderr);
}

*value = (int) strtol ((char *) strtok (char_buf, delimiters), NULL, 10);/*takes what is before #*/
info = (char *) strtok (NULL, delimiters);/*takes what is after #*/

//fprintf(stdout, "\n%s: %d", info, *value);/*prints for check*/
//fflush (stdout);

delete [] char_buf;
}

/*-----------------------------------------------------------------*/

void read_line_2get_string (FILE *in, char **value)
{
/*reads current line and stores value before # () as a string value*/
size_t nchar = 1024;
char *char_buf = new char[nchar];
const char delimiters[4] = {'#','\n','\0'};
char *info;
char *tmp;

if (getline (&char_buf, &nchar, in) < 2) /*reads line*/
{
    fprintf (stderr, "\nread_line_2get_string: file reading error!\a\n");
    fflush (stderr);
}

tmp = (char *) strtok (char_buf, delimiters);/*takes what is before #*/
info = (char *) strtok (NULL, delimiters);/*takes what is after #*/
strcpy (*value, trim(tmp));

//fprintf(stdout, "\n%s: %s", info, *value);/*prints for check*/
//fflush (stdout);

delete [] char_buf;
}

/*-----------------------------------------------------------------*/

void read_line_2skip_it (FILE *in, char **value)
{
/*reads current line*/
size_t nchar = 1024;
char *char_buf = new char[nchar];

getline (&char_buf, &nchar, in); /*reads line*/

//fprintf(stdout, "\n%s", char_buf);/*prints for check*/
//fflush (stdout);

delete [] char_buf;
}

/*-----------------------------------------------------------------*/

int load_and_alloc_profile (char *name, int *dim, double **rgrid, double **qgrid)
{
FILE *in;

size_t nchar = 1024;

char *file_name = new char[nchar];

strcpy (file_name, name);
strcat (file_name, ".dim");

if ((in=fopen(file_name, "r")) == NULL)
{
    fprintf (stderr, "\nfailed to open file %s\a\a\a\n", file_name);
    exit(0);
}

char *char_buf = new char[nchar];

getline (&char_buf, &nchar, in); /*reads line*/

fclose(in);

*dim = (int) strtol (char_buf, NULL, 0); /*converts to an integer*/

*rgrid = new double[*dim];
*qgrid = new double[*dim];

strcpy (file_name, name);
strcat (file_name, ".dat");

/*fprintf(stdout, "\nProfile %s has %d points.", file_name, *dim);*/

if ((in=fopen(file_name, "r")) == NULL)
{
    fprintf (stderr, "\nfailed to open file %s\a\a\a\n", file_name);
    exit(0);
}

char *tail_ptr1, *tail_ptr2;

int i;
for (i=0; i<*dim; i++)
{
    if (feof(in) || ferror(in)) break;
    if (getline (&char_buf, &nchar, in) < 1) /*reads line*/
    {
        fprintf (stderr, "\nload_and_alloc_profile: file reading error!\a\n");
    }
    (*rgrid)[i] = strtod (char_buf, &tail_ptr1);
    (*qgrid)[i] = strtod (tail_ptr1, &tail_ptr2);
}

if(i != *dim)
{
    fprintf(stderr, "\nerror: load_and_alloc_profile: read error or false array dimension!\a\n");
    exit(0);
}

fclose (in);

delete [] char_buf;
delete [] file_name;

return 0;
}

/*-----------------------------------------------------------------*/

int load_profile (char *name, int dim, double *rgrid, double *qgrid)
{
FILE *in;

size_t nchar = 1024;

char *file_name = new char[nchar];

strcpy (file_name, name);
strcat (file_name, ".dim");

if ((in=fopen(file_name, "r")) == NULL)
{
    fprintf (stderr, "\nfailed to open file %s\a\a\a\n", file_name);
    exit (0);
}

char *char_buf = new char[nchar];

getline (&char_buf, &nchar, in); /*reads line*/

fclose(in);

int dimt = (int) strtol (char_buf, NULL, 0); /*converts to an integer*/

if (dim != dimt)
{
    fprintf (stderr, "\ndimension=%d of the profile %s is different from the value=%d in set.in!", dimt,  file_name, dim);
    exit (0);
}

strcpy (file_name, name);
strcat (file_name, ".dat");

/*fprintf(stdout, "\nProfile %s has %d points.", file_name, *dim);*/

if ((in=fopen(file_name, "r")) == NULL)
{
    fprintf (stderr, "\nfailed to open file %s\a\a\a\n", file_name);
    exit(0);
}

char *tail_ptr1, *tail_ptr2;

int err = 0;

int i;
for (i=0; i<dim; i++)
{
    if (feof(in) || ferror(in)) break;
    if (getline (&char_buf, &nchar, in) < 1) /*reads line*/
    {
        fprintf (stderr, "\nload_profile: file reading error!\a\n");
    }
    rgrid[i] = strtod (char_buf, &tail_ptr1);
    if (tail_ptr1 == char_buf) err++;

    qgrid[i] = strtod (tail_ptr1, &tail_ptr2);
    if (tail_ptr2 == tail_ptr1) err++;
}

if (i != dim)
{
    fprintf(stderr, "\nerror: load_profile: read error or false array dimension!\a\n");
    exit(0);
}

if (err)
{
    fprintf(stderr, "\nerror: load_profile: possible read error: count = %d!\n", err);
    exit(0);
}

fclose (in);

delete [] char_buf;
delete [] file_name;

return err;
}

/*-----------------------------------------------------------------*/

int count_lines_in_file (char *filename, int flag_print)
{
FILE *in;

if ((in = fopen (filename, "r")) == NULL)
{
    fprintf (stderr, "\ncount_lines_in_file: failed to open file %s\n", filename);
}

size_t nchar = 65536;
char *char_buf = new char[nchar];

int num_lines = 0;

while (getline (&char_buf, &nchar, in) > 0) num_lines++;

fclose (in);

delete [] char_buf;

if (flag_print)
{
    fprintf (stdout, "\nfile %s contains %d non-empty lines\n", filename, num_lines);
    fflush (stdout);
}

return num_lines;
}

/*-----------------------------------------------------------------*/

int load_complex_profile (char *name, int dim, double *rgrid, double *qgrid)
{
FILE *in;

size_t nchar = 1024;

char *file_name = new char[nchar];
char *char_buf = new char[nchar];

strcpy (file_name, name);
//strcat (file_name, ".dat");

/*fprintf(stdout, "\nProfile %s has %d points.", file_name, *dim);*/

if ((in=fopen(file_name, "r")) == NULL)
{
    fprintf (stderr, "\nfailed to open file %s\a\a\a\n", file_name);
    exit(0);
}

char *tail_ptr1, *tail_ptr2, *tail_ptr3;

int i;
int err = 0;

for (i=0; i<dim; i++)
{
    if (feof(in) || ferror(in)) break;
    if (getline (&char_buf, &nchar, in) < 1) /*reads line*/
    {
        fprintf (stderr, "\nload_profile: file reading error!\a\n");
    }

    rgrid[i] = strtod (char_buf, &tail_ptr1);
    if (tail_ptr1 == char_buf) err++;

    qgrid[i] = strtod (tail_ptr1, &tail_ptr2);
    if (tail_ptr2 == tail_ptr1) err++;

    qgrid[i+dim] = strtod (tail_ptr2, &tail_ptr3);
    if (tail_ptr3 == tail_ptr2) err++;
}

if (i != dim)
{
    fprintf(stderr, "\nerror: load_complex_profile: read error or false array dimension!\a\n");
    exit(0);
}

if (err)
{
    fprintf(stderr, "\nerror: load_complex_profile: possible read error: count = %d!\n", err);
    exit(0);
}

fclose (in);

delete [] char_buf;
delete [] file_name;

return err;
}

/*-----------------------------------------------------------------*/

int load_data_file (char *file_name, int dim, int ncols, double *rgrid, double *qgrid)
{
//dim and ncols to be checked inside; ncols means number of y data columns!
//(x, y1, y2,...,yncols)

if (dim != count_lines_in_file (file_name, 0))
{
    fprintf(stderr, "\nerror: load_data_file: read error or false input data dimension: %s\n", file_name);
    exit(1);
}

FILE *in;

size_t nchar = ncols*1024;

char *char_buf = new char[nchar];

if ((in=fopen(file_name, "r")) == NULL)
{
    fprintf (stderr, "\nload_data_file: failed to open file %s\n", file_name);
    exit(1);
}

char *str1, *str2;

int i, k;
int err = 0;

for (i=0; i<dim; i++)
{
    if (feof(in) || ferror(in)) break;

    if (getline (&char_buf, &nchar, in) < 1) /*reads line*/
    {
        fprintf (stderr, "\nload_data_file: file %s reading error!\n", file_name);
        exit(1);
    }
    rgrid[i] = strtod (char_buf, &str2);
    if (str2 == char_buf) err++;
    str1 = str2;

    for (k=0; k<ncols; k++)
    {
        qgrid[i+k*dim] = strtod (str1, &str2);
        if (str2 == str1) err++;
        str1 = str2;
    }
}

if (i != dim)
{
    fprintf(stderr, "\nerror: load_data_file: %s: read error or false data dimension!\n", file_name);
    exit(1);
}

if (err)
{
    fprintf(stderr, "\nerror: load_data_file: %s: %d read errors detected!\n", file_name, err);
    exit(1);
}

fclose (in);

delete [] char_buf;

return err;
}

/*-----------------------------------------------------------------*/

int count_lines_in_file_with_comments (char *filename, int flag_print)
{
FILE *in;

if ((in = fopen (filename, "r")) == NULL)
{
    fprintf (stderr, "\ncount_lines_in_file: failed to open file %s\n", filename);
}

size_t nchar = 65536;
char *char_buf = new char[nchar];

int num_lines = 0;

while (getline (&char_buf, &nchar, in) > 0)
{
    if (!(strchr (char_buf, '%') || strchr (char_buf, '#') || strchr (char_buf, '!')))
    {
        num_lines++;
    }
}

fclose (in);

delete [] char_buf;

if (flag_print)
{
    fprintf (stdout, "\nfile %s contains %d non-empty lines\n", filename, num_lines);
    fflush (stdout);
}

return num_lines;
}

/*-----------------------------------------------------------------*/

int load_data_file_with_comments (char *file_name, int dim, int ncols, double *rgrid, double *qgrid)
{
//dim and ncols to be checked inside; ncols means number of y data columns!
//(x, y1, y2,...,yncols)

if (dim != count_lines_in_file_with_comments (file_name, 0))
{
    fprintf(stderr, "\nerror: load_data_file: read error or false input data dimension: %s\n", file_name);
    exit(1);
}

FILE *in;

size_t nchar = ncols*1024;

char *char_buf = new char[nchar];

if ((in=fopen(file_name, "r")) == NULL)
{
    fprintf (stderr, "\nload_data_file: failed to open file %s\n", file_name);
    exit(1);
}

char *str1, *str2;

int i, k;
int err = 0;

int ind = 0;

int dimf = count_lines_in_file (file_name, 0);

for (i=0; i<dimf; i++)
{
    if (feof(in) || ferror(in)) break;

    if (getline (&char_buf, &nchar, in) < 1) /*reads line*/
    {
        fprintf (stderr, "\nload_data_file: file %s reading error!\n", file_name);
        exit(1);
    }

    if (strchr (char_buf, '%') || strchr (char_buf, '#') || strchr (char_buf, '!')) continue;

    rgrid[ind] = strtod (char_buf, &str2);
    if (str2 == char_buf) err++;
    str1 = str2;

    for (k=0; k<ncols; k++)
    {
        qgrid[ind+k*dim] = strtod (str1, &str2);
        if (str2 == str1) err++;
        str1 = str2;
    }

    ind++;
}

if (ind != dim)
{
    fprintf(stderr, "\nerror: load_data_file: %s: read error or false data dimension!\n", file_name);
    exit(1);
}

if (err)
{
    fprintf(stderr, "\nerror: load_data_file: %s: %d read errors detected!\n", file_name, err);
    exit(1);
}

fclose (in);

delete [] char_buf;

return err;
}

/*-----------------------------------------------------------------*/

void count_lines_in_file_ (char *filename, int *flag_print, int *dim)
{
*dim = count_lines_in_file (filename, *flag_print);
}

/*-----------------------------------------------------------------*/
