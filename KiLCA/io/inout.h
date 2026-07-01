/*! \file
    \brief The declarations of functions for input/output routines.
*/

#ifndef INOUT_INCLUDE

#define INOUT_INCLUDE

#include <complex>
#include <stdio.h>

using namespace std;

// All routines are implemented in Fortran (kilca_inout_m, io/inout_m.f90) and
// exported with C linkage via bind(C). The FILE*-based readers keep the handle
// opaque and read it through libc, so C++ callers pass their fopen'd FILE*
// through unchanged.
extern "C"
{
int save_cmplx_matrix (int Nrows, int Ncols, int Npoints, const double *xgrid, const double *arr, const char *path_name);

int save_cmplx_matrix_to_one_file (int Nrows, int Ncols, int Npoints, const double *xgrid, const double *arr, const char *full_name);

int save_real_array (int dim, const double *xgrid, const double *arr, const char *full_name);

int load_data_file (char *file_name, int dim, int ncols, double *rgrid, double *qgrid);

int count_lines_in_file (char *filename, int flag_print);

void read_line_2get_double (FILE *in, double *value);

void read_line_2get_complex (FILE *in, complex<double> *value);

void read_line_2get_int (FILE *in, int *value);

void read_line_2get_string (FILE *in, char **value);

void read_line_2skip_it (FILE *in, char **value);

int count_lines_in_file_with_comments (char *filename, int flag_print);

int load_data_file_with_comments (char *file_name, int dim, int ncols, double *rgrid, double *qgrid);
}

#endif
