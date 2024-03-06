/*! \file
    \brief The declarations of functions for input/output routines.
*/

#ifndef INOUT_INCLUDE

#define INOUT_INCLUDE

#include <complex>

using namespace std;

int save_cmplx_matrix (int Nrows, int Ncols, int Npoints, const double *xgrid, const double *arr, const char *path_name);

int save_cmplx_matrix_to_one_file (int Nrows, int Ncols, int Npoints, const double *xgrid, const double *arr, const char *full_name);

int save_real_matrix_to_one_file (int order, int Nrows, int Ncols, int Npoints, const double *xgrid, const double *arr, const char *full_name);

int save_real_array (int dim, const double *xgrid, const double *arr, const char *full_name);

int save_complex_array (int dim, const double *xgrid, const double *arr, const char *full_name);

char * trim (char *str);

void read_line_2get_double (FILE *in, double *value);

void read_line_2get_complex (FILE *in, complex<double> *value);

void read_line_2get_int (FILE *in, int *value);

void read_line_2get_string (FILE *in, char **value);

void read_line_2skip_it (FILE *in, char **value);

int load_profile (char *name, int dim, double *rgrid, double *qgrid);

int load_and_alloc_profile (char *name, int *dim, double **rgrid, double **qgrid);

int count_lines_in_file (char *filename, int flag_print);

int load_complex_profile (char *name, int dim, double *rgrid, double *qgrid);

int load_data_file (char *file_name, int dim, int ncols, double *rgrid, double *qgrid);

int count_lines_in_file_with_comments (char *filename, int flag_print);

int load_data_file_with_comments (char *file_name, int dim, int ncols, double *rgrid, double *qgrid);

#if defined(__cplusplus)
extern "C"
{
#endif

#if defined(__cplusplus)
}
#endif

extern "C"
{
void count_lines_in_file_ (char *filename, int *flag_print, int *dim);
}

#endif
