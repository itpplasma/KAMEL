/*! \file eigmode_sett.cpp
    \brief The implementation of eigmode_sett class.
*/

#include <stdlib.h>
#include <stdio.h>
#include <complex>
#include <cstring>

#include "eigmode_sett.h"
#include "inout.h"

/*****************************************************************************/

void eigmode_sett::read_settings (char *path)
{
char *file_set = new char[1024];

sprintf (file_set, "%s/eigmode.in", path);

FILE *in;

if ((in=fopen (file_set, "r"))==NULL)
{
    fprintf(stderr, "\nerror: read_settings: failed to open file %s\a\n", file_set);
    exit(0);
}

char *str_buf = new char[1024]; //str buffer

//Output:
read_line_2skip_it (in, &str_buf);
fname = new char[1024];
read_line_2get_string (in, &(fname));
read_line_2skip_it (in, &str_buf);

//frequency scan or root search:
read_line_2skip_it (in, &str_buf);
read_line_2get_int (in, &(search_flag));
read_line_2skip_it (in, &str_buf);

//frequency grid:
read_line_2skip_it (in, &str_buf);
read_line_2get_int (in, &(rdim));
read_line_2get_double (in, &(rfmin));
read_line_2get_double (in, &(rfmax));
read_line_2get_int (in, &(idim));
read_line_2get_double (in, &(ifmin));
read_line_2get_double (in, &(ifmax));
read_line_2skip_it (in, &str_buf);

//Stopping criteria:
read_line_2skip_it (in, &str_buf);
read_line_2get_int (in, &(stop_flag));
read_line_2get_double (in, &(eps_res));
read_line_2get_double (in, &(eps_abs));
read_line_2get_double (in, &(eps_rel));
read_line_2skip_it (in, &str_buf);

//For derivative:
read_line_2skip_it (in, &str_buf);
read_line_2get_double (in, &(delta));
read_line_2skip_it (in, &str_buf);

//For testing:
read_line_2skip_it (in, &str_buf);
read_line_2get_int (in, &(test_roots));
read_line_2get_int (in, &(flag_debug));
read_line_2skip_it (in, &str_buf);

// ZerSol parameters:
read_line_2skip_it (in, &str_buf);
read_line_2get_int (in, &(n_zeros));
read_line_2get_int (in, &(use_winding));
read_line_2skip_it (in, &str_buf);

//Starting points:
read_line_2skip_it (in, &str_buf);
read_line_2get_int (in, &(Nguess));
read_line_2get_int (in, &(kmin));
read_line_2get_int (in, &(kmax));
read_line_2skip_it (in, &str_buf);

fstart = new complex<double>[Nguess];

int k;

read_line_2skip_it (in, &str_buf);

for (k=0; k<Nguess; k++)
{
    read_line_2get_complex (in, &(fstart[k]));
}

fclose (in);

delete [] file_set;
delete [] str_buf;

if (flag_debug) print_settings ();
}

/*****************************************************************************/

void eigmode_sett::print_settings (void)
{
fprintf(stdout, "\nCheck for eigmode settings below:\n");

fprintf (stdout, "\nfile name: %s", fname);
fprintf (stdout, "\nsearch flag: %d", search_flag);
fprintf (stdout, "\nreal freq mesh dim: %d", rdim);
fprintf (stdout, "\nreal freq mesh minimum: %le", rfmin);
fprintf (stdout, "\nreal freq mesh maximum: %le", rfmax);
fprintf (stdout, "\nimag freq mesh dim: %d", idim);
fprintf (stdout, "\nimag freq mesh minimum: %le", ifmin);
fprintf (stdout, "\nimag freq mesh maximum: %le", ifmax);
fprintf (stdout, "\nstopping criteria: %d", stop_flag);
fprintf (stdout, "\nresidual error parameter: %le", eps_res);
fprintf (stdout, "\nabs error parameter: %le", eps_abs);
fprintf (stdout, "\nrel error parameter: %le", eps_rel);
fprintf (stdout, "\ndelta for derivative: %le", delta);
fprintf (stdout, "\ntest roots flag: %d", test_roots);
fprintf (stdout, "\nflag_debug: %d", flag_debug);
fprintf (stdout, "\nn_zeros: %d", n_zeros);
fprintf (stdout, "\nuse_winding: %d", use_winding);
fprintf (stdout, "\nguess array dimension: %d", Nguess);
fprintf (stdout, "\nkmin: %d", kmin);
fprintf (stdout, "\nkmax: %d", kmax);

int k;
fprintf (stdout, "\nguess array:");
for (k=0; k<Nguess; k++)
{
    printf ("\nk=%d\tf=(%le, %le)", k, real(fstart[k]), imag(fstart[k]));
}

fprintf(stdout, "\n");
}

/*****************************************************************************/
