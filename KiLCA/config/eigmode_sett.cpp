#include "eigmode_sett.h"

#include <stdlib.h>
#include <stdio.h>
#include <complex>
#include <cstring>

void eigmode_sett::print_settings()
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
