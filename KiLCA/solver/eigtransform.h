#include "lapack.h"
#include "shared.h"
#include "mat_interp.h"

int coeff_start_vals (int Nw, int Nfs, double r, double *zstart);
int coeffs2solution (int Nw, int Nfs, int dim, double *rgrid, double *sol);
