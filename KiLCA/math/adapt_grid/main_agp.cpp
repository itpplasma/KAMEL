#include <stdlib.h>
#include <stdio.h>
#include <cmath>
#include <complex>
#include <cstring>
#include <climits>

#include "adaptive_grid_pol.h"

using namespace std;

/*****************************************************************************/

void test1(double *x, double *y, void *p);

/*****************************************************************************/

int main(void)
{
double a = 0.0;
double b = 60.0;

int dimy = 3;
int xdim = 1000;

int deg = 9;
double eps = 1.0e-12;

double *x = new double[xdim];
double *y = new double[dimy*xdim];

adaptive_grid_polynom (&test1, NULL, a, b, dimy, deg, &xdim, &eps, x, y);

FILE *out;
if (!(out = fopen ("resgrid.dat", "w")))
{
    fprintf (stderr, "\nFailed to open file\a\n");
}

for (int k=0; k<xdim; k++)
{
    fprintf (out, "\n%24.16le\t%24.16le\t%24.16le\t%24.16le", x[k], y[dimy*k+0], y[dimy*k+1], y[dimy*k+2]);
    //fprintf (out, "\n%24.16le\t%24.16le", x[k], y[dimy*k+0]);

}

fclose(out);

if (!(out = fopen ("checkgrid.dat", "w")))
{
    fprintf (stderr, "\nFailed to open file\a\n");
}

//check interpolation:
int dimt = xdim*10, ind = 0;
double xt[dimt], ye[dimy*dimt], yi[dimy*dimt];

for (int k=0; k<dimt; k++)
{
    xt[k] = a + k*(b-a)/(dimt-1);

    test1 (xt+k, ye+k*dimy, NULL);

    //interpolation:
    find_index_for_interp (deg, xt[k], xdim, x, &ind);

    eval_interp_polynom (deg, x+ind, y+ind*dimy, dimy, xt[k], yi+k*dimy);

    //output:
  fprintf (out, "\n%24.16le\t%24.16le\t%24.16le\t%24.16le\t%24.16le\t%24.16le\t%24.16le\t%24.16le\t%24.16le\t%24.16le", xt[k], ye[dimy*k+0], yi[dimy*k+0], abs(ye[dimy*k+0]-yi[dimy*k+0])/abs(ye[dimy*k+0]), ye[dimy*k+1], yi[dimy*k+1], abs(ye[dimy*k+1]-yi[dimy*k+1])/abs(ye[dimy*k+1]), ye[dimy*k+2], yi[dimy*k+2], abs(ye[dimy*k+2]-yi[dimy*k+2])/abs(ye[dimy*k+2]));

  //  fprintf (out, "\n%24.16le\t%24.16le\t%24.16le\t%24.16le", xt[k], ye[dimy*k+0], yi[dimy*k+0], abs(ye[dimy*k+0]-yi[dimy*k+0])/abs(ye[dimy*k+0]));
}
fclose (out);

delete [] x;
delete [] y;

return 0;
}

/*****************************************************************************/

void test1(double *x, double *y, void *p)
{
double x1 = 0.1*(*x);
double x2 = 0.1*(*x);

//y[0] = sin(*x)/(*x+0.001);
//y[0] = (*x)*sin(1/(*x+0.1));

y[0] = (x1)*sin(x1);
y[1] = sin(x1)/((x1)+0.001);
y[2] = 1+(*x)-(x1)*(x1)+(x2)*(x2)*(x2)-pow(x2,7);
}

/*****************************************************************************/
