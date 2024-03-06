/*! \file shared.cpp
    \brief The implementation of some common and frequently used functions.
*/

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <string.h>

#include "constants.h"
#include "shared.h"

/*******************************************************************/

extern "C"
{
void dgesv_ (int *, int *, double *, int *, int *, double *, int *, int *);
}

/*******************************************************************/

void * xmalloc (size_t size)
{
register void *value = (void *) malloc (size);
if (!value)
{
	fprintf(stderr,"\nxmalloc: memory allocation failed: exit!");
	exit(1);
}
return value;
}

/*******************************************************************/

void * xrealloc (void *ptr, size_t size)
{
register void *value = (void *) realloc (ptr, size);
if (!value)
{
	fprintf(stderr,"\nxrealloc: memory reallocation failed: exit!");
	exit(1);
}
return value;
}

/*******************************************************************/

int signum(double x)
{
return (x<0.0) ? (-1):((x==0) ? 0:1);
}

/*******************************************************************/

int compare_doubles (const void *a, const void *b)
{
const double *da = (const double *) a;
const double *db = (const double *) b;

return (*da > *db) - (*da < *db);
}

/*******************************************************************/

void binomial_coefficients (int N, double *BC)
{
//computes C^k_n = n!/k!/(n-k)! coefficients for n=0..N, k=0..n
//fortran: C(n,k) = C^k_n; C: C[k][n] = C^k_n;
int k, n;
double tmp;

for (n=0; n<=N; n++)
{
    tmp = 1.0;
    BC[n] = tmp; //C^0_n
    for (k=1; k<=n; k++)
    {
        tmp *= double(n-k+1)/double(k);
        BC[n+k*(N+1)] = tmp;
    }
}
}

/*******************************************************************/

void binomial_coefficients_ (int N, double *BC)
{
binomial_coefficients (N, BC);
}

/*******************************************************************/

void Ckn (int N, double *BC)
{
//computes C[k][n] = C^k_n = n!/k!/(n-k)! coefficients for n=0..N, k=0..n
int k, n;
double tmp;

typedef double bincoeffs[N+1][N+1];
bincoeffs &C = *((bincoeffs *)(BC));

for (n=0; n<=N; n++)
{
    tmp = 1.0;
    C[0][n] = tmp; //C^0_n
    for (k=1; k<=n; k++)
    {
        tmp *= double(n-k+1)/double(k);
        C[k][n] = tmp;
    }
}
}

/*******************************************************************/

int calc_interp_polynom (int N, const double *x, const double *y, double *C)
{
//calculates an interpolating polinom of degree N through the points:
//(x[0],y[0]),...,(x[N], y[N])
//strores coefficients in array C for the later evaluation

//the code is not checked!!!

typedef double matrix[N+1][N+1];

matrix &A = *((matrix *)(new double[(N+1)*(N+1)]));

int D = N+1, NRHS = 1, LDA = D, LDB = D, INFO;

int *IPIV = new int[D];

int p, i;

for (i=0; i<=N; i++)
{
    A[0][i] = 1.0;
    A[1][i] = x[i]-x[0];
    for (p=2; p<=N; p++) //fortran matrix ordering
    {
        A[p][i] = A[p-1][i]*A[1][i];
    }
    C[i] = y[i];
}

//linear system:
dgesv_ (&D, &NRHS, (double *)A, &LDA, IPIV, C, &LDB, &INFO);
if (INFO)
{
    fprintf (stderr, "\nerror: calc_interp_polynom: system failed: INFO=%d.", INFO);
    return INFO;
}

delete [] IPIV;
delete [] A;
return 0;
}

/*******************************************************************/

int eval_interp_polynom (int N, const double *x, const double *C, double r, int Dmin, int Dmax, double *R)
{
if (r<x[0] || r>x[N])
{
    fprintf (stderr, "\nerror: eval_interp_polynom: argument is outside the range: r=%le", r);
    return 1;
}

//the code is not finished!!!
fprintf (stderr, "\nerror: eval_interp_polynom: not implemented!");

return 0;
}

/*******************************************************************/

void localizator (double x, double x0, double L, double *W)
{
double dir, x1, x2, t, fac;

dir = signum (x-x0);

if (dir <= 0.0)
{
    x1 = x0 - L;
    x2 = x1 + L/2.0;
}
else
{
    x2 = x0 + L;
    x1 = x2 - L/2.0;
}

if (dir > 0.0)
{
    t = (x-x1)/(x2-x1);
    fac = 1.0;
}
else
{
    t = (x2-x)/(x2-x1);
    fac = -1.0;
}

double c1 = 2.0*pi, c2 = 1.414213562373095;

if (t <= 0)
{
    W[0] = 1.0;
    W[1] = 0.0;
}
else if (t >= 1.0)
{
    W[0] = 0.0;
    W[1] = 0.0;
}
else
{
    W[0] = exp(-c1/(1.0-t)*exp(-c2/t));
    W[1] = - W[0]*c1/(1.0-t)*exp(-c2/t)*(1.0/(1.0-t)+c2/t/t)*fac/(x2-x1);
}
}

/*******************************************************************/

void localizator_4_derivs (double x, double x0, double L, double *W)
{
double dir, x1, x2, t, fac;

dir = signum (x-x0);

if (dir <= 0.0)
{
    x1 = x0 - L;
    x2 = x1 + L/2.0;
}
else
{
    x2 = x0 + L;
    x1 = x2 - L/2.0;
}

if (dir > 0.0)
{
    t = (x-x1)/(x2-x1);
    fac = 1.0;
}
else
{
    t = (x2-x)/(x2-x1);
    fac = -1.0;
}

double c1 = 2.0*pi, c2 = 1.414213562373095;
double E = 2.718281828459045235360287;

if (t <= 1.0e-2) //0.0 - needed to avoid nans
{
    W[0] = 1.0;
    W[1] = 0.0;
    W[2] = 0.0;
    W[3] = 0.0;
    W[4] = 0.0;
}
else if (t >= 1.0 - 1.0e-2) //1.0 - needed to avoid nans
{
    W[0] = 0.0;
    W[1] = 0.0;
    W[2] = 0.0;
    W[3] = 0.0;
    W[4] = 0.0;
}
else
{
    W[0] = exp(-c1/(1.0-t)*exp(-c2/t));

    W[1] = - c1/(1.0-t)*exp(-c2/t)*(1.0/(1.0-t)+c2/t/t);
    W[1] *= W[0]*pow(fac/(x2-x1),1);

    W[2] = (c1*(c1*pow(c2 - c2*t + pow(t,2),2) + pow(E,c2/t)*(-1 + t)*
           (pow(c2,2)*pow(-1 + t,2) + 2*pow(t,4) - 2*c2*t*(1 - 3*t + 2*pow(t,2)))))/(pow(E,(2*c2)/t)*pow(-1 + t,4)*pow(t,4));
    W[2] *= W[0]*pow(fac/(x2-x1),2);

    W[3] = (c1*(pow(c1,2)*pow(c2*(-1 + t) - pow(t,2),3) - 3*c1*pow(E,c2/t)*
           (-1 + t)*(-(pow(c2,3)*pow(-1 + t,3)) + 2*pow(t,6) +
           pow(c2,2)*pow(-1 + t,2)*t*(-2 + 5*t) - 2*c2*pow(t,3)*(1 - 4*t + 3*pow(t,2)))
           - pow(E,(2*c2)/t)*pow(-1 + t,2)*(-(pow(c2,3)*pow(-1 + t,3)) + 6*pow(t,6) + 
           3*pow(c2,2)*pow(-1 + t,2)*t*(-2 + 3*t) - 6*c2*pow(t,2)*(-1 + 4*t - 6*pow(t,2)
           +3*pow(t,3)))))/(pow(E,(3*c2)/t)*pow(-1 + t,6)*pow(t,6));
    W[3] *= W[0]*pow(fac/(x2-x1),3);

    W[4] = (c1*(pow(c1,3)*pow(c2 - c2*t + pow(t,2), 4) + 6*pow(c1,2)*pow(E,c2/t)*
           (-1 + t)*pow(c2 - c2*t + pow(t,2), 2)*(pow(c2,2)*pow(-1 + t,2) + 2*pow(t,4) -
           2*c2*t*(1 - 3*t + 2*pow(t,2))) + c1*pow(E,(2*c2)/t)*pow(-1 + t,2)*
           (7*pow(c2,4)*pow(-1 + t,4) + 36*pow(t,8) - 4*pow(c2,3)*pow(-1 + t,3)*t*
           (-9 + 16*t) + 12*pow(c2,2)*pow(-1 + t,2)*pow(t,2)*(3 - 12*t + 14*pow(t,2))
           - 24*c2*pow(t,4)*(-1 + 5*t - 10*pow(t,2) + 6*pow(t,3))) + pow(E,(3*c2)/t)*
           pow(-1 + t,3)*(pow(c2,4)*pow(-1 + t,4) + 24*pow(t,8) - 4*pow(c2,3)*
           pow(-1 + t,3)*t*(-3 + 4*t) + 12*pow(c2,2)*pow(-1 + t,2)*
           pow(t,2)*(3 - 8*t + 6*pow(t,2)) - 24*c2*pow(t,3)*(1 - 5*t + 10*pow(t,2) - 
           10*pow(t,3) + 4*pow(t,4)))))/(pow(E,(4*c2)/t)*pow(-1 + t,8)*pow(t,8));
    W[4] *= W[0]*pow(fac/(x2-x1),4);
}
}

/*******************************************************************/
