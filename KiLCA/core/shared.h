/*! \file shared.h
    \brief Some common and frequently used functions declarations.
*/

#ifndef SHARED_INCLUDE

#define SHARED_INCLUDE

#include <stdio.h>
#include <limits.h>
#include <math.h>
#include <stdlib.h>

//#if defined(__cplusplus)
//extern "C" {
//#endif

void * xmalloc (size_t size);

void * xrealloc (void *ptr, size_t size);

int signum(double x);

int compare_doubles (const void *a, const void *b);

void binomial_coefficients (int N, double *BC);

struct cmplx_number
{
  double re;
  double im;
};

typedef struct cmplx_number cmplx;

// #if defined(__cplusplus)
// }
// #endif

extern "C"
{
void binomial_coefficients_ (int N, double *BC);
void Ckn (int N, double *C);
void localizator (double x, double x0, double L, double *W);
void localizator_4_derivs (double x, double x0, double L, double *W);
}

#endif
