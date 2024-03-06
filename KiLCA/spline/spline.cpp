/*! \file
    \brief The implementation of spline_data class.
*/

#include <cmath>
#include <cstring>
#include <stdio.h>
#include <stdlib.h>

#include "spline.h"

/*
To do:
2. basic checks in init subr
3. optimize further.
4. add more comments.
5. implement other boundary conditions.
*/

/*******************************************************************/

void spline_alloc_(int N, int type, int dimx, const double *x, double *C,
                   uintptr_t *sid) {
  /*
  N (input)- degree of the spline polinoms must be odd!
  type (input) - type of the boundary: natural, zero or periodic
  dimx (input) - abscissa grid size
  x (input) - abscissa aray

  C - (output) array of spline coefficients stores as: C[0][0][0...N],
  C[0][1][0...N], ..., C[dimy-1][...]

  sid (output)- splines id
  */

  // some checks:
  if (!(N > 0 && (N % 2) == 1 && dimx >= N && type == 1)) {
    fprintf(stderr,
            "\nerror: check of input parameters failed in the function '%s' at "
            "line %d of the file '%s'.\n",
            __FUNCTION__, __LINE__, __FILE__);
    return;
  }

  // allocates spline structure:
  spline_data *sd = new spline_data;
  *sid = (uintptr_t)sd;

  // set parameters:
  sd->N = N;

  sd->type = type;

  sd->dimx = dimx;
  sd->x = x;

  sd->ind = 0.5 * dimx;

  sd->C = C;

  // calcs binomial coefficients:
  sd->BC = new double[(N + 1) * (N + 1)];
  set_bc_array(N, sd->BC);

  // calcs fac coeffs used in spline evaluation:
  sd->fac = new double[(N + 1) * (N + 1)];
  set_fac_array(N, sd->fac);
}

/*******************************************************************/

void spline_calc_(uintptr_t sid, const double *y, int Imin, int Imax, double *W,
                  int *ierr) {
  /*
  sid - spline id

  y (input) - array of interpolated functions stored as: y[0][0], y[0][1], ...,
  y[0][dimx-1], y[1][0], y[1][1], ..., y[1][dimx-1], y[dimy-1][0], y[dimy-1][1],
  ..., y[dimy-1][dimx-1].

  W (input) - working array of the minimum dimension
  (N+1)(2dimy+dimx(4+3((N-1)/2))) - not needed after return from this function.
  If W=NULL W is allocated inside.

  Imin, ..., Imax - indices of y data in C array

  ierr (output) - error code: zero if successful.
  */

  // checks: jmin, jmax
  //  if (Dmax > N || Dmin<0 || Dmax<Dmin || dimz<1)
  //  {
  //      fprintf (stderr, "\nerror: spline_eval: check input parameters
  //      failed!");
  //  }

  spline_data *sd = (spline_data *)(sid);

  int dimy = Imax - Imin + 1;

  int flag_alloc;
  if (W == NULL) {
    W = new double[(sd->N + 1) *
                   (2 * dimy + (sd->dimx) * (4 + 3 * ((sd->N - 1) / 2)))];
    flag_alloc = 1;
  }

  *ierr = 0;

  // for boundary conditions:
  *ierr = sd->calc_spline_boundaries(dimy, y, W);
  if (*ierr) {
    fprintf(stderr, "\nerror: spline_calc: calc_spline_boundaries failed!");
  }

  // coefficients:
  *ierr = sd->calc_spline_coefficients(Imin, dimy, y, W);
  if (*ierr) {
    fprintf(stderr, "\nerror: spline_calc: calc_spline_coefficients failed!");
  }

  if (flag_alloc)
    delete[] W;
}

/*******************************************************************/

int spline_data::calc_spline_boundaries(int dimy, const double *y, double *W) {
  double *Cbnd = W; // for boundaries: 2(N+1)dimy

  double *A = W + 2 * (N + 1) * dimy; //[(N+1)*(N+1)];

  double *b = NULL;

  int D = N + 1, NRHS = dimy, LDA = D, LDB = D, INFO;

  int *IPIV = new int[D];

  int ind, k, p, i, j;

  // type = 0 - zero boundary conditions: not implemented

  // type = 1 - natural boundary conditions:

  for (k = 0; k < dimx; k += dimx - 1) // two boundaries
  {
    b = Cbnd + (k / (dimx - 1)) * (N + 1) * dimy;
    for (i = 0; i <= N; i++) {
      ind = k + i * sign(1 - k);
      A[i] = 1.0;
      A[i + N + 1] = x[ind] - x[k];
      for (p = 2; p <= N; p++) // fortran matrix ordering
      {
        A[i + p * (N + 1)] = A[i + (p - 1) * (N + 1)] * A[i + N + 1];
      }
      for (j = 0; j < dimy; j++)
        b[i + j * (N + 1)] = y[ind + j * dimx];
    }

    // linear system:
    dgesv_(&D, &NRHS, A, &LDA, IPIV, b, &LDB, &INFO);
    if (INFO) {
      fprintf(stderr,
              "\nerror: calc_spline_boundaries: failed to interpolate near "
              "boundaries: INFO=%d dimx=%d dimy=%d k=%d.",
              INFO, dimx, dimy, k);
    }
  }

  // type = 1 - periodic boundary conditions: not implemented

  delete[] IPIV;
  return INFO;
}

/*******************************************************************/

int spline_data::calc_spline_coefficients(int Imin, int dimy, const double *y,
                                          double *W) {
  int s = (N - 1) / 2;

  int KL = s + 1, KU = s + 1, LDAB = 2 * KL + KU + 1;

  int len = (N + 1) * dimx, LEN = LDAB * len;

  double *Cbnd = W;

  double *M =
      W + 2 * (N + 1) * dimy; // new double[LEN]; //main band system matrix

  int j;
  for (j = 0; j < LEN; j++)
    M[j] = 0;

  int ieqn = 0; // an equation index
  int iunk;     // an unknown index

  double *SC =
      C + (N + 1) * dimx *
              Imin; // beginning of the spline coeffs array for Imin...Imax

  // first point k=0, left boundary values: are known from
  // calc_spline_boundaries()
  int n, KLKU = KL + KU;

  // system of matching equations:

  for (n = 0; n <= s; n++) {
    iunk = n;
    M[KLKU + ieqn + iunk * (LDAB - 1)] = 1.0;
    for (j = 0; j < dimy; j++)
      SC[ieqn + j * len] = Cbnd[n + j * (N + 1)];
    ieqn++;
  }

  // check if band system:
  //  if (!(fmax(0,iunk-KU) <= ieqn && ieqn <= fmin(len-1,iunk+KL)))
  //  {
  //      fprintf (stderr,"\nband system faeiled: ieqn=%d\tiunk=%d", ieqn,
  //      iunk);
  //  }

  // further points:
  int k, p, ind;

  double dx;

  for (k = 1; k < dimx - 1; k++) // over points
  {
    for (n = 0; n < N; n++) // over derivatives: matching equations
    {
      p = n;
      iunk = (k - 1) * (N + 1) + p;
      dx = 1.0;
      ind = n * (N + 1);
      M[KLKU + ieqn + iunk * (LDAB - 1)] = BC[p + ind] * dx;

      for (p = n + 1; p <= N; p++) // over powers
      {
        iunk = (k - 1) * (N + 1) + p;
        dx *= x[k] - x[k - 1];
        M[KLKU + ieqn + iunk * (LDAB - 1)] = BC[p + ind] * dx;
        // M[...] = BC[p+n*(N+1)]*pow(x[k]-x[k-1], p-n); p=n...N
      }

      iunk = k * (N + 1) + n;
      M[KLKU + ieqn + iunk * (LDAB - 1)] = -1.0;
      for (j = 0; j < dimy; j++)
        SC[ieqn + j * len] = 0.0; // rhs vector component
      ieqn++;
    }

    // C(k,0)=f(k): //interpolation condition
    iunk = k * (N + 1);
    M[KLKU + ieqn + iunk * (LDAB - 1)] = 1.0;
    for (j = 0; j < dimy; j++)
      SC[ieqn + j * len] = y[k + j * dimx];
    ieqn++;
  }

  // last point:
  k = dimx - 1;
  for (n = 0; n <= s; n++) {
    p = n;
    iunk = (k - 1) * (N + 1) + p;
    dx = 1.0;
    ind = n * (N + 1);
    M[KLKU + ieqn + iunk * (LDAB - 1)] = BC[p + ind] * dx;

    for (p = n + 1; p <= N; p++) {
      iunk = (k - 1) * (N + 1) + p;
      dx *= x[k] - x[k - 1];
      M[KLKU + ieqn + iunk * (LDAB - 1)] = BC[p + ind] * dx;
      // M[...] = BC[p+n*(N+1)]*pow(x[k]-x[k-1], p-n); p=n...N
    }

    iunk = k * (N + 1) + n;
    M[KLKU + ieqn + iunk * (LDAB - 1)] = -1.0;
    for (j = 0; j < dimy; j++)
      SC[ieqn + j * len] = 0.0;
    ieqn++;
  }

  // right boundary values: are known from calc_spline_boundaries()
  for (n = 0; n <= N; n++) {
    iunk = k * (N + 1) + n;
    M[KLKU + ieqn + iunk * (LDAB - 1)] = 1.0;
    for (j = 0; j < dimy; j++)
      SC[ieqn + j * len] = Cbnd[(N + 1) * dimy + n + j * (N + 1)];
    ieqn++;
  }

  // for linear system solver:
  int D = len, NRHS = dimy, LDB = D, INFO;

  int *IPIV = new int[D];

  dgbsv_(&D, &KL, &KU, &NRHS, M, &LDAB, IPIV, SC, &LDB, &INFO);
  if (INFO) {
    fprintf(stderr,
            "\nerror: calc_spline_coefficients: failed to solve the band "
            "spline system: INFO=%d.",
            INFO);
  }

  delete[] IPIV;

  return INFO;
}

/*******************************************************************/

void spline_free_(uintptr_t sid) {
  spline_data *sd = (spline_data *)(sid);

  if (sd && sd->BC) {
    delete[] sd->BC;
    sd->BC = NULL;
  }
  if (sd && sd->fac) {
    delete[] sd->fac;
    sd->fac = NULL;
  }

  delete sd;
}

/*******************************************************************/

inline void search_array(double x, int dimx, const double *xa, int *ind) {
  // 0 <= ind <= dimx-2
  // warning: returns dim-2 for x>=xa[dim-1] - Ok for splines!
  // be careful to change something here:
  if (x < xa[*ind]) {
    *ind = binary_search(x, xa, 0, *ind);
  } else if (x >= xa[*ind + 1]) {
    *ind = binary_search(x, xa, *ind, dimx - 1);
  } else {
  }
}

/*******************************************************************/

inline int binary_search(double x, const double *xa, int ilo, int ihi) {
  // warning: returns dim-2 for x>=xa[dim-1] - Ok for splines!
  int i;
  while (ihi > ilo + 1) {
    i = (ihi + ilo) / 2;
    if (xa[i] > x)
      ihi = i;
    else
      ilo = i;
  }
  return ilo;
}

/*******************************************************************/

inline int sign(double x) { return (x < 0.0) ? (-1) : ((x == 0) ? 0 : 1); }

/*******************************************************************/

inline void set_bc_array(int N, double *BC) {
  // computes C^k_n = n!/k!/(n-k)! coefficients for n=0..N, k=0..n
  int k, n;
  double tmp;

  for (n = 0; n <= N; n++) {
    tmp = 1.0;
    BC[n] = tmp; // C^0_n
    for (k = 1; k <= n; k++) {
      tmp *= double(n - k + 1) / double(k);
      BC[n + k * (N + 1)] = tmp;
    }
  }
}

/*******************************************************************/

inline void set_fac_array(int N, double *fac) {
  // computes fac^p_n = p(p-1)...(p-n+1) for p=0..N, n=0..p coefficients
  int p, n;
  double tmp;

  for (p = 0; p <= N; p++) {
    tmp = 1.0;
    fac[p] = tmp;
    for (n = 1; n <= p; n++) {
      tmp *= (p - n + 1);
      fac[p + n * (N + 1)] = tmp;
    }
  }
}

/*******************************************************************/

void spline_eval_(uintptr_t sid, int dimz, double *z, int Dmin, int Dmax,
                  int Imin, int Imax, double *R) {
  /*
  sid - spline id
  dimz - z array dimension
  z - array where spline values are needed
  Dmin - lowest derivative to compute
  Dmax - highest derivative to compute
  Imin - lowest function index to compute
  Imax - highest function index to compute
  R - interpolated functions and derivatives values at z grid
  */

  spline_data *sd = (spline_data *)(sid);

  int N = sd->N;

  int p, n, j, k, ic0, ic1, ir0, ir1, ind;

  int len = (N + 1) * (sd->dimx);

  int D1 = Dmax - Dmin + 1, D2 = D1 * (Imax - Imin + 1);

  double tmp;

  for (k = 0; k < dimz; k++) // over grid points z
  {
    // finds nearest left grid point:
    search_array(z[k], sd->dimx, sd->x, &(sd->ind));

    ic0 = (sd->ind) * (N + 1);
    ir0 = k * D2 - Dmin;

    for (j = Imin; j <= Imax; j++) // over y data arrays
    {
      ic1 = ic0 + j * len;         // index of a jth y at kth x in C
      ir1 = ir0 + (j - Imin) * D1; // index of Dmin deriv of jth y at kth z in R

      for (n = Dmin; n <= Dmax; n++) // over spline derivs
      {
        // horner evaluation scheme:
        ind = n * (N + 1);
        tmp = (sd->fac[N + ind]) * (sd->C[N + ic1]);

        for (p = N - n; p > 0; p--) {
          tmp = (sd->fac[p + n - 1 + ind]) * (sd->C[p + n - 1 + ic1]) +
                (z[k] - (sd->x[sd->ind])) * tmp;
        }

        R[n + ir1] = tmp; // R[n-Dmin+(j-Imin)*D1+k*D2]
      }
    }
  }
}

/*******************************************************************/

void spline_eval_d_(uintptr_t sid, int dimz, double *z, int Dmin, int Dmax,
                    int Imin, int Imax, double *R) {
  /*
  sid - spline id
  dimz - z array dimension
  z - array where spline values are needed
  Dmin - lowest derivative to compute
  Dmax - highest derivative to compute
  Imin - lowest function index to compute
  Imax - highest function index to compute
  R - interpolated functions and derivatives values at z grid - another order of
  output

  check for consistency of the parameters:
  */

  spline_data *sd = (spline_data *)(sid);

  int N = sd->N;

  int p, n, j, k, ind, ir0, ir1, ic0, ic1;

  int len = (N + 1) * (sd->dimx);

  int D1 = Imax - Imin + 1, D2 = D1 * (Dmax - Dmin + 1);

  double tmp;

  for (k = 0; k < dimz; k++) // over grid points z
  {
    // finds nearest left grid point:
    search_array(z[k], sd->dimx, sd->x, &(sd->ind));

    ic0 = (sd->ind) * (N + 1);
    ir0 = k * D2 - (Imin + D1 * Dmin);

    for (j = Imin; j <= Imax; j++) // over y data arrays
    {
      ic1 = ic0 + j * len; // index of a jth y at kth x in C
      ir1 = ir0 + j;

      for (n = Dmin; n <= Dmax; n++) // over spline derivs
      {
        // horner evaluation scheme:
        ind = n * (N + 1);
        tmp = (sd->fac[N + ind]) * (sd->C[N + ic1]);

        for (p = N - n; p > 0; p--) {
          tmp = (sd->fac[p + n - 1 + ind]) * (sd->C[p + n - 1 + ic1]) +
                (z[k] - (sd->x[sd->ind])) * tmp;
        }
        // the derivs order is changed after quants:
        R[ir1 + D1 * n] = tmp; // R[j-Imin+D1*(n-Dmin)+k*D2]
      }
    }
  }
}

/*******************************************************************/
