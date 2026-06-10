/*! \file quadpack_qagi.cpp
    \brief Adaptive quadrature over infinite ranges (QUADPACK dqagie).

    Translation of the public-domain Netlib QUADPACK routines dqagie,
    dqk15i, dqelg and dqpsrt (Piessens, de Doncker-Kapenga, Ueberhuber,
    Kahaner, "QUADPACK", Springer 1983).  The arithmetic order of the
    original is kept, so results agree with other faithful ports of
    QUADPACK to floating-point accuracy.  Thread-safe: all state lives
    in caller-owned stack or heap storage.
*/

//no C++ runtime dependencies here: parts of KAMEL link this file into
//executables driven by gfortran without the C++ standard library

#include <math.h>
#include <float.h>
#include <stdio.h>
#include <stdlib.h>

#include "adaptive_quad.h"

/*******************************************************************/

static const double epmach = DBL_EPSILON;
static const double uflow = DBL_MIN;
static const double oflow = DBL_MAX;

struct qagi_fun
{
    quad_integrand f;
    void *params;
};

/*******************************************************************/

/*dqk15i: 15-point Gauss-Kronrod rule on a subrange (a, b) of (0, 1),
  after mapping of the infinite range onto (0, 1); arrays are 1-based
  as in the Fortran original, entry 0 is unused*/

static void qk15i (const qagi_fun *fn, double boun, int inf, double a, double b,
                   double *result, double *abserr, double *resabs, double *resasc)
{
static const double wg[9] = {0.0,
    0.0, 0.129484966168869693270611432679082,
    0.0, 0.279705391489276667901467771423780,
    0.0, 0.381830050505118944950369775488975,
    0.0, 0.417959183673469387755102040816327};
static const double xgk[9] = {0.0,
    0.991455371120812639206854697526329, 0.949107912342758524526189684047851,
    0.864864423359769072789712788640926, 0.741531185599394439863864773280788,
    0.586087235467691130294144838258730, 0.405845151377397166906606412076961,
    0.207784955007898467600689403773245, 0.000000000000000000000000000000000};
static const double wgk[9] = {0.0,
    0.022935322010529224963732008058970, 0.063092092629978553290700663189204,
    0.104790010322250183839876322541518, 0.140653259715525918745189590510238,
    0.169004726639267902826583426598550, 0.190350578064785409913256402421014,
    0.204432940075298892414161999234649, 0.209482141084727828012999174891714};

double fv1[8], fv2[8];

double dinf = (inf < 1) ? (double) inf : 1.0;

double centr = 0.5*(a + b);
double hlgth = 0.5*(b - a);
double tabsc1 = boun + dinf*(1.0 - centr)/centr;
double fval1 = fn->f (tabsc1, fn->params);
if (inf == 2) fval1 = fval1 + fn->f (-tabsc1, fn->params);
double fc = (fval1/centr)/centr;

double resg = wg[8]*fc;
double resk = wgk[8]*fc;
*resabs = fabs (resk);

for (int j = 1; j <= 7; j++)
{
    double absc = hlgth*xgk[j];
    double absc1 = centr - absc;
    double absc2 = centr + absc;
    double t1 = boun + dinf*(1.0 - absc1)/absc1;
    double t2 = boun + dinf*(1.0 - absc2)/absc2;
    double f1 = fn->f (t1, fn->params);
    double f2 = fn->f (t2, fn->params);
    if (inf == 2) f1 = f1 + fn->f (-t1, fn->params);
    if (inf == 2) f2 = f2 + fn->f (-t2, fn->params);
    f1 = (f1/absc1)/absc1;
    f2 = (f2/absc2)/absc2;
    fv1[j] = f1;
    fv2[j] = f2;
    double fsum = f1 + f2;
    resg = resg + wg[j]*fsum;
    resk = resk + wgk[j]*fsum;
    *resabs = *resabs + wgk[j]*(fabs (f1) + fabs (f2));
}

double reskh = resk*0.5;
*resasc = wgk[8]*fabs (fc - reskh);

for (int j = 1; j <= 7; j++)
{
    *resasc = *resasc + wgk[j]*(fabs (fv1[j] - reskh) + fabs (fv2[j] - reskh));
}

*result = resk*hlgth;
*resasc = (*resasc)*hlgth;
*resabs = (*resabs)*hlgth;
*abserr = fabs ((resk - resg)*hlgth);

if (*resasc != 0.0 && *abserr != 0.0)
{
    *abserr = (*resasc)*fmin (1.0, pow (200.0*(*abserr)/(*resasc), 1.5));
}
if (*resabs > uflow/(50.0*epmach))
{
    *abserr = fmax ((epmach*50.0)*(*resabs), *abserr);
}
}

/*******************************************************************/

/*dqelg: Wynn epsilon algorithm; epstab is 1-based with 52 usable slots,
  res3la 1-based with 3 slots*/

static void qelg (int *n, double epstab[53], double *result, double *abserr,
                  double res3la[4], int *nres)
{
*nres = *nres + 1;
*abserr = oflow;
*result = epstab[*n];

if (*n < 3)
{
    *abserr = fmax (*abserr, 5.0*epmach*fabs (*result));
    return;
}

const int limexp = 50;

epstab[*n + 2] = epstab[*n];
int newelm = (*n - 1)/2;
epstab[*n] = oflow;
int num = *n;
int k1 = *n;
int converged = 0;

for (int i = 1; i <= newelm; i++)
{
    int k2 = k1 - 1;
    int k3 = k1 - 2;
    double res = epstab[k1 + 2];
    double e0 = epstab[k3];
    double e1 = epstab[k2];
    double e2 = res;
    double e1abs = fabs (e1);
    double delta2 = e2 - e1;
    double err2 = fabs (delta2);
    double tol2 = fmax (fabs (e2), e1abs)*epmach;
    double delta3 = e1 - e0;
    double err3 = fabs (delta3);
    double tol3 = fmax (e1abs, fabs (e0))*epmach;

    if (err2 <= tol2 && err3 <= tol3)
    {
        //e0, e1 and e2 are equal to within machine accuracy
        *result = res;
        *abserr = err2 + err3;
        converged = 1;
        break;
    }

    double e3 = epstab[k1];
    epstab[k1] = e1;
    double delta1 = e1 - e3;
    double err1 = fabs (delta1);
    double tol1 = fmax (e1abs, fabs (e3))*epmach;

    if (err1 <= tol1 || err2 <= tol2 || err3 <= tol3)
    {
        //two elements very close to each other: omit part of the table
        *n = i + i - 1;
        break;
    }

    double ss = 1.0/delta1 + 1.0/delta2 - 1.0/delta3;
    double epsinf = fabs (ss*e1);

    if (epsinf <= 1.0e-4)
    {
        //irregular behaviour in the table: omit part of the table
        *n = i + i - 1;
        break;
    }

    res = e1 + 1.0/ss;
    epstab[k1] = res;
    k1 = k1 - 2;
    double error = err2 + fabs (res - e2) + err3;
    if (error <= *abserr)
    {
        *abserr = error;
        *result = res;
    }
}

if (!converged)
{
    //shift the table
    if (*n == limexp) *n = 2*(limexp/2) - 1;
    int ib = ((num/2)*2 == num) ? 2 : 1;
    int ie = newelm + 1;
    for (int i = 1; i <= ie; i++)
    {
        int ib2 = ib + 2;
        epstab[ib] = epstab[ib2];
        ib = ib2;
    }
    if (num != *n)
    {
        int indx = num - *n + 1;
        for (int i = 1; i <= *n; i++)
        {
            epstab[i] = epstab[indx];
            indx = indx + 1;
        }
    }
    if (*nres < 4)
    {
        res3la[*nres] = *result;
        *abserr = oflow;
    }
    else
    {
        *abserr = fabs (*result - res3la[3]) + fabs (*result - res3la[2])
                  + fabs (*result - res3la[1]);
        res3la[1] = res3la[2];
        res3la[2] = res3la[3];
        res3la[3] = *result;
    }
}

*abserr = fmax (*abserr, 5.0*epmach*fabs (*result));
}

/*******************************************************************/

/*dqpsrt: maintain the descending ordering in the list of error
  estimates; all index arrays and counters are 1-based*/

static void qpsrt (int limit, int last, int *maxerr, double *ermax,
                   const double *elist, int *iord, int *nrmax)
{
if (last <= 2)
{
    iord[1] = 1;
    iord[2] = 2;
    *maxerr = iord[*nrmax];
    *ermax = elist[*maxerr];
    return;
}

double errmax = elist[*maxerr];

if (*nrmax != 1)
{
    int ido = *nrmax - 1;
    for (int i = 1; i <= ido; i++)
    {
        int isucc = iord[*nrmax - 1];
        if (errmax <= elist[isucc]) break;
        iord[*nrmax] = isucc;
        *nrmax = *nrmax - 1;
    }
}

int jupbn = last;
if (last > limit/2 + 2) jupbn = limit + 3 - last;
double errmin = elist[last];

int jbnd = jupbn - 1;
int ibeg = *nrmax + 1;
int i = ibeg;
int found = 0;

for (i = ibeg; i <= jbnd; i++)
{
    int isucc = iord[i];
    if (errmax >= elist[isucc])
    {
        found = 1;
        break;
    }
    iord[i - 1] = isucc;
}

if (!found)
{
    iord[jbnd] = *maxerr;
    iord[jupbn] = last;
}
else
{
    //insert errmin by traversing the list bottom-up
    iord[i - 1] = *maxerr;
    int k = jbnd;
    int placed = 0;
    for (int j = i; j <= jbnd; j++)
    {
        int isucc = iord[k];
        if (errmin < elist[isucc])
        {
            iord[k + 1] = last;
            placed = 1;
            break;
        }
        iord[k + 1] = isucc;
        k = k - 1;
    }
    if (!placed) iord[i] = last;
}

*maxerr = iord[*nrmax];
*ermax = elist[*maxerr];
}

/*******************************************************************/

/*dqagie: globally adaptive integration over an infinite range with
  extrapolation by the epsilon algorithm*/

static int qagie (const qagi_fun *fn, double bound, int inf,
                  double epsabs, double epsrel, size_t limit_in,
                  double *result, double *abserr)
{
int limit = (int) limit_in;

int ier, ierro, iroff1, iroff2, iroff3, ksgn, ktmin;
int last, maxerr, nres, nrmax, numrl2, id, jupbnd, k;
int extrap, noext;
double boun, dres, errbnd, errmax, erlast, errsum, erlarg, ertest;
double area, area1, area2, area12, a1, a2, b1, b2;
double error1, error2, erro12, defabs, defab1, defab2, resabs;
double reseps, abseps, correc, small;
double rlist2[53], res3la[4];

double *work = (double *) malloc (4*(limit + 1)*sizeof (double));
int *iord = (int *) malloc ((limit + 1)*sizeof (int));
if (!work || !iord)
{
    fprintf (stderr, "\nqagie: memory allocation failed: exit!");
    exit (1);
}

double *alist = work;
double *blist = work + (limit + 1);
double *rlist = work + 2*(limit + 1);
double *elist = work + 3*(limit + 1);

ier = 0;
*result = 0.0;
*abserr = 0.0;
alist[1] = 0.0;
blist[1] = 1.0;
rlist[1] = 0.0;
elist[1] = 0.0;
iord[1] = 0;
correc = 0.0;
small = 0.0;
erlarg = 0.0;
ertest = 0.0;

if (limit < 1 || (epsabs <= 0.0 && epsrel < fmax (50.0*epmach, 0.5e-28)))
{
    ier = 6;
    goto L999;
}

//first approximation to the integral: map onto (0, 1); if inf = 2 the
//integral is computed as the sum over (-infinity, 0) and (0, +infinity)
boun = bound;
if (inf == 2) boun = 0.0;
qk15i (fn, boun, inf, 0.0, 1.0, result, abserr, &defabs, &resabs);

last = 1;
rlist[1] = *result;
elist[1] = *abserr;
iord[1] = 1;
dres = fabs (*result);
errbnd = fmax (epsabs, epsrel*dres);
if (*abserr <= 100.0*epmach*defabs && *abserr > errbnd) ier = 2;
if (limit == 1) ier = 1;
if (ier != 0 || (*abserr <= errbnd && *abserr != resabs)
    || *abserr == 0.0) goto L130;

rlist2[1] = *result;
errmax = *abserr;
maxerr = 1;
area = *result;
errsum = *abserr;
*abserr = oflow;
nrmax = 1;
nres = 0;
ktmin = 0;
numrl2 = 2;
extrap = 0;
noext = 0;
ierro = 0;
iroff1 = 0;
iroff2 = 0;
iroff3 = 0;
ksgn = -1;
if (dres >= (1.0 - 50.0*epmach)*defabs) ksgn = 1;

//main do-loop
for (last = 2; last <= limit; last++)
{
    //bisect the subinterval with the nrmax-th largest error estimate
    a1 = alist[maxerr];
    b1 = 0.5*(alist[maxerr] + blist[maxerr]);
    a2 = b1;
    b2 = blist[maxerr];
    erlast = errmax;
    qk15i (fn, boun, inf, a1, b1, &area1, &error1, &resabs, &defab1);
    qk15i (fn, boun, inf, a2, b2, &area2, &error2, &resabs, &defab2);

    area12 = area1 + area2;
    erro12 = error1 + error2;
    errsum = errsum + erro12 - errmax;
    area = area + area12 - rlist[maxerr];
    if (defab1 == error1 || defab2 == error2) goto L15;
    if (fabs (rlist[maxerr] - area12) > 1.0e-5*fabs (area12)
        || erro12 < 0.99*errmax) goto L10;
    if (extrap) iroff2 = iroff2 + 1;
    if (!extrap) iroff1 = iroff1 + 1;
L10:
    if (last > 10 && erro12 > errmax) iroff3 = iroff3 + 1;
L15:
    rlist[maxerr] = area1;
    rlist[last] = area2;
    errbnd = fmax (epsabs, epsrel*fabs (area));

    //test for roundoff error and eventually set error flag
    if (iroff1 + iroff2 >= 10 || iroff3 >= 20) ier = 2;
    if (iroff2 >= 5) ierro = 3;

    //number of subintervals equals limit
    if (last == limit) ier = 1;

    //bad integrand behaviour at some points of the integration range
    if (fmax (fabs (a1), fabs (b2))
        <= (1.0 + 100.0*epmach)*(fabs (a2) + 1000.0*uflow)) ier = 4;

    //append the newly-created intervals to the list
    if (error2 > error1) goto L20;
    alist[last] = a2;
    blist[maxerr] = b1;
    blist[last] = b2;
    elist[maxerr] = error1;
    elist[last] = error2;
    goto L30;
L20:
    alist[maxerr] = a2;
    alist[last] = a1;
    blist[last] = b1;
    rlist[maxerr] = area2;
    rlist[last] = area1;
    elist[maxerr] = error2;
    elist[last] = error1;
L30:
    qpsrt (limit, last, &maxerr, &errmax, elist, iord, &nrmax);
    if (errsum <= errbnd) goto L115;
    if (ier != 0) goto L100;
    if (last == 2) goto L80;
    if (noext) goto L90;
    erlarg = erlarg - erlast;
    if (fabs (b1 - a1) > small) erlarg = erlarg + erro12;
    if (extrap) goto L40;

    //test whether the interval to be bisected next is the smallest one
    if (fabs (blist[maxerr] - alist[maxerr]) > small) goto L90;
    extrap = 1;
    nrmax = 2;
L40:
    if (ierro == 3 || erlarg <= ertest) goto L60;

    //the smallest interval has the largest error: before bisecting,
    //decrease the sum of the errors over the larger intervals
    id = nrmax;
    jupbnd = last;
    if (last > 2 + limit/2) jupbnd = limit + 3 - last;
    for (k = id; k <= jupbnd; k++)
    {
        maxerr = iord[nrmax];
        errmax = elist[maxerr];
        if (fabs (blist[maxerr] - alist[maxerr]) > small) goto L90;
        nrmax = nrmax + 1;
    }
L60:
    //perform extrapolation
    numrl2 = numrl2 + 1;
    rlist2[numrl2] = area;
    qelg (&numrl2, rlist2, &reseps, &abseps, res3la, &nres);
    ktmin = ktmin + 1;
    if (ktmin > 5 && *abserr < 1.0e-3*errsum) ier = 5;
    if (abseps >= *abserr) goto L70;
    ktmin = 0;
    *abserr = abseps;
    *result = reseps;
    correc = erlarg;
    ertest = fmax (epsabs, epsrel*fabs (reseps));
    if (*abserr <= ertest) goto L100;
L70:
    //prepare bisection of the smallest interval
    if (numrl2 == 1) noext = 1;
    if (ier == 5) goto L100;
    maxerr = iord[1];
    errmax = elist[maxerr];
    nrmax = 1;
    extrap = 0;
    small = small*0.5;
    erlarg = errsum;
    goto L90;
L80:
    small = 0.375;
    erlarg = errsum;
    ertest = errbnd;
    rlist2[2] = area;
L90:;
}

L100:
//set final result and error estimate
if (*abserr == oflow) goto L115;
if (ier + ierro == 0) goto L110;
if (ierro == 3) *abserr = *abserr + correc;
if (ier == 0) ier = 3;
if (*result != 0.0 && area != 0.0) goto L105;
if (*abserr > errsum) goto L115;
if (area == 0.0) goto L130;
goto L110;
L105:
if (*abserr/fabs (*result) > errsum/fabs (area)) goto L115;
L110:
//test on divergence
if (ksgn == -1
    && fmax (fabs (*result), fabs (area)) <= defabs*0.01) goto L130;
if (0.01 > (*result/area) || (*result/area) > 100.0
    || errsum > fabs (area)) ier = 6;
goto L130;
L115:
//compute global integral sum
*result = 0.0;
for (k = 1; k <= last; k++)
{
    *result = *result + rlist[k];
}
*abserr = errsum;
L130:
if (ier > 2) ier = ier - 1;
L999:
free (work);
free (iord);

return ier;
}

/*******************************************************************/

int quad_qagi (quad_integrand f, void *params,
               double epsabs, double epsrel, size_t limit,
               double *result, double *abserr)
{
qagi_fun fn = {f, params};

return qagie (&fn, 0.0, 2, epsabs, epsrel, limit, result, abserr);
}

/*******************************************************************/

int quad_qagiu (quad_integrand f, void *params, double a,
                double epsabs, double epsrel, size_t limit,
                double *result, double *abserr)
{
qagi_fun fn = {f, params};

return qagie (&fn, a, 1, epsabs, epsrel, limit, result, abserr);
}
