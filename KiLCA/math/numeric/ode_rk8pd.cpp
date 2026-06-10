/*! \file ode_rk8pd.cpp
    \brief The implementation of the functions declared in ode_rk8pd.h.

    Tableau: Prince & Dormand (1981), RK8(7)13M. Step-size control: the
    standard error control of adaptive embedded Runge-Kutta methods with
    scale D_i = eps_abs + eps_rel*|y_i|, safety factor 0.9 and growth
    limits [0.2, 5] (see e.g. Hairer, Norsett, Wanner, "Solving Ordinary
    Differential Equations I", sec. II.4).
*/

#include <math.h>
#include <float.h>
#include <stdio.h>
#include <stdlib.h>

#include "ode_rk8pd.h"

/*******************************************************************/

static const int nstages = 13;

static const double c[13] = {
    0.0,
    1.0/18.0,
    1.0/12.0,
    1.0/8.0,
    5.0/16.0,
    3.0/8.0,
    59.0/400.0,
    93.0/200.0,
    5490023248.0/9719169821.0,
    13.0/20.0,
    1201146811.0/1299019798.0,
    1.0,
    1.0};

static const double a2[1] = {1.0/18.0};
static const double a3[2] = {1.0/48.0, 1.0/16.0};
static const double a4[3] = {1.0/32.0, 0.0, 3.0/32.0};
static const double a5[4] = {5.0/16.0, 0.0, -75.0/64.0, 75.0/64.0};
static const double a6[5] = {3.0/80.0, 0.0, 0.0, 3.0/16.0, 3.0/20.0};
static const double a7[6] = {
    29443841.0/614563906.0, 0.0, 0.0, 77736538.0/692538347.0,
    -28693883.0/1125000000.0, 23124283.0/1800000000.0};
static const double a8[7] = {
    16016141.0/946692911.0, 0.0, 0.0, 61564180.0/158732637.0,
    22789713.0/633445777.0, 545815736.0/2771057229.0, -180193667.0/1043307555.0};
static const double a9[8] = {
    39632708.0/573591083.0, 0.0, 0.0, -433636366.0/683701615.0,
    -421739975.0/2616292301.0, 100302831.0/723423059.0, 790204164.0/839813087.0,
    800635310.0/3783071287.0};
static const double a10[9] = {
    246121993.0/1340847787.0, 0.0, 0.0, -37695042795.0/15268766246.0,
    -309121744.0/1061227803.0, -12992083.0/490766935.0, 6005943493.0/2108947869.0,
    393006217.0/1396673457.0, 123872331.0/1001029789.0};
static const double a11[10] = {
    -1028468189.0/846180014.0, 0.0, 0.0, 8478235783.0/508512852.0,
    1311729495.0/1432422823.0, -10304129995.0/1701304382.0,
    -48777925059.0/3047939560.0, 15336726248.0/1032824649.0,
    -45442868181.0/3398467696.0, 3065993473.0/597172653.0};
static const double a12[11] = {
    185892177.0/718116043.0, 0.0, 0.0, -3185094517.0/667107341.0,
    -477755414.0/1098053517.0, -703635378.0/230739211.0,
    5731566787.0/1027545527.0, 5232866602.0/850066563.0,
    -4093664535.0/808688257.0, 3962137247.0/1805957418.0, 65686358.0/487910083.0};
static const double a13[12] = {
    403863854.0/491063109.0, 0.0, 0.0, -5068492393.0/434740067.0,
    -411421997.0/543043805.0, 652783627.0/914296604.0,
    11173962825.0/925320556.0, -13158990841.0/6184727034.0,
    3936647629.0/1978049680.0, -160528059.0/685178525.0,
    248638103.0/1413531060.0, 0.0};

static const double *atab[13] = {
    0, a2, a3, a4, a5, a6, a7, a8, a9, a10, a11, a12, a13};

//8th order weights
static const double b8[13] = {
    14005451.0/335480064.0, 0.0, 0.0, 0.0, 0.0,
    -59238493.0/1068277825.0, 181606767.0/758867731.0, 561292985.0/797845732.0,
    -1041891430.0/1371343529.0, 760417239.0/1151165299.0,
    118820643.0/751138087.0, -528747749.0/2220607170.0, 1.0/4.0};

//embedded 7th order weights
static const double b7[13] = {
    13451932.0/455176623.0, 0.0, 0.0, 0.0, 0.0,
    -808719846.0/976000145.0, 1757004468.0/5645159321.0, 656045339.0/265891186.0,
    -3867574721.0/1518517206.0, 465885868.0/322736535.0,
    53011238.0/667516719.0, 2.0/45.0, 0.0};


/*******************************************************************/

ode_rk8pd_state::ode_rk8pd_state (size_t n, double epsabs, double epsrel)
    : dim (n), eps_abs (epsabs), eps_rel (epsrel)
{
double *work = (double *) malloc (16*n*sizeof (double));
if (!work)
{
    fprintf (stderr, "\node_rk8pd_state: memory allocation failed: exit!");
    exit (1);
}

k = work;          //13 stage derivative blocks
ytmp = work + 13*n;
y0 = work + 14*n;
yerr = work + 15*n;
}

/*******************************************************************/

ode_rk8pd_state::~ode_rk8pd_state ()
{
free (k);
}

/*******************************************************************/

//one step of size h from (t, y0) stored in s; writes y and s.yerr;
//s.k[0..dim) must hold f(t, y0) on entry
static int rk8pd_step (ode_rk8pd_state &s, ode_rhs f, void *params,
                       double t, double h, double y[])
{
size_t n = s.dim;

for (int i = 1; i < nstages; i++)
{
    const double *ai = atab[i];

    for (size_t j = 0; j < n; j++)
    {
        double sum = 0.0;
        for (int l = 0; l < i; l++)
        {
            if (ai[l] != 0.0) sum += ai[l]*s.k[l*n + j];
        }
        s.ytmp[j] = s.y0[j] + h*sum;
    }

    int status = f (t + c[i]*h, &s.ytmp[0], &s.k[i*n], params);
    if (status != ODE_SUCCESS) return status;
}

for (size_t j = 0; j < n; j++)
{
    double sum8 = 0.0, sum7 = 0.0;
    for (int i = 0; i < nstages; i++)
    {
        if (b8[i] != 0.0) sum8 += b8[i]*s.k[i*n + j];
        if (b7[i] != 0.0) sum7 += b7[i]*s.k[i*n + j];
    }
    y[j] = s.y0[j] + h*sum8;
    s.yerr[j] = h*(sum7 - sum8);
}

return ODE_SUCCESS;
}

/*******************************************************************/

//standard step-size control; returns -1 (decrease, reject), 0 (keep) or
//+1 (increase) and updates *h accordingly
static int control_hadjust (const ode_rk8pd_state &s, const double y[], double *h)
{
const double safety = 0.9;
const int ord = 8;

double rmax = 0.0;

for (size_t j = 0; j < s.dim; j++)
{
    double d0 = s.eps_abs + s.eps_rel*fabs (y[j]);
    double r = fabs (s.yerr[j])/fabs (d0);
    if (r > rmax) rmax = r;
}

if (rmax > 1.1)
{
    double r = safety/pow (rmax, 1.0/ord);
    if (r < 0.2) r = 0.2;
    *h *= r;
    return -1;
}

if (rmax < 0.5)
{
    double r = safety/pow (rmax, 1.0/(ord + 1.0));
    if (r > 5.0) r = 5.0;
    if (r < 1.0) r = 1.0;
    *h *= r;
    return 1;
}

return 0;
}

/*******************************************************************/

int ode_rk8pd_step_once (ode_rk8pd_state &s, ode_rhs f, void *params,
                         double t, double h, double y[], double yerr[])
{
size_t n = s.dim;

for (size_t j = 0; j < n; j++) s.y0[j] = y[j];

int status = f (t, y, &s.k[0], params);
if (status != ODE_SUCCESS) return status;

status = rk8pd_step (s, f, params, t, h, y);
if (status != ODE_SUCCESS) return status;

for (size_t j = 0; j < n; j++) yerr[j] = s.yerr[j];

return ODE_SUCCESS;
}

/*******************************************************************/

int ode_rk8pd_evolve (ode_rk8pd_state &s, ode_rhs f, void *params,
                      double *t, double t1, double *h, double y[])
{
size_t n = s.dim;

double t0 = *t;
double h0 = *h;

int final_step = 0;
double dt = t1 - t0;

for (size_t j = 0; j < n; j++) s.y0[j] = y[j];

//k1 is reused across step-size retries (t and y are unchanged)
int status = f (t0, y, &s.k[0], params);
if (status != ODE_SUCCESS) return status;

while (1)
{
    if ((dt >= 0.0 && h0 > dt) || (dt < 0.0 && h0 < dt))
    {
        h0 = dt;
        final_step = 1;
    }
    else
    {
        final_step = 0;
    }

    double h_used = h0;

    status = rk8pd_step (s, f, params, t0, h_used, y);
    if (status != ODE_SUCCESS)
    {
        for (size_t j = 0; j < n; j++) y[j] = s.y0[j];
        return status;
    }

    if (final_step)
    {
        *t = t1;
    }
    else
    {
        *t = t0 + h_used;
    }

    int adjust = control_hadjust (s, y, &h0);

    if (adjust < 0)
    {
        //retry only if the decrease is genuine and still advances time
        double t_curr = *t;
        double t_next = *t + h0;

        if (fabs (h0) < fabs (h_used) && t_next != t_curr)
        {
            for (size_t j = 0; j < n; j++) y[j] = s.y0[j];
            continue;
        }

        h0 = h_used; //no progress possible: accept and keep the step size
    }

    *h = h0; //suggestion for the next call

    return ODE_SUCCESS;
}
}
