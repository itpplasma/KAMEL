/*! \file adaptive_quad.cpp
    \brief The implementation of the functions declared in adaptive_quad.h.

    Nodes and weights are the public-domain Netlib QUADPACK dqk15/dqk21/dqk31
    data; the adaptive bisection and the error tests follow QUADPACK dqage.
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

//Kronrod nodes of the (2n+1)-point rule interlace the n Gauss nodes, so
//even-index entries (1-based) of the half array xgk are Gauss nodes
//(QUADPACK layout)

static const double xgk15[8] = {
    0.9914553711208126392068546975263, 0.9491079123427585245261896840479,
    0.8648644233597690727897127886409, 0.7415311855993944398638647732808,
    0.5860872354676911302941448382587, 0.4058451513773971669066064120770,
    0.2077849550078984676006894037732, 0.0};
static const double wgk15[8] = {
    0.02293532201052922496373200805897, 0.06309209262997855329070066318920,
    0.1047900103222501838398763225415, 0.1406532597155259187451895905102,
    0.1690047266392679028265834265986, 0.1903505780647854099132564024210,
    0.2044329400752988924141619992346, 0.2094821410847278280129991748917};
static const double wg15[4] = {
    0.1294849661688696932706114326791, 0.2797053914892766679014677714238,
    0.3818300505051189449503697754890, 0.4179591836734693877551020408163};

static const double xgk21[11] = {
    0.9956571630258080807355272806890, 0.9739065285171717200779640120845,
    0.9301574913557082260012071800595, 0.8650633666889845107320966884235,
    0.7808177265864168970637175783450, 0.6794095682990244062343273651149,
    0.5627571346686046833390000992727, 0.4333953941292471907992659431658,
    0.2943928627014601981311266031039, 0.1488743389816312108848260011297,
    0.0};
static const double wgk21[11] = {
    0.01169463886737187427806439606219, 0.03255816230796472747881897245939,
    0.05475589657435199603138130024458, 0.07503967481091995276704314091619,
    0.09312545458369760553506546508337, 0.1093871588022976418992105903258,
    0.1234919762620658510779581098311, 0.1347092173114733259280540017717,
    0.1427759385770600807970942731387, 0.1477391049013384913748415159721,
    0.1494455540029169056649364683898};
static const double wg21[5] = {
    0.06667134430868813759356880989333, 0.1494513491505805931457763396577,
    0.2190863625159820439955349342282, 0.2692667193099963550912269215695,
    0.2955242247147528701738929946513};

static const double xgk31[16] = {
    0.9980022986933970602851728401523, 0.9879925180204854284895657185866,
    0.9677390756791391342573479787843, 0.9372733924007059043077589477102,
    0.8972645323440819008825096564545, 0.8482065834104272162006483207742,
    0.7904185014424659329676492948179, 0.7244177313601700474161860546139,
    0.6509967412974169705337358953133, 0.5709721726085388475372267372539,
    0.4850818636402396806936557402324, 0.3941513470775633698972073709810,
    0.2991800071531688121667800242664, 0.2011940939974345223006283033946,
    0.1011420669187174990270742314474, 0.0};
static const double wgk31[16] = {
    0.005377479872923348987792051430128, 0.01500794732931612253837476307581,
    0.02546084732671532018687400101965, 0.03534636079137584622203794847836,
    0.04458975132476487660822729937328, 0.05348152469092808726534314723943,
    0.06200956780067064028513923096080, 0.06985412131872825870952007709915,
    0.07684968075772037889443277748266, 0.08308050282313302103828924728610,
    0.08856444305621177064727544369377, 0.09312659817082532122548687274735,
    0.09664272698362367850517990762759, 0.09917359872179195933239317348460,
    0.1007698455238755950449466626176, 0.1013300070147915490173747927675};
static const double wg31[8] = {
    0.03075324199611726835462839357720, 0.07036604748810812470926741645067,
    0.1071592204671719350118695466859, 0.1395706779261543144478047945110,
    0.1662692058169939335532008604812, 0.1861610000155622110268005618664,
    0.1984314853271115764561183264438, 0.2025782419255612728806201999675};

/*******************************************************************/

struct gk_rule_data
{
    int nh;            //pairs in the half node array
    const double *xgk; //Kronrod nodes in (0, 1], descending, plus center 0
    const double *wgk; //Kronrod weights
    int ngw;           //number of Gauss weights
    const double *wg;  //Gauss weights
    int center_gauss;  //center node carries the last Gauss weight
};

static const gk_rule_data rule15 = {8, xgk15, wgk15, 4, wg15, 1};
static const gk_rule_data rule21 = {11, xgk21, wgk21, 5, wg21, 0};
static const gk_rule_data rule31 = {16, xgk31, wgk31, 8, wg31, 1};

/*******************************************************************/

static void gk_apply (const gk_rule_data *rule, quad_integrand f, void *params,
                      double a, double b, double *result, double *abserr,
                      double *resabs, double *resasc)
{
double fv1[16], fv2[16];

double centr = 0.5*(a + b);
double hlgth = 0.5*(b - a);
double dhlgth = fabs (hlgth);

int nh = rule->nh;

double fc = f (centr, params);

double resg;
int npair_gauss;

if (rule->center_gauss)
{
    resg = rule->wg[rule->ngw - 1]*fc;
    npair_gauss = rule->ngw - 1;
}
else
{
    resg = 0.0;
    npair_gauss = rule->ngw;
}

double resk = rule->wgk[nh - 1]*fc;
*resabs = fabs (resk);

for (int j = 1; j <= npair_gauss; j++)
{
    int jtw = 2*j - 1; //0-based even Kronrod index = Gauss node
    double absc = hlgth*rule->xgk[jtw];
    double fval1 = f (centr - absc, params);
    double fval2 = f (centr + absc, params);
    fv1[jtw] = fval1;
    fv2[jtw] = fval2;
    double fsum = fval1 + fval2;
    resg += rule->wg[j - 1]*fsum;
    resk += rule->wgk[jtw]*fsum;
    *resabs += rule->wgk[jtw]*(fabs (fval1) + fabs (fval2));
}

for (int j = 1; j <= nh/2; j++)
{
    int jtwm1 = 2*j - 2; //0-based odd Kronrod index
    double absc = hlgth*rule->xgk[jtwm1];
    double fval1 = f (centr - absc, params);
    double fval2 = f (centr + absc, params);
    fv1[jtwm1] = fval1;
    fv2[jtwm1] = fval2;
    double fsum = fval1 + fval2;
    resk += rule->wgk[jtwm1]*fsum;
    *resabs += rule->wgk[jtwm1]*(fabs (fval1) + fabs (fval2));
}

double reskh = 0.5*resk;
*resasc = rule->wgk[nh - 1]*fabs (fc - reskh);

for (int j = 0; j < nh - 1; j++)
{
    *resasc += rule->wgk[j]*(fabs (fv1[j] - reskh) + fabs (fv2[j] - reskh));
}

*result = resk*hlgth;
*resabs *= dhlgth;
*resasc *= dhlgth;
*abserr = fabs ((resk - resg)*hlgth);

if (*resasc != 0.0 && *abserr != 0.0)
{
    *abserr = (*resasc)*fmin (1.0, pow (200.0*(*abserr)/(*resasc), 1.5));
}
if (*resabs > uflow/(50.0*epmach))
{
    *abserr = fmax (*abserr, 50.0*epmach*(*resabs));
}
}

/*******************************************************************/

static int qag_bisect (const gk_rule_data *rule, quad_integrand f, void *params,
                       double a, double b, double epsabs, double epsrel,
                       size_t limit, double *result, double *abserr)
{
int ierr = 0;

double *work = (double *) malloc (4*limit*sizeof (double));
if (!work)
{
    fprintf (stderr, "\nqag_bisect: memory allocation failed: exit!");
    exit (1);
}

double *alist = work;
double *blist = work + limit;
double *rlist = work + 2*limit;
double *elist = work + 3*limit;

alist[0] = a;
blist[0] = b;
rlist[0] = *result;
elist[0] = *abserr;

double area = *result;
double errsum = *abserr;

int iroff1 = 0, iroff2 = 0;

size_t last;

for (last = 2; last <= limit; last++)
{
    //bisect the interval with the largest error estimate
    size_t maxerr = 0;
    for (size_t i = 1; i < last - 1; i++)
    {
        if (elist[i] > elist[maxerr]) maxerr = i;
    }
    double errmax = elist[maxerr];

    double a1 = alist[maxerr];
    double b1 = 0.5*(alist[maxerr] + blist[maxerr]);
    double a2 = b1;
    double b2 = blist[maxerr];

    double area1, error1, resabs1, defab1;
    double area2, error2, resabs2, defab2;

    gk_apply (rule, f, params, a1, b1, &area1, &error1, &resabs1, &defab1);
    gk_apply (rule, f, params, a2, b2, &area2, &error2, &resabs2, &defab2);

    double area12 = area1 + area2;
    double erro12 = error1 + error2;

    errsum += erro12 - errmax;
    area += area12 - rlist[maxerr];

    if (defab1 != error1 && defab2 != error2)
    {
        if (fabs (rlist[maxerr] - area12) <= 1.0e-5*fabs (area12)
            && erro12 >= 0.99*errmax) iroff1++;
        if (last > 10 && erro12 > errmax) iroff2++;
    }

    rlist[maxerr] = area1;
    elist[maxerr] = error1;
    blist[maxerr] = b1;
    alist[last - 1] = a2;
    blist[last - 1] = b2;
    rlist[last - 1] = area2;
    elist[last - 1] = error2;

    double errbnd = fmax (epsabs, epsrel*fabs (area));
    if (errsum <= errbnd) break;

    if (iroff1 >= 6 || iroff2 >= 20) ierr = QUAD_EROUND;
    if (last == limit) ierr = QUAD_EMAXITER;
    if (fmax (fabs (a1), fabs (b2))
        <= (1.0 + 100.0*epmach)*(fabs (a2) + 1000.0*uflow)) ierr = QUAD_ESINGULAR;
    if (ierr != 0) break;
}

if (last > limit) last = limit;

//fresh sums avoid cancellation in the running errsum update
*result = 0.0;
*abserr = 0.0;
for (size_t i = 0; i < last; i++)
{
    *result += rlist[i];
    *abserr += elist[i];
}

free (work);

return ierr;
}

/*******************************************************************/

static int qag_adapt (const gk_rule_data *rule, quad_integrand f, void *params,
                      double a, double b, double epsabs, double epsrel,
                      size_t limit, double *result, double *abserr)
{
*result = 0.0;
*abserr = 0.0;

if (limit < 1 || (epsabs <= 0.0 && epsrel < fmax (50.0*epmach, 0.5e-28)))
{
    return QUAD_EINVAL;
}

double resabs, resasc;

gk_apply (rule, f, params, a, b, result, abserr, &resabs, &resasc);

double errbnd = fmax (epsabs, epsrel*fabs (*result));

if (*abserr <= 50.0*epmach*resabs && *abserr > errbnd) return QUAD_EROUND;
if (limit == 1 && *abserr > errbnd) return QUAD_EMAXITER;
if ((*abserr <= errbnd && *abserr != resasc) || *abserr == 0.0) return QUAD_SUCCESS;

return qag_bisect (rule, f, params, a, b, epsabs, epsrel, limit, result, abserr);
}

/*******************************************************************/

int quad_qag (quad_integrand f, void *params, double a, double b,
              double epsabs, double epsrel, size_t limit, int key,
              double *result, double *abserr)
{
const gk_rule_data *rule;

switch (key)
{
    case 15: rule = &rule15; break;
    case 21: rule = &rule21; break;
    case 31: rule = &rule31; break;
    default:
        *result = 0.0;
        *abserr = 0.0;
        return QUAD_EINVAL;
}

return qag_adapt (rule, f, params, a, b, epsabs, epsrel, limit, result, abserr);
}
