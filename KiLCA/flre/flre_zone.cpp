/*! \file flre_zone.cpp
    \brief The implementation of flre_zone class.
*/

#include "flre_zone.h"
#include "mode.h"
#include "shared.h"
#include "inout.h"
#include "eval_back.h"
#include "adaptive_grid_pol.h"
#include "rhs_func.h"
#include "solver.h"
#include "transforms.h"
#include "typedefs.h"
#include "calc_flre_quants.h"

/*****************************************************************************/

void flre_zone::read_settings (char * file)
{
read (file); //base class read function

//derived class specific:
FILE *in;

if ((in=fopen (file, "r"))==NULL)
{
    fprintf(stderr, "\nerror: hmedium_zone: read_settings: failed to open file %s\a\n", file);
    exit(0);
}

char *str_buf = new char[1024];

//skip lines (already readed and values are set):
read_line_2skip_it (in, &str_buf);
read_line_2skip_it (in, &str_buf);
read_line_2skip_it (in, &str_buf);
read_line_2skip_it (in, &str_buf);
read_line_2skip_it (in, &str_buf);
read_line_2skip_it (in, &str_buf);
read_line_2skip_it (in, &str_buf);
read_line_2skip_it (in, &str_buf);

read_line_2skip_it (in, &str_buf);
read_line_2get_int (in, &(flre_order));
read_line_2get_int (in, &(Nmax));
read_line_2get_int (in, &(gal_corr));
read_line_2get_int (in, &(N));
read_line_2get_int (in, &(max_dim_c));
read_line_2get_double (in, &(D));
read_line_2get_double (in, &(eps_out));
read_line_2get_double (in, &(eps_res));
read_line_2get_int (in, &(hom_sys));
read_line_2skip_it (in, &str_buf);

//ODE solver settings:
read_line_2skip_it (in, &str_buf);
read_line_2get_int (in, &(max_dim));
read_line_2get_double (in, &(eps_rel));
read_line_2get_double (in, &(eps_abs));
read_line_2get_int (in, &(Nort));
read_line_2get_double (in, &(norm_fac));
read_line_2get_double (in, &(dr_out));
read_line_2get_double (in, &(dr_res));
read_line_2get_double (in, &(del));
read_line_2skip_it (in, &str_buf);

//Space out settings:
read_line_2skip_it (in, &str_buf);
read_line_2get_int (in, &(deg));
read_line_2get_double (in, &(reps));
read_line_2get_double (in, &(aeps));
read_line_2get_double (in, &(step));
read_line_2skip_it (in, &str_buf);

//Debugging settings:
read_line_2skip_it (in, &str_buf);
read_line_2get_int (in, &(flag_debug));
read_line_2skip_it (in, &str_buf);

//Collisions model flag:
read_line_2skip_it (in, &str_buf);
read_line_2get_int (in, &(collmod[0]));
read_line_2get_int (in, &(collmod[1]));
read_line_2skip_it (in, &str_buf);

fclose (in);

delete [] str_buf;

if (hom_sys != 0)
{
    fprintf (stdout, "warning: homogenious limit should not be normally used.");
    exit (1);
}

Nwaves = 6*flre_order - 2;
Ncomps = 6*flre_order + 7;

if (bc1 == BOUNDARY_CENTER   || bc1 == BOUNDARY_IDEALWALL ||
    bc2 == BOUNDARY_INFINITY || bc2 == BOUNDARY_IDEALWALL)
{
    Nfs = Nwaves/2;
}
else
{
    Nfs = Nwaves;
}

if (flag_debug) print_settings ();
}

/*****************************************************************************/

void flre_zone::print_settings (void)
{
print (); //base class print function

fprintf(stdout, "\nCheck for flre specific settings below:\n");

fprintf (stdout, "\norder of FLR expansion: %d", flre_order);
fprintf (stdout, "\nhighest cyclotron harmonic: %d", Nmax);
fprintf (stdout, "\nflag if to use correction term in conductivity: %d", gal_corr);
fprintf (stdout, "\nsplines degree: %d", N);
fprintf (stdout, "\nmaximum dimension of the grid for conductivity martices: %d", max_dim_c);
fprintf (stdout, "\nresonant layer width: %le", D);
fprintf (stdout, "\nerror parameter used for adaptive radial grid outside the resonant layer: %le", eps_out);
fprintf (stdout, "\nerror parameter used for the adaptive radial grid in the resonant layer: %le", eps_res);

fprintf (stdout, "\nflag for homogenious system: %d", hom_sys);
fprintf (stdout, "\nmax number of orthonormalization steps (ONS) for the solver: %d", Nort);
fprintf (stdout, "\ncontrolling factor for ONS by QR: %lg", norm_fac);
fprintf (stdout, "\noutput grid step outside the resonance region for the ME solutions: %lg", dr_out);
fprintf (stdout, "\noutput grid step inside the resonance region for the ME solutions: %lg", dr_res);
fprintf (stdout, "\nwidth of the resonance region: %lg", del);

fprintf (stdout, "\ncollisions model flags: %d %d", collmod[0], collmod[1]);

fprintf(stdout, "\n");
}

/*****************************************************************************/

void flre_zone::calc_basis_fields (int flag)
{
flre_zone *pointer = this;

if (flag) rsp = 0; else rsp = 1;

setup_flre_data_module_ (&pointer);

//calculates various data related to maxwell equations:
calc_and_set_maxwell_system_parameters_module_ ();

me = new maxwell_eqs_data (Nwaves);

copy_module_data_to_maxwell_eqs_data_struct (me);

if (flag_debug) me->print ();

allocate_and_set_conductivity_arrays_ ();

cp = new cond_profiles;

set_cond_profiles_in_mode_data_module_ (&cp);

cp->calc_and_spline_main_conductivity_profiles (this, flag);

if (flag)
{
    deallocate_conductivity_arrays_ ();

    clean_maxwell_system_parameters_module_ ();

    clean_flre_data_module_ ();

    return;
}

//allocates and calculates system matrix:
sp = new sysmat_profiles;

set_sysmat_profiles_in_mode_data_module_ (&sp);

sp->calc_and_spline_sysmatrix_profiles (this);

//calc dispersion if needed:
if (sd->os->flag_dispersion > 1)
{
    calc_dispersion ();
    save_dispersion ();
}

calculate_field_profiles_orth ();

//cleaning:
deallocate_conductivity_arrays_ ();

clean_maxwell_system_parameters_module_ ();

clean_flre_data_module_ ();

//transform basis to the lab frame:
//this is done later after spacing out, together with the evaluation of the full system
//vector, and the transformation to the lab frame.
//for stitching equations, fields (state vector) in the moving frame are used.
//check for consistency flre.f90!
}

/*****************************************************************************/

void flre_zone::copy_E_and_B_fields (double *EB_p)
{
complex<double> EBrsp[6], EBcyl[6];

for (int node=0; node<dim; node++)
{
    for (int comp=0; comp<6; comp++) //Er, Es, Ep, Br, Bs, Bp
    {
        EBrsp[comp] = EB[iF(node, comp, 0)] + I*EB[iF(node, comp, 1)];
    }

    transform_EB_from_rsp_to_cyl (bp, r[node], EBrsp, EBcyl);

    for (int comp=0; comp<6; comp++) //Er, Et, Ez, Br, Bt, Bz
    {
        EB_p[iFFM(node, comp, 0)] = real (EBcyl[comp]);
        EB_p[iFFM(node, comp, 1)] = imag (EBcyl[comp]);
    }
}
}

/*****************************************************************************/

void set_equations_settings_c_ (flre_zone **ptr, int *hom_sys, int *Nwaves, int *Nfs, int *Nphys, int *flag_debug)
{
flre_zone *zone = (flre_zone *)(*ptr);

*hom_sys = zone->hom_sys;

*Nwaves = zone->Nwaves;

*Nfs = zone->Nfs;

*Nphys = 0;

*flag_debug = zone->flag_debug;
}

/*****************************************************************************/

void set_conductivity_settings_c_ (flre_zone **ptr, int *flre_order, int *Nmax, int *gal_corr, int *rsp)
{
flre_zone *zone = (flre_zone *)(*ptr);

*flre_order = zone->flre_order;

*Nmax = zone->Nmax;

*gal_corr = zone->gal_corr;

*rsp = zone->rsp;
}

/*****************************************************************************/

void set_collisions_settings_c_ (flre_zone **ptr, int *collmod)
{
flre_zone *zone = (flre_zone *)(*ptr);

collmod[0] =  zone->collmod[0];
collmod[1] =  zone->collmod[1];
}

/*****************************************************************************/

void flre_zone::calculate_field_profiles_orth (void)
{
double *y = new double[2*Nwaves*Nwaves]; //state vectors at a point

int ind1 = 0, ind2 = 0;

double ri, rf; //initial and final r values for integration

if (bc1 == BOUNDARY_CENTER)
{
    ri = r1;
    rf = r2;

    //calc_start_values_anywhere_low_derivs_ (&ri, y);
    calc_start_values_center_with_correct_asymptotic_ (&ri, y);

    Nfs = Nwaves/2;
    ind1 = 0;     //min index of a fundamental solution which is regular at the boundary
    ind2 = Nfs-1; //max index of a fundamental solution which is regular at the boundary
}
else if (bc1 == BOUNDARY_IDEALWALL)
{
    ri = r1;
    rf = r2;

    calc_start_values_anywhere_low_derivs_ (&ri, y);
    //calc_start_values_ideal_wall_without_fake_modes_ (&ri, y);

    Nfs = Nwaves/2;
    ind1 = 0;     //min index of a fundamental solution which is regular at the boundary
    ind2 = Nfs-1; //max index of a fundamental solution which is regular at the boundary
}
else if (bc2 == BOUNDARY_INFINITY)
{
    ri = r2;
    rf = r1;

    //calc_start_values_anywhere_low_derivs_ (&ri, y);
    //calc_start_values_infinity_without_fake_modes_ (&ri, y);

    Nfs = Nwaves/2;
    ind1 = Nfs;      //min index of a fundamental solution which is regular at the boundary
    ind2 = Nwaves-1; //max index of a fundamental solution which is regular at the boundary

    fprintf (stdout, "\nflre_zone::calculate_field_profiles_orth: not implemented!");
    exit (1);
}
else if (bc2 == BOUNDARY_IDEALWALL)
{
    ri = r2;
    rf = r1;

    calc_start_values_anywhere_low_derivs_ (&ri, y);
    //calc_start_values_ideal_wall_without_fake_modes_ (&ri, y);

    Nfs = Nwaves/2;
    ind1 = 0;     //min index of a fundamental solution which is regular at the boundary
    ind2 = Nfs-1; //max index of a fundamental solution which is regular at the boundary
}
else
{
    //general case:
    ri = r1;
    rf = r2;

    //calc_start_values_general_without_fake_modes_ (&ri, y);

    Nfs = Nwaves;
    ind1 = 0;        //min index of a fundamental solution which is regular at the boundary
    ind2 = Nwaves-1; //max index of a fundamental solution which is regular at the boundary

    fprintf (stdout, "\nflre_zone::calculate_field_profiles_orth: unknown BCs!");
    exit (1);
}

//number of equations without an approximation for short modes:
int Neq = 2*Nwaves*Nfs;

//radial grid generation for a zone with resonance:
double *grid;

double sgn = signum(rf-ri);

if (r1 <= wd->r_res && r2 >= wd->r_res)
{
    double r_res = wd->r_res;

    double r = ri;

    int node = 0;

    while ((r-rf)*sgn < 0)
    {
        r += ((dr_out-dr_res)*(1.0 - exp(-(r-r_res)*(r-r_res)/del/del)) + dr_res)*sgn;
        node++;
    }

    dim = node+1; //last point should be replaced by rf

    grid = new double[dim];

    r = ri;

    node = 0;
    while ((r-rf)*sgn < 0)
    {
        grid[node++] = r;

        r += ((dr_out-dr_res)*(1.0 - exp(-(r-r_res)*(r-r_res)/del/del)) + dr_res)*sgn;
    }

    grid[node] = rf;
}
else
{
    dim = (int) abs(rf-ri)/dr_out;

    grid = new double[dim];

    double dr = (rf - ri)/(dim-1);

    for (int node=0; node<dim; node++) grid[node] = ri + dr*node;

    grid[dim-1] = rf; //last point (to have exactly rf):
}

double *state = new double[dim*Neq];

for (int j=0; j<Nfs; j++) //over starting vectors
{
    galilean_transform_of_flre_state_vector (this, sd->bs->V_gal_sys, wd->olab, ri,
                                             &y[2*Nwaves*j], &state[2*Nwaves*j]);
}

//settings for solver:
double *Dmat = new double[2*Nwaves*(Nwaves + 2*Nfs)];

rhs_func_params params = {Nwaves, Nwaves, Nfs, Dmat, sp};

solver_settings ss = {Nort, eps_rel, eps_abs, norm_fac, flag_debug};

integrate_basis_vecs (&rhs_func, Nfs, Nwaves, dim, grid, state, &ss, (void *)&params);

delete [] Dmat;

//base zone class data allocation:
r = new double[dim];

int len = dim*Nwaves*Ncomps*2;

basis = new double[len];

for (int i=0; i<len; i++) basis[i] = 0.0e0;

int j1 = 0, j2 = 0, dj = 0;

if (grid[0] == r1 && grid[dim-1] == r2)
{
    j1 = 0;
    j2 = dim-1;
    dj = 1;
}
else if (grid[0] == r2 && grid[dim-1] == r1) //sort array in inverse order
{
    j1 = dim-1;
    j2 = 0;
    dj = -1;
}
else
{
    fprintf (stdout, "\nflre_zone::calculate_field_profiles_orth: error!");
    exit (1);
}

int k = 0;
for (int j=j1; (j-j2)*dj<=0; j+=dj)
{
    r[k] = grid[j];

    for (int ind=ind1; ind<=ind2; ind++)
    {
        //warning: basis array is filled (partly) by state vectors, not system vectors!
        state_to_system_copy (this, state + is(j,ind-ind1,0,0), basis + ib(k,ind,0,0));
    }
    k++;
}

delete [] grid;
delete [] state;
delete [] y;

//normalization of the basis functions:
normalize_flre_basis_ (&Ncomps, &Nwaves, &dim, basis, me->iErsp_sys, &ind1, &ind2, &j2);
}

/*******************************************************************/

void flre_zone::calc_final_fields (void)
{
//makes spacing out of the grid with the solution;
//evaluates the whole system vector from the existing state vector stored in EB array;
//transforms the solution to the lab frame;

//for spacing out:
int dimnew;
int   *ind = new int[dim];
double *rnew = new double[dim];
double *snew = new double[dim*Nwaves*2];
double *sold = new double[dim*Nwaves*2];

//here EB is in a mov frame
for (int k=0; k<dim; k++) system_to_state_copy (this, &EB[2*Ncomps*k], &sold[2*Nwaves*k]);

sparse_grid_polynom (r, dim, sold, 2*Nwaves, deg, &aeps, &reps, step, &dimnew, rnew, snew, ind);

int num_vars = me->num_vars;
int num_eqs = me->num_eqs;

double *rhs = new double[2*num_eqs];

for (int j=0; j<2*num_eqs; j++) rhs[j] = 0.0e0;

//store spaced out values:
flre_zone *ptr = this;

activate_fortran_modules_for_zone_ (&ptr);

dim = dimnew;

EB_mov = new double[dim*Ncomps*2]; //system vector in a mov frame

for (int k=0; k<dim; k++)
{
    r[k] = rnew[k];
    state2sys_ (&r[k], sd->bs->flag_back, &snew[2*Nwaves*k], &EB_mov[2*Ncomps*k], rhs, 1);
}

deactivate_fortran_modules_for_zone_ (&ptr);

//galilean transformation of the solution system vector EB_mov to lab frame:
for (int k=0; k<dim; k++)
{
    galilean_transform_of_flre_system_vector (this, -sd->bs->V_gal_sys, wd->omov,
                                              r[k], &EB_mov[2*Ncomps*k], &EB[2*Ncomps*k]);
}

if (DEBUG_FLAG) save_system_vector_in_mov_frame ();

//clean up:
delete [] rnew;
delete [] snew;
delete [] sold;
delete [] ind;
delete [] rhs;
}

/*****************************************************************************/

template <typename T1, typename T2, typename T3> inline void
calc_derive_of_func_product (int N, double *C, int n, const T1 *f, const T2 *g, T3 *fg)
{
//N - max n index in C^k_n array
//C - binomial coefficients (C^k_n = C[n+k*(N+1)]) array
//n - order of the der, n<=N
//f - array of the f derivs
//g - array of the g derivs
//fg - the value of the fg deriv of n-th order

*fg = 0; for (int k=0; k<=n; k++) *fg += C[n+k*(N+1)]*f[k]*g[n-k];
}

/*****************************************************************************/

inline void galilean_transform_of_flre_state_vector (const flre_zone *zone, double V, complex<double> omega, double r, const double *E1, double *E2)
{
int *dim_Ersp_state = (int *) zone->me->dim_Ersp_state;
int *iErsp_state    = (int *) zone->me->iErsp_state;

//max orders of the derivative for Er, Es, Ep:
int mo[3] = {dim_Ersp_state[0]-1, dim_Ersp_state[1]-1, dim_Ersp_state[2]-1};

//wave data:
int m = zone->cp->wd->m;
double kz = (zone->wd->n)/(zone->sd->bs->rtor);

int ho = mo[2]; //max der order (2N-1);

//binomial coefficients:
double *C = new double[(ho+1)*(ho+1)];
binomial_coefficients (ho, C);

double *htz = new double[2*(ho+1)];
double *ht = htz, *hz = htz + ho + 1;

eval_hthz (r, 0, ho, zone->bp, htz); //ht & hz derivs

double *mor = new double[ho+1]; //(m/r) derivatives

mor[0] = m/r; for (int o=1; o<=ho; o++) mor[o] = (-o/r)*mor[o-1];

double *ks = new double[ho+1];
double *kp = new double[ho+1];

double htmor, hzmor;

for (int o=0; o<=ho; o++)
{
    calc_derive_of_func_product (ho, C, o, ht, mor, &htmor); //ht*(m/r) derivatives
    calc_derive_of_func_product (ho, C, o, hz, mor, &hzmor); //hz*(m/r) derivatives

    ks[o] = hzmor - ht[o]*kz;
    kp[o] = htmor + hz[o]*kz;
}

complex<double> *E1c[3], *E2c[3]; //Er, Es, Ep derivative arrays in 2 frames

for (int k=0; k<3; k++)
{
    E1c[k] = new complex<double>[mo[k]+1];
    E2c[k] = new complex<double>[mo[k]+1];
}

//getting complex fields from input double array in frame 1:
for (int k=0; k<3; k++) //over (rsp)-components
{
    for (int o=0; o<=mo[k]; o++)
    {
        int ind = 2*(iErsp_state[k] + o);
        E1c[k][o] = E1[ind + 0] + I*E1[ind + 1];
    }
}

complex<double> *Br   = new complex<double>[ho+1];
complex<double> *Ez   = new complex<double>[ho+1];
complex<double> *htBr = new complex<double>[ho+1];
complex<double> *hzBr = new complex<double>[ho+1];

complex<double> kpEs, ksEp, hzEp, htEs;

//aux quants in frame 1:
for (int o=0; o<=ho; o++)
{
    calc_derive_of_func_product (ho, C, o, kp, E1c[1], &kpEs);
    calc_derive_of_func_product (ho, C, o, ks, E1c[2], &ksEp);

    Br[o] = (c/omega)*(ksEp - kpEs);

    calc_derive_of_func_product (ho, C, o, ht, Br, &htBr[o]);
    calc_derive_of_func_product (ho, C, o, hz, Br, &hzBr[o]);

    calc_derive_of_func_product (ho, C, o, ht, E1c[1], &htEs);
    calc_derive_of_func_product (ho, C, o, hz, E1c[2], &hzEp);

    Ez[o] = hzEp - htEs;
}

//Er, Es, Ep transformation:
for (int o=0; o<=mo[0]; o++)
{
    E2c[0][o] = E1c[0][o] - V/(omega)*(kz*E1c[0][o] + I*Ez[o+1]);
}

for (int o=0; o<=mo[1]; o++)
{
    E2c[1][o] = E1c[1][o] + (V/c)*hzBr[o];
}

for (int o=0; o<=mo[2]; o++)
{
    E2c[2][o] = E1c[2][o] + (V/c)*htBr[o];
}

//Store new values to state vector:
for (int k=0; k<3; k++) //(r,s,p)
{
    for (int o=0; o<=mo[k]; o++) //der orders
    {
        int ind = 2*(iErsp_state[k] + o);

        E2[ind + 0] = real (E2c[k][o]);
        E2[ind + 1] = imag (E2c[k][o]);
    }
}

for (int k=0; k<3; k++)
{
    delete [] E1c[k];
    delete [] E2c[k];
}

delete [] C;
delete [] htz;
delete [] mor;
delete [] ks;
delete [] kp;
delete [] Br;
delete [] Ez;
delete [] htBr;
delete [] hzBr;
}

/*****************************************************************************/

inline void galilean_transform_of_flre_system_vector (const flre_zone *zone, double V, complex<double> omega, double r, const double *EB1, double *EB2)
{
int *dim_Ersp_sys = (int *) zone->me->dim_Ersp_sys;
int *iErsp_sys    = (int *) zone->me->iErsp_sys;

int *dim_Brsp_sys = (int *) zone->me->dim_Brsp_sys;
int *iBrsp_sys    = (int *) zone->me->iBrsp_sys;

//max orders of the derivative for Er, Es, Ep:
int moe[3] = {dim_Ersp_sys[0]-1, dim_Ersp_sys[1]-1, dim_Ersp_sys[2]-1};
int mob[3] = {dim_Brsp_sys[0]-1, dim_Brsp_sys[1]-1, dim_Brsp_sys[2]-1};

//wave data:
int m = zone->cp->wd->m;
double kz = (zone->wd->n)/(zone->sd->bs->rtor);

int ho = moe[2]; //max der order (2N);

//binomial coefficients:
double *C = new double[(ho+1)*(ho+1)];
binomial_coefficients (ho, C);

double *htz = new double[2*(ho+1)];
double *ht = htz, *hz = htz + ho + 1;

eval_hthz (r, 0, ho, zone->bp, htz); //ht & hz derivs

double *mor = new double[ho+1]; //(m/r) derivatives

mor[0] = m/r; for (int o=1; o<=ho; o++) mor[o] = (-o/r)*mor[o-1];

double *ks = new double[ho+1];
double *kp = new double[ho+1];

double htmor, hzmor;

for (int o=0; o<=ho; o++)
{
    calc_derive_of_func_product (ho, C, o, ht, mor, &htmor); //ht*(m/r) derivatives
    calc_derive_of_func_product (ho, C, o, hz, mor, &hzmor); //hz*(m/r) derivatives

    ks[o] = hzmor - ht[o]*kz;
    kp[o] = htmor + hz[o]*kz;
}

complex<double> *E1c[3], *E2c[3]; //Er, Es, Ep derivative arrays in 2 frames
complex<double> *B1c[3], *B2c[3]; //Br, Bs, Bp derivative arrays in 2 frames

for (int k=0; k<3; k++)
{
    E1c[k] = new complex<double>[moe[k]+1];
    E2c[k] = new complex<double>[moe[k]+1];

    B1c[k] = new complex<double>[mob[k]+1];
    B2c[k] = new complex<double>[mob[k]+1];
}

//getting complex fields from input double array in frame 1:
for (int k=0; k<3; k++) //over (rsp)-components
{
    for (int o=0; o<=moe[k]; o++)
    {
        int ind = 2*(iErsp_sys[k] + o);
        E1c[k][o] = EB1[ind + 0] + I*EB1[ind + 1];
    }

    for (int o=0; o<=mob[k]; o++)
    {
        int ind = 2*(iBrsp_sys[k] + o);
        B1c[k][o] = EB1[ind + 0] + I*EB1[ind + 1];
    }
}

//aux quants in frame 1:
complex<double> *Br   = new complex<double>[ho+1];
complex<double> *Ez   = new complex<double>[ho+1];
complex<double> *Et   = new complex<double>[ho+1];
complex<double> *htBr = new complex<double>[ho+1];
complex<double> *hzBr = new complex<double>[ho+1];
complex<double> *htEr = new complex<double>[ho+1];
complex<double> *hzEr = new complex<double>[ho+1];

complex<double> kpEs, ksEp, hzEp, htEs, hzEs, htEp;

for (int o=0; o<=ho; o++)
{
    calc_derive_of_func_product (ho, C, o, kp, E1c[1], &kpEs);
    calc_derive_of_func_product (ho, C, o, ks, E1c[2], &ksEp);

    Br[o] = (c/omega)*(ksEp - kpEs);

    calc_derive_of_func_product (ho, C, o, ht, Br, &htBr[o]);
    calc_derive_of_func_product (ho, C, o, hz, Br, &hzBr[o]);

    calc_derive_of_func_product (ho, C, o, ht, E1c[1], &htEs);
    calc_derive_of_func_product (ho, C, o, hz, E1c[2], &hzEp);

    calc_derive_of_func_product (ho, C, o, ht, E1c[2], &htEp);
    calc_derive_of_func_product (ho, C, o, hz, E1c[1], &hzEs);

    Et[o] = hzEs + htEp;
    Ez[o] = hzEp - htEs;
}

for (int o=0; o<=mob[2]; o++)
{
    calc_derive_of_func_product (ho, C, o, ht, E1c[0], &htEr[o]);
    calc_derive_of_func_product (ho, C, o, hz, E1c[0], &hzEr[o]);
}

//Er, Es, Ep transformation:
for (int o=0; o<=moe[0]; o++)
{
    E2c[0][o] = E1c[0][o] - V/(omega)*(kz*E1c[0][o] + I*Ez[o+1]);
}

for (int o=0; o<=moe[1]; o++)
{
    E2c[1][o] = E1c[1][o] + (V/c)*hzBr[o];
}

for (int o=0; o<=moe[2]; o++)
{
    E2c[2][o] = E1c[2][o] + (V/c)*htBr[o];
}

//Br, Bs, Bp transformation:
for (int o=0; o<=mob[0]; o++)
{
    B2c[0][o] = B1c[0][o];// + (V/c)*Et[o];
}

for (int o=0; o<=mob[1]; o++)
{
    B2c[1][o] = B1c[1][o];// - (V/c)*hzEr[o];
}

for (int o=0; o<=mob[2]; o++)
{
    B2c[2][o] = B1c[2][o];// - (V/c)*htEr[o];
}

//Store new values to system vector:
for (int k=0; k<3; k++) //(r,s,p)
{
    for (int o=0; o<=moe[k]; o++) //der orders
    {
        int ind = 2*(iErsp_sys[k] + o);

        EB2[ind + 0] = real (E2c[k][o]);
        EB2[ind + 1] = imag (E2c[k][o]);
    }

    for (int o=0; o<=mob[k]; o++) //der orders
    {
        int ind = 2*(iBrsp_sys[k] + o);

        EB2[ind + 0] = real (B2c[k][o]);
        EB2[ind + 1] = imag (B2c[k][o]);
    }
}

for (int k=0; k<3; k++)
{
    delete [] E1c[k];
    delete [] E2c[k];
    delete [] B1c[k];
    delete [] B2c[k];
}

delete [] C;
delete [] htz;
delete [] mor;
delete [] ks;
delete [] kp;
delete [] Br;
delete [] Ez;
delete [] Et;
delete [] htBr;
delete [] hzBr;
delete [] htEr;
delete [] hzEr;
}

/*****************************************************************************/

inline void system_to_state_copy (const flre_zone *zone, double *system, double *state)
{
for (int k=0; k<3; k++)
{
    for (int i=0; i<zone->me->dim_Ersp_state[k]; i++)
    {
        state[2*(zone->me->iErsp_state[k]+i)+0] = system[2*(zone->me->iErsp_sys[k]+i)+0];
        state[2*(zone->me->iErsp_state[k]+i)+1] = system[2*(zone->me->iErsp_sys[k]+i)+1];
    }
}
}

/*****************************************************************************/

inline void state_to_system_copy (const flre_zone *zone, double *state, double *system)
{
for (int k=0; k<3; k++)
{
    for (int i=0; i<zone->me->dim_Ersp_state[k]; i++)
    {
        system[2*(zone->me->iErsp_sys[k]+i)+0] = state[2*(zone->me->iErsp_state[k]+i)+0];
        system[2*(zone->me->iErsp_sys[k]+i)+1] = state[2*(zone->me->iErsp_state[k]+i)+1];
    }
}
}

/*****************************************************************************/

void get_sys_ind_array_ (flre_zone **ptr, int *sys_ind)
{
flre_zone *zone = (flre_zone *)(*ptr);

for (int k=0; k<zone->Nwaves; k++)
{
    sys_ind[k] = zone->me->sys_ind[k] + 1; //in Fortran index starts from 1
}
}

/*******************************************************************/

void activate_fortran_modules_for_zone_ (flre_zone **ptr)
{
flre_zone *zone = (flre_zone *)(*ptr);

setup_flre_data_module_ (&zone);

calc_and_set_maxwell_system_parameters_module_ ();

allocate_and_set_conductivity_arrays_ ();

set_cond_profiles_in_mode_data_module_ (&(zone->cp));

set_sysmat_profiles_in_mode_data_module_ (&(zone->sp));
}

/*******************************************************************/

void deactivate_fortran_modules_for_zone_ (flre_zone **ptr)
{
flre_zone *zone = (flre_zone *)(*ptr);

deallocate_conductivity_arrays_ ();

clean_maxwell_system_parameters_module_ ();

clean_flre_data_module_ ();
}

/*******************************************************************/

void flre_zone::save_system_vector_in_mov_frame (void)
{
char *path2linear = new char[1024];

eval_path_to_linear_data (sd->path2project, wd->m, wd->n, wd->olab, path2linear);

char *fname = new char[1024];

sprintf (fname, "%szone_%d_EB_mov.dat", path2linear, index);

FILE *out;
if (!(out = fopen (fname, "w")))
{
    fprintf (stderr, "\nFailed to open file %s\a\n", fname);
}

for (int i=0; i<dim; i++)
{
    fprintf (out, "%.16le", r[i]);

    for (int comp=0; comp<Ncomps; comp++)
    {
       fprintf (out, "\t%.16le\t%.16le", EB_mov[iEB(i, comp, 0)], EB_mov[iEB(i, comp, 1)]);
    }
     fprintf (out, "\n");
}

fclose (out);

delete [] path2linear;
delete [] fname;
}

/*****************************************************************************/

void calc_flre_basis_in_lab_cyl_frame_with_full_system_vectors_ (flre_zone **ptr, double *r, double *EB1, double *EB2)
{
flre_zone *zone = (flre_zone *)(*ptr);

double *rhs = new double[2*(zone->me->num_eqs)];

for (int j=0; j<2*(zone->me->num_eqs); j++) rhs[j] = 0.0e0;

double *state   = new double[2*(zone->Nwaves)];
double *sys_mov = new double[2*(zone->Ncomps)];
double *sys_lab = new double[2*(zone->Ncomps)];

activate_fortran_modules_for_zone_ (ptr);

for (int k=0; k<zone->Nwaves; k++) //over basis system vectors:
{
    int ind = zone->ib(0, k, 0, 0);

    //calc full flre system vector in mov frame:
    system_to_state_copy (zone, EB1 + ind, state);

    state2sys_ (r, zone->sd->bs->flag_back, state, sys_mov, rhs, 1);

    galilean_transform_of_flre_system_vector (zone, - zone->sd->bs->V_gal_sys,
    zone->wd->omov, *r, sys_mov, sys_lab);

    transform_of_flre_system_vector_to_cyl_coordinates (zone, *r, sys_lab, EB2 + ind);
}

deactivate_fortran_modules_for_zone_ (ptr);

delete [] rhs;
delete [] state;
delete [] sys_mov;
delete [] sys_lab;
}

/*****************************************************************************/

inline void transform_of_flre_system_vector_to_cyl_coordinates (const flre_zone *zone, double r, const double *EB1, double *EB2)
{
int *dim_Ersp_sys = (int *) zone->me->dim_Ersp_sys;
int *iErsp_sys    = (int *) zone->me->iErsp_sys;

int *dim_Brsp_sys = (int *) zone->me->dim_Brsp_sys;
int *iBrsp_sys    = (int *) zone->me->iBrsp_sys;

//max orders of the derivative for Er, Es, Ep:
int moe[3] = {dim_Ersp_sys[0]-1, dim_Ersp_sys[1]-1, dim_Ersp_sys[2]-1};
int mob[3] = {dim_Brsp_sys[0]-1, dim_Brsp_sys[1]-1, dim_Brsp_sys[2]-1};

int ho = moe[2]; //max der order (2N) of (s, p or theta, z) components;

//binomial coefficients:
double *C = new double[(ho+1)*(ho+1)];
binomial_coefficients (ho, C);

double *htz = new double[2*(ho+1)];
double *ht = htz, *hz = htz + ho + 1;

eval_hthz (r, 0, ho, zone->bp, htz); //ht & hz derivs

complex<double> *E1c[3], *E2c[3]; //Er, Es, Ep and Er, Et, Ez derivative arrays
complex<double> *B1c[3], *B2c[3]; //Br, Bs, Bp and Br, Br, Bz derivative arrays

for (int k=0; k<3; k++)
{
    E1c[k] = new complex<double>[moe[k]+1];
    E2c[k] = new complex<double>[moe[k]+1];

    B1c[k] = new complex<double>[mob[k]+1];
    B2c[k] = new complex<double>[mob[k]+1];
}

//getting complex fields from input double array in rsp:
for (int k=0; k<3; k++) //over (rsp)-components
{
    for (int o=0; o<=moe[k]; o++)
    {
        int ind = 2*(iErsp_sys[k] + o);
        E1c[k][o] = EB1[ind + 0] + I*EB1[ind + 1];
    }

    for (int o=0; o<=mob[k]; o++)
    {
        int ind = 2*(iBrsp_sys[k] + o);
        B1c[k][o] = EB1[ind + 0] + I*EB1[ind + 1];
    }
}

//aux quants:
complex<double> hzAp, htAs, hzAs, htAp;

//electric field transformation:
for (int o=0; o<=moe[0]; o++) E2c[0][o] = E1c[0][o]; //r components are the same

for (int o=0; o<=moe[1]; o++) //theta
{
    calc_derive_of_func_product (ho, C, o, hz, E1c[1], &hzAs);
    calc_derive_of_func_product (ho, C, o, ht, E1c[2], &htAp);

    E2c[1][o] = hzAs + htAp;
}

for (int o=0; o<=moe[2]; o++) //z
{
    calc_derive_of_func_product (ho, C, o, hz, E1c[2], &hzAp);
    calc_derive_of_func_product (ho, C, o, ht, E1c[1], &htAs);

    E2c[2][o] = hzAp - htAs;
}

//magnetic field transformation:
for (int o=0; o<=mob[0]; o++) B2c[0][o] = B1c[0][o]; //r components are the same

for (int o=0; o<=mob[1]; o++) //theta
{
    calc_derive_of_func_product (ho, C, o, hz, B1c[1], &hzAs);
    calc_derive_of_func_product (ho, C, o, ht, B1c[2], &htAp);

    B2c[1][o] = hzAs + htAp;
}

for (int o=0; o<=mob[2]; o++) //z
{
    calc_derive_of_func_product (ho, C, o, hz, B1c[2], &hzAp);
    calc_derive_of_func_product (ho, C, o, ht, B1c[1], &htAs);

    B2c[2][o] = hzAp - htAs;
}

//Store new values to system vector:
for (int k=0; k<3; k++) //(r,s,p)
{
    for (int o=0; o<=moe[k]; o++) //der orders
    {
        int ind = 2*(iErsp_sys[k] + o);

        EB2[ind + 0] = real (E2c[k][o]);
        EB2[ind + 1] = imag (E2c[k][o]);
    }

    for (int o=0; o<=mob[k]; o++) //der orders
    {
        int ind = 2*(iBrsp_sys[k] + o);

        EB2[ind + 0] = real (B2c[k][o]);
        EB2[ind + 1] = imag (B2c[k][o]);
    }
}

for (int k=0; k<3; k++)
{
    delete [] E1c[k];
    delete [] E2c[k];
    delete [] B1c[k];
    delete [] B2c[k];
}

delete [] C;
delete [] htz;
}

/*****************************************************************************/

void get_iersp_sys_array_ (flre_zone **ptr, int *iErsp_sys)
{
flre_zone *Z = (flre_zone *)(*ptr);

for (uchar k=0; k<3; k++) iErsp_sys[k] = Z->me->iErsp_sys[k] + 1;
}

/*****************************************************************************/

void get_ibrsp_sys_array_ (flre_zone **ptr, int *iBrsp_sys)
{
flre_zone *Z = (flre_zone *)(*ptr);

for (uchar k=0; k<3; k++) iBrsp_sys[k] = Z->me->iBrsp_sys[k] + 1;
}

/*****************************************************************************/

void flre_zone::calc_all_quants (void)
{
qp = new flre_quants (this);

calculate_JaE (qp);

qp->calculate_local_profiles ();

qp->calculate_integrated_profiles ();

transform_quants_to_lab_cyl_frame (qp); //transforms and saves
}

/*****************************************************************************/

void flre_zone::save_all_quants (void)
{
qp->save_profiles ();
}

/*****************************************************************************/

void flre_zone::eval_diss_power_density (double x, int type, int spec, double * dpd)
{
interp_diss_power_density (qp, x, type, spec, dpd);
}

/*****************************************************************************/

void flre_zone::eval_current_density (double x, int type, int spec, int comp, double * J)
{
interp_current_density (qp, x, type, spec, comp, J);
}

/*****************************************************************************/

void flre_zone::calc_dispersion (void)
{
dp = new disp_profiles (Nwaves, sp->dimx, sp->x, sd->bs->flag_back);

dp->calculate_dispersion_profiles ();
}

/*****************************************************************************/

void flre_zone::save_dispersion (void)
{
char * path2disp = new char[512];
eval_path_to_dispersion_data (sd->path2project, wd->m, wd->n, wd->olab, path2disp);

char * filename = new char[1024];
sprintf (filename, "%szone_%d_%s", path2disp, index, "kr.dat");

dp->save_dispersion_profiles (filename);

delete [] path2disp;
delete [] filename;
}

/*****************************************************************************/
