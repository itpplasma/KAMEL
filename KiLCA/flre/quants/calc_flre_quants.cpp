/*! \file calc_flre_quants.cpp
    \brief The implementation of functions declared in calc_flre_quants.h.
*/

#include "flre_zone.h"
#include "flre_quants.h"
#include "calc_flre_quants.h"
#include "inout.h"
#include "eval_back.h"
#include "eval_cond.h"
#include "shared.h"
#include "spline.h"
#include "interp.h"

/*******************************************************************/

void calculate_JaE (flre_quants *qp)
{
for (int k=0; k<qp->dimx; k++)
{
    qp->jaE[k]  = 0.0e0;
    qp->jaEi[k] = 0.0e0;
}

if (qp->zone->sd->as->wa == 0.0)
{
    calculate_JaE_delta (qp);
}
else
{
    calculate_JaE_distributed (qp);
}
}

/*******************************************************************/

void calculate_JaE_delta (flre_quants *qp)
{
int ia;

if (qp->zone->bc1 == BOUNDARY_ANTENNA)
{
    ia = 0;
}
else if (qp->zone->bc2 == BOUNDARY_ANTENNA)
{
    ia = qp->dimx-1;
}
else
{
    return;
}

//computes total absorbed power as a work of electric field on antenna current
double jsurf[4], jsurft[4];

current_density_ (jsurf);
cyl2rsp_ (&(qp->zone->sd->as->ra), jsurf, jsurf+2, jsurft, jsurft+2);

complex<double> ja[2] = {jsurft[0]+jsurft[1]*I, jsurft[2]+jsurft[3]*I}; //ja_s, ja_p

int Ncomps = qp->zone->Ncomps;

int const * iErsp_sys = qp->zone->me->iErsp_sys;

//electric field at the antenna location:
double Es_re = qp->zone->EB_mov[2*Ncomps*ia + 2*iErsp_sys[1] + 0];
double Es_im = qp->zone->EB_mov[2*Ncomps*ia + 2*iErsp_sys[1] + 1];
double Ep_re = qp->zone->EB_mov[2*Ncomps*ia + 2*iErsp_sys[2] + 0];
double Ep_im = qp->zone->EB_mov[2*Ncomps*ia + 2*iErsp_sys[2] + 1];

complex<double> Ef[2] = {Es_re+Es_im*I, Ep_re+Ep_im*I}; //Es_a, Ep_a

//total energy absorbed in the plasma: = - (work of Ef on Ja)
qp->JaE = 0.5*(qp->vol_fac)*(qp->x[ia])*real(ja[0]*conj(Ef[0]) + ja[1]*conj(Ef[1]));

#if DEBUG_FLAG
fprintf(stdout, "\nzone %d: - JaE: %le.\n", qp->zone->index, -qp->JaE);
#endif

if (qp->zone->sd->os->flag_emfield > 1)
{
    //save total absorbed energy:
    char *full_name = new char[1024];

    sprintf (full_name, "%szone_%d_JaE.dat", qp->path2linear, qp->zone->index);
    save_real_array (1, &(qp->x[ia]), &(qp->JaE), full_name);

    delete [] full_name;
}
}

/*******************************************************************/

void calculate_JaE_distributed (flre_quants *qp)
{
//!computes total absorbed power as a work of electric field on antenna current

double Ja_rsp[4];

complex<double> ja[2];
complex<double> Ef[2];

int Ncomps = qp->zone->Ncomps;

int const * iErsp_sys = qp->zone->me->iErsp_sys;

for (int k=0; k<qp->dimx; k++) //over r grid
{
    calc_current_density_r_s_p_ (&(qp->x[k]), Ja_rsp);

    ja[0] = Ja_rsp[0] + Ja_rsp[1]*I; //ja_s(r)
    ja[1] = Ja_rsp[2] + Ja_rsp[3]*I; //ja_p(r)

    //electric field:
    double Es_re = qp->zone->EB_mov[2*Ncomps*k + 2*iErsp_sys[1] + 0];
    double Es_im = qp->zone->EB_mov[2*Ncomps*k + 2*iErsp_sys[1] + 1];
    double Ep_re = qp->zone->EB_mov[2*Ncomps*k + 2*iErsp_sys[2] + 0];
    double Ep_im = qp->zone->EB_mov[2*Ncomps*k + 2*iErsp_sys[2] + 1];

    Ef[0] = Es_re + Es_im*I; //Es_a
    Ef[1] = Ep_re + Ep_im*I; //Ep_a

    qp->jaE[k] = 0.5*real(ja[0]*conj(Ef[0]) + ja[1]*conj(Ef[1]));
}

integrate_over_cylinder (qp->dimx, qp->x, qp->jaE, qp->vol_fac, qp->jaEi);

qp->JaE = qp->jaEi[qp->dimx-1];

#if DEBUG_FLAG
fprintf(stdout, "\nzone %d: - JaE: %le.\n", qp->zone->index, -qp->JaE);
#endif

if (qp->zone->sd->os->flag_emfield > 1)
{
    //save total absorbed energy:
    char *full_name = new char[1024];

    sprintf (full_name, "%szone_%d_JaE.dat", qp->path2linear, qp->zone->index);
    save_real_array (1, &(qp->x[qp->dimx-1]), &(qp->JaE), full_name);

    sprintf (full_name, "%szone_%d_jaE_dens.dat", qp->path2linear, qp->zone->index);
    save_real_array (qp->dimx, qp->x, qp->jaE, full_name);

    sprintf (full_name, "%szone_%d_jaE_int.dat", qp->path2linear, qp->zone->index);
    save_real_array (qp->dimx, qp->x, qp->jaEi, full_name);

    delete [] full_name;
}
}

/*******************************************************************/

void calc_current_density (flre_quants *qp)
{
//!computes amplitudes of the perturbed current density for the given r-node value

//check yourself that all needed quants are already computed:
if (!(qp->flagC))
{
    fprintf (stderr, "\n\aerror: consistency check failed in the function '%s' at line %d of the file '%s'.",  __FUNCTION__, __LINE__, __FILE__);
    return;
}

//it is fun to consider 1D array as ND one: the same trick is used everythere

//C matrix: spec=0:1, [type=0:1, s=0..2flreo, i=0:2, j=0:2], (re, im)*/
//here deriv order=0 - must match to call eval_all_C_matrices()
typedef double cond_C_matrix[2][2][2*(qp->flreo)+1][3][3][2];
cond_C_matrix &C = *((cond_C_matrix *)qp->C);

//EM field profiles: {0...dimx}, {0...num_vars}, {re,im}
typedef double em_field_data[qp->dimx][qp->zone->Ncomps][2];
em_field_data &F = *((em_field_data *)qp->zone->EB_mov);

complex<double> cd, Cm, Ef;

int i, j, spec, type, order;

//in arrays it is asumed that one has 2 types of the current density type=0:1

//current density for each r value: {{{re,im},comp={r,s,p},type={0,1},spec={i,e,t}}
typedef double curr_dens[3][2][3][2][qp->dimx];
curr_dens &CD = *((curr_dens *)(qp->qloc[qp->CURRENT_DENS])); //address of cd is treated as 5D array

for (spec=0; spec<2; spec++) //over species
{
    for (type=0; type<2; type++) //over current types
    {
        for (i=0; i<3; i++) //over current density components
        {
            cd = O; //complex zero
            for (j=0; j<3; j++) //over electric field (ef) components
            {
                for (order=0; order<=qp->zone->me->der_order[i][j]; order++) //derivatives
                {
                    //conductivity component for the given r node:
                    Cm = C[spec][type][order][i][j][0]+C[spec][type][order][i][j][1]*I;

                    //derivative of ef for the given r node:
                    Ef = F[qp->node][qp->zone->me->iErsp_sys[j]+order][0] +
                         F[qp->node][qp->zone->me->iErsp_sys[j]+order][1]*I;

                    //contribution to the curent:
                    cd += Cm*Ef;
                }
            }
            CD[spec][type][i][0][qp->node] = real(cd); //stores value to qloc array
            CD[spec][type][i][1][qp->node] = imag(cd); //stores value to qloc array
        }
    }
}

//total i+e current density: spec=2
for (type=0; type<2; type++) //over current types
{
    for (i=0; i<3; i++) //over current density components
    {
        CD[2][type][i][0][qp->node] = CD[0][type][i][0][qp->node] +
                                      CD[1][type][i][0][qp->node];

        CD[2][type][i][1][qp->node] = CD[0][type][i][1][qp->node] +
                                      CD[1][type][i][1][qp->node];
    }
}
qp->flag[qp->CURRENT_DENS] = 1; //cd is computed
}

/*******************************************************************/

void save_current_density (const flre_quants *qp)
{
char *filename = new char[1024];

char sort[3] = {'i','e','t'};
char comp[3] = {'r','s','p'};

FILE *outfile;

//current density for each r value: {{{re,im},comp={r,s,p},type={0,1},spec={i,e,t}}
typedef double curr_dens[3][2][3][2][qp->dimx];
curr_dens &CD = *((curr_dens *)(qp->qloc[qp->CURRENT_DENS])); //address of cd is treated as 5D array

int i, k, spec, type;

for (spec=0; spec<3; spec++) //over species
{
    for (type=0; type<2; type++) //over current types
    {
        for (i=0; i<3; i++) //over current density components
        {
            //file name:
            sprintf (filename, "%szone_%d_%s_dens_%c_%d_%c.dat", qp->path2linear, qp->zone->index, qp->name[qp->CURRENT_DENS], comp[i], type, sort[spec]);

            if (!(outfile = fopen (filename, "w")))
            {
                fprintf (stderr, "\nFailed to open file %s\a\n", filename);
            }

            for (k=0; k<qp->dimx; k++) //over r grid
            {
                fprintf (outfile, "%.16le\t%.16le\t%.16le\n",
                         qp->x[k], CD[spec][type][i][0][k], CD[spec][type][i][1][k]);
            }
            fclose (outfile);
        }
    }
}

delete [] filename;

#if DEBUG_FLAG
fprintf(stdout, "\n%s is saved.", qp->name[qp->CURRENT_DENS]);
#endif
}

/*******************************************************************/

void calc_absorbed_power_density (flre_quants *qp)
{
//check yourself that all needed quants are already computed:
if (!(qp->flag[qp->CURRENT_DENS]))
{
    fprintf (stderr, "\n\aerror: consistency check failed in the function '%s' at line %d of the file '%s'.",  __FUNCTION__, __LINE__, __FILE__);
    return;
}

//it is fun to consider 1D array as ND one: the same trick is used everythere

//EM field profiles: {0...dimx}, {0...num_vars}, {re,im}
typedef double em_field_data[qp->dimx][qp->zone->Ncomps][2];
em_field_data &F = *((em_field_data *)qp->zone->EB_mov);

//current density for each r value: {{{re,im},comp={r,s,p},type={0,1},spec={i,e,t}}
typedef double curr_dens[3][2][3][2][qp->dimx];
curr_dens &CD = *((curr_dens *)(qp->qloc[qp->CURRENT_DENS])); //address of cd is treated as 5D array

//absorbed power density each r value: //{{re},type={0,1},spec={i,e,t}}
typedef double abs_pow_dens[3][2][qp->dimx];

//address of abs pow dens in qloc array treated as 3D array:
abs_pow_dens &APD = *((abs_pow_dens *)(qp->qloc[qp->ABS_POWER_DENS]));

complex<double> cd, Ef;
double apd;

for (int spec=0; spec<3; spec++) //over species
{
    for (int type=0; type<2; type++) //over current types
    {
        apd = 0.0;
        for (int i=0; i<3; i++) //over components
        {
            cd = CD[spec][type][i][0][qp->node] + CD[spec][type][i][1][qp->node]*I;

            Ef = F[qp->node][qp->zone->me->iErsp_sys[i]][0] +
                 F[qp->node][qp->zone->me->iErsp_sys[i]][1]*I;

            apd += 0.5*real(cd*conj(Ef));
        }
        APD[spec][type][qp->node] = apd;
    }
}
qp->flag[qp->ABS_POWER_DENS] = 1; //apd is computed
}

/*******************************************************************/

void calc_absorbed_power_in_cylinder (flre_quants *qp)
{
//check yourself that all needed quants are already computed:
if (!(qp->zone->sd->os->flag_quants[qp->ABS_POWER_DENS]))
{
    fprintf (stderr, "\n\aerror: consistency check failed in the function '%s' at line %d of the file '%s'.",  __FUNCTION__, __LINE__, __FILE__);
    return;
}

//it is fun to consider 1D array as ND one: the same trick is used everythere

//absorbed power at each r value: //{{re},type={0,1},spec={i,e,t}}
typedef double abs_pow[3][2][qp->dimx];

//address of abs pow dens in qloc array treated as 3D array:
abs_pow &APD = *((abs_pow *)(qp->qloc[qp->ABS_POWER_DENS]));

//address of integrated abs pow in qint array treated as 3D array:
abs_pow &API = *((abs_pow *)(qp->qint[qp->ABS_POWER_DENS]));

for (int spec=0; spec<3; spec++) //over species
{
    for (int type=0; type<2; type++) //over current types
    {
        integrate_over_cylinder (qp->dimx, qp->x, APD[spec][type], qp->vol_fac, API[spec][type]);
    }
}
}

/*******************************************************************/

void save_absorbed_power (const flre_quants *qp)
{
char *filename = new char[1024];

char sort[3] = {'i','e','t'};

//absorbed power at each r value: //{{re},type={0,1},spec={i,e,t}}
typedef double abs_pow[3][2][qp->dimx];

//address of abs pow dens in qloc array treated as 3D array:
abs_pow &APD = *((abs_pow *)(qp->qloc[qp->ABS_POWER_DENS]));

//address of integrated abs pow in qint array treated as 3D array:
abs_pow &API = *((abs_pow *)(qp->qint[qp->ABS_POWER_DENS]));

for (int spec=0; spec<3; spec++) //over species
{
    for (int type=0; type<2; type++) //over current types
    {
        sprintf (filename, "%szone_%d_%s_dens_%d_%c.dat", qp->path2linear, qp->zone->index, qp->name[qp->ABS_POWER_DENS], type, sort[spec]);

        save_real_array (qp->dimx, qp->x, APD[spec][type], filename);

        sprintf (filename, "%szone_%d_%s_int_%d_%c.dat", qp->path2linear, qp->zone->index, qp->name[qp->ABS_POWER_DENS], type, sort[spec]);

        save_real_array (qp->dimx, qp->x, API[spec][type], filename);
    }
}

delete [] filename;

#if DEBUG_FLAG
fprintf(stdout, "\n%s is saved.", qp->name[qp->ABS_POWER_DENS]);
#endif
}

/*******************************************************************/

void calc_dissipated_power_density (flre_quants *qp)
{
//calculates dissipated power density, see thesis p. 53
//check yourself that all needed quants are already computed:
if (!(qp->flagK))
{
    fprintf (stderr, "\n\aerror: consistency check failed in the function '%s' at line %d of the file '%s'.",  __FUNCTION__, __LINE__, __FILE__);
    return;
}

int flreo = qp->flreo;

//it is fun to consider 1D array as ND one: the same trick is used everythere

/*K: deriv=0:Dmax, spec=0:1, [type=0:1, p=0:flreo, q=0:flreo, i=0:2, j=0:2], (re, im)*/
//max derivatives order (1 dimension) must match to real call of eval_all_K_matrices!
typedef double cond_K_matrix[qp->Dmax+1][2][2][flreo+1][flreo+1][3][3][2];
cond_K_matrix &K = *((cond_K_matrix *)qp->K);

//EM field profiles: {0...dimx}, {0...num_vars}, {re,im}
typedef double em_field_data[qp->dimx][qp->zone->Ncomps][2];
em_field_data &F = *((em_field_data *)qp->zone->EB_mov);

//dissipated power density for each r value: ////{{re},type={0,1,2=1-0},spec={i,e,t}}
typedef double diss_pow_dens[3][3][qp->dimx];

//address of diss pow dens in qloc array treated as 3D array:
diss_pow_dens &DPD = *((diss_pow_dens *)(qp->qloc[qp->DISS_POWER_DENS]));

complex<double> dpd, Km, Ef_1, Ef_2, fac;

double *charge = qp->zone->sd->bs->charge;

int i, j, n1, n2, spec, type;

for (spec=0; spec<2; spec++) //over species
{
    fac = pi*(charge[spec])*(charge[spec])/(qp->zone->wd->omov)/(qp->r);

    for (type=0; type<2; type++) //over current types
    {
        dpd = O;
        for (n1=0; n1<=flreo; n1++)
        {
            for (n2=0; n2<=flreo; n2++)
            {
                for (i=0; i<3; i++)
                {
                    for (j=0; j<3; j++)
                    {
                        //conductivity K matrix for the given r node:
                        Km = K[0][spec][type][n1][n2][i][j][0] +
                             K[0][spec][type][n1][n2][i][j][1]*I;

                        //derivative of ef for the given r node:
                        Ef_1 = F[qp->node][qp->zone->me->iErsp_sys[i]+n1][0] +
                               F[qp->node][qp->zone->me->iErsp_sys[i]+n1][1]*I;

                        Ef_2 = F[qp->node][qp->zone->me->iErsp_sys[j]+n2][0] +
                               F[qp->node][qp->zone->me->iErsp_sys[j]+n2][1]*I;

                        dpd += Km*conj(Ef_1)*Ef_2;
                    }
                }
            }
        }
        DPD[spec][type][qp->node] = real(fac*I*dpd);
    }
}

//total i+e dissipated power density: spec=2
for (type=0; type<2; type++) //over current types
{
    DPD[2][type][qp->node] = DPD[0][type][qp->node] + DPD[1][type][qp->node];
}

for (spec=0; spec<3; spec++) //over species to remove contribution of j_0
{
    DPD[spec][2][qp->node] = DPD[spec][1][qp->node] - DPD[spec][0][qp->node];
}

qp->flag[qp->DISS_POWER_DENS] = 1; //dpd is computed
}

/*******************************************************************/

void calc_dissipated_power_in_cylinder (flre_quants *qp)
{
//check yourself that all needed quants are already computed:
if (!(qp->zone->sd->os->flag_quants[qp->DISS_POWER_DENS]))
{
    fprintf (stderr, "\n\aerror: consistency check failed in the function '%s' at line %d of the file '%s'.",  __FUNCTION__, __LINE__, __FILE__);
    return;
}

//it is fun to consider 1D array as ND one: the same trick is used everythere

//dissipated power for each r value: ////{{re},type={0,1,2=1-0},spec={i,e,t}}
typedef double diss_pow[3][3][qp->dimx];

//address of diss pow density in qloc array treated as 3D array:
diss_pow &DPD = *((diss_pow *)(qp->qloc[qp->DISS_POWER_DENS]));

//address of integrated abs pow in qint array treated as 3D array:
diss_pow &DPI = *((diss_pow *)(qp->qint[qp->DISS_POWER_DENS]));

for (int spec=0; spec<3; spec++) //over species
{
    for (int type=0; type<3; type++) //over current types
    {
        integrate_over_cylinder (qp->dimx, qp->x, DPD[spec][type], qp->vol_fac, DPI[spec][type]);
    }
}
}

/*******************************************************************/

void save_dissipated_power (const flre_quants *qp)
{
char *filename = new char[1024];

char sort[3] = {'i','e','t'};

//dissipated power for each r value: ////{{re},type={0,1,2=1-0},spec={i,e,t}}
typedef double diss_pow[3][3][qp->dimx];

//address of diss pow density in qloc array treated as 3D array:
diss_pow &DPD = *((diss_pow *)(qp->qloc[qp->DISS_POWER_DENS]));

//address of integrated abs pow in qint array treated as 3D array:
diss_pow &DPI = *((diss_pow *)(qp->qint[qp->DISS_POWER_DENS]));

for (int spec=0; spec<3; spec++) //over species
{
    for (int type=0; type<3; type++) //over current types
    {
        sprintf (filename, "%szone_%d_%s_dens_%d_%c.dat", qp->path2linear, qp->zone->index, qp->name[qp->DISS_POWER_DENS], type, sort[spec]);

        save_real_array (qp->dimx, qp->x, DPD[spec][type], filename);

        sprintf (filename, "%szone_%d_%s_int_%d_%c.dat", qp->path2linear, qp->zone->index, qp->name[qp->DISS_POWER_DENS], type, sort[spec]);

        save_real_array (qp->dimx, qp->x, DPI[spec][type], filename);
    }
}

delete [] filename;

#if DEBUG_FLAG
fprintf(stdout, "\n%s is saved.", qp->name[qp->DISS_POWER_DENS]);
#endif
}

/*******************************************************************/

void calc_kinetic_flux (flre_quants *qp)
{
//calculates radial component of kinetic flux out of the cylinder of radius r,
//see thesis p. 53 for an expression for kinetic flux density

//check yourself that all needed quants are already computed:
if (!(qp->flagK))
{
    fprintf (stderr, "\n\aerror: consistency check failed in the function '%s' at line %d of the file '%s'.",  __FUNCTION__, __LINE__, __FILE__);
    return;
}

int const & flreo = qp->flreo;

//it is fun to consider 1D array as ND one: the same trick is used everythere

/*K: deriv=0:Dmax, spec=0:1, [type=0:1, p=0:flreo, q=0:flreo, i=0:2, j=0:2], (re, im)*/
//max derivatives order (1 dimension) must match to real call of eval_all_K_matrices!
typedef double cond_K_matrix[qp->Dmax+1][2][2][flreo+1][flreo+1][3][3][2];
cond_K_matrix &K = *((cond_K_matrix *)qp->K);

//EM field profiles: {0...dimx}, {0...num_vars}, {re,im}
typedef double em_field_data[qp->dimx][qp->zone->Ncomps][2];
em_field_data &F = *((em_field_data *)qp->zone->EB_mov);

//kinetic flux density for each r value: ////{{re},type={0,1,2=1-0},spec={i,e,t}}
typedef double kin_flux[3][3][qp->dimx];

//address of kinetic flux dens in qloc array treated as 3D array:
kin_flux &KF = *((kin_flux *)(qp->qloc[qp->KIN_FLUX]));

//binomial coefficients: //C^k_n = n!/k!/(n-k)! coefficients for n=0..flreo, k=0..n
typedef double binom_coeffs[flreo+1][flreo+1]; //must match to allocation and calc!
binom_coeffs &bico = *((binom_coeffs *)(qp->bico));

complex<double> kf, Km, Ef_1, Ef_2, fac;

double coeff;

double *charge = qp->zone->sd->bs->charge;

int i, j, n1, n2, spec, type, p, s;

for (spec=0; spec<2; spec++) //over species
{
    fac = pi*(charge[spec])*(charge[spec])/(qp->zone->wd->omov)/(qp->r);

    for (type=0; type<2; type++) //over current types
    {
        kf = O;
        for (n1=0; n1<=flreo; n1++)
        {
            for (n2=0; n2<=flreo; n2++)
            {
                for (p=0; p<=n1-1; p++)
                {
                    for (s=0; s<=n1-p-1; s++)
                    {
                        coeff = pow(-1.0, n1+p)*bico[s][n1-p-1]; //C^k_n->bico[k][n]

                        for (i=0; i<3; i++)
                        {
                            for (j=0; j<3; j++)
                            {
                                //conductivity K matrix for the given r node:
                                Km = K[n1-p-s-1][spec][type][n1][n2][i][j][0] +
                                     K[n1-p-s-1][spec][type][n1][n2][i][j][1]*I;

                                //derivative of ef for the given r node:
                                Ef_1 = F[qp->node][qp->zone->me->iErsp_sys[i]+p][0] +
                                       F[qp->node][qp->zone->me->iErsp_sys[i]+p][1]*I;

                                Ef_2 = F[qp->node][qp->zone->me->iErsp_sys[j]+n2+s][0] +
                                       F[qp->node][qp->zone->me->iErsp_sys[j]+n2+s][1]*I;

                                kf += coeff*Km*conj(Ef_1)*Ef_2;
                            }
                        }
                    }
                }
            }
        }
        KF[spec][type][qp->node] = real(fac*I*kf)*(qp->vol_fac)*(qp->r);
    }
}

//total i+e kinetic flux: spec=2
for (type=0; type<2; type++) //over current types
{
    KF[2][type][qp->node] = KF[0][type][qp->node] + KF[1][type][qp->node];
}

//removes contribution of j_0:
for (spec=0; spec<3; spec++)
{
    KF[spec][2][qp->node] = KF[spec][1][qp->node] - KF[spec][0][qp->node];
}

qp->flag[qp->KIN_FLUX] = 1; //dpd is computed
}

/*******************************************************************/

void save_kinetic_flux (const flre_quants *qp)
{
char *filename = new char[1024];

char sort[3] = {'i','e','t'};

//kinetic flux density for each r value: ////{{re},type={0,1,2=1-0},spec={i,e,t}}
typedef double kin_flux[3][3][qp->dimx];

//address of kinetic flux dens in qloc array treated as 3D array:
kin_flux &KF = *((kin_flux *)(qp->qloc[qp->KIN_FLUX]));

for (int spec=0; spec<3; spec++) //over species
{
    for (int type=0; type<3; type++) //over current types
    {
        //file name:
        sprintf (filename, "%szone_%d_%s_%d_%c.dat", qp->path2linear, qp->zone->index, qp->name[qp->KIN_FLUX], type, sort[spec]);

        save_real_array (qp->dimx, qp->x, KF[spec][type], filename);
    }
}

delete [] filename;

#if DEBUG_FLAG
fprintf(stdout, "\n%s is saved.", qp->name[qp->KIN_FLUX]);
#endif
}

/*******************************************************************/

void calc_poynting_flux (flre_quants *qp)
{
//Poynting flux along r-direction out of the cylinder surface: depends on a radius

//it is fun to consider 1D array as ND one: the same trick is used everythere

//EM field profiles: {0...dimx}, {0...num_vars}, {re,im}
typedef double em_field_data[qp->dimx][qp->zone->Ncomps][2];
em_field_data &F = *((em_field_data *)qp->zone->EB_mov);

//poynting flux for each r value of the cylinder: //{{re}}
typedef double poy_flux[qp->dimx];

//address of poynting flux in qloc array treated as 1D array:
poy_flux &PF = *((poy_flux *)(qp->qloc[qp->POY_FLUX]));

complex<double> Es, Ep, Bs, Bp; //electric and magnetic fields

int const * iErsp_sys = qp->zone->me->iErsp_sys;
int const * iBrsp_sys = qp->zone->me->iBrsp_sys;

Es = F[qp->node][iErsp_sys[1]][0] + F[qp->node][iErsp_sys[1]][1]*I;
Ep = F[qp->node][iErsp_sys[2]][0] + F[qp->node][iErsp_sys[2]][1]*I;
Bs = F[qp->node][iBrsp_sys[1]][0] + F[qp->node][iBrsp_sys[1]][1]*I;
Bp = F[qp->node][iBrsp_sys[2]][0] + F[qp->node][iBrsp_sys[2]][1]*I;

PF[qp->node] = (qp->vol_fac)*(qp->r)*c/(4.0*pi)*0.5*real(conj(Es)*Bp - conj(Ep)*Bs);

qp->flag[qp->POY_FLUX] = 1; //PFD is computed
}

/*******************************************************************/

void save_poynting_flux (const flre_quants *qp)
{
char *filename = new char[1024];

//file name:
sprintf (filename, "%szone_%d_%s.dat", qp->path2linear, qp->zone->index, qp->name[qp->POY_FLUX]);

//poynting flux for each r value of the cylinder: //{{re}}
typedef double poy_flux[qp->dimx];

//address of poynting flux in qloc array treated as 1D array:
poy_flux &PF = *((poy_flux *)(qp->qloc[qp->POY_FLUX]));

save_real_array (qp->dimx, qp->x, PF, filename);

delete [] filename;

#if DEBUG_FLAG
fprintf(stdout, "\n%s is saved.", qp->name[qp->POY_FLUX]);
#endif
}

/*******************************************************************/

void calculate_field_profiles_poy_test (flre_quants *qp)
{
// if (qp->flag[qp->ABS_POWER_DENS] == 0) calc_absorbed_power_density (qp);
// if (qp->flag[qp->POY_FLUX] == 0)       calc_poynting_flux (qp);

if (qp->flag[qp->ABS_POWER_DENS] == 0 || qp->flag[qp->POY_FLUX] == 0)
{
    fprintf (stderr, "\n\aerror: consistency check failed in the function '%s' at line %d of the file '%s'.",  __FUNCTION__, __LINE__, __FILE__);
    return;
}

//absorbed power density each r value: //{{re},type={0,1},spec={i,e,t}}
typedef double abs_pow_dens[3][2][qp->dimx];

//address of abs pow dens in qloc array treated as 3D array:
abs_pow_dens &APD = *((abs_pow_dens *)(qp->qloc[qp->ABS_POWER_DENS]));

//poynting flux for each r value of the cylinder: //{{re}}
typedef double poy_flux[qp->dimx];

//address of poynting flux in qloc array treated as 1D array:
poy_flux &PF = *((poy_flux *)(qp->qloc[qp->POY_FLUX]));

//splining PF:
double* SPF = new double[(qp->N+1)*(qp->dimx)*1];

uintptr_t sidPF;
spline_alloc_ (qp->N, 1, qp->dimx, qp->x, SPF, &sidPF);

int ierr;
spline_calc_ (sidPF, PF, 0, 0, NULL, &ierr);

double divPFD;

double avrg_err = 0.0, max_err = 0.0;

double *err = new double[qp->dimx];

int k, ind = 0;

for (k=0; k<qp->dimx-1; k++) //the last point is excluded
{
    qp->node = k;
    qp->r = qp->x[k];

    spline_eval_d_ (sidPF, 1, &(qp->r), 1, 1, 0, 0, &divPFD);

    divPFD /= -(qp->r)*(qp->vol_fac);

    err[k] = fabs(divPFD - (APD[2][0][k] + qp->jaE[k]))/fmax(1.0, fabs(divPFD));

    if(err[k] > max_err)
    {
        max_err = err[k];
        ind = k;
    }

    avrg_err += err[k];
}

//average error:
avrg_err /= (qp->dimx-1);

#if DEBUG_FLAG
fprintf(stdout, "\naverage error of the solution: %le.\n", avrg_err);
#endif

if (qp->zone->sd->os->flag_additional > 1)
{
    //save:
    char *full_name = new char[1024];
    sprintf (full_name, "%szone_%d_poy_test_err.dat", qp->path2linear, qp->zone->index);
    save_real_array (qp->dimx-1, qp->x, err, full_name);
    delete [] full_name;
}

spline_free_ (sidPF);

delete [] SPF;

delete [] err;
}

/*******************************************************************/

void calc_total_flux (flre_quants *qp)
{
//total flux out of the cylinder: must be equal to
//(-1)*(dissipated power density integrated over cylinder of a radius r)

if (!(qp->flag[qp->KIN_FLUX] && qp->flag[qp->POY_FLUX]))
{
    fprintf (stderr, "\n\aerror: consistency check failed in the function '%s' at line %d of the file '%s'.",  __FUNCTION__, __LINE__, __FILE__);
    return;
}

//kinetic flux density for each r value: ////{{re},type={0,1,2=1-0},spec={i,e,t}}
typedef double kin_flux[3][3][qp->dimx];

//address of kinetic flux dens in qloc array treated as 3D array:
kin_flux &KF = *((kin_flux *)(qp->qloc[qp->KIN_FLUX]));

//poynting flux for each r value of the cylinder: //{{re}}
typedef double poy_flux[qp->dimx];

//address of poynting flux in qloc array treated as 1D array:
poy_flux &PF = *((poy_flux *)(qp->qloc[qp->POY_FLUX]));

//total flux for each r value of the cylinder: //{{re}}
typedef double tot_flux[qp->dimx];

//address of total flux in qloc array treated as 1D array:
tot_flux &TF = *((tot_flux *)(qp->qloc[qp->TOT_FLUX]));

TF[qp->node] = PF[qp->node] + KF[2][0][qp->node];

qp->flag[qp->TOT_FLUX] = 1; //TF is computed
}

/*******************************************************************/

void save_total_flux (const flre_quants *qp)
{
char *filename = new char[1024];

//file name:
sprintf (filename, "%szone_%d_%s.dat", qp->path2linear, qp->zone->index, qp->name[qp->TOT_FLUX]);

//total flux for each r value of the cylinder: //{{re}}
typedef double tot_flux[qp->dimx];

//address of total flux in qloc array treated as 1D array:
tot_flux &TF = *((tot_flux *)(qp->qloc[qp->TOT_FLUX]));

save_real_array (qp->dimx, qp->x, TF, filename);

delete [] filename;

#if DEBUG_FLAG
fprintf(stdout, "\n%s is saved.", qp->name[qp->TOT_FLUX]);
#endif
}

/*******************************************************************/

void calc_splines_for_current_density (flre_quants *qp)
{
//current density must be already computed:
if (!(qp->zone->sd->os->flag_quants[qp->CURRENT_DENS]))
{
    fprintf (stderr, "\n\aerror: consistency check failed in the function '%s' at line %d of the file '%s'.",  __FUNCTION__, __LINE__, __FILE__);
    return;
}

typedef double curr_dens[3][2][3][2][qp->dimx];

curr_dens &CD = *((curr_dens *)(qp->qloc[qp->CURRENT_DENS]));

int type = 0; //basic current density
int comp = 0; //r - component

for (int spec=0; spec<3; spec++) //over species
{
    for (int part=0; part<2; part++) //(re, im)
    {
        int ind = (qp->dimx)*(part + 2*(spec));
        for (int k=0; k<qp->dimx; k++) //over r grid
        {
            //[spec][re,im][k]:
            qp->Y[ind+k] = CD[spec][type][comp][part][k];
        }
    }
}

int ierr;
spline_calc_ (qp->sidY, qp->Y, 0, qp->nY-1, NULL, &ierr);

qp->flagS = 1;
}

/*******************************************************************/

void calc_number_density (flre_quants *qp)
{
//computes complex amplitudes of the density perturbations

if (!(qp->flagS))
{
    fprintf (stderr, "\n\aerror: consistency check failed in the function '%s' at line %d of the file '%s'.",  __FUNCTION__, __LINE__, __FILE__);
    return;
}

//current density for each r value: {{{re,im},comp={r,s,p},type={0,1},spec={i,e,t}}
typedef double curr_dens[3][2][3][2][qp->dimx];

curr_dens &CD = *((curr_dens *)(qp->qloc[qp->CURRENT_DENS]));

//number density each r value: //{{re,im},spec={i,e,t}}
typedef double number_dens[3][2][qp->dimx];

//address of number density in qloc array treated as 3D array:
number_dens &ND = *((number_dens *)(qp->qloc[qp->NUMBER_DENS]));

//eval wave numbers for given r point:
double kvals[3];

char * flag_back = qp->zone->sd->bs->flag_back;

eval_and_set_background_parameters_spec_independent_ (&(qp->r), flag_back);
eval_and_set_wave_parameters_ (&(qp->r), flag_back);
get_wave_parameters_ (kvals);

complex<double> j[3], dj[1], nd;

double djr[2];

for (int spec=0; spec<2; spec++) //over species
{
    int type = 0;
    for (int i=0; i<3; i++) //over current density components
    {
        j[i] = CD[spec][type][i][0][qp->node] + CD[spec][type][i][1][qp->node]*I;
    }

    int Imin = 2*spec; //index of real part
    int Imax = Imin+1; //index of imag part

    spline_eval_d_ (qp->sidY, 1, &(qp->r), 1, 1, Imin, Imax, djr); //[spec][re,im]

    dj[0] = djr[0] + djr[1]*I;

    nd = -I/(qp->zone->wd->omov)/(qp->zone->sd->bs->charge[spec])*
         (j[0]/(qp->r) + dj[0] + I*(kvals[1]*j[1] + kvals[2]*j[2]));

    ND[spec][0][qp->node] = real(nd);
    ND[spec][1][qp->node] = imag(nd);
}

//'total' density: enters to the total (j_t = j_i + j_e) current transformation:
//j_t' = j_t - V*(ei*n_i + ee*n_e).
for (int part=0; part<2; part++) //over {re, im}
{
    ND[2][part][qp->node] = ND[0][part][qp->node]*(qp->zone->sd->bs->charge[0]/e) + //i
                            ND[1][part][qp->node]*(qp->zone->sd->bs->charge[1]/e);  //e
}

qp->flag[qp->NUMBER_DENS] = 1;
}

/*******************************************************************/

void save_number_density (const flre_quants *qp)
{
char *filename = new char[1024];

FILE *outfile;

char sort[3] = {'i','e','t'};

//number density each r value: //{{re,im},spec={i,e,t}}
typedef double number_dens[3][2][qp->dimx];

//address of number density in qloc array treated as 3D array:
number_dens &ND = *((number_dens *)(qp->qloc[qp->NUMBER_DENS]));

for (int spec=0; spec<3; spec++) //over species
{
    //file name:
    sprintf (filename, "%szone_%d_%s_%c.dat", qp->path2linear, qp->zone->index, qp->name[qp->NUMBER_DENS], sort[spec]);

    if (!(outfile = fopen (filename, "w")))
    {
        fprintf (stderr, "\nFailed to open file %s\a\n", filename);
    }

    for (int k=0; k<qp->dimx; k++) //over r grid
    {
        fprintf (outfile, "%.16le\t%.16le\t%.16le\n", qp->x[k], ND[spec][0][k], ND[spec][1][k]);
    }
    fclose (outfile);
}

delete [] filename;

#if DEBUG_FLAG
fprintf(stdout, "\n%s is saved.", qp->name[qp->NUMBER_DENS]);
#endif
}

/*******************************************************************/

void calc_lorentz_torque_density (flre_quants *qp)
{
if (!(qp->flag[qp->NUMBER_DENS]))
{
    fprintf (stderr, "\n\aerror: consistency check failed in the function '%s' at line %d of the file '%s'.",  __FUNCTION__, __LINE__, __FILE__);
    return;
}

//current density for each r value: {{{re,im},comp={r,s,p},type={0,1},spec={i,e,t}}
typedef double curr_dens[3][2][3][2][qp->dimx];

curr_dens &CD = *((curr_dens *)(qp->qloc[qp->CURRENT_DENS]));

//EM field profiles: {0...dimx}, {0...num_vars}, {re,im}
typedef double em_field_data[qp->dimx][qp->zone->Ncomps][2];
em_field_data &F = *((em_field_data *)qp->zone->EB_mov);

//number density each r value: //{{re,im},spec={i,e,t}}
typedef double number_dens[3][2][qp->dimx];

//address of number density in qloc array treated as 3D array:
number_dens &ND = *((number_dens *)(qp->qloc[qp->NUMBER_DENS]));

//Lorentz torque density: each r value: //{{{re},comp={r,th,z}},spec={i,e,t}}
typedef double lor_torq_dens[3][3][qp->dimx];

//address of Lorentz torque density in qloc array treated as 3D array:
lor_torq_dens &LTD = *((lor_torq_dens *)(qp->qloc[qp->LOR_TORQUE_DENS]));

double h[3];
eval_and_set_background_parameters_spec_independent_ (&(qp->r),qp->zone->sd->bs->flag_back);
get_magnetic_field_parameters_ (h);

complex<double> j_rsp[3], E_rsp[3], B_rsp[3];  //rsp sys
complex<double> j_cyl[3], E_cyl[3], Bc_cyl[3];  //cyl sys

complex<double> jxBc[3]; //j x B*
complex<double> nd;      //density perturbation

for (int spec=0; spec<2; spec++) //over species (i,e)
{
    int type = 0;

    for (int i=0; i<3; i++) //over components (r,s,p)
    {
        j_rsp[i] = CD[spec][type][i][0][qp->node] +
                   CD[spec][type][i][1][qp->node]*I;

        E_rsp[i] = F[qp->node][qp->zone->me->iErsp_sys[i]][0] +
                   F[qp->node][qp->zone->me->iErsp_sys[i]][1]*I;

        B_rsp[i] = F[qp->node][qp->zone->me->iBrsp_sys[i]][0] +
                   F[qp->node][qp->zone->me->iBrsp_sys[i]][1]*I;
    }

    //transformation to cyl system (r,th,z):
    //current density:
    j_cyl[0] = j_rsp[0];
    j_cyl[1] = h[2]*j_rsp[1] + h[1]*j_rsp[2];
    j_cyl[2] = h[2]*j_rsp[2] - h[1]*j_rsp[1];

    //electric field:
    E_cyl[0] = E_rsp[0];
    E_cyl[1] = h[2]*E_rsp[1] + h[1]*E_rsp[2];
    E_cyl[2] = h[2]*E_rsp[2] - h[1]*E_rsp[1];

    //conjugate magnetic field:
    Bc_cyl[0] = conj(B_rsp[0]);
    Bc_cyl[1] = conj(h[2]*B_rsp[1] + h[1]*B_rsp[2]);
    Bc_cyl[2] = conj(h[2]*B_rsp[2] - h[1]*B_rsp[1]);

    vec_product_3D (j_cyl, Bc_cyl, jxBc);

    nd = ND[spec][0][qp->node] + ND[spec][1][qp->node]*I; //density perturbation

    //force density:
    for (int i=0; i<3; i++) //over components (r,th,z)
    {
        LTD[spec][i][qp->node] = 0.5*real((qp->zone->sd->bs->charge[spec])*
                                          nd*conj(E_cyl[i]) + E/c*jxBc[i]);
    }

    //torque density:
    LTD[spec][1][qp->node] *= qp->r;                  //T_theta = F_theta*r
    LTD[spec][2][qp->node] *= qp->zone->sd->bs->rtor; //T_z = F_z*R
}

//total torque t = i+e:
for (int i=0; i<3; i++) //over components (r,th,z)
{
    LTD[2][i][qp->node] = LTD[0][i][qp->node] + LTD[1][i][qp->node];
}

qp->flag[qp->LOR_TORQUE_DENS] = 1; //LTD is computed
}

/*******************************************************************/

void calc_lorentz_torque_on_cylinder (flre_quants *qp)
{
if (!(qp->zone->sd->os->flag_quants[qp->LOR_TORQUE_DENS]))
{
    fprintf (stderr, "\n\aerror: consistency check failed in the function '%s' at line %d of the file '%s'.",  __FUNCTION__, __LINE__, __FILE__);
    return;
}

//Lorentz torque: each r value: //{{{re},comp={r,th,z}},spec={i,e,t}}
typedef double lor_torq[3][3][qp->dimx];

//address of Lorentz torque density in qloc array treated as 3D array:
lor_torq &LTD = *((lor_torq *)(qp->qloc[qp->LOR_TORQUE_DENS]));

//address of integrated Lorentz torque in qloc array treated as 3D array:
lor_torq &LTI = *((lor_torq *)(qp->qint[qp->LOR_TORQUE_DENS]));

double factor; //must be zero for r component - no total force in r direction

for (int spec=0; spec<3; spec++) //over species (i,e,t)
{
    for (int i=0; i<3; i++) //over components (r,th,z)
    {
        if (i == 0) factor = 0.0; else factor = qp->vol_fac;
        integrate_over_cylinder (qp->dimx, qp->x, LTD[spec][i], factor, LTI[spec][i]);
    }
}
}

/*******************************************************************/

void save_lorentz_torque (const flre_quants *qp)
{
char *filename = new char[1024];

char sort[3] = {'i','e','t'};
char comp[3] = {'r','t','z'};

//Lorentz torque: each r value: //{{{re},comp={r,th,z}},spec={i,e,t}}
typedef double lor_torq[3][3][qp->dimx];

//address of Lorentz torque density in qloc array treated as 3D array:
lor_torq &LTD = *((lor_torq *)(qp->qloc[qp->LOR_TORQUE_DENS]));

//address of integrated Lorentz torque in qloc array treated as 3D array:
lor_torq &LTI = *((lor_torq *)(qp->qint[qp->LOR_TORQUE_DENS]));

for (int spec=0; spec<3; spec++) //over species
{
    for (int i=0; i<3; i++) //over (r,theta,z)
    {
        //file name:
        sprintf (filename, "%szone_%d_%s_dens_%c_%c.dat", qp->path2linear, qp->zone->index, qp->name[qp->LOR_TORQUE_DENS], comp[i], sort[spec]);

        save_real_array (qp->dimx, qp->x, LTD[spec][i], filename);

        sprintf (filename, "%szone_%d_%s_int_%c_%c.dat", qp->path2linear, qp->zone->index, qp->name[qp->LOR_TORQUE_DENS], comp[i], sort[spec]);

        save_real_array (qp->dimx, qp->x, LTI[spec][i], filename);
    }
}

delete [] filename;

#if DEBUG_FLAG
fprintf(stdout, "\n%s is saved.", qp->name[qp->LOR_TORQUE_DENS]);
#endif
}

/*******************************************************************/

void vec_product_3D (complex<double> *a, complex<double> *b, complex<double> *res)
{
res[0] =  (a[1]*b[2] - a[2]*b[1]);
res[1] = -(a[0]*b[2] - a[2]*b[0]);
res[2] =  (a[0]*b[1] - a[1]*b[0]);
}

/*******************************************************************/

void integrate_over_cylinder (int dim, double *x, double *q, double vol_fac, double *qi)
{
//simple trapezoidal formula: must be improved!
qi[0] = 0.0;

for (int i=1; i<dim; i++)
{
    qi[i] = qi[i-1] + vol_fac*0.5*(x[i]-x[i-1])*(x[i-1]*q[i-1] + x[i]*q[i]);
}
}

/*******************************************************************/

void transform_quants_to_lab_cyl_frame (const flre_quants *qp)
{
//!transforms all relevant quants to the laboratory frame and cylindrical coordinates

if (qp->zone->sd->os->flag_quants[qp->CURRENT_DENS] > 0)
{
    eval_current_dens_in_lab_frame (qp, qp->cdlab);
}

//current density:
if (qp->zone->sd->os->flag_quants[qp->CURRENT_DENS] > 1)
{
    char lframe[]  = "lab";
    char cylcomp[] = "rtz";
    char rspcomp[] = "rsp";

    double *cd = new double[(qp->dimq[qp->CURRENT_DENS])*(qp->dimx)];

    eval_current_dens_in_lab_frame (qp, cd);
    save_current_density (qp, cd, (const char *)lframe, (const char *)rspcomp);

    eval_current_dens_in_cyl_sys (qp, cd);
    save_current_density (qp, cd, (const char *)lframe, (const char *)cylcomp);

    delete [] cd;
}
}

/*******************************************************************/

void eval_current_dens_in_lab_frame (const flre_quants *qp, double *cd)
{
//!transforms current density to lab frame
//current density for each r value: {{{re,im},comp={r,s,p},type={0,1},spec={i,e,t}}

typedef double curr_dens[3][2][3][2][qp->dimx];

curr_dens &CD = *((curr_dens *)(qp->qloc[qp->CURRENT_DENS]));

//number density each r value: //{{re,im},spec={i,e,t}}
typedef double number_dens[3][2][qp->dimx];

number_dens &ND = *((number_dens *)(qp->qloc[qp->NUMBER_DENS]));

curr_dens &cude = *((curr_dens *)(cd));

complex<double> J;

double htz[2];
double vel[3];

for (int k=0; k<qp->dimx; k++) //over r grid
{
    eval_hthz (qp->x[k], 0, 0, qp->zone->bp, htz);

    vel[0] =   0.0; //r component
    vel[1] = - htz[0]*(qp->zone->sd->bs->V_gal_sys); //s component
    vel[2] =   htz[1]*(qp->zone->sd->bs->V_gal_sys); //p component

    for (int type=0; type<2; type++) //over current types
    {
        for (int i=0; i<3; i++) //over current density components
        {
            for (int spec=0; spec<2; spec++) //over species (i,e)
            {
                //j = j' + e n' V
                J = (CD[spec][type][i][0][k] + I*CD[spec][type][i][1][k]) +
                    (qp->zone->sd->bs->charge[spec])*(vel[i]) *
                    (ND[spec][0][k] +I*ND[spec][1][k]);

                cude[spec][type][i][0][k] = real(J);
                cude[spec][type][i][1][k] = imag(J);
            }
            //total current density:
            cude[2][type][i][0][k] = cude[0][type][i][0][k] + cude[1][type][i][0][k];
            cude[2][type][i][1][k] = cude[0][type][i][1][k] + cude[1][type][i][1][k];
        }
    }
}
}

/*******************************************************************/

void eval_current_dens_in_cyl_sys (const flre_quants *qp, double *cd)
{
//!transforms current density to cylindrical system and overwrites it to cd

//current density for each r value: {{{re,im},comp={r,s,p},type={0,1},spec={i,e,t}}
typedef double curr_dens[3][2][3][2][qp->dimx];
curr_dens &cude = *((curr_dens *)(cd));

complex<double> Jrsp[3], Jcyl[3];

double htz[2];

for (int k=0; k<qp->dimx; k++) //over r grid
{
    eval_hthz (qp->x[k], 0, 0, qp->zone->bp, htz);

    for (int type=0; type<2; type++) //over current types
    {
        for (int spec=0; spec<3; spec++) //over species (i,e,t)
        {
            for (int i=0; i<3; i++) //over current density components
            {
                Jrsp[i] = cude[spec][type][i][0][k] + I*cude[spec][type][i][1][k];
            }

            //transform to cylyndrical system:
            Jcyl[0] =   Jrsp[0];
            Jcyl[1] =   htz[1]*Jrsp[1] + htz[0]*Jrsp[2];
            Jcyl[2] = - htz[0]*Jrsp[1] + htz[1]*Jrsp[2];

            for (int i=0; i<3; i++) //over current density components
            {
                cude[spec][type][i][0][k] = real(Jcyl[i]);
                cude[spec][type][i][1][k] = imag(Jcyl[i]);
            }
        }
    }
}
}

/*******************************************************************/

void save_current_density (const flre_quants *qp, double *cd, const char *frame, const char *comp)
{
char *filename = new char[1024];

char sort[3] = {'i','e','t'};

FILE *outfile;

//current density for each r value: {{{re,im},comp={r,s,p},type={0,1},spec={i,e,t}}
typedef double curr_dens[3][2][3][2][qp->dimx];
curr_dens &cude = *((curr_dens *)(cd));

for (int spec=0; spec<3; spec++) //over species
{
    for (int type=0; type<2; type++) //over current types
    {
        for (int i=0; i<3; i++) //over current density components
        {
            //file name:
            sprintf (filename, "%szone_%d_%s_dens_%c_%d_%c_%s.dat", qp->path2linear, qp->zone->index, qp->name[qp->CURRENT_DENS], comp[i], type, sort[spec], frame);

            if (!(outfile = fopen (filename, "w")))
            {
                fprintf (stderr, "\nFailed to open file %s\a\n", filename);
            }

            for (int k=0; k<qp->dimx; k++) //over r grid
            {
                fprintf (outfile, "%24.16le\t%24.16le\t%24.16le\n",
                         qp->x[k], cude[spec][type][i][0][k], cude[spec][type][i][1][k]);
            }
            fclose (outfile);
        }
    }
}

delete [] filename;

#if DEBUG_FLAG
fprintf(stdout, "\n%s in %s frame in %s coordinates is saved.",
        qp->name[qp->CURRENT_DENS], frame, comp);
#endif
}

/*******************************************************************/

void interp_diss_power_density (const flre_quants *qp, double x, int type, int spec, double * dpd)
{
//dissipated power density for each r value: ////{{re},type={0,1,2=1-0},spec={i,e,t}}
typedef double diss_pow_dens[3][3][qp->dimx];

//address of diss pow dens in qloc array treated as 3D array:
diss_pow_dens & DPD = *((diss_pow_dens *)(qp->qloc[qp->DISS_POWER_DENS]));

double * xg = qp->x; //radial grid

double * yg = &DPD[spec][type][0]; //1D array of DPD of a given type and spec

int deg = 5, ind = 0.5*(qp->dimx);

eval_neville_polynom (qp->dimx, xg, yg, deg, x, 0, 0, &ind, dpd);
}

/*******************************************************************/

void interp_current_density (const flre_quants *qp, double x, int type, int spec, int comp, double * J)
{
//current density for each r value: {{{re,im},comp={r,s,p},type={0,1},spec={i,e,t}}
typedef double curr_dens[3][2][3][2][qp->dimx];
curr_dens & CD = *((curr_dens *)(qp->cdlab));

double * xg = qp->x; //radial grid

double * ygr = CD[spec][type][comp][0]; //1D array of CD of a given type, spec and comp: real
double * ygi = CD[spec][type][comp][1]; //1D array of CD of a given type, spec and comp: imag

int deg = 5, ind = 0.5*(qp->dimx);

eval_neville_polynom (qp->dimx, xg, ygr, deg, x, 0, 0, &ind, &J[0]);
eval_neville_polynom (qp->dimx, xg, ygi, deg, x, 0, 0, &ind, &J[1]);
}

/*******************************************************************/
