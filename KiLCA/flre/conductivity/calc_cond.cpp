/*! \file calc_cond.cpp
    \brief The implementation of functions declared in calc_cond.h.
*/

#include <stdlib.h>
#include <stdio.h>
#include <cmath>
#include <cstring>
#include <complex>

#include <gsl/gsl_heapsort.h>

#include "calc_cond.h"
#include "cond_profs.h"
#include "shared.h"

/*******************************************************************/

void sample_cond_func (double *r, double *f, void *p)
{
cond_profiles *cp = (cond_profiles *)(p);

cp->xt[cp->ind] = *r;

char *flag_back = cp->flag_back;

eval_and_set_background_parameters_spec_independent_ (r, flag_back, 1);
eval_and_set_wave_parameters_ (r, flag_back, 1);

eval_a_matrix_ ();
calc_dem_djmi_arrays_ (r);

int type;

for (int spec=0; spec<2; spec++)
{
    eval_and_set_background_parameters_spec_dependent_ (r, &spec, flag_back, 1);
    eval_and_set_f0_parameters_nu_and_derivs_ (r, &spec, flag_back, 1);
    eval_electric_drift_velocities_ ();
    eval_fgi_arrays_ ();
    calc_w2_array_ (&spec);
    calc_d_array_ ();

    //print_conduct_params_ ();
    //print_conduct_arrays_ ();

    type = 0;
    calc_k_matrices_ (cp->yt + cp->iKa(spec, type, 0, 0, 0, 0, 0, cp->ind));

    type = 1;
    calc_k1_matrices_ (cp->yt + cp->iKa(spec, type, 0, 0, 0, 0, 0, cp->ind));
}

//target function: imag(k0033): weigthed sum of imag parts.
*f  = 2.0e3*(cp->yt[cp->iKa(0, 0, 0, 0, 2, 2, 1, cp->ind)]); //ions
*f += 1.0e0*(cp->yt[cp->iKa(1, 0, 0, 0, 2, 2, 1, cp->ind)]); //electrons:

(cp->ind)++;
}

/*******************************************************************/

void sample_cond_func_polynom (double *r, double *f, void *p)
{
cond_profiles *cp = (cond_profiles *)(p);

char *flag_back = cp->flag_back;

eval_and_set_background_parameters_spec_independent_ (r, flag_back, 1);
eval_and_set_wave_parameters_ (r, flag_back, 1);

eval_a_matrix_ ();
calc_dem_djmi_arrays_ (r);

int type;

for (int spec=0; spec<2; spec++)
{
    eval_and_set_background_parameters_spec_dependent_ (r, &spec, flag_back, 1);
    eval_and_set_f0_parameters_nu_and_derivs_ (r, &spec, flag_back, 1);
    eval_electric_drift_velocities_ ();
    eval_fgi_arrays_ ();
    calc_w2_array_ (&spec);
    calc_d_array_ ();

    type = 0;
    calc_k_matrices_ (f + cp->iKa(spec, type, 0, 0, 0, 0, 0, 0));

    type = 1;
    calc_k1_matrices_ (f + cp->iKa(spec, type, 0, 0, 0, 0, 0, 0));
}
}

/*******************************************************************/

void calc_splines_for_K (cond_profiles *cp)
{
/*K profiles memory aligment:
{spec=0:1, {[{p=0:flreo, {q=0:flreo, {i=0:2, {j=0:2], {(re, im), {i=0:dimx-1}}}}}}}}
*/
/*sorting grid points*/
size_t *perm = new size_t[cp->dimx];

gsl_heapsort_index (perm, cp->xt, cp->dimx, sizeof(double), compare_doubles);

//rearanging K array for splining:
for (int node=0; node<cp->dimx; node++)
{
    cp->x[node] = cp->xt[perm[node]];
    for (int spec=0; spec<2; spec++)
    {
        for (int type=0; type<cp->dimt; type++)
        {
           for (int p=0; p<=cp->flreo; p++)
           {
                for (int q=0; q<=cp->flreo; q++)
                {
                    for (int i=0; i<3; i++)
                    {
                        for (int j=0; j<3; j++)
                        {
                            for (int part=0; part<2; part++)
                            {
                                cp->K [cp->iKs(spec, type, p, q, i, j, part, node)] =
                                cp->yt[cp->iKa(spec, type, p, q, i, j, part, perm[node])];
                            }
                        }
                    }
                }
            }
        }
    }
}

delete [] perm;

spline_alloc_ (cp->NK, 1, cp->dimx, cp->x, cp->CK, &(cp->sidK));

int ierr;
spline_calc_ (cp->sidK, cp->K, 0, cp->dimK-1, NULL, &ierr);

#if DEBUG_FLAG
fprintf(stdout, "\nK profiles are splined...\n");
#endif
}

/*******************************************************************/

void calc_splines_for_K_polynom (cond_profiles *cp)
{
/*K profiles memory aligment:
{spec=0:1, {[{p=0:flreo, {q=0:flreo, {i=0:2, {j=0:2], {(re, im), {i=0:dimx-1}}}}}}}}
*/

//rearanging K array for splining:
for (int node=0; node<cp->dimx; node++)
{
    cp->x[node] = cp->xt[node];
    for (int spec=0; spec<2; spec++)
    {
        for (int type=0; type<cp->dimt; type++)
        {
            for (int p=0; p<=cp->flreo; p++)
           {
                for (int q=0; q<=cp->flreo; q++)
                {
                    for (int i=0; i<3; i++)
                    {
                        for (int j=0; j<3; j++)
                        {
                            for (int part=0; part<2; part++)
                            {
                                cp->K [cp->iKs(spec, type, p, q, i, j, part, node)] =
                                cp->yt[cp->iKa(spec, type, p, q, i, j, part, node)];
                            }
                        }
                    }
                }
            }
        }
    }
}

spline_alloc_ (cp->NK, 1, cp->dimx, cp->x, cp->CK, &(cp->sidK));

int ierr;
spline_calc_ (cp->sidK, cp->K, 0, cp->dimK-1, NULL, &ierr);

#if DEBUG_FLAG
fprintf(stdout, "\nK profiles are splined...\n");
#endif
}

/*******************************************************************/

void calc_splines_for_C (cond_profiles *cp)
{
for (int spec=0; spec<2; spec++) //over species
{
    for (int type=0; type<cp->dimt; type++) //over types
    {
        for (int node=0; node<cp->dimx; node++) //over r grid
        {
            calc_C_matrices (cp, spec, type, cp->x[node], cp->RC);

            for (int s=0; s<=2*(cp->flreo); s++)
            {
                for (int i=0; i<3; i++)
                {
                    for (int j=0; j<3; j++)
                    {
                        for (int part=0; part<2; part++)
                        {

                            cp->C[cp->iCs(spec, type, s, i, j, part, node)] =
                            cp->RC[cp->iCa(s, i, j, part)];
                        }
                    }
                }
            }
        }
    }
}

spline_alloc_ (cp->NC, 1, cp->dimx, cp->x, cp->CC, &(cp->sidC));

int ierr;
spline_calc_ (cp->sidC, cp->C, 0, cp->dimC-1, NULL, &ierr);

#if DEBUG_FLAG
fprintf(stdout, "\nC profiles are splined...\n");
#endif
}

/*******************************************************************/

void calc_C_matrices (const cond_profiles *cp, int spec, int type, double r, double *C)
{
int flreo = cp->flreo;
int dimc = 9*(2*flreo+1);

//K-matrices and derivs up to flreo: huge_factor is included!
eval_K_matrices (cp, spec, type, 0, flreo, r, cp->RK);

complex<double> *Cm = new complex<double>[dimc];

for (int k=0; k<dimc; k++) Cm[k] = O;

int dimk = 2*9*(flreo+1)*(flreo+1);

for (int p=0; p<=2*flreo; p++)
{
    int nmin = fmax ((double)0, (double)(p-flreo));
    int nmax = fmin ((double)p, (double)flreo);

    for (int n=nmin; n<=nmax; n++)
    {
        for (int m=0; m<=flreo-(p-n); m++)
        {
            double coeff = pow(-1.0, m+p-n)*(cp->bico[m+p-n+(p-n)*(flreo+1)]);
            for (int i=0; i<3; i++)
            {
                for (int j=0; j<3; j++)
                {
                    //K[p,q,i,j]: (re, im)
                    int ind = 2*(j+3*(i+3*(n+(flreo+1)*(m+p-n))));

                    Cm[j+3*(i+3*p)] += coeff*((cp->RK[ind+0+m*dimk])+
                                              (cp->RK[ind+1+m*dimk])*I);
                }
            }
        }
    }
}

double scale_fac;

if (cp->flag_back[0] != 'f') scale_fac = cp->sd->bs->huge_factor;
else                         scale_fac = 1.0e0;

complex<double> cft = 2.0*pi*I*(cp->sd->bs->charge[spec])*(cp->sd->bs->charge[spec])/
                      (cp->wd->omov)/(r*scale_fac);

for (int k=0; k<dimc; k++)
{
    C[2*k+0] = real(cft*Cm[k]);
    C[2*k+1] = imag(cft*Cm[k]);
}

//transformation of C matrices: not needed if K is transformed already.
//transform_c_matrices_to_rsp_ (&r, C);

//adds a current for better invariance of the total curent: already in rsp:
if (cp->gal_corr==1 && cp->flag_back[0]=='f' && type==0 && flreo==1)
{
    calc_and_add_galilelian_correction_ (&r, &spec, cp->flag_back, C, 1);
}

//if use conduct factor: insert here!

delete [] Cm;
}

/*******************************************************************/

void save_K_matrices (const cond_profiles *cp, int spec, int type)
{
//spec=-1 means total = i+e
char reim[2][3] = {{"re"}, {"im"}};
char *full_name = new char[1024];

FILE *outfile;

int flreo = cp->flreo;

double k_m;

for (int p=0; p<=flreo; p++)
{
    for (int q=0; q<=flreo; q++)
    {
        for (int part=0; part<2; part++)
        {
            sprintf (full_name, "%s%s%d%d%s%s%s", cp->path2linear, "debug-data/kt_", p, q, "_", reim[part], ".dat");

            if (!(outfile = fopen (full_name, "w")))
            {
                fprintf (stderr, "\nFailed to open file %s\a\n", full_name);
            }

            for (int k=0; k<cp->dimx; k++)
            {
                //fprintf (outfile, "%.16le", cp->x[k]);
                for (int i=0; i<3; i++)
                {
                    for (int j=0; j<3; j++)
                    {
                        if (spec == -1)
                        {
                            k_m = cp->K[cp->iKs(0, type, p, q, i, j, part, k)] +
                                  cp->K[cp->iKs(1, type, p, q, i, j, part, k)];
                        }
                        else
                        {
                            k_m = cp->K[cp->iKs(spec, type, p, q, i, j, part, k)];
                        }
                        fprintf (outfile, "%.16le\t", k_m);
                    }
                }
                fprintf (outfile, "%.16le", cp->x[k]);
                fprintf (outfile, "\n");
            }
            fclose (outfile);
        }
    }
}
delete [] full_name;
}

/*******************************************************************/

void save_C_matrices (const cond_profiles *cp, int spec, int type)
{
//spec=-1 means total = i+e
char reim[2][3] = {{"re"}, {"im"}};
char *full_name = new char[1024];

FILE *outfile;

//int flreo = cp->sd->cs->flre_order;
//int dimc = 9*(2*flreo+1);
//int D = 2*(cp->dimx)*dimc*2; //type=0:1

double c_m;

for (int p=0; p<=2*(cp->flreo); p++)
{
    for (int part=0; part<2; part++)
    {
        sprintf (full_name, "%s%s%d%s%s%s", cp->path2linear, "debug-data/ct_", p, "_", reim[part], ".dat");

        if (!(outfile = fopen (full_name, "w")))
        {
            fprintf (stderr, "\nFailed to open file %s\a\n", full_name);
        }

        for (int k=0; k<cp->dimx; k++)
        {
            //fprintf (outfile, "%.16le", cp->x[k]);
            for (int i=0; i<3; i++)
            {
                for (int j=0; j<3; j++)
                {
                    //int ind_c = j + 3*i + 9*p + 9*(2*flreo+1)*type;
                    //int ind_a = k + part*(cp->dimx) + 2*(cp->dimx)*ind_c;

                    if (spec == -1)
                    {
                        //c_m = cp->C[ind_a+0*D] + cp->C[ind_a+1*D];
                        c_m = cp->C[cp->iCs(0, type, p, i, j, part, k)] +
                              cp->C[cp->iCs(1, type, p, i, j, part, k)];
                    }
                    else
                    {
                        //c_m = cp->C[ind_a + spec*D];
                        c_m = cp->C[cp->iCs(spec, type, p, i, j, part, k)];
                    }
                    //fprintf (outfile, "\t%.16le", k_tot);
                    fprintf (outfile, "%.16le\t", c_m);
                }
            }
            fprintf (outfile, "%.16le", cp->x[k]);
            fprintf (outfile, "\n");
        }
        fclose (outfile);
    }
}
delete [] full_name;
}

/*******************************************************************/

void save_C_matrices_fine (const cond_profiles *cp, int spec, int type, int dimf)
{
//spec=-1 means total = i+e
char *full_name = new char[1024];

FILE *outfile;

int dimC = 2*9*(2*(cp->flreo)+1);

double *Ci = new double[dimC];
double *Ce = new double[dimC];
double *C  = new double[dimC];

sprintf (full_name, "%s%s", cp->path2linear, "debug-data/ct.dat");

if (!(outfile = fopen (full_name, "w")))
{
    fprintf (stderr, "\nFailed to open file %s\a\n", full_name);
}

for (int k=0; k<cp->dimx-1; k++)
{
    for (int s=0; s<dimf; s++)
    {
        double r = cp->x[k] + s*(cp->x[k+1] - cp->x[k])/dimf;

        if (spec == -1)
        {
            eval_C_matrices (cp, 0, type, 0, 0, r, Ci);
            eval_C_matrices (cp, 1, type, 0, 0, r, Ce);

            for (int l=0; l<dimC; l++) C[l] = Ci[l] + Ce[l];
        }
        else
        {
            eval_C_matrices (cp, spec, type, 0, 0, r, C);
        }

        fprintf (outfile, "%24.16le", r);
        for (int l=0; l<dimC; l++) fprintf (outfile, "\t%24.16le", C[l]);
        fprintf (outfile, "\n");
    }
}

fclose (outfile);

delete [] Ci;
delete [] Ce;
delete [] C;

delete [] full_name;
}

/*******************************************************************/

void save_K_matrices_fine (const cond_profiles *cp, int spec, int type, int dimf)
{
//spec=-1 means total = i+e
char *full_name = new char[1024];

FILE *outfile;

int dimK = (cp->flreo+1)*(cp->flreo+1)*3*3*2;

double *Ki = new double[dimK];
double *Ke = new double[dimK];
double *K  = new double[dimK];

sprintf (full_name, "%s%s", cp->path2linear, "debug-data/kt.dat");

if (!(outfile = fopen (full_name, "w")))
{
    fprintf (stderr, "\nFailed to open file %s\a\n", full_name);
}

for (int k=0; k<cp->dimx-1; k++)
{
    for (int s=0; s<dimf; s++)
    {
        double r = cp->x[k] + s*(cp->x[k+1] - cp->x[k])/dimf;

        if (spec == -1)
        {
            eval_K_matrices (cp, 0, type, 0, 0, r, Ki);
            eval_K_matrices (cp, 1, type, 0, 0, r, Ke);

            for (int l=0; l<dimK; l++) K[l] = Ki[l] + Ke[l];
        }
        else
        {
            eval_K_matrices (cp, spec, type, 0, 0, r, K);
        }

        fprintf (outfile, "%24.16le", r);
        for (int l=0; l<dimK; l++) fprintf (outfile, "\t%24.16le", K[l]);
        fprintf (outfile, "\n");
    }
}

fclose (outfile);

delete [] Ki;
delete [] Ke;
delete [] K;

delete [] full_name;
}

/*******************************************************************/

void set_arrays_for_K (cond_profiles *cp)
{
/*K profiles memory aligment:
{spec=0:1, {[{p=0:flreo, {q=0:flreo, {i=0:2, {j=0:2], {(re, im), {i=0:dimx-1}}}}}}}}
*/
/*sorting grid points*/
size_t *perm = new size_t[cp->dimx];

gsl_heapsort_index (perm, cp->xt, cp->dimx, sizeof(double), compare_doubles);

//rearanging K array for splining:
for (int node=0; node<cp->dimx; node++)
{
    cp->x[node] = cp->xt[perm[node]];
    for (int spec=0; spec<2; spec++)
    {
        for (int type=0; type<cp->dimt; type++)
        {
            for (int p=0; p<=cp->flreo; p++)
           {
                for (int q=0; q<=cp->flreo; q++)
                {
                    for (int i=0; i<3; i++)
                    {
                        for (int j=0; j<3; j++)
                        {
                            for (int part=0; part<2; part++)
                            {
                                cp->K [cp->iKs(spec, type, p, q, i, j, part, node)] =
                                cp->yt[cp->iKa(spec, type, p, q, i, j, part, perm[node])];
                            }
                        }
                    }
                }
            }
        }
    }
}

delete [] perm;
}

/*******************************************************************/

void set_arrays_for_K_polynom (cond_profiles *cp)
{
/*K profiles memory aligment:
{spec=0:1, {[{p=0:flreo, {q=0:flreo, {i=0:2, {j=0:2], {(re, im), {i=0:dimx-1}}}}}}}}
*/

//rearanging K array for splining:
for (int node=0; node<cp->dimx; node++)
{
    cp->x[node] = cp->xt[node];
    for (int spec=0; spec<2; spec++)
    {
        for (int type=0; type<cp->dimt; type++)
        {
            for (int p=0; p<=cp->flreo; p++)
           {
                for (int q=0; q<=cp->flreo; q++)
                {
                    for (int i=0; i<3; i++)
                    {
                        for (int j=0; j<3; j++)
                        {
                            for (int part=0; part<2; part++)
                            {
                                cp->K [cp->iKs(spec, type, p, q, i, j, part, node)] =
                                cp->yt[cp->iKa(spec, type, p, q, i, j, part, node)];
                            }
                        }
                    }
                }
            }
        }
    }
}
}

/*******************************************************************/

void smooth_arrays_for_K (cond_profiles *cp)
{
/*K profiles memory aligment:
{spec=0:1, {[{p=0:flreo, {q=0:flreo, {i=0:2, {j=0:2], {(re, im), {i=0:dimx-1}}}}}}}}
*/
/*sorting grid points*/
size_t *perm = new size_t[cp->dimx];

gsl_heapsort_index (perm, cp->xt, cp->dimx, sizeof(double), compare_doubles);

for (int node=0; node<cp->dimx; node++)
{
    cp->x[node] = cp->xt[perm[node]];
}

delete [] perm;

//smooth the grid:
int ngauss = 50;  smooth_array_gauss_(&(cp->dimx), &ngauss, cp->x);

//resampling K array:
char * flag_back = cp->flag_back;

for (int node=0; node<cp->dimx; node++)
{
    double * r = &(cp->x[node]);

    eval_and_set_background_parameters_spec_independent_ (r, flag_back, 1);
    eval_and_set_wave_parameters_ (r, flag_back, 1);

    eval_a_matrix_ ();
    calc_dem_djmi_arrays_ (r);

    for (int spec=0; spec<2; spec++)
    {
        eval_and_set_background_parameters_spec_dependent_ (r, &spec, flag_back, 1);
        eval_and_set_f0_parameters_nu_and_derivs_ (r, &spec, flag_back, 1);
        eval_electric_drift_velocities_ ();
        eval_fgi_arrays_ ();
        calc_w2_array_ (&spec);
        calc_d_array_ ();

        //print_conduct_params_ ();
        //print_conduct_arrays_ ();

        int type = 0;
        calc_k_matrices_ (cp->yt + cp->iKa(spec, type, 0, 0, 0, 0, 0, 0));
    }

    for (int spec=0; spec<2; spec++)
    {
        for (int type=0; type<1; type++)
        {
           for (int p=0; p<=cp->flreo; p++)
           {
                for (int q=0; q<=cp->flreo; q++)
                {
                    for (int i=0; i<3; i++)
                    {
                        for (int j=0; j<3; j++)
                        {
                            for (int part=0; part<2; part++)
                            {
                                cp->K [cp->iKs(spec, type, p, q, i, j, part, node)] =
                                cp->yt[cp->iKa(spec, type, p, q, i, j, part, 0)];
                            }
                        }
                    }
                }
            }
        }
    }
}
}

/*******************************************************************/
