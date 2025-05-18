/*! \file disp_profs.cpp
    \brief The implementation of disp_profiles class.
*/

#include "disp_profs.h"
#include "inout.h"

#include <cfloat>
#include <climits>
#include <cmath>
#include <cstring>
#include <stdio.h>
#include <stdlib.h>

/*******************************************************************/

disp_profiles::disp_profiles (int Nw, int dimx_p, double * x_p, const char * flag_back_p)
{
Nwaves = Nw;

dimx = dimx_p;
x = x_p;

flag_back = new char[8];
strcpy (flag_back, flag_back_p);

dimk = 2*Nwaves;
k = new double[(dimk)*(dimx)];

dimp = 2*(Nwaves)*(Nwaves);
p = new double[(dimp)*(dimx)];

S = 0;
R = 0;
sid = 0;
}

/*******************************************************************/

void disp_profiles::calculate_dispersion_profiles (void)
{
int f = 0;

for (int i=0; i<dimx; i++)
{
    calc_dispersion_ (x+i, flag_back, &f, k+dimk*i, p+dimp*i, 1);
}

#if SORT_DISPERSION_PROFILES == 1

sort_dispersion_profiles ();

#endif

#if DEBUG_FLAG
    fprintf (stdout, "\n\ndispersion profiles are computed.\n");
#endif
}

/*******************************************************************/

void disp_profiles::sort_dispersion_profiles (void)
{
//sorts kvs at first point over abs vals, then tries to keep the order of kvs
//at the end sorts over average modulus of kvs
int i, j;

int Nw = Nwaves, Nw2 = 2*Nw, NwNw = Nw*Nw, NwNw2 = 2*NwNw;

double *kv = k;
double *um = p;
double *r = x;

double *tmp = new double[Nw];
double *tmpk = new double[Nw2];
double *tmpu = new double[NwNw2];

int *ind = new int[Nw];

//sorts at first r point over abs values of kvs:
for (i=0; i<Nw; i++) tmp[i] = sqrt(kv[2*i]*kv[2*i] + kv[2*i+1]*kv[2*i+1]);

int ier, flag = 1;

dpsort_ (tmp, &Nw, ind, &flag, &ier);
if (ier != 0)
{
    fprintf (stderr, "\nsort_dispersion_profiles: sorting failed: ier = %d", ier);
}

for (j=0; j<Nw; j++) ind[j]--; /*in fortran indices start from 1!*/

//rearranging k's at first point:
for (j=0; j<Nw2; j++) tmpk[j] = kv[j];
for (j=0; j<Nw; j++) /* over eigvals */
{
    if (ind[j] == j) continue;
    kv[2*j + 0] = tmpk[2*ind[j] + 0];
    kv[2*j + 1] = tmpk[2*ind[j] + 1];
}

//rearranging the whole eigvec matrix at first point: eigvec(:,i) = tmp2(:,ind(i,k))
for (i=0; i<NwNw2; i++) tmpu[i] = um[i];
for (i=0; i<Nw; i++) //over eigvectors
{
    if (ind[i] == i) continue;
    for (j=0; j<Nw; j++) //over components
    {
        um[Nw2*i + 2*j + 0] = tmpu[Nw2*ind[i] + 2*j + 0];
        um[Nw2*i + 2*j + 1] = tmpu[Nw2*ind[i] + 2*j + 1];
    }
}

double xm, ym, dx, dy;

int I;

double dst, min1, min2;

//tries to keep the order of eigenvalues and get continious curves:
for (int l=1; l<dimx; l++) //over rgrid
{
    I = 2*Nw*l; //offset of l's values for lth point

    for (j=0; j<Nw; j++) //over sorted vals at previous r-point with index l-1
    {
        if (l == 1)
        {
            xm = kv[2*j + 0];
            ym = kv[2*j + 1];
        }
        else
        {
            //extrapolated position of jth k's:
            xm = kv[I - Nw2 + 2*j + 0] +
                (kv[I - Nw2 + 2*j + 0] - kv[I - 2*Nw2 + 2*j + 0])/
                (r[l-1] - r[l-2])*(r[l] - r[l-1]);

            ym = kv[I - Nw2 + 2*j + 1] +
                (kv[I - Nw2 + 2*j + 1] - kv[I - 2*Nw2 + 2*j + 1])/
                (r[l-1] - r[l-2])*(r[l] - r[l-1]);
        }

        i = 0;
        dx = kv[I+2*i+0] - xm;
        dy = kv[I+2*i+1] - ym ;
        dst = sqrt(dx*dx+dy*dy);
        min1 = dst;
        min2 = DBL_MAX;
        ind[j] = i;

        for (i=1; i<Nw; i++) //over eigvals at current r point to find closest to j-th
        {
            dx = kv[I+2*i+0] - xm;
            dy = kv[I+2*i+1] - ym ;

            dst = sqrt(dx*dx+dy*dy);

            if (dst < min1)
            {
                min2 = min1;
                min1 = dst;
                ind[j] = i;
            }
            if (dst < min2 && dst > min1) min2 = dst;
        }

        if (min2 < 1.0e1*min1)
        {
            fprintf (stdout, "\nsort_dispersion_profiles: warning: min1 = %lg min2 = %lg j = %d ind = %d", min1, min2, j, ind[j]);
        }
    }

    //rearranging eigenvalues and vectors to form continious curves:
    for (j=0; j<Nw2; j++) tmpk[j] = kv[I+j];

    for (j=0; j<Nw; j++) //over eigvals
    {
        if (ind[j] == j) continue;
        kv[I + 2*j + 0] = tmpk[2*ind[j] + 0];
        kv[I + 2*j + 1] = tmpk[2*ind[j] + 1];
    }

    //the whole eigvec matrix at l-th r point: eigvec(:,i) = tmp2(:,ind(i,l))
    for (i=0; i<NwNw2; i++) tmpu[i] = um[NwNw2*l + i];
    for (i=0; i<Nw; i++) //over eigvectors
    {
        if (ind[i] == i) continue;
        for (j=0; j<Nw; j++) /* over components */
        {
            um[NwNw2*l + Nw2*i + 2*j + 0] = tmpu[Nw2*ind[i] + 2*j + 0];
            um[NwNw2*l + Nw2*i + 2*j + 1] = tmpu[Nw2*ind[i] + 2*j + 1];
        }
    }
}

//sorts the eigencurves by average value of their modulus:
for (i=0; i<Nw; i++) tmp[i] = 0.0;

for (int l=0; l<dimx; l++) //over rgrid
{
    for (i=0; i<Nw; i++)
    {
        tmp[i] += sqrt(kv[Nw2*l+2*i]*kv[Nw2*l+2*i] + kv[Nw2*l+2*i+1]*kv[Nw2*l+2*i+1]);
    }
}

flag = 1;
dpsort_ (tmp, &Nw, ind, &flag, &ier);
if (ier != 0)
{
    fprintf (stderr, "\nsort_dispersion_profiles: sorting failed: ier = %d", ier);
}

for (j=0; j<Nw; j++) ind[j]--; //in fortran indices start from 1!

fprintf (stdout, "\naverage |k| values:");
for (j=0; j<Nw; j++)
{
    fprintf (stdout, "\nj = %d: |k| = %le\tind = %d", j, tmp[j]/(dimx), ind[j]);
}
fprintf (stdout, "\n");

//rearranges eigcurves: the last are the fake modes:
for (int l=0; l<dimx; l++)
{
    I = 2*Nw*l; //offset of k's values for lth point

    for (j=0; j<Nw2; j++) tmpk[j] = kv[I+j];

    for (j=0; j<Nw; j++) //over eigvals
    {
        if (ind[j] == j) continue;
        kv[I + 2*j + 0] = tmpk[2*ind[j] + 0];
        kv[I + 2*j + 1] = tmpk[2*ind[j] + 1];
    }

    //the whole eigvec matrix at l-th r point: eigvec(:,i) = tmp2(:,ind(i,l))
    for (i=0; i<NwNw2; i++) tmpu[i] = um[NwNw2*l + i];
    for (i=0; i<Nw; i++) //over eigvectors
    {
        if (ind[i] == i) continue;
        for (j=0; j<Nw; j++) /* over components */
        {
            um[NwNw2*l + Nw2*i + 2*j + 0] = tmpu[Nw2*ind[i] + 2*j + 0];
            um[NwNw2*l + Nw2*i + 2*j + 1] = tmpu[Nw2*ind[i] + 2*j + 1];
        }
    }
}

free (ind);
free (tmp);
free (tmpk);
free (tmpu);
}

/*******************************************************************/

void disp_profiles::save_dispersion_profiles (char * filename)
{
save_cmplx_matrix_to_one_file (Nwaves, 1, dimx, x, k, filename);
}

/*******************************************************************/
