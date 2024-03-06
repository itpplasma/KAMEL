#include "eigtransform.h"

/*-----------------------------------------------------------------*/

int coeff_start_vals (int Nw, int Nfs, double r, double *zstart)
{
double *eigmat = (double *) xmalloc (2*Nw*Nw*sizeof(double));

eig_mat_eval (r, eigmat);

int k,j;

for (k=0; k<Nfs; k++)
	{
		for (j=0; j<2*Nw; j++)
		{
			fprintf(stdout, "\ncoeffs: k=%d j=%d: zstart=%25.15e", k, j, zstart[2*Nw*k+j]);
		}
}

for (k=0; k<Nw; k++)
{
	for (j=0; j<Nw; j++)
	{
		fprintf(stdout, "\neigmat: k=%d j=%d: val=(%25.15e, %25.15e)", k, j, eigmat[2*Nw*k+2*j], eigmat[2*Nw*k+2*j+1]);
	}
}

int *ipiv = (int*) xmalloc (Nw*sizeof(int));
int info;

zgesv_ (&Nw, &Nfs, eigmat, &Nw, ipiv, zstart, &Nw, &info);
if (info)
{
	fprintf(stderr, "\nerror: coeff_start_vals: zgesv_ failed!: info=%d", info);
}

	for (k=0; k<Nfs; k++)
	{
		for (j=0; j<2*Nw; j++)
		{
			fprintf(stdout, "\ncoeffs: k=%d j=%d: zstart=%25.15e", k, j, zstart[2*Nw*k+j]);
		}
	}

free(eigmat);
free(ipiv);

return info;
}

/*-----------------------------------------------------------------*/

int coeffs2solution (int Nw, int Nfs, int dim, double *rgrid, double *sol)
{
/* computes a solution sol = eig*coeff on a r grid; on entrance sol = coeff */
/* sol and coeffs are Nw x Nfs complex matrices at dim points stored in 2*Nw*Nfs*dim real array */

int Neq = 2*Nw*Nfs;

double *eigmat = (double *) xmalloc (2*Nw*Nw*sizeof(double));
double *tmp = (double *) xmalloc (Neq*sizeof(double));

double alpha[2] = {1.0,0.0}, beta[2] = {0.0,0.0};
char trans = 'N';
int info = 0;

int i, k;

for(i=0; i<dim; i++)
{
eig_mat_eval (rgrid[i], eigmat);

zgemm_ (&trans, &trans, &Nw, &Nfs, &Nw, alpha, eigmat, &Nw, sol+Neq*i, &Nw, beta, tmp, &Nw);

for(k=0; k<Neq; k++) *(sol+Neq*i+k) = tmp[k];
}

free(eigmat);
free(tmp);

return info;
}

/*-----------------------------------------------------------------*/
