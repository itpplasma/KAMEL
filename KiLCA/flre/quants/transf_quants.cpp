/*! \file transf_quants.cpp
    \brief The implementation of functions declared in transf_quants.h.
*/

#include "transf_quants.h"
#include "quants_profs.h"
#include "constants.h"
#include "inout.h"
#include "eval_back.h"

/*******************************************************************/

void transform_quants_to_lab_sys (const quants_profiles *qp)
{
//transforms all relevant quants to the laboratory frame:

//fields:
double *field = new double[2*6*(qp->fp->dimx)];
eval_and_save_EB_fields_in_lab_frame (qp, field);

delete [] field;

//current density:
if (qp->sd->os->flag_quants[qp->CURRENT_DENS] == 2)
{
    double *cd = new double[3*2*3*2*(qp->dimx)];

    char lframe[] = "lab";
    char cylcomp[] = "rtz";
    char rspcomp[] = "rsp";

    eval_current_dens_in_lab_frame (qp, cd);
    save_current_density (qp, cd, (const char *)lframe, (const char *)rspcomp);

    eval_current_dens_in_cyl_sys (qp, cd);
    save_current_density (qp, cd, (const char *)lframe, (const char *)cylcomp);

    delete [] cd;
}
}

/*******************************************************************/

void eval_current_dens_in_lab_frame (const quants_profiles *qp, double *cd)
{
//transforms current density to lab frame
//current density for each r value: {{{re,im},comp={r,s,p},type={0,1},spec={i,e,t}}
typedef double curr_dens[3][2][3][2][qp->dimx];
uintptr_t ind_cd = (qp->dimx)*(qp->offset[qp->ind[qp->CURRENT_DENS]]);
curr_dens &CD = *((curr_dens *)(qp->qloc+ind_cd));

//number density each r value: //{{re,im},spec={i,e,t}}
typedef double number_dens[3][2][qp->dimx];
uintptr_t ind_nd = (qp->dimx)*(qp->offset[qp->ind[qp->NUMBER_DENS]]); //index of number density
number_dens &ND = *((number_dens *)(qp->qloc+ind_nd)); //address of number density

curr_dens &cude = *((curr_dens *)(cd));

complex<double> J;

double htz[2];
double vel[3];

for (int k=0; k<qp->dimx; k++) //over r grid
{
    eval_hthz (qp->x[k], 0, 0, qp->bp, htz);

    vel[0] =   0.0; //r component
    vel[1] = - htz[0]*(qp->bp->V_gal_sys); //s component
    vel[2] =   htz[1]*(qp->bp->V_gal_sys); //p component

    for (int type=0; type<2; type++) //over current types
    {
        for (int i=0; i<3; i++) //over current density components
        {
            for (int spec=0; spec<2; spec++) //over species (i,e)
            {
                //j = j' + e n' V
                J = (CD[spec][type][i][0][k] + I*CD[spec][type][i][1][k]) +
                    (qp->bp->charge[spec])*(vel[i]) *
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

void eval_current_dens_in_cyl_sys (const quants_profiles *qp, double *cd)
{
//transforms current density to cylindrical system and overwrite it to cd
//current density for each r value: {{{re,im},comp={r,s,p},type={0,1},spec={i,e,t}}
typedef double curr_dens[3][2][3][2][qp->dimx];
curr_dens &cude = *((curr_dens *)(cd));

complex<double> Jrsp[3], Jcyl[3];

double htz[2];

for (int k=0; k<qp->dimx; k++) //over r grid
{
    eval_hthz (qp->x[k], 0, 0, qp->bp, htz);

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

void save_current_density (const quants_profiles *qp, double *cd, const char *frame, const char *comp)
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
            sprintf (filename, "%s%s_dens_%c_%d_%c_%s.dat", qp->path2linear, qp->name[qp->CURRENT_DENS], comp[i], type, sort[spec], frame);

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

fprintf (stdout, "\n%s in %s frame in %s coordinates is saved.", qp->name[qp->CURRENT_DENS], frame, comp);
}

/*******************************************************************/

void eval_and_save_EB_fields_in_lab_frame (const quants_profiles *qp, double *field)
{
//transforms fields to lab frame
//EM field profiles: {0...dimx}, {0...num_vars}, {re,im}
typedef double em_field_data[qp->fp->dimx][qp->fp->dimEB][2];
em_field_data &F = *((em_field_data *)qp->fp->EB);

typedef double eb_field[6][2][qp->fp->dimx];
eb_field &eb = *((eb_field *)field);

double *fieldcm = new double[2*6*(qp->fp->dimx)]; //cyl mov
eb_field &ebcm = *((eb_field *)fieldcm);

double *fieldcl = new double[2*6*(qp->fp->dimx)]; //cyl lab
eb_field &ebcl = *((eb_field *)fieldcl);

complex<double> Erspmov[3], Ecylmov[3], Ecyllab[3], Ersplab[3];
complex<double> Brspmov[3], Bcylmov[3], Bcyllab[3], Brsplab[3];

double htz[2];

for (int k=0; k<qp->fp->dimx; k++) //over r grid
{
    for (int i=0; i<3; i++) //over field rsp components
    {
        Erspmov[i] = F[k][qp->fp->iEB_sys[0][i]][0] + F[k][qp->fp->iEB_sys[0][i]][1]*I;
        Brspmov[i] = F[k][qp->fp->iEB_sys[1][i]][0] + F[k][qp->fp->iEB_sys[1][i]][1]*I;
    }

    eval_hthz (qp->fp->x[k], 0, 0, qp->bp, htz);

    //transform to cylyndrical system:
    Ecylmov[0] =   Erspmov[0];
    Ecylmov[1] =   htz[1]*Erspmov[1] + htz[0]*Erspmov[2];
    Ecylmov[2] = - htz[0]*Erspmov[1] + htz[1]*Erspmov[2];

    Bcylmov[0] =   Brspmov[0];
    Bcylmov[1] =   htz[1]*Brspmov[1] + htz[0]*Brspmov[2];
    Bcylmov[2] = - htz[0]*Brspmov[1] + htz[1]*Brspmov[2];

   //store:
    for (int i=0; i<3; i++)
    {
        ebcm[i][0][k]   = real (Ecylmov[i]);
        ebcm[i][1][k]   = imag (Ecylmov[i]);
        ebcm[3+i][0][k] = real (Bcylmov[i]);
        ebcm[3+i][1][k] = imag (Bcylmov[i]);
    }

    //transform to laboratory frame: Landau II, p. 91
    //E' = E + 1/c V x B,  B' = B - 1/c V x E
    //E = E' - 1/c V x B', B = B' + 1/c V x E'
    //V x A = e_r (-V A_th) + e_th (V A_r) + e_z (0);
    Ecyllab[0] = Ecylmov[0] + (qp->bp->V_gal_sys)/c*Bcylmov[1];
    Ecyllab[1] = Ecylmov[1] - (qp->bp->V_gal_sys)/c*Bcylmov[0];
    Ecyllab[2] = Ecylmov[2];

    Bcyllab[0] = Bcylmov[0] - (qp->bp->V_gal_sys)/c*Ecylmov[1];
    Bcyllab[1] = Bcylmov[1] + (qp->bp->V_gal_sys)/c*Ecylmov[0];
    Bcyllab[2] = Bcylmov[2];

    //store:
    for (int i=0; i<3; i++)
    {
        ebcl[i][0][k]   = real (Ecyllab[i]);
        ebcl[i][1][k]   = imag (Ecyllab[i]);
        ebcl[3+i][0][k] = real (Bcyllab[i]);
        ebcl[3+i][1][k] = imag (Bcyllab[i]);
    }

    //transform to rsp system:
    Ersplab[0] = Ecyllab[0];
    Ersplab[1] = htz[1]*Ecyllab[1] - htz[0]*Ecyllab[2];
    Ersplab[2] = htz[0]*Ecyllab[1] + htz[1]*Ecyllab[2];

    Brsplab[0] = Bcyllab[0];
    Brsplab[1] = htz[1]*Bcyllab[1] - htz[0]*Bcyllab[2];
    Brsplab[2] = htz[0]*Bcyllab[1] + htz[1]*Bcyllab[2];

    //store:
    for (int i=0; i<3; i++)
    {
        eb[i][0][k]   = real (Ersplab[i]);
        eb[i][1][k]   = imag (Ersplab[i]);
        eb[3+i][0][k] = real (Brsplab[i]);
        eb[3+i][1][k] = imag (Brsplab[i]);
    }
}

char mframe[] = "mov";
char lframe[] = "lab";
char cylcomp[] = "rtz";
char rspcomp[] = "rsp";

save_EB_fields (qp, fieldcm, (const char *)mframe, (const char *)cylcomp);
save_EB_fields (qp, fieldcl, (const char *)lframe, (const char *)cylcomp);
save_EB_fields (qp, field,   (const char *)lframe, (const char *)rspcomp);

delete [] fieldcm;
delete [] fieldcl;
}

/*******************************************************************/

void save_EB_fields (const quants_profiles *qp, const double *fields, const char *frame, const char *comp)
{
char *filename = new char[1024];

char field[2] = {'e','b'};
FILE *outfile;

typedef double eb_field[6][2][qp->fp->dimx];
eb_field &eb = *((eb_field *)fields);

for (int i=0; i<2; i++) //e, b
{
    for (int j=0; j<3; j++) //r,s,p
    {
        sprintf (filename, "%s%c%c_%s.dat", qp->path2linear, field[i], comp[j], frame);
        if (!(outfile = fopen (filename, "w")))
        {
            fprintf (stderr, "\nFailed to open file %s\a\n", filename);
        }

        for (int k=0; k<qp->fp->dimx; k++)
        {
            fprintf (outfile, "%24.16le", qp->fp->x[k]);
            fprintf (outfile, "\t%24.16le", eb[j+3*i][0][k]);
            fprintf (outfile, "\t%24.16le", eb[j+3*i][1][k]);
            fprintf (outfile, "\n");
        }
        fclose (outfile);
    }
}
fprintf (stdout, "\nE & B fields in %s frame in %s coordinates are saved.", frame, comp);
delete [] filename;
}

/*******************************************************************/

void transform_EB_fields_to_mov_frame (const background *bp, int dim, double *x, double *field)
{
//transforms fields to mov frame
//EM field profiles: {0...dimx}, {0...num_vars}, {re,im}
typedef double eb_field[6][2][dim];
eb_field &eb = *((eb_field *)field);

complex<double> Erspmov[3], Ecylmov[3], Ecyllab[3], Ersplab[3];
complex<double> Brspmov[3], Bcylmov[3], Bcyllab[3], Brsplab[3];

double htz[2];

for (int k=0; k<dim; k++) //over r grid
{
    for (int i=0; i<3; i++) //over field rsp components
    {
        Ersplab[i] = eb[i][0][k] + eb[i][1][k]*I;
        Brsplab[i] = eb[i+3][0][k] + eb[i+3][1][k]*I;
    }

    eval_hthz (x[k], 0, 0, bp, htz);

    //transform to cylyndrical system:
    Ecyllab[0] =   Ersplab[0];
    Ecyllab[1] =   htz[1]*Ersplab[1] + htz[0]*Ersplab[2];
    Ecyllab[2] = - htz[0]*Ersplab[1] + htz[1]*Ersplab[2];

    Bcyllab[0] =   Brsplab[0];
    Bcyllab[1] =   htz[1]*Brsplab[1] + htz[0]*Brsplab[2];
    Bcyllab[2] = - htz[0]*Brsplab[1] + htz[1]*Brsplab[2];

    //transform to moving frame: Landau II, p. 91
    //E' = E + 1/c V x B,  B' = B - 1/c V x E
    //E = E' - 1/c V x B', B = B' + 1/c V x E'
    //V x A = e_r(-V A_th) + e_th (V A_r) + e_z 0;
    Ecylmov[0] = Ecyllab[0] - (bp->V_gal_sys)/c*Bcyllab[1];
    Ecylmov[1] = Ecyllab[1] + (bp->V_gal_sys)/c*Bcyllab[0];
    Ecylmov[2] = Ecyllab[2];

    Bcylmov[0] = Bcyllab[0] + (bp->V_gal_sys)/c*Ecyllab[1];
    Bcylmov[1] = Bcyllab[1] - (bp->V_gal_sys)/c*Ecyllab[0];
    Bcylmov[2] = Bcyllab[2];

    //transform to rsp system:
    Erspmov[0] = Ecylmov[0];
    Erspmov[1] = htz[1]*Ecylmov[1] - htz[0]*Ecylmov[2];
    Erspmov[2] = htz[0]*Ecylmov[1] + htz[1]*Ecylmov[2];

    Brspmov[0] = Bcylmov[0];
    Brspmov[1] = htz[1]*Bcylmov[1] - htz[0]*Bcylmov[2];
    Brspmov[2] = htz[0]*Bcylmov[1] + htz[1]*Bcylmov[2];

    //store:
    for (int i=0; i<3; i++)
    {
        eb[i][0][k]   = real (Erspmov[i]);
        eb[i][1][k]   = imag (Erspmov[i]);
        eb[3+i][0][k] = real (Brspmov[i]);
        eb[3+i][1][k] = imag (Brspmov[i]);
    }
}
}

/*******************************************************************/

template <typename T1, typename T2, typename T3> void
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

/*******************************************************************/

void transform_start_values_to_mov_frame (const field_profiles *fp, double r, double *Smat)
{
int *dim_Ersp_state = (int *) fp->me->dim_Ersp_state;
int *iErsp_state    = (int *) fp->me->iErsp_state;

//max orders of the derivative for Er, Es, Ep:
int mo[3] = {dim_Ersp_state[0]-1, dim_Ersp_state[1]-1, dim_Ersp_state[2]-1};

//wave data:
int m = fp->wd->m;
double kz = (fp->wd->n)/(fp->bp->rtor);
complex<double> omlab = fp->wd->omov + kz*(fp->bp->V_gal_sys);

int ho = mo[2]; //max der order (2N-1);

//binomial coefficients:
double *C = new double[(ho+1)*(ho+1)];
binomial_coefficients (ho, C);

double *htz = new double[2*(ho+1)];
double *ht = htz, *hz = htz + ho + 1;

eval_hthz (r, 0, ho, fp->bp, htz); //ht & hz derivs

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

complex<double> *Elab[3], *Emov[3]; //Er, Es, Ep

for (int k=0; k<3; k++)
{
    Elab[k] = new complex<double>[mo[k]+1];
    Emov[k] = new complex<double>[mo[k]+1];
}

complex<double> *Br   = new complex<double>[ho+1];
complex<double> *Ez   = new complex<double>[ho+1];
complex<double> *htBr = new complex<double>[ho+1];
complex<double> *hzBr = new complex<double>[ho+1];

complex<double> kpEs, ksEp, hzEp, htEs;

int Nfs = fp->sd->es->Nfs;
int Nwaves = fp->sd->es->Nwaves;

double V = fp->bp->V_gal_sys;

for (int j=0; j<Nfs; j++) //over starting vectors
{
    //getting fields from start vectors:
    for (int k=0; k<3; k++) //over (rsp)-components
    {
        //complex(dpc), dimension(Nwaves,Nfs), intent(out) :: zstart
        for (int o=0; o<=mo[k]; o++)
        {
            int ind = 2*(j*Nwaves + iErsp_state[k] + o);
            Elab[k][o] = Smat[ind + 0] + I*Smat[ind + 1];
        }
    }

    //aux quants:
    for (int o=0; o<=ho; o++)
    {
        calc_derive_of_func_product (ho, C, o, kp, Elab[1], &kpEs);
        calc_derive_of_func_product (ho, C, o, ks, Elab[2], &ksEp);

        Br[o] = (c/omlab)*(ksEp - kpEs);

        calc_derive_of_func_product (ho, C, o, ht, Br, &htBr[o]);
        calc_derive_of_func_product (ho, C, o, hz, Br, &hzBr[o]);

        calc_derive_of_func_product (ho, C, o, ht, Elab[1], &htEs);
        calc_derive_of_func_product (ho, C, o, hz, Elab[2], &hzEp);

        Ez[o] = hzEp - htEs;
    }

    //Er, Es, Ep transformation:
    for (int o=0; o<=mo[0]; o++)
    {
        Emov[0][o] = Elab[0][o] - V/(omlab)*(kz*Elab[0][o] + I*Ez[o+1]);
    }

    for (int o=0; o<=mo[1]; o++)
    {
        Emov[1][o] = Elab[1][o] + (V/c)*hzBr[o];
    }

    for (int o=0; o<=mo[2]; o++)
    {
        Emov[2][o] = Elab[2][o] + (V/c)*htBr[o];
    }

    //Store new values to starting vectors:
    for (int k=0; k<3; k++) //(r,s,p)
    {
        for (int o=0; o<=mo[k]; o++) //der orders
        {
            int ind = 2*(j*Nwaves + iErsp_state[k] + o);

            Smat[ind + 0] = real (Emov[k][o]);
            Smat[ind + 1] = imag (Emov[k][o]);

            //fprintf (stdout, "\nj=%d k=%d o=%d: ind=%d, Elab=(%lg, %lg), Emov=(%lg, %lg)", j, k, o, ind, real(Elab[k][o]), imag(Elab[k][o]), real(Emov[k][o]), imag(Emov[k][o]));
            //fflush (stdout);
        }
    }
}

for (int k=0; k<3; k++)
{
    delete [] Elab[k];
    delete [] Emov[k];
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

/*******************************************************************/
