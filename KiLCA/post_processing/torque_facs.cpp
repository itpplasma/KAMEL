/*! \file
    \brief The post processing program to evaluate form-factors for torques acting on a plasma.
*/

#include <stdlib.h>
#include <stdio.h>
#include <cmath>
#include <complex>
#include <cstring>
#include <climits>

#include "inout.h"
#include "spline.h"
#include "constants.h"
#include "background.h"
#include "settings.h"
#include "adaptive_grid_pol.h"

/*******************************************************************/

extern "C"
{
void save_form_facs_fortran_ (int *modes_dim, int *modes_mn, int *m_min, int *m_max,
                              int *n_min, int *n_max, int *modes_index,
                              int *dim, double *sqr_psi_grid, double *ffpsi, int *format_flag);
}

/*******************************************************************/

int main (int argc, char **argv)
{
//read file:
char file_ff[] = "torque_facs.in";
FILE *in;

if ((in=fopen (file_ff, "r"))==NULL)
{
    fprintf(stderr, "\nerror: torque_facs: failed to open file %s\a\n", file_ff);
    exit(0);
}

char *str_buf = new char[1024]; //str buffer

//Projects settings:
read_line_2skip_it (in, &str_buf);
char *fname = new char [1024];
read_line_2get_string (in, &fname);

char **path2projects = new char *[2];

path2projects[0] = new char[1024];
read_line_2get_string (in, &(path2projects[0]));

path2projects[1] = new char[1024];
read_line_2get_string (in, &(path2projects[1]));

read_line_2skip_it (in, &str_buf);

//grid settings:
int dim;
double rmin = 0.0, rmini, rmax;
read_line_2skip_it (in, &str_buf);
read_line_2get_int (in, &dim);
read_line_2get_double (in, &rmini);
read_line_2get_double (in, &rmax);
read_line_2skip_it (in, &str_buf);

double *r = new double[dim];
for (int k=0; k<dim; k++) r[k] = rmin + k*(rmax-rmin)/(dim-1);

//modes dimension:
int dma, imin, imax;
read_line_2skip_it (in, &str_buf);
read_line_2get_int (in, &dma);
read_line_2get_int (in, &imin);
read_line_2get_int (in, &imax);
read_line_2skip_it (in, &str_buf);

//flagsmodes dimension:
int ffflag;
read_line_2skip_it (in, &str_buf);
read_line_2get_int (in, &ffflag);
read_line_2skip_it (in, &str_buf);

char *quant_file = new char [1024];
read_line_2skip_it (in, &str_buf);
read_line_2get_string (in, &quant_file);
read_line_2skip_it (in, &str_buf);

int spl_deg;
read_line_2skip_it (in, &str_buf);
read_line_2get_int (in, &spl_deg);
read_line_2skip_it (in, &str_buf);

int grid_flag;
read_line_2skip_it (in, &str_buf);
read_line_2get_int (in, &grid_flag);
read_line_2skip_it (in, &str_buf);

read_line_2skip_it (in, &str_buf);
char *eq_file = new char[1024];
read_line_2get_string (in, &eq_file);
read_line_2skip_it (in, &str_buf);

int pol_deg;
read_line_2skip_it (in, &str_buf);
read_line_2get_int (in, &pol_deg);
read_line_2skip_it (in, &str_buf);

read_line_2skip_it (in, &str_buf);
char *fname_psi = new char [1024];
read_line_2get_string (in, &fname_psi);
read_line_2skip_it (in, &str_buf);

int m_min, m_max, n_min, n_max;
read_line_2skip_it (in, &str_buf);
read_line_2get_int (in, &m_min);
read_line_2get_int (in, &m_max);
read_line_2get_int (in, &n_min);
read_line_2get_int (in, &n_max);
read_line_2skip_it (in, &str_buf);

int format_flag;
read_line_2skip_it (in, &str_buf);
read_line_2get_int (in, &format_flag);
read_line_2skip_it (in, &str_buf);

fclose (in);

//reading of the equlibrium file:
int dim_eq = 0;
double *x_eq = NULL, *eq_prof = NULL, *sqr_psi = NULL;
double *sqr_psi_grid = NULL;
double *rval = NULL;

if (grid_flag > 0)
{
    dim_eq = count_lines_in_file_with_comments (eq_file, 1);

    x_eq = new double[dim_eq];
    eq_prof = new double[5*dim_eq];

    load_data_file_with_comments (eq_file, dim_eq, 5, x_eq, eq_prof);

    sqr_psi = new double[dim_eq];

    for (int k=0; k<dim_eq; k++)
    {
        sqr_psi[k] = sqrt(eq_prof[k + 1*dim_eq]);
    }

    delete [] eq_prof;

    sqr_psi_grid = new double[dim];
    rval = new double[dim];

    int ind = 0;

    for (int k=0; k<dim; k++)
    {
        sqr_psi_grid[k] = sqr_psi[0] + k*(sqr_psi[dim_eq-1]-sqr_psi[0])/(dim-1);

        //interpolation of r at sqr_psi point:
        find_index_for_interp (pol_deg, sqr_psi_grid[k], dim_eq, sqr_psi, &ind);

        eval_interp_polynom (pol_deg, sqr_psi+ind, x_eq+ind, 1, sqr_psi_grid[k], rval+k);
    }

    FILE *outfile;

    if (!(outfile = fopen ("sqrt_psi_grid.dat", "w")))
    {
        fprintf (stderr, "\nFailed to open file %s\a\n", "sqrt_psi_grid.dat");
    }

    for (int k=0; k<dim; k++)
    {
        fprintf (outfile, "%24.16le\t%24.16le\n", sqr_psi_grid[k], rval[k]);
    }

    fclose (outfile);
}

//Reading modes.in file:
char file_modes[] = "modes.in";

if ((in=fopen (file_modes, "r"))==NULL)
{
    fprintf(stderr, "\nerror: read_input_settings: failed to open file %s\a\n", file_modes);
    exit(0);
}

int *modes = new int[2*dma]; //(m,n)

int i;
for (i=0; i<dma; i++)
{
    if (feof(in) || ferror(in)) break;
    fscanf (in, "(%d, %d)\n", &(modes[2*i]), &(modes[2*i+1]));
}

if(i != dma)
{
    fprintf(stderr, "\nerror: read_input_settings: read error or false modes array dimension!\a\n");
    exit(0);
}

fclose (in);

//wave field data:
int dim_prof;
double *tmp_prof;

int dimx[2];
double *x[2], *y[2], *s[2];
uintptr_t sid_wav[2];

complex<double> FP[2];
double fp[2];

//complex<double> form_facs[imax-imin+1][dim]; //does not work for big values!!!
complex<double> *form_facs = new complex<double>[(imax-imin+1)*dim];

//complex<double> form_facs_psi[imax-imin+1][dim];
double *form_facs_psi = NULL;
double *ffgrid = NULL;

if (grid_flag > 0)
{
    form_facs_psi = new double[2*(imax-imin+1)*dim];
    ffgrid = new double[2*dim];
}

settings *sd = new settings;
read_and_set_settings (sd, path2projects[0]);

background *bp = new background;
read_background_settings (bp, path2projects[0]);

int l;

char *filename = new char[1024];
int ierr;

int Br_norm_point = 0; //1

for (i=imin; i<=imax; i++) //over modes
{
    fprintf (stdout, "\n(m=%d, n=%d) mode...", modes[2*i], modes[2*i+1]); fflush(stdout);

    int m = sd->as->modes[2*i];
    int n = sd->as->modes[2*i+1];

    //set wave frequency in a moving frame for the mode:
    complex<double> flab = sd->as->flab;
    complex<double> omov = 2.0*pi*flab - n*(bp->V_gal_sys)/(bp->rtor);

    fprintf (stdout, "\nomov = (%le, %le)", real(omov), imag(omov)); fflush(stdout);

    complex<double> coeff = ((double) n)/omov; //n/omega'

    for (int p=0; p<1; p++) //over paths
    {
        sprintf (filename, "%slinear-data/n_%d_m_%d/%s", path2projects[p], modes[2*i+1], modes[2*i], quant_file);

        dim_prof = count_lines_in_file (filename, 0);

        x[p] = new double[dim_prof];
        tmp_prof = new double[2*dim_prof]; //re and im

        load_data_file (filename, dim_prof, 1, x[p], tmp_prof);

        //take a part before antenna:
        for (l=0; l<dim_prof; l++) if (x[p][l]>=rmax) break;
        dimx[p] = l;

        y[p] = new double[dimx[p]*2];

        for (int k=0; k<dimx[p]; k++)
        {
            y[p][k] = tmp_prof[k];
            //y[p][k+dimx[p]] = tmp_prof[k+dim_prof];
        }

        s[p] = new double[dimx[p]*(spl_deg+1)*2];

        //spline:
        spline_alloc_ (spl_deg, 1, dimx[p], x[p], s[p], &sid_wav[p]);
        spline_calc_ (sid_wav[p], y[p], 0, 0, NULL, &ierr);

        delete [] tmp_prof;
    }

        int p = 1;

        sprintf (filename, "%slinear-data/n_%d_m_%d/br_lab.dat", path2projects[p], modes[2*i+1], modes[2*i]);

        dim_prof = count_lines_in_file (filename, 0);

        x[p] = new double[dim_prof];
        tmp_prof = new double[2*dim_prof]; //re and im

        load_data_file (filename, dim_prof, 2, x[p], tmp_prof);

        //take a part before antenna:
        for (l=0; l<dim_prof; l++) if (x[p][l]>=rmax) break;
        dimx[p] = l;

        y[p] = new double[dimx[p]*2];

        for (int k=0; k<dimx[p]; k++)
        {
            y[p][k] = tmp_prof[k];
            y[p][k+dimx[p]] = tmp_prof[k+dim_prof];
        }

        s[p] = new double[dimx[p]*(spl_deg+1)*2];

        //spline:
        spline_alloc_ (spl_deg, 1, dimx[p], x[p], s[p], &sid_wav[p]);
        spline_calc_ (sid_wav[p], y[p], 0, 1, NULL, &ierr);

        delete [] tmp_prof;

        if (Br_norm_point == 1)
        {
            p = 1;
            double rmid = 0.5*x_eq[dim_eq-1];
            spline_eval_ (sid_wav[p], 1, &rmid, 0, 0, 0, 1, fp);

            FP[p] = fp[0]*fp[0] + fp[1]*fp[1];

            fprintf (stdout, "\nrnorm = %le\t|Br|^2 = %le", rmid, real(FP[1]));
        }

    //interpolation:
    for (int k=0; k<dim; k++) //over grid
    {
        for (int p=0; p<1; p++) //over paths
        {
            spline_eval_ (sid_wav[p], 1, r+k, 0, 0, 0, 0, fp);
            //FP[p] = fp[0] + I*fp[1];
            FP[p] = fp[0];
        }

        if (Br_norm_point == 0)
        {
            p = 1;
            spline_eval_ (sid_wav[p], 1, r+k, 0, 0, 0, 1, fp);

            FP[p] = fp[0]*fp[0] + fp[1]*fp[1];
        }

        //form factors:
        //if (ffflag == 0)
        //{
            //form_facs[i-imin][k] = FP[0]/FP[1];
        //    form_facs[k + (i-imin)*dim] = FP[0]/FP[1];
        //}
        //else if (ffflag == 1)
        //{
            //form_facs[k + (i-imin)*dim] = FP[0];
            //form_facs[k + (i-imin)*dim] = FP[0]/FP[1];
            //form_facs[k + (i-imin)*dim] = FP[1];
            form_facs[k + (i-imin)*dim] = coeff * FP[0]/FP[1];
        //}
    }

    //removes noise at the center:
    for (l=0; l<dim; l++) if (r[l]>rmini) break;
    for (int k=0; k<l; k++)
    {
        //form_facs[i-imin][k] = form_facs[i-imin][l] +
        //(r[k]-r[l])/(r[l+1]-r[l])*(form_facs[i-imin][l+1]-form_facs[i-imin][l]);

        form_facs[k + (i-imin)*dim] = form_facs[l+0 + (i-imin)*dim] +
           (r[k]-r[l])/(r[l+1]-r[l])*(form_facs[l+1 + (i-imin)*dim] -
                                      form_facs[l+0 + (i-imin)*dim]);

    }

    //calculation ff on a sqrt(psi) grid:
    if (grid_flag > 0)
    {
        for (int k=0; k<dim; k++)
        {
            //ffgrid[2*k+0] = real(form_facs[i-imin][k]);
            //ffgrid[2*k+1] = imag(form_facs[i-imin][k]);
            ffgrid[2*k+0] = real(form_facs[k + (i-imin)*dim]);
            ffgrid[2*k+1] = imag(form_facs[k + (i-imin)*dim]);
        }

        int ind = 0;

        for (int k=0; k<dim; k++)
        {
            //interpolation of ff at rval[k] point:
            find_index_for_interp (pol_deg, rval[k], dim, r, &ind);

            double ffval[2];
            eval_interp_polynom (pol_deg, r+ind, ffgrid+2*ind, 2, rval[k], ffval);

            //form_facs_psi[i-imin][k] = ffval[0] + I*ffval[1];
            form_facs_psi[0 + 2*k + (i-imin)*dim*2] = ffval[0];
            form_facs_psi[1 + 2*k + (i-imin)*dim*2] = ffval[1];
        }
    }

    //clean up:
    for (int p=0; p<2; p++)
    {
        delete [] x[p];
        delete [] y[p];
        delete [] s[p];
        spline_free_ (sid_wav[p]);
    }

    fprintf (stdout, " Ok."); fflush (stdout);
}

//save:
fprintf (stdout, "\nSaving..."); fflush (stdout);

FILE *outfile;

if (grid_flag != 1)
{
    if (!(outfile = fopen (fname, "w")))
    {
        fprintf (stderr, "\nFailed to open file %s\a\n", fname);
    }

    fprintf (outfile, "%%modes:");
    for (i=imin; i<=imax; i++)
    {
        fprintf (outfile, "\t(%d, %d)", modes[2*i+0], modes[2*i+1]);
    }
    fprintf (outfile, "\n");

    fprintf (outfile, "%%r");
    for (i=imin; i<=imax; i++)
    {
        fprintf (outfile, "\treal(T_{mn}) imag(T_{mn})");
    }
    fprintf (outfile, "\n");

    for (int k=0; k<dim; k++)
    {
        fprintf (outfile, "%24.16le", r[k]);
        for (i=imin; i<=imax; i++)
        {
            //fprintf (outfile, "\t%24.16le\t%24.16le",
            //         real(form_facs[i-imin][k]), imag(form_facs[i-imin][k]));
            fprintf (outfile, "\t%24.16le\t%24.16le",
                     real(form_facs[k + (i-imin)*dim]), imag(form_facs[k + (i-imin)*dim]));
        }
        fprintf (outfile, "\n");
    }

fclose (outfile);
}

if (grid_flag >= 1)
{
    if (!(outfile = fopen (fname_psi, "w")))
    {
        fprintf (stderr, "\nFailed to open file %s\a\n", fname_psi);
    }

    fprintf (outfile, "%%modes:");
    for (i=imin; i<=imax; i++)
    {
        fprintf (outfile, "\t(%d, %d)", modes[2*i+0], modes[2*i+1]);
    }
    fprintf (outfile, "\n");

    fprintf (outfile, "%%r");
    for (i=imin; i<=imax; i++)
    {
        fprintf (outfile, "\treal(T_{mn}) imag(T_{mn})");
    }
    fprintf (outfile, "\n");

    for (int k=0; k<dim; k++)
    {
        fprintf (outfile, "%24.16le", sqr_psi_grid[k]);
        for (i=imin; i<=imax; i++)
        {
            //fprintf (outfile, "\t%24.16le\t%24.16le",
            //         real(form_facs_psi[i-imin][k]), imag(form_facs_psi[i-imin][k]));
            fprintf (outfile, "\t%24.16le\t%24.16le",
                              form_facs_psi[0 + 2*k + (i-imin)*dim*2],
                              form_facs_psi[1 + 2*k + (i-imin)*dim*2]);
        }
        fprintf (outfile, "\n");
    }

    fclose (outfile);

    //for field line tracing:
    int modes_dim = imax-imin+1;
    int *modes_mn = new int[2*modes_dim];

    int ind = 0;
    for (int i=imin; i<=imax; i++)
    {
        modes_mn[2*ind+0] =  modes[2*i+0];
        modes_mn[2*ind+1] =  modes[2*i+1];
        ind++;
    }

    int *modes_index = new int[(m_max-m_min+1)*(n_max-n_min+1)];

    for (int m=m_min; m<=m_max; m++)
    {
        for (int n=n_min; n<=n_max; n++)
        {
            modes_index[m-m_min + (m_max-m_min+1)*(n-n_min)] = -1;
        }
    }

    for (int i=0; i<modes_dim; i++)
    {
        int m = modes_mn[2*i+0];
        int n = modes_mn[2*i+1];

        modes_index[m-m_min + (m_max-m_min+1)*(n-n_min)] = i + 1;
    }

    save_form_facs_fortran_ (&modes_dim, modes_mn, &m_min, &m_max, &n_min, &n_max, modes_index, &dim, sqr_psi_grid, form_facs_psi, &format_flag);
}

fprintf (stdout, " Ok.\n"); fflush (stdout);

//clean up;
delete [] str_buf;
delete [] fname;
delete [] path2projects[0];
delete [] path2projects[1];
delete [] path2projects;
delete [] quant_file;
delete [] eq_file;
delete [] fname_psi;

delete [] r;
delete [] modes;
delete [] filename;
delete [] form_facs;

if (grid_flag > 0)
{
delete [] x_eq;
delete [] sqr_psi;
delete [] sqr_psi_grid;
delete [] rval;
delete [] ffgrid;
delete [] form_facs_psi;
}

return 0;
}

/*****************************************************************************/
