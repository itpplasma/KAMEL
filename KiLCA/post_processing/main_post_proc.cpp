/*! \file
    \brief The post processing program to evaluate form-factors in various variables (see input file).
*/

#include <stdlib.h>
#include <stdio.h>
#include <cmath>
#include <complex>
#include <cstring>
#include <climits>

#include "inout.h"
#include "constants.h"
#include "background.h"
#include "adaptive_grid_pol.h"

/*******************************************************************/

extern "C"
{
void save_form_facs_fortran_ (char *fname, int *modes_dim, int *modes_mn, int *m_min, int *m_max,
                              int *n_min, int *n_max, int *modes_index,
                              int *dim, double *sqr_psi_grid, double *ffpsi, int *format_flag);
}

/*******************************************************************/

void transform_cyl_radius_to_surface_label (int dim, double * r, int pol_deg,
                                            int label_flag, double * label);

void transform_surface_label_to_cyl_radius (int dim, double * r, int pol_deg,
                                            int label_flag, double * label);

inline void evaluate_Br_formfactors (int dim1, double * r_EB1, double * prof_EB1,
                                     int dim2, double * r_EB2, double * prof_EB2,
                                     int pol_deg, double r, complex<double> *form_facs);

inline void evaluate_jp_over_Br (int dim_EB, double * r_EB, double * prof_EB,
                                 int dim_jp, double * r_jp, double * prof_jp,
                                 int pol_deg, double r, complex<double> *form_facs);

/*******************************************************************/

int main (int argc, char **argv)
{
char name_out[][32] = {"Br_ff", "Jp_over_Br_over_r", "Tphi_ff"};

char name_label[][32] = {"cyl_radius", "sqrt_psi_pol", "psi_pol", "norm_psi_pol", "sqrt_psi_tor", "psi_tor", "norm_psi_tor"};

//read file:
char file_ff[] = "post_proc.in";
FILE *in;

if ((in=fopen (file_ff, "r"))==NULL)
{
    fprintf(stderr, "\nerror: post_proc: failed to open file %s\a\n", file_ff);
    exit(0);
}

char *str_buf = new char[1024]; //str buffer

//configuration settings:
read_line_2skip_it (in, &str_buf);
char *conf = new char [1024];
read_line_2get_string (in, &conf);

char **path2projects = new char *[2];

path2projects[0] = new char[1024];
read_line_2get_string (in, &(path2projects[0]));

path2projects[1] = new char[1024];
read_line_2get_string (in, &(path2projects[1]));

read_line_2skip_it (in, &str_buf);

//grid settings:
int dim;
double label_min, label_max, rmini;
read_line_2skip_it (in, &str_buf);
read_line_2get_int (in, &dim);
read_line_2get_double (in, &label_min);
read_line_2get_double (in, &label_max);
read_line_2get_double (in, &rmini);
read_line_2skip_it (in, &str_buf);

//modes settings:
complex<double> flab;
read_line_2skip_it (in, &str_buf);
read_line_2get_complex (in, &flab);

int dma, imin, imax;
read_line_2get_int (in, &dma);
read_line_2get_int (in, &imin);
read_line_2get_int (in, &imax);

int m_min, m_max, n_min, n_max;
read_line_2get_int (in, &m_min);
read_line_2get_int (in, &m_max);
read_line_2get_int (in, &n_min);
read_line_2get_int (in, &n_max);
read_line_2skip_it (in, &str_buf);

int pol_deg;
read_line_2skip_it (in, &str_buf);
read_line_2get_int (in, &pol_deg);
read_line_2skip_it (in, &str_buf);

int output_flag;
read_line_2skip_it (in, &str_buf);
read_line_2get_int (in, &output_flag);
read_line_2skip_it (in, &str_buf);

int label_flag;
read_line_2skip_it (in, &str_buf);
read_line_2get_int (in, &label_flag);
read_line_2skip_it (in, &str_buf);

int format_flag;
read_line_2skip_it (in, &str_buf);
read_line_2get_int (in, &format_flag);
read_line_2skip_it (in, &str_buf);

fclose (in);

//output label grid:
double * label = new double[dim]; for (int k=0; k<dim; k++) label[k] = label_min + k*(label_max-label_min)/(dim-1);

//output radial grid:
double * rout = new double[dim];

transform_surface_label_to_cyl_radius (dim, rout, pol_deg, label_flag, label);

//Reading modes.in file:
char file_modes[] = "modes.in";

if ((in=fopen (file_modes, "r"))==NULL)
{
    fprintf(stderr, "\nerror: failed to open file %s\a\n", file_modes);
    exit(0);
}

int *modes = new int[2*dma]; //(m,n)

int i;
for (i=0; i<dma; i++)
{
    if (feof(in) || ferror(in)) break;
    fscanf (in, "(%d, %d)\n", &(modes[2*i]), &(modes[2*i+1]));
}

if (i != dma)
{
    fprintf(stderr, "\nerror: read error or false modes array dimension!\a\n");
    exit(0);
}

fclose (in);

//EB data:
int ncols_EB = 12;
int dim_EB[2] = {NULL, NULL};
double * r_EB[2] = {NULL, NULL};
double * prof_EB[2] = {NULL, NULL};

//jp over Br:
int ncols_jp = 2;
int dim_jp[2] = {NULL, NULL};
double * r_jp[2] = {NULL, NULL};
double * prof_jp[2] = {NULL, NULL};

//torques:
//...

//form_factors:
complex<double> *form_facs = new complex<double>[(imax-imin+1)*dim];

int l;

char *filename = new char[1024];

for (i=imin; i<=imax; i++) //over modes
{
    fprintf (stdout, "\n(m=%d, n=%d) mode...", modes[2*i], modes[2*i+1]); fflush(stdout);

    //EB is read always:
    for (int p=0; p<2; p++) //over paths
    {
        //EB arrays:
        sprintf (filename, "%slinear-data/m_%d_n_%d_flab_[%g,%g]/EB.dat", path2projects[p], modes[2*i], modes[2*i+1], real(flab), imag(flab));

        dim_EB[p] = count_lines_in_file (filename, 0);

        r_EB[p]    = new double[dim_EB[p]];
        prof_EB[p] = new double[ncols_EB*dim_EB[p]];

        load_data_file (filename, dim_EB[p], ncols_EB, r_EB[p], prof_EB[p]);
    }

    //current density: loads apropriate files
    if (output_flag == 1)
    {
        int p = 0; //plasma path

        sprintf (filename,"%slinear-data/m_%d_n_%d_flab_[%g,%g]/zone_0_current_dens_p_0_t_lab.dat", path2projects[p], modes[2*i], modes[2*i+1], real(flab), imag(flab));

        dim_jp[p] = count_lines_in_file (filename, 0);

        r_jp[p]    = new double[dim_jp[p]];
        prof_jp[p] = new double[ncols_jp*dim_jp[p]];

        load_data_file (filename, dim_jp[p], ncols_jp, r_jp[p], prof_jp[p]);
    }

    //torques: loads apropriate files
    if (output_flag == 2)
    {

    }

    //interpolation:
    for (int k=0; k<dim; k++) //over grid
    {
        switch (output_flag)
        {
            case 0: //Br formfactors
                evaluate_Br_formfactors (dim_EB[0], r_EB[0], prof_EB[0], dim_EB[1], r_EB[1], prof_EB[1], pol_deg, rout[k], &form_facs[k + (i-imin)*dim]);
            break;

            case 1: //current density over Br
                evaluate_jp_over_Br (dim_EB[0], r_EB[0], prof_EB[0], dim_jp[0], r_jp[0], prof_jp[0], pol_deg, rout[k], &form_facs[k + (i-imin)*dim]);
            break;

            case 2: //torques
                fprintf (stdout, "\nerror: not implemented!\n");
            break;

            default:
                fprintf (stdout, "\nerror: unknown output flag = %d\n", output_flag);
                exit(1);
            break;
        }
    }

    //removes noise at the center by linear interpolation:
    for (l=0; l<dim; l++) if (rout[l]>rmini) break;

    for (int k=0; k<l; k++)
    {
        //form_facs[i-imin][k] = form_facs[i-imin][l] +
        //(r[k]-r[l])/(r[l+1]-r[l])*(form_facs[i-imin][l+1]-form_facs[i-imin][l]);

        form_facs[k + (i-imin)*dim] = form_facs[l+0 + (i-imin)*dim] +
                                      (rout[k]-rout[l])/(rout[l+1]-rout[l])*
        (form_facs[l+1 + (i-imin)*dim] - form_facs[l+0 + (i-imin)*dim]);
    }

    //clean up:
    for (int p=0; p<2; p++)
    {
        if (r_EB[p])    delete [] r_EB[p];
        if (prof_EB[p]) delete [] prof_EB[p];
        if (r_jp[p])    delete [] r_jp[p];
        if (prof_jp[p]) delete [] prof_jp[p];
    }

    fprintf (stdout, " Ok."); fflush (stdout);
}

//save:
fprintf (stdout, "\nSaving..."); fflush (stdout);

//double * label = new double[dim];

transform_cyl_radius_to_surface_label (dim, rout, pol_deg, label_flag, label);

FILE *outfile;

sprintf (filename, "%s.%s.%s", conf, name_out[output_flag], name_label[label_flag]);

if (!(outfile = fopen (filename, "w")))
{
    fprintf (stderr, "\nFailed to open file %s\a\n", filename);
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
    fprintf (outfile, "%24.16le", label[k]);
    for (i=imin; i<=imax; i++)
    {
        fprintf (outfile, "\t%24.16le\t%24.16le",
                 real(form_facs[k + (i-imin)*dim]), imag(form_facs[k + (i-imin)*dim]));
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

if (format_flag == 0) sprintf (filename, "%s.uff", filename); //unformatted fortran file
else                  sprintf (filename, "%s.fff", filename); //formatted fortran file

double *form_facs_fort = new double[2*(imax-imin+1)*dim];

for (i=imin; i<=imax; i++)
{
    for (int k=0; k<dim; k++)
    {
        form_facs_fort[0 + 2*k + (i-imin)*dim*2] = real(form_facs[k + (i-imin)*dim]);
        form_facs_fort[1 + 2*k + (i-imin)*dim*2] = imag(form_facs[k + (i-imin)*dim]);
    }
}

save_form_facs_fortran_ (filename, &modes_dim, modes_mn, &m_min, &m_max, &n_min, &n_max,
                         modes_index, &dim, label, form_facs_fort, &format_flag);

fprintf (stdout, " Ok.\n"); fflush (stdout);

//clean up;
delete [] str_buf;
delete [] conf;
delete [] path2projects[0];
delete [] path2projects[1];

delete [] modes;
delete [] filename;
delete [] form_facs;
delete [] form_facs_fort;

delete [] rout;
delete [] label;

delete [] modes_mn;
delete [] modes_index;

return 0;
}

/*****************************************************************************/

void transform_cyl_radius_to_surface_label (int dim, double * r, int pol_deg,
                                            int label_flag, double * label)
{
if (label_flag == 0) //cylindrical radius
{
    for (int k=0; k<dim; k++)
    {
        label[k] = r[k];
    }
    return;
}

//reading of the equlibrium file:
int ncols = 5;
double * x_eq = NULL, * prof_eq = NULL, * label_eq = NULL;

char file_eq[] = "equil_r_q_psi.dat";

int dim_eq = count_lines_in_file_with_comments (file_eq, 1);

x_eq = new double[dim_eq];
prof_eq = new double[ncols*dim_eq];

load_data_file_with_comments (file_eq, dim_eq, ncols, x_eq, prof_eq);

label_eq = new double[dim_eq];

double norm_coeff;

switch (label_flag)
{
    case 1: //sqrt(|psi_pol|)
    for (int k=0; k<dim_eq; k++) label_eq[k] = sqrt(abs(prof_eq[k + 1*dim_eq]));
    break;

    case 2: //psi_pol
    for (int k=0; k<dim_eq; k++) label_eq[k] = prof_eq[k + 1*dim_eq];
    break;

    case 3: //s_pol = psi_pol/psi_pol_sep
    norm_coeff = prof_eq[dim_eq-1 + 1*dim_eq];
    for (int k=0; k<dim_eq; k++) label_eq[k] = prof_eq[k + 1*dim_eq]/norm_coeff;
    break;

    case 4: //sqrt(|psi_tor|)
    for (int k=0; k<dim_eq; k++) label_eq[k] = sqrt(abs(prof_eq[k + 2*dim_eq]));
    break;

    case 5: //psi_tor
    for (int k=0; k<dim_eq; k++) label_eq[k] = prof_eq[k + 2*dim_eq];
    break;

    case 6: //s_tor = psi_tor/psi_tor_sep
    norm_coeff = prof_eq[dim_eq-1 + 2*dim_eq];
    for (int k=0; k<dim_eq; k++) label_eq[k] = prof_eq[k + 2*dim_eq]/norm_coeff;
    break;

    default:
    fprintf (stdout, "\nerror: unknown grid flag = %d\n", label_flag);
    exit(1);
}

delete [] prof_eq;

int ind = 0;

for (int k=0; k<dim; k++)
{
    find_index_for_interp (pol_deg, r[k], dim_eq, x_eq, &ind);

    eval_neville_polynom (x_eq+ind, label_eq+ind, 1, pol_deg, r[k], label+k);
}

FILE *outfile;

if (!(outfile = fopen ("label_grid.dat", "w")))
{
    fprintf (stderr, "\nFailed to open file %s\a\n", "label_grid.dat");
}

for (int k=0; k<dim; k++)
{
    fprintf (outfile, "%24.16le\t%24.16le\n", r[k], label[k]);
}

fclose (outfile);

delete [] x_eq;
delete [] label_eq;
}

/*****************************************************************************/

void transform_surface_label_to_cyl_radius (int dim, double * r, int pol_deg,
                                            int label_flag, double * label)
{
if (label_flag == 0) //cylindrical radius
{
    for (int k=0; k<dim; k++)
    {
        r[k] = label[k];
    }
    return;
}

//reading of the equlibrium file:
int ncols = 5;
double * x_eq = NULL, * prof_eq = NULL, * label_eq = NULL;

char file_eq[] = "equil_r_q_psi.dat";

int dim_eq = count_lines_in_file_with_comments (file_eq, 1);

x_eq = new double[dim_eq];
prof_eq = new double[ncols*dim_eq];

load_data_file_with_comments (file_eq, dim_eq, ncols, x_eq, prof_eq);

label_eq = new double[dim_eq];

double norm_coeff;

switch (label_flag)
{
    case 1: //sqrt(|psi_pol|)
    for (int k=0; k<dim_eq; k++) label_eq[k] = sqrt(abs(prof_eq[k + 1*dim_eq]));
    break;

    case 2: //psi_pol
    for (int k=0; k<dim_eq; k++) label_eq[k] = prof_eq[k + 1*dim_eq];
    break;

    case 3: //s_pol = psi_pol/psi_pol_sep
    norm_coeff = prof_eq[dim_eq-1 + 1*dim_eq];
    for (int k=0; k<dim_eq; k++) label_eq[k] = prof_eq[k + 1*dim_eq]/norm_coeff;
    break;

    case 4: //sqrt(|psi_tor|)
    for (int k=0; k<dim_eq; k++) label_eq[k] = sqrt(abs(prof_eq[k + 2*dim_eq]));
    break;

    case 5: //psi_tor
    for (int k=0; k<dim_eq; k++) label_eq[k] = prof_eq[k + 2*dim_eq];
    break;

    case 6: //s_tor = psi_tor/psi_tor_sep
    norm_coeff = prof_eq[dim_eq-1 + 2*dim_eq];
    for (int k=0; k<dim_eq; k++) label_eq[k] = prof_eq[k + 2*dim_eq]/norm_coeff;
    break;

    default:
    fprintf (stdout, "\nerror: unknown grid flag = %d\n", label_flag);
    exit(1);
}

delete [] prof_eq;

int ind = 0;

for (int k=0; k<dim; k++)
{
    find_index_for_interp (pol_deg, label[k], dim_eq, label_eq, &ind);

    eval_neville_polynom (label_eq+ind, x_eq+ind, 1, pol_deg, label[k], r+k);
}

FILE *outfile;

if (!(outfile = fopen ("radial_grid.dat", "w")))
{
    fprintf (stderr, "\nFailed to open file %s\a\n", "radial_grid.dat");
}

for (int k=0; k<dim; k++)
{
    fprintf (outfile, "%24.16le\t%24.16le\n", r[k], label[k]);
}

fclose (outfile);

delete [] x_eq;
delete [] label_eq;
}

/*****************************************************************************/

inline void evaluate_Br_formfactors (int dim1, double * r_EB1, double * prof_EB1,
                                     int dim2, double * r_EB2, double * prof_EB2,
                                     int pol_deg, double r, complex<double> *form_facs)
{
complex<double> Q[2];

double re_Q, im_Q;

int re_Br_ind = 7-1; //index starts from zero
int im_Br_ind = 8-1; //index starts from zero

static int ind1 = dim1/2; find_index_for_interp (pol_deg, r, dim1, r_EB1, &ind1);

eval_neville_polynom (r_EB1+ind1, prof_EB1+dim1*re_Br_ind+ind1, 1, pol_deg, r, &re_Q);
eval_neville_polynom (r_EB1+ind1, prof_EB1+dim1*im_Br_ind+ind1, 1, pol_deg, r, &im_Q);

Q[0] = re_Q + I*im_Q; //Br

static int ind2 = dim2/2; find_index_for_interp (pol_deg, r, dim2, r_EB2, &ind2);

eval_neville_polynom (r_EB2+ind2, prof_EB2+dim2*re_Br_ind+ind2, 1, pol_deg, r, &re_Q);
eval_neville_polynom (r_EB2+ind2, prof_EB2+dim2*im_Br_ind+ind2, 1, pol_deg, r, &im_Q);

Q[1] = re_Q + I*im_Q; //Br

//*form_facs = Q[0];
*form_facs = Q[0]/Q[1];
}

/*****************************************************************************/

inline void evaluate_jp_over_Br (int dim_EB, double * r_EB, double * prof_EB,
                                 int dim_jp, double * r_jp, double * prof_jp,
                                 int pol_deg, double r, complex<double> *form_facs)
{
complex<double> Q[2];

double re_Q, im_Q;

int re_Br_ind = 7-1; //index starts from zero
int im_Br_ind = 8-1; //index starts from zero

static int ind1 = dim_EB/2; find_index_for_interp (pol_deg, r, dim_EB, r_EB, &ind1);

eval_neville_polynom (r_EB+ind1, prof_EB+dim_EB*re_Br_ind+ind1, 1, pol_deg, r, &re_Q);
eval_neville_polynom (r_EB+ind1, prof_EB+dim_EB*im_Br_ind+ind1, 1, pol_deg, r, &im_Q);

Q[0] = re_Q + I*im_Q; //Br

int re_jp_ind = 0; //index starts from zero
int im_jp_ind = 1; //index starts from zero

static int ind2 = dim_jp/2; find_index_for_interp (pol_deg, r, dim_jp, r_jp, &ind2);

eval_neville_polynom (r_jp+ind2, prof_jp+dim_jp*re_jp_ind+ind2, 1, pol_deg, r, &re_Q);
eval_neville_polynom (r_jp+ind2, prof_jp+dim_jp*im_jp_ind+ind2, 1, pol_deg, r, &im_Q);

Q[1] = re_Q + I*im_Q; //jp

//*form_facs = Q[1];
//*form_facs = Q[1]/Q[0];
*form_facs = Q[1]/Q[0]/r;
}

/*****************************************************************************/
