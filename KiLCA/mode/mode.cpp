/*! \file
    \brief The implementation of mode_data class.
*/

#include <stdlib.h>
#include <stdio.h>
#include <string.h>

#include "settings.h"
#include "background.h"
#include "mode.h"

/*******************************************************************/

mode_data::mode_data (int m, int n, complex<double> olab, const settings *sd_p, const background *bp_p)
{
//set const pointers to common structures:
sd = sd_p;
bp = bp_p;

set_settings_in_mode_data_module_ (&sd);
set_back_profiles_in_mode_data_module_ (&bp);

complex<double> omov = olab - n*(get_background_V_gal_sys_())/(get_background_rtor_());

wd = wave_data_create_ (m, n, real(olab), imag(olab), real(omov), imag(omov));

set_wave_data_in_mode_data_module_ (&wd);

//set_wave_parameters_in_mode_data_module_ (&m, &n, &(real(olab)), &(imag(olab)),
//                                                  &(real(omov)), &(imag(omov)));

double olab_re = real(olab), olab_im = imag(olab);
double omov_re = real(omov), omov_im = imag(omov);

set_wave_parameters_in_mode_data_module_ (&m, &n, &olab_re, &olab_im, &omov_re, &omov_im);

if (DEBUG_FLAG)
{
    fprintf (stdout, "\n(m=%d, n=%d): f_lab=(%le, %le) f_mov=(%le, %le)\n",
                     m, n, wave_data_get_olab_re_(wd)/2.0/pi, wave_data_get_olab_im_(wd)/2.0/pi,
                           get_wave_data_obj_omov_re_(wd)/2.0/pi, get_wave_data_obj_omov_im_(wd)/2.0/pi);
}

if (get_output_flag_background_() > 0)
{
    find_resonance_location ();
    double r_res_local = wave_data_get_r_res_ (wd);
    set_resonance_location_in_mode_data_module_ (&r_res_local);
}

allocate_and_setup_zones ();

//Directories for the given mode:
path2linear = new char[1024];
path2dispersion = new char[1024];
path2poincare = new char[1024];

eval_path_to_linear_data (sd->path2project, m, n, olab, path2linear);

eval_path_to_dispersion_data (sd->path2project, m, n, olab, path2dispersion);

eval_path_to_poincare_data (sd->path2project, m, n, olab, path2poincare);

//makes directories for the mode data:
set_and_make_mode_data_directories ();

mode_data *tmp = this;
copy_mode_paths_to_mode_data_module_ (&tmp);

// initialize pointers:
r = 0;
EB = 0;
EB_int = 0;
index = 0;
A = 0;
B = 0;
S = 0;
}

/*******************************************************************/

void eval_path_to_linear_data (const char *path2project, int m, int n, complex<double> olab, char *path2linear)
{
complex<double> flab = olab/(2.0*pi);

//linear data directory:
sprintf (path2linear, "%s%s%d%s%d%s%.15lg%s%.15lg%s", path2project, "linear-data/m_", m, "_n_", n, "_flab_[", real(flab), ",", imag(flab), "]/");
}

/*******************************************************************/

void eval_path_to_dispersion_data (const char *path2project, int m, int n, complex<double> olab, char *path2dispersion)
{
complex<double> flab = olab/(2.0*pi);

//dispersion data directory:
sprintf (path2dispersion, "%s%s%d%s%d%s%.15lg%s%.15lg%s", path2project, "dispersion-data/m_", m, "_n_", n, "_flab_[", real(flab), ",", imag(flab), "]/");
}

/*******************************************************************/

void eval_path_to_poincare_data (const char *path2project, int m, int n, complex<double> olab, char *path2poincare)
{
complex<double> flab = olab/(2.0*pi);

//dispersion data directory:
sprintf (path2poincare, "%s%s%d%s%d%s%.15lg%s%.15lg%s", path2project, "poincare-data/m_", m, "_n_", n, "_flab_[", real(flab), ",", imag(flab), "]/");
}

/*******************************************************************/

void mode_data::set_and_make_mode_data_directories (void)
{
char *sys_command = new char[1024];

if (get_output_flag_emfield_() > 1)
{
    sprintf (sys_command, "%s%s", "mkdir -p ", path2linear);
    if (system (sys_command)==-1)
    {
        fprintf (stderr, "\nerror: make_mode_dirs: system()!");
    }

    //makes debug data directory:
    sprintf (sys_command, "%s%s%s", "mkdir -p ", path2linear, "debug-data/");
    if (system (sys_command)==-1)
    {
        fprintf (stderr, "\nerror: make_mode_dirs: system()!");
    }

    //file with the mode data:
    FILE *outfile;

    sprintf (sys_command, "%s%s", path2linear, "mode_data.dat");

    if (!(outfile = fopen (sys_command, "w")))
    {
        fprintf (stderr, "\nFailed to open file %s\a\n", sys_command);
    }

    fprintf (outfile, "%%m n\n%%Re(flab) Im(flab)\n%%Re(fmov) Im(fmov)\n%%r_res 0.0\n");
    fprintf (outfile, "%d %d\n", wave_data_get_m_(wd), wave_data_get_n_(wd));
    fprintf (outfile, "%.16lg %.16lg\n", wave_data_get_olab_re_(wd)/2.0/pi, wave_data_get_olab_im_(wd)/2.0/pi);
    fprintf (outfile, "%.16lg %.16lg\n", get_wave_data_obj_omov_re_(wd)/2.0/pi, get_wave_data_obj_omov_im_(wd)/2.0/pi);
    fprintf (outfile, "%.16lg %.16lg\n", wave_data_get_r_res_(wd), 0.0);

    fclose (outfile);
}

if (get_output_flag_dispersion_() > 1)
{
    //makes dispersion data directory:
    sprintf (sys_command, "%s%s", "mkdir -p ", path2dispersion);
    if (system (sys_command)==-1)
    {
        fprintf (stderr, "\nerror: make_mode_dirs: system()!");
    }
}

//makes poincare data directory:
//sprintf (sys_command, "%s%s", "mkdir -p ", path2poincare);
//if (system (sys_command)==-1)
//{
//    fprintf (stderr, "\nerror: make_mode_dirs: system()!");
//}

delete [] sys_command;
}

/*****************************************************************************/

void copy_mode_paths_from_mode_data_struct_ (mode_data **md, char *path2linear, char *path2dispersion, char *path2poincare)
{
strcpy (path2linear,     (*md)->path2linear);
strcpy (path2dispersion, (*md)->path2dispersion);
strcpy (path2poincare,   (*md)->path2poincare);
}

/*****************************************************************************/

void mode_data::save_mode_det_data (void)
{
char *filename = new char[1024];
FILE *outfile;

sprintf (filename, "%smode_%d_%d_[%.16le,%.16le].dat", path2linear, wave_data_get_m_(wd), wave_data_get_n_(wd),
         wave_data_get_olab_re_(wd), wave_data_get_olab_im_(wd));

if (!(outfile = fopen (filename, "w")))
{
    fprintf (stderr, "\nFailed to open file %s\a\n", filename);
}

fprintf (outfile,"%.16le  %.16le\t%.16le  %.16le\n", wave_data_get_olab_re_(wd), wave_data_get_olab_im_(wd),
          wave_data_get_det_re_(wd), wave_data_get_det_im_(wd));

fclose (outfile);

delete [] filename;
}

/*****************************************************************************/
