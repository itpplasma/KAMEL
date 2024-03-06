/*! \file
    \brief The additional driver program to make parameter studies of eigenmodes.
*/

#include <stdlib.h>
#include <stdio.h>
#include <cmath>
#include <cstring>
#include <climits>
#include <inttypes.h>
#include <cstring>
#include <dirent.h>
#include <fnmatch.h>
#include <iostream>
#include <fstream>
#include <unistd.h>

#include "inout.h"
#include "core.h"
#include "shared.h"
#include "main_linear.h"

using namespace std;

/*******************************************************************/

void make_run (char * fullpath, char * param_name, double p_curr, complex<double> f_est);

void edit_background_file_Vscale (char * fullpath, double p_curr, complex<double> f_est);

void edit_eigmode_file_Vscale (char * fullpath, double p_curr, complex<double> f_est);

void edit_input_files_Vscale (char * fullpath, double p_curr, complex<double> f_est);

void run_kilca (char * path);

void remove_run (char * fullpath);

void clean_run (char * fullpath, complex<double> funct_first, complex<double> funct_last);

void save_result (char * filename, double p_curr, complex<double> f_curr);

int get_root_value (char * fullpath, double eps_det, double eps_root, complex<double> f_est, complex<double> * f_curr);

/*******************************************************************/

//!The main program
/*!
Calls KilCA library functions
*/

int main (int argc, char **argv)
{
std::ofstream file ("kilca_version", std::ofstream::out); file << "KiLCA version used for the last run: " << argv[0]; file.close();

//gets path to the project
char * path = new char[1024];
char * fullpath = new char[1024];

if (argc == 1)
{
   getcwd (path, 1024); //current directory
}
else if (argc >= 2)
{
   strcpy (path, argv[1]); //command line argument
}

if (path[strlen(path)-1] != '/') strcat(path, "/");

//read input settings:
double p1, p2;
double dp_est, dp_min, dp_max;
double eps_det, eps_root;
complex<double> f_start;

double del = 0.1;

int max_iter_num = 10000;

char *file_set = new char[1024];

sprintf (file_set, "%s/search.in", path);

FILE *in;

if ((in=fopen (file_set, "r"))==NULL)
{
    fprintf(stderr, "\nerror: read_settings: failed to open file %s\a\n", file_set);
    exit(0);
}

//reading: check for consistence with file!
char *str_buf = new char[1024]; //str buffer

read_line_2skip_it (in, &str_buf);

char * config = new char[1024];
read_line_2get_string (in, &config);

char * param_name = new char[1024];
read_line_2get_string (in, &param_name);

char * func_name = new char[1024];
read_line_2get_string (in, &func_name);

read_line_2get_double (in, &p1);
read_line_2get_double (in, &p2);
read_line_2get_double (in, &dp_min);
read_line_2get_double (in, &dp_est);
read_line_2get_double (in, &dp_max);
read_line_2get_double (in, &eps_det);
read_line_2get_double (in, &eps_root);
read_line_2get_complex(in, &f_start);
read_line_2skip_it (in, &str_buf);

fclose (in);

delete [] str_buf;
delete [] file_set;

char * filename = new char[1024];
sprintf (filename, "%s%s_%s_%s.dat", path, func_name, param_name, config);

double * param = new double[max_iter_num];
complex<double> * funct = new complex<double>[max_iter_num];

//for loop:
double p_curr, p_prev, dp;
complex<double> f_prev, f_curr, f_est;

p_prev = p1;
p_curr = p1;
dp = dp_est;

f_prev = f_start;
f_curr = f_start;
f_est  = f_start;

int step = 0, ind = -1, num_tries = 0;

while (abs(p_curr) < abs(p2))
{
    step++;

    if (ind == max_iter_num)
    {
        fprintf (stderr, "\nerror: maximum number of steps is reached!");
        exit (1);
    }

    sprintf (fullpath, "%s%s_%lg/", path, param_name, p_curr);

    make_run (fullpath, param_name, p_curr, f_start);

    int status = get_root_value (fullpath, eps_det, eps_root, f_est, &f_curr);

    if (status > 0) //code failed
    {
//        num_tries++;

//        fprintf (stdout, "\nwarning: failed to find a root: attempt = %d.", num_tries);

        remove_run (fullpath);

        if (step == 1)
        {
            fprintf (stdout, "\nerror: root approval failed for the first value: %le, try another starting guess.", p_curr);
            if (status != 1)
            {
                fprintf (stdout, "\nroot is found: %le%+lei, use it as initial guess.", real(f_curr), imag(f_curr));
            }

            exit (1);
        }

        dp /= 2.0;

        dp = signum (dp) * min (abs(dp), abs(dp_max));

        if (abs(dp) < abs(dp_min))
        {
            fprintf (stdout, "\nwarning: too small dp = %lg.", dp);
        }

        p_curr = p_prev + dp;

        if (ind > 1)
        {
            f_est = funct[ind] + ((funct[ind] - funct[ind-1])/(param[ind] - param[ind-1]))*
                                 (p_curr - param[ind]);
        }
        else
        {
            f_est = f_prev;
        }

        fprintf (stdout, "\nf_est = %lg %lgi", real(f_est), imag(f_est));

        f_start = f_est;
    }
    else
    {
        ind++;

        fprintf (stdout, "\nmessage: succeeded to find a root.");

        fprintf (stdout, "\nind = %d: at %s = %lg found a root = %lg%+lgi", ind, param_name, p_curr, real(f_curr), imag(f_curr));

        clean_run (fullpath, f_start, f_curr);

        param[ind] = p_curr;
        funct[ind] = f_curr;

        save_result (filename, p_curr, f_curr);

        p_prev = p_curr;
        f_prev = f_curr;

        dp *= 2.0;

        dp = signum (dp) * min (abs(dp), abs(dp_max));

        p_curr = p_prev + dp;

        if (ind > 1)
        {
            f_est = funct[ind] + ((funct[ind] - funct[ind-1])/(param[ind] - param[ind-1]))*
                                 (p_curr - param[ind]);
        }
        else
        {
            f_est = f_prev;
        }

        fprintf (stdout, "\nf_est = %lg%+lgi", real(f_est), imag(f_est));

        f_start = f_est;
    }
}

delete [] filename;
delete [] config;
delete [] path;
delete [] fullpath;
delete [] param;
delete [] funct;

return 0;
}

/*******************************************************************/

void make_run (char * fullpath, char * param_name, double p_curr, complex<double> f_start)
{
fprintf (stdout, "\n\nThe next run is: %s = %lg\tstart = %lg%+lgi", param_name, p_curr, real(f_start), imag(f_start)); fflush (stdout);

char * sys_command = new char[1024];

//make run directory:
sprintf (sys_command, "%s%s", "mkdir -p ", fullpath);
if (system (sys_command)==-1)
{
    fprintf (stderr, "\nerror: system()!");
    exit (1);
}

//copy input files:
sprintf (sys_command, "%s%s", "cp -L -p *.in ", fullpath);
if (system (sys_command)==-1)
{
    fprintf (stderr, "\nerror: system()!");
    exit (1);
}

//modify files:
edit_input_files_Vscale (fullpath, p_curr, f_start);

//run code:
run_kilca (fullpath);

delete [] sys_command;
}

/*******************************************************************/

void run_kilca (char * path)
{
//!Allocates core data structure containing pointers to all important code data
core_data *cd = new core_data (path);
set_core_data_in_core_module_ (&cd);

cd->calc_and_set_mode_independent_core_data ();

if (cd->sd->as->flag_eigmode == 0)
{
    cd->calc_and_set_mode_dependent_core_data_antenna ();
}
else
{
    cd->calc_and_set_mode_dependent_core_data_eigmode ();
}

delete cd;
}

/*******************************************************************/

void edit_input_files_Vscale (char * fullpath, double p_curr, complex<double> f_prev)
{
edit_background_file_Vscale (fullpath, p_curr, f_prev);
edit_eigmode_file_Vscale (fullpath, p_curr, f_prev);
}

/*******************************************************************/

int get_root_value (char * fullpath, double eps_det, double eps_root, complex<double> f_est, complex<double> * f_curr)
{
char * filename = new char[1024];

sprintf (filename, "%sroots.dat", fullpath);

FILE *in;
if (!(in = fopen (filename, "r")))
{
    fprintf (stderr, "\nFailed to open file %s\a\n", filename);
}

size_t nchar = 1024;

char * str1 = new char[nchar];
char * str2 = new char[nchar];
char * buff = new char[nchar];

int status = 1;

int iter, ind = 0;
double freq_re, freq_im, det_re, det_im;
double freq_re_g, freq_im_g, det_re_g, det_im_g;

while (getline (&buff, &nchar, in) > 0)
{
    ind++;

    strcpy (str1, str2);
    strcpy (str2, buff);

    if (ind == 2) //second line in the file
    {
        sscanf (str2, "%5u\t%le  %le\t%le  %le", &iter, &freq_re_g, &freq_im_g, &det_re_g, &det_im_g);
        //fprintf (stdout, "\ndet at first point: %le %le", det_re_g, det_im_g);
    }

    char * p1 = strstr (buff, "status");

    if (p1 == NULL) continue; //not finished

    sscanf (str1, "%5u\t%le  %le\t%le  %le", &iter, &freq_re, &freq_im, &det_re, &det_im);

    //determine if one gets a root:

    char * p2 = strstr (buff, "success");

    if (p2 != NULL) //root has found
    {
        fprintf (stdout, "\nsolver returned success.");
    }
    else
    {
        fprintf (stdout, "\nwarning: solver failed to locate a root.");
    }

    //compare determinant values:
    double ratio = abs(det_re+I*det_im)/abs(det_re_g+I*det_im_g);

    fprintf (stdout, "\ncheck: ratio = %lg", ratio);

    //check if the root is close to previous roots:
    //double err = abs(freq_re - freq_re_g + I*(freq_im - freq_im_g))/abs(freq_re+I*freq_im);

    //double err_re = abs(freq_re - real(f_est)); if (abs(freq_re) > 1.0) err_re /= abs(freq_re);
    double err_re = 0.0;
    double err_im = abs(freq_im - imag(f_est)); if (abs(freq_im) > 1.0) err_im /= abs(freq_im);
    double err = max (err_re, err_im);

    fprintf (stdout, "\ncheck: accuracy of the initial guess = %lg", err);

    if (ratio < eps_det && err < eps_root)
    {
        fprintf (stdout, "\nThe ratio of determinants and accuracy of the guess: PASSED.");
        status = 0;
    }
    else
    {
        fprintf (stdout, "\nThe ratio of determinants and accuracy of the guess: FAILED.");
        if (ratio > eps_det) status = 1;
        else if (ratio < eps_det && err > eps_root)  status = 2;
    }

    //more checks...

    if (status == 0)
    {
        fprintf (stdout, "\nroot: freq = %le%+lei\tdet = %le%+lei", freq_re, freq_im, det_re, det_im);
    }
    else
    {
        fprintf (stdout, "\nfailed to find a proper root.");
    }
}

*f_curr = freq_re + I*freq_im;

delete [] filename;
delete [] str1;
delete [] str2;
delete [] buff;

return status;
}

/*******************************************************************/

void remove_run (char * fullpath)
{
char * sys_command = new char[1024];

sprintf (sys_command, "%s%s", "rm -R -f ", fullpath);

if (system (sys_command)==-1)
{
    fprintf (stderr, "\nerror: system()!");
    exit (1);
}

delete [] sys_command;
}

/*******************************************************************/

void clean_run (char * fullpath, complex<double> funct_first, complex<double> funct_last)
{
char * buff = new char[1024];

sprintf (buff, "%slinear-data/", fullpath);

DIR *dp;
struct dirent *ep;

char * str[5]; for (int i=0; i<5; i++) str[i] = new char[1024];

sprintf (str[0], "*\\[%.15lg%s%.15lg\\]*", real(funct_first), ",", imag(funct_first));
sprintf (str[1], "*\\[%.15lg%s%.15lg\\]*", real(funct_last), ",", imag(funct_last));
sprintf (str[2], ".");
sprintf (str[3], "..");
sprintf (str[4], "...");

if (dp = opendir (buff))
{
    while (ep = readdir (dp))
    {
        int flag = 0;

        for (int i=0; i<5; i++)
        {
            if (!fnmatch (str[i], ep->d_name, 0))
            {
                flag = 1;
                break;
            }
        }

        if (flag) continue;

        sprintf (buff, "rm -R -f %slinear-data/%s", fullpath, ep->d_name);

        if (system (buff) == -1)
        {
            fprintf (stderr, "\nerror: system()!");
        }

    }

    closedir (dp);
}
else
{
    fprintf (stderr, "\nclean_run: faled to open the directory %s.", buff);
}

sprintf (buff, "rm -R -f %sdispersion-data", fullpath); system (buff);
sprintf (buff, "rm -R -f %spoincare-data", fullpath);   system (buff);

delete [] buff;

for (int i=0; i<5; i++) delete [] str[i];
}

/*******************************************************************/

void edit_background_file_Vscale (char * fullpath, double p_curr, complex<double> f_prev)
{
char * filename_in = new char[1024];

sprintf (filename_in, "%sbackground.in", fullpath);

FILE *in;
if (!(in = fopen (filename_in, "r")))
{
    fprintf (stderr, "\nFailed to open file %s\a\n", filename_in);
}

char * filename_out = new char[1024];

sprintf (filename_out, "%sbackground.out", fullpath);

FILE *out;
if (!(out = fopen (filename_out, "w")))
{
    fprintf (stderr, "\nFailed to open file %s\a\n", filename_out);
}

size_t nchar = 1024;

char * buff = new char[nchar];

int count = 0;

while (getline (&buff, &nchar, in) > 0)
{
    count++;

    if (count == 12)
    {
        fprintf (out, "%lg\t\t#V_scale: scale factor for the Vz velocity profile: Vz = V_scale*Vz - V_gal_sys\n", p_curr);
    }
    else
    {
        fprintf (out, "%s", buff);
    }
}

fclose (in);
fclose (out);

//copy:
char * sys_command = new char[1024];

sprintf (sys_command, "mv -f %s %s", filename_out, filename_in);

if (system (sys_command)==-1)
{
    fprintf (stderr, "\nerror: system()!");
    exit (1);
}

delete [] sys_command;
delete [] buff;
delete [] filename_in;
delete [] filename_out;
}

/*******************************************************************/

void edit_eigmode_file_Vscale (char * fullpath, double p_curr, complex<double> f_prev)
{
char * filename_in = new char[1024];

sprintf (filename_in, "%seigmode.in", fullpath);

FILE *in;
if (!(in = fopen (filename_in, "r")))
{
    fprintf (stderr, "\nFailed to open file %s\a\n", filename_in);
}

char * filename_out = new char[1024];

sprintf (filename_out, "%seigmode.out", fullpath);

FILE *out;
if (!(out = fopen (filename_out, "w")))
{
    fprintf (stderr, "\nFailed to open file %s\a\n", filename_out);
}

size_t nchar = 1024;

char * buff = new char[nchar];

int count = 0;

while (getline (&buff, &nchar, in) > 0)
{
    count++;

    if (count == 34)
    {
        fprintf (out, "(%.15lg, %.15lg)\t\t#starting point\n", real(f_prev), imag(f_prev));
    }
    else
    {
        fprintf (out, "%s", buff);
    }
}

fclose (in);
fclose (out);

//copy:
char * sys_command = new char[1024];

sprintf (sys_command, "mv -f %s %s", filename_out, filename_in);

if (system (sys_command)==-1)
{
    fprintf (stderr, "\nerror: system()!");
    exit (1);
}

delete [] sys_command;
delete [] buff;
delete [] filename_in;
delete [] filename_out;
}

/*******************************************************************/

void save_result (char * filename, double p_curr, complex<double> f_curr)
{
FILE *out;
if (!(out = fopen (filename, "a")))
{
    fprintf (stderr, "\nFailed to open file %s\a\n", filename);
}

fprintf (out, "%.20le\t%.20le\t%.20le\n", p_curr, real(f_curr), imag(f_curr));

fclose (out);
}

/*******************************************************************/
