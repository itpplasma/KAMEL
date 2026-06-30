/*! \file antenna.h
    \brief C entry points for the antenna settings, now owned by the Fortran
           antenna_data module (read_antenna_settings_ parses antenna.in and
           modes.in). The former C++ antenna class has been translated away.
*/

#ifndef ANTENNA_INCLUDE

#define ANTENNA_INCLUDE

extern "C"
{
void read_antenna_settings_ (char *path);

int get_antenna_dma_ (void);

void get_antenna_flab_ (double *re, double *im);

void get_antenna_mode_ (int ind, int *m, int *n);

double get_antenna_ra_ (void);

double get_antenna_wa_ (void);

int get_antenna_flag_eigmode_ (void);
}

#endif
