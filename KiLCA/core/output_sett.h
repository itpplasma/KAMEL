/*! \file output_sett.h
    \brief C entry points for output settings, now owned by the Fortran
           output_data module (read_output_settings_ parses output.in). The
           former C++ output_sett class has been translated away.
*/

#ifndef OUTPUT_SETTINGS_INCLUDE

#define OUTPUT_SETTINGS_INCLUDE

#include "constants.h"

extern "C"
{
void read_output_settings_ (char *path);

int get_output_flag_background_ (void);

int get_output_flag_emfield_ (void);

int get_output_flag_additional_ (void);

int get_output_flag_dispersion_ (void);

int get_output_num_quants_ (void);

int get_output_flag_quants_ (int i);
}

#endif
