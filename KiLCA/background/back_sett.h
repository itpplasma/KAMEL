/*! \file back_sett.h
    \brief C entry points for background settings, now owned by the Fortran
           background module (read_background_settings_ parses background.in).
           The former C++ back_sett class has been translated away.
*/

#ifndef BACK_SETT_INCLUDE

#define BACK_SETT_INCLUDE

extern "C"
{
void read_background_settings_ (char *path);

double get_background_rtor_ (void);

double get_background_rp_ (void);

double get_background_B0_ (void);

double get_background_V_gal_sys_ (void);

double get_background_V_scale_ (void);

double get_background_zele_ (void);

double get_background_zion_ (void);

int get_background_flag_debug_ (void);

double get_background_huge_factor_ (void);

int get_background_calc_back_ (void);

int get_background_N_ (void);

double get_background_mass_ (int i);

double get_background_charge_ (int i);

char get_background_flag_back_ (void);

void get_background_path2profiles_ (char *out);
}

#endif
