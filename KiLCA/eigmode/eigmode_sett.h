/*! \file eigmode_sett.h
    \brief C entry points for eigenmode search settings, now owned by the
           Fortran eigmode_sett_data module (read_eigmode_settings_ parses
           eigmode.in). The former C++ eigmode_sett class has been translated
           away.
*/

#ifndef EIGMODE_SETT_INCLUDE

#define EIGMODE_SETT_INCLUDE

extern "C"
{
void read_eigmode_settings_ (char *path);

int get_eigmode_search_flag_ (void);

double get_eigmode_delta_ (void);

double get_eigmode_rfmin_ (void);

double get_eigmode_rfmax_ (void);

double get_eigmode_ifmin_ (void);

double get_eigmode_ifmax_ (void);

int get_eigmode_rdim_ (void);

int get_eigmode_idim_ (void);

int get_eigmode_n_zeros_ (void);

int get_eigmode_use_winding_ (void);

int get_eigmode_Nguess_ (void);

int get_eigmode_kmin_ (void);

int get_eigmode_kmax_ (void);

int get_eigmode_test_roots_ (void);

double get_eigmode_eps_abs_ (void);

double get_eigmode_eps_rel_ (void);

double get_eigmode_eps_res_ (void);

void get_eigmode_fname_ (char *out);

double get_eigmode_fstart_re_ (int k);

double get_eigmode_fstart_im_ (int k);
}

#endif
