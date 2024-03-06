/*! \file antenna.h
    \brief The declaration of antenna class representing the antenna settings.
*/

#ifndef ANTENNA_INCLUDE

#define ANTENNA_INCLUDE

#include <inttypes.h>

#include "constants.h"

/*! \class antenna
    \brief The class contains the antenna parameters.
*/
class antenna
{
public:
    double ra; //!<small radius (cm) of antenna location

    double wa; //!<current density layer width

    double I0; //!<current in antenna coils (statamp)

    complex<double> flab; //!<frequency (Hz) in the laboratory frame

    int dma; //!<dimension of modes array

    int *modes; //!<array of modes (m,n)

    int flag_debug; //!>debug flag

    int flag_eigmode; //!<flag for eigmode search

public:
    antenna (void) {}

    ~antenna (void)
    {
        delete [] modes;
    }

    void read_settings (char *path);
    void print_settings ();
};

extern "C"
{
void set_antenna_settings_c_ (antenna **ptr, double *ra, double *wa, double *I0, double *flab_re, double *flab_im, int *dma, int *flag_debug);
}

#endif
