#ifndef ANTENNA_INCLUDE
#define ANTENNA_INCLUDE

#include <complex>

struct antenna_sett
{
    double ra; // small radius (cm) of antenna location
    double wa; // current density layer width
    double I0; // current in antenna coils (statamp)
    std::complex<double> flab; // frequency (Hz) in the laboratory frame
    int dma; // dimension of modes array
    int *modes; // array of modes (m,n)
    bool flag_eigmode; // flag for eigmode search
    bool flag_debug; // debug flag

    explicit antenna_sett() = default;

    void print_settings();
};

extern "C"
{
void set_antenna_settings_c_(antenna_sett **ptr, double *ra, double *wa,
                             double *I0, double *flab_re, double *flab_im,
                             int *dma, bool *flag_debug);
}

#endif
