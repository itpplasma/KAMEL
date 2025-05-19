#ifndef EIGMODE_SETT_INCLUDE
#define EIGMODE_SETT_INCLUDE

#include <complex>

struct eigmode_sett
{
    char fname[256]; // name of output file

    bool search_flag; // flag specifies the search option

    ///parameters of a grid: real and imag parts
    int rdim;
    double rfmin;
    double rfmax;
    int idim;
    double ifmin;
    double ifmax;

    bool stop_flag;

    ///accuracies
    double eps_res;
    double eps_abs;
    double eps_rel;

    double delta; // a delta for a numerical derivative

    bool test_roots;

    bool flag_debug;

    ///starting values for root search
    int Nguess;
    int kmin;
    int kmax;

    std::complex<double> *fstart;

    ///number of zeros to be found
    int n_zeros;

    ///flag specifies whether to use winding number evaluation or not
    bool use_winding;

    explicit eigmode_sett() = default;

    void print_settings();
};

#endif
