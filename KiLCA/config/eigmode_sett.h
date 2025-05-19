/*! \file eigmode_sett.h
    \brief The declaration of eigmode_sett class.
*/

#ifndef EIGMODE_SETT_INCLUDE

#define EIGMODE_SETT_INCLUDE

#include <complex>

using namespace std;

/*****************************************************************************/

/*! \class eigmode_sett
    \brief Class represents set of settings for eigenmode search.
*/
class eigmode_sett
{
public:
    char *fname; //!<name of output file

    int search_flag; //!<flag specifies the search option

    ///parameters of a grid: real and imag parts
    int rdim;
    double rfmin;
    double rfmax;
    int idim;
    double ifmin;
    double ifmax;

    int stop_flag;

    ///accuracies
    double eps_res;
    double eps_abs;
    double eps_rel;

    double delta; //!<a delta for a numerical derivative

    int test_roots;

    int flag_debug;

    ///starting values for root search
    int Nguess;
    int kmin;
    int kmax;

    complex<double> *fstart;

    ///number of zeros to be found
    int n_zeros;

    ///flag specifies whether to use winding number evaluation or not
    int use_winding;

public:
    //eigmode_sett (void);

    ~eigmode_sett (void)
    {
        delete [] fname;
        delete [] fstart;
    }

    void read_settings (char *path);
    void print_settings ();
};

/*****************************************************************************/

#endif
