/*! \file output_sett.h
    \brief The declaration of output_sett class.
*/

#ifndef OUTPUT_SETTINGS_INCLUDE

#define OUTPUT_SETTINGS_INCLUDE

#include "constants.h"

/*****************************************************************************/

/*! \class output_sett
    \brief Class for output settings.
*/
class output_sett
{
public:
    //Output settings:
    int flag_background;   //!<1 if compute background data, 2 - store
    int flag_emfield;      //!<1 if compute em field data, 2 - store
    int flag_additional;   //!<1 if compute additional quants, 2 - store
    int flag_dispersion;   //!<1 if compute dispersion, 2 - store

    int num_quants;        //!<number of flags
    int *flag_quants;      //!<flags for each quantitity if compute it

    int flag_debug;        //!<flag for debugging mode

public:
    output_sett (void) {}

    ~output_sett (void)
    {
        delete [] flag_quants;
    }

    void read_settings (char *path);
    void print_settings (void);
};

/*****************************************************************************/

#endif
