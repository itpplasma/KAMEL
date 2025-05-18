/*! \file settings.h
    \brief The declaration of settings class.
*/

#ifndef SETTINGS_INCLUDE
#define SETTINGS_INCLUDE

#include "antenna.h"
#include "back_sett.h"
#include "output_sett.h"
#include "eigmode_sett.h"

#include <cstring>

using namespace :: std;

/*! \class settings
    \brief The class contains all KiLCA settings.
*/
class settings
{
public:
    char *path2project; //!<project path

    antenna      *as; //!<antenna settings

    back_sett    *bs; //!<background settings

    output_sett  *os; //!<output settings

    eigmode_sett *es; //!<settings for eigmode search

public:
    settings (char *path)
    {
        path2project = new char[1024];
        strcpy (path2project, path);
    }

    ~settings (void)
    {
        delete [] path2project;

        delete as;
        delete bs;
        delete os;
        delete es;
    }

    void read_settings (void); //!<read settings from KiLCA input files
};

/*******************************************************************/

extern "C"
{
void copy_antenna_data_to_antenna_module_ (antenna **);

void copy_background_data_to_background_module_ (back_sett **);
}

/*******************************************************************/

#endif
