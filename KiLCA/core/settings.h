/*! \file settings.h
    \brief The declaration of settings class.
*/

#ifndef SETTINGS_INCLUDE

#define SETTINGS_INCLUDE

#include <cstring>

#include "code_settings.h"

#include "antenna.h"
#include "back_sett.h"
#include "output_sett.h"
#include "eigmode_sett.h"

using namespace :: std;

/*! \class settings
    \brief The class contains all KiLCA settings.
*/
class settings
{
public:
    char *path2project; //!<project path

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

        delete es;
    }

    void read_settings (void); //!<read settings from KiLCA input files
};

/*******************************************************************/

#endif
