/*! \file settings.cpp
    \brief The implementation of settings class.
*/

#include <stdlib.h>
#include <iostream>
#include <cmath>
#include <cstring>

#include "settings.h"

/*******************************************************************/

void settings::read_settings (void)
{
    std::cout << ">> KiLCA: Reading settings from " << path2project << '\n';

    read_antenna_settings_ (path2project);

    read_background_settings_ (path2project);

    read_output_settings_ (path2project);

    es = new eigmode_sett;
    es->read_settings (path2project);
    std::cout << ">> KiLCA: Settings read successfully.\n";
}

/*******************************************************************/
