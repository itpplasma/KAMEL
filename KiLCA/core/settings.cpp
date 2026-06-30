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

    bs = new back_sett;
    bs->read_settings (path2project);
    copy_background_data_to_background_module_ (&bs);

    os = new output_sett;
    os->read_settings (path2project);

    es = new eigmode_sett;
    es->read_settings (path2project);
    std::cout << ">> KiLCA: Settings read successfully.\n";
}

/*******************************************************************/
