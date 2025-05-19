/*! \file settings.cpp
    \brief The implementation of settings class.
*/

#include <stdlib.h>
#include <stdio.h>
#include <cmath>
#include <cstring>

#include "settings.h"

/*******************************************************************/

void settings_t::read_settings (void)
{
as = new antenna;
as->read_settings (path2project);
copy_antenna_data_to_antenna_module_ (&as);

bs = new back_sett;
bs->read_settings (path2project);
copy_background_data_to_background_module_ (&bs);

os = new output_sett;
os->read_settings (path2project);

es = new eigmode_sett;
es->read_settings (path2project);
}

/*******************************************************************/
