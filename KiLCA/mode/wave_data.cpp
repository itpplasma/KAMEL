/*! \file wave_data.cpp
    \brief The implementation of the wave_data accessor functions declared in wave_data.h.
*/

#include "constants.h"
#include "wave_data.h"

/*-----------------------------------------------------------------*/

double get_wave_data_obj_omov_re_ (const wave_data *wd)
{
return real (wd->omov);
}

double get_wave_data_obj_omov_im_ (const wave_data *wd)
{
return imag (wd->omov);
}

/*-----------------------------------------------------------------*/
