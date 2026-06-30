/*! \file settings.h
    \brief C entry points for the aggregated KiLCA settings, now owned by
           the Fortran kilca_settings_m module (KiLCA/core/settings_m.f90).
           The former C++ settings class has been translated away; still-C++
           callers (core.cpp, eigmode/*.cpp) hold an opaque intptr_t handle
           instead of a `settings*`.
*/

#ifndef SETTINGS_INCLUDE

#define SETTINGS_INCLUDE

#include <cstdint>

#include "code_settings.h"

#include "antenna.h"
#include "back_sett.h"
#include "output_sett.h"
#include "eigmode_sett.h"

extern "C"
{
intptr_t settings_create_ (const char *path);

void settings_destroy_ (intptr_t handle);

void settings_read_settings_ (intptr_t handle);

void settings_get_path2project_ (intptr_t handle, char *buf);
}

/*******************************************************************/

#endif
