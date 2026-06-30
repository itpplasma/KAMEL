/*! \file
    \brief Declarations for flre_zone-specific entry points (kilca_flre_zone_m,
    KiLCA/flre/flre_zone_m.f90) not covered by the generic zone_t dispatch
    shim (zone_dispatch.h). Still-C++ callers (wave_code_interface.cpp) that
    used to `static_cast<flre_zone*>` a `zone*` now just pass the same
    opaque intptr_t handle straight through - the Fortran side does its own
    `select type` check.
*/

#ifndef FLRE_ZONE_DISPATCH_INCLUDE
#define FLRE_ZONE_DISPATCH_INCLUDE

#include <cstdint>

extern "C"
{
void activate_fortran_modules_for_zone_ (intptr_t *ptr);
void deactivate_fortran_modules_for_zone_ (intptr_t *ptr);

int flre_zone_get_flre_order_ (intptr_t handle);
intptr_t flre_zone_get_cp_ (intptr_t handle);
}

#endif
