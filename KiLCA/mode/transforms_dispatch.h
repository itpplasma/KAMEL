/*! \file
    \brief Declarations for the Fortran field/coordinate transforms
    (kilca_transforms_m, KiLCA/mode/transforms_m.f90), replacing
    transforms.h now that flre_zone (the other caller) is Fortran too.
    bind(C) under the oracle's original (unmangled) names.
*/

#ifndef TRANSFORMS_DISPATCH_INCLUDE
#define TRANSFORMS_DISPATCH_INCLUDE

#include <complex>

#include "background.h"

extern "C"
{
void galilean_transform_of_EB_fields (std::complex<double> *EB1, std::complex<double> *EB2, double V);
void transform_EB_from_cyl_to_rsp (const background *bp, double r, std::complex<double> *EBcyl, std::complex<double> *EBrsp);
void transform_EB_from_rsp_to_cyl (const background *bp, double r, std::complex<double> *EBrsp, std::complex<double> *EBcyl);
}

#endif
