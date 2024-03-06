/*! \file
    \brief The declarations of functions for transformations of fields.
*/

#ifndef TRANSFORMS_INCLUDE

#define TRANSFORMS_INCLUDE

#include <complex>

#include "background.h"

using namespace std;

/*******************************************************************/

void galilean_transform_of_EB_fields (complex<double> *EB1, complex<double> *EB2, double V);

void transform_EB_from_cyl_to_rsp (const background *bp, double r, complex<double> *EBcyl, complex<double> *EBrsp);

void transform_EB_from_rsp_to_cyl (const background *bp, double r, complex<double> *EBrsp, complex<double> *EBcyl);

/*******************************************************************/

#endif
