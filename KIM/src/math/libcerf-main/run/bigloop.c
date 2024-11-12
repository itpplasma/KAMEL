/* Library libcerf:
 *   Compute complex error functions, based on a new implementation of
 *   Faddeeva's w_of_z. Also provide Dawson and Voigt functions.
 *
 * File bigloop.c:
 *   Compute many function values, to measure timing.
 *
 * Copyright:
 *   (C) 2022 Forschungszentrum Jülich GmbH
 *
 * Licence:
 *   Public domain.
 *
 * Author:
 *   Joachim Wuttke, Forschungszentrum Jülich, 2022
 *
 * Website:
 *   http://apps.jcns.fz-juelich.de/libcerf
 */

/*
   Prevent temperature dependent variation of processor frequency.
   Under Windows:
       set SpeedStep and Turbo(boost).
   Under Linux:
       echo "1" > /sys/devices/system/cpu/intel_pstate/no_turbo
*/

#include "cerf.h"
#include <stdio.h>
#define size_t long int

int main()
{
    const size_t N = 1<<29; // number of function calls
    const size_t n = N>>7;  // number of sweeps through the x range
    double sum = 0;

    /* We run n times through the range x=0..1, each time with slightly different x values.
       This is intended to prevent the function under test from keeping code for nearby x values
       in the L1 cache.
       To scan one single time in tiny steps through the range x=0..1, just exchange the order
       of the following two for loops.
    */

    for (size_t j=0; j<n; ++j) {
        for (size_t i=0; i<N; i += n) {
            const double x = (i+j) * (1.0/N);
            sum += im_w_of_x(x*2.);
            sum += im_w_of_x(x*10.);
            sum += im_w_of_x(x*50.);
        }
    }
    return (int)sum;
}
