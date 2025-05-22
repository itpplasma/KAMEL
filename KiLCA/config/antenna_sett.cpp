#include "antenna_sett.h"

#include <cmath>
#include <cstring>
#include <stdio.h>
#include <stdlib.h>

/*-----------------------------------------------------------------*/

void antenna_sett::print_settings (void)
{
fprintf(stdout, "\nCheck for antenna parameters below:\n");

//Antenna settings:
fprintf (stdout, "\nantenna radius: %lg cm", ra);
fprintf (stdout, "\nantenna current layer width: %lg cm", wa);
fprintf (stdout, "\nantenna coils current: %lg statamps", I0);
fprintf (stdout, "\nantenna lab frequency: (%lg, %lg) 1/s", real(flab), imag(flab));
fprintf (stdout, "\ndimension of modes array: %d", dma);
fprintf (stdout, "\nflag for debugging mode: %d", flag_debug);
fprintf (stdout, "\nflag for eigmode search: %d", flag_eigmode);

fprintf (stdout, "\narray of mode numbers (m, n): ");

int i;
for (i=0; i<dma; i++)
{
    fprintf (stdout, "(%d, %d) ", modes[2*i], modes[2*i+1]);
}

fprintf(stdout, "\n");
}

/*-----------------------------------------------------------------*/

void set_antenna_settings_c_ (antenna_sett **ptr, double *ra, double *wa, double *I0, double *flab_re, double *flab_im, int *dma, bool *flag_debug)
{
/*! \fn set_antenna_settings_c_ (antenna **ptr, double *ra, double *wa, double *I0, double *flab_re, double *flab_im, int *dma, int *flag_debug)
    \brief The function copies class variables to the parameters which are in fact members of the antenna module. The parameter ptr is an address of an address of the antenna object passed.

    \param
*/

antenna_sett *ap = *ptr;

*ra = ap->ra;
*wa = ap->wa;
*I0 = ap->I0;
*flab_re = real(ap->flab);
*flab_im = imag(ap->flab);
*dma = ap->dma;
*flag_debug = ap->flag_debug;
}

/*-----------------------------------------------------------------*/
