/*! \file
    \brief The implementation of functions declared in transforms.h.
*/

#include "constants.h"
#include "transforms.h"
#include "eval_back.h"

/*******************************************************************/

void galilean_transform_of_EB_fields (complex<double> *EB1, complex<double> *EB2, double V)
{
//transforms (E, B) fields in cylindrical coordinates
//from a frame 1 to a frame 2 moving with relative velocity V
//Landau II, p. 91
//E' = E + 1/c V x B,  B' = B - 1/c V x E
//E = E' - 1/c V x B', B = B' + 1/c V x E'
//V x A = e_r (-V A_th) + e_th (V A_r) + e_z (0);

EB2[0] = EB1[0] - (V/c)*EB1[4]; //Er
EB2[1] = EB1[1] + (V/c)*EB1[3]; //Etheta
EB2[2] = EB1[2];                //Ez

//for consistence, magnetic fields should not transform
EB2[3] = EB1[3]; // + (V/c)*EB1[1]; //Br
EB2[4] = EB1[4]; // - (V/c)*EB1[0]; //Btheta
EB2[5] = EB1[5];                    //Bz
}

/*******************************************************************/

void transform_EB_from_cyl_to_rsp (const background *bp, double r, complex<double> *EBcyl, complex<double> *EBrsp)
{
double htz[2];

eval_hthz (r, 0, 0, bp, htz);

//transform EB from cyl to rsp system:
EBrsp[0] = EBcyl[0];
EBrsp[1] = htz[1]*EBcyl[1] - htz[0]*EBcyl[2];
EBrsp[2] = htz[0]*EBcyl[1] + htz[1]*EBcyl[2];

EBrsp[3] = EBcyl[3];
EBrsp[4] = htz[1]*EBcyl[4] - htz[0]*EBcyl[5];
EBrsp[5] = htz[0]*EBcyl[4] + htz[1]*EBcyl[5];
}

/*****************************************************************************/

void transform_EB_from_rsp_to_cyl (const background *bp, double r, complex<double> *EBrsp, complex<double> *EBcyl)
{
double htz[2];

eval_hthz (r, 0, 0, bp, htz);

//transform EB from rsp to cyl system:
EBcyl[0] =   EBrsp[0];
EBcyl[1] =   htz[1]*EBrsp[1] + htz[0]*EBrsp[2];
EBcyl[2] = - htz[0]*EBrsp[1] + htz[1]*EBrsp[2];

EBcyl[3] =    EBrsp[3];
EBcyl[4] =    htz[1]*EBrsp[4] + htz[0]*EBrsp[5];
EBcyl[5] =  - htz[0]*EBrsp[4] + htz[1]*EBrsp[5];
}

/*****************************************************************************/
