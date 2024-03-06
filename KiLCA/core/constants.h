/*! \file constants.h
    \brief The definitions of some constants used throughout the library.
*/

#ifndef CONSTANTS_INCLUDE

#define CONSTANTS_INCLUDE

#include <complex>

using namespace std;

const double pi = 3.141592653589793238462643383279502884197;
const double eul = 0.5772156649015328606065120900824024310422;

const double boltz = 1.60216428e-12; /* erg/eV */

const double c = 29979245800.0;

const double m_p = 1.67262158e-24; /*proton mass*/
const double m_e = 9.10938185917485e-28; /*electron mass*/

const double e = 4.8032e-10; /*elementary charge*/

const double GAMMA = 5.0/3.0; /*adiabatic constant*/

const complex<double> O(0.0e0, 0.0e0), E(1.0e0, 0.0e0), I(0.0e0, 1.0e0);

const double sqrt2pi = sqrt(2.0e0*pi);

#endif
