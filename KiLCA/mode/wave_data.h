/*! \file
    \brief The declaration of wave_data class.
*/

#ifndef WAVE_DATA_INCLUDE

#define WAVE_DATA_INCLUDE

#include <complex>

/*-----------------------------------------------------------------*/

/*! \class wave_data
    \brief Class describes the properties of a perturbation mode.
*/
class wave_data
{
public:
    int m;                  //!<m harmonic number
    int n;                  //!<n harmonic number

    complex<double> olab;   //!<wave frequency in a lab frame
    complex<double> omov;   //!<wave frequency in a moving frame

    double r_res;           //!<resonance location for the given mode (if present)

    complex<double> det;    //!<determinant

public:
    wave_data (int m1, int n1, complex<double> olab1, complex<double> omov1)
    {
        m = m1;
        n = n1;
        olab = olab1;
        omov = omov1;
    }

    ~wave_data (void) {}
};

/*-----------------------------------------------------------------*/

#endif
