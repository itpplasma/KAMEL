#ifndef BACK_SETT_INCLUDE
#define BACK_SETT_INCLUDE

#include "constants.h"

#include <cmath>

struct back_sett {
  /// Machine settings:
  double rtor; // big torus radius (cm) of the machine
  double rp;   // plasma radius (cm)
  double B0;   // toroidal magnetic field (G) at the center

  /// Backround field and plasma settings:
  char path2profiles[256]; // path to input background profiles
  int calc_back; // the sign of the flag shows whether profiles must be loaded
                 // from files or the interface, the value shows how profiles
                 // have to be recalculated
  char flag_back[8]; // flag for background: normal (full) or homogenious
  int N; // splines degree:  >= NC + 2N+1, where N - order of flr expansion,
         // NC - spl degree for C matrices, must be odd
  double V_gal_sys; // velocity (cm/c) of a moving frame
  double V_scale;   // scale of the Vz velocity profile: Vz = V_scale*Vz -
                  // V_gal_sys
  double m_i;  // ions mass in units of proton mass
  double zele; // collision coefficient for electrons: = 1.0 for realistic
               // collision frequency
  double zion; // collision coefficient for ions: = 1.0 for realistic
               // collision frequency
  bool flag_debug; // flag for debugging mode (additional checks are performed)

  /// Particles settings:
  double mass[2] = {NAN, m_e};   //(ions, electrons) masses
  double charge[2] = {e, -e}; //(ions, electrons) charges

  /// Other misc parameters:
  double huge_factor{1e20}; // big factor used in special cases

  back_sett() = default;

  void print_settings();
};

/*-----------------------------------------------------------------*/

extern "C" // C interface for Fortran
{
void set_background_settings_c_ (back_sett **isett, double *rtor, double *rp, double *B0, char *flag_back, double *V_gal_sys, double *V_scale, double *zele, double *zion, int *flag_debug);

void set_particles_settings_c_ (back_sett **isett, double *mass, double *charge);

void set_huge_factor_c_ (back_sett **isett, double *fac);
}

/*-----------------------------------------------------------------*/

#endif
