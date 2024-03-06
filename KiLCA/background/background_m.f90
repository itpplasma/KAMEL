!<The module for some of the background data

module background

use constants, only: dp, dpc;

implicit none;

!Machine settings:
real (dp) :: rtor;  !big torus radius (cm) of the machine
real (dp) :: rp;    !plasma radius (cm)
real (dp) :: B0;    !toroidal magnetic field (G) at the center

character (1) :: flag_back; !flag for background

real (dp) :: V_gal_sys; !velocity (cm/c) of a moving frame
real (dp) :: V_scale;   !scale of the Vz velocity profile: Vz = V_scale*Vz - V_gal_sys

real (dp) :: zele;  !collision coefficient for electrons
real (dp) :: zion;  !collision coefficient for ions

integer :: flag_debug;  !flag for debugging mode (additional checks are performed in the

!Particles settings:
real (dp), dimension (0:1) :: mass;
real (dp), dimension (0:1) :: charge;

!Other misc parameters:
real (dp) :: huge_factor;

end module

!------------------------------------------------------------------------------

subroutine copy_background_data_to_background_module (bp)

use constants, only: pp;
use background;

implicit none;

integer(pp), intent(in) :: bp;

call set_background_settings_c (bp, rtor, rp, B0, flag_back, V_gal_sys, V_scale, zele, zion, flag_debug);

call set_particles_settings_c (bp, mass, charge);

call set_huge_factor_c (bp, huge_factor);

if (flag_debug > 0) then
    call print_background_module_data ();
end if

end subroutine

!--------------------------------------------------------------------

subroutine print_background_module_data ()

use background;

implicit none;

print *;

!Machine vessel settings:
print *, "torus big radius: ", rtor;
print *, "plasma radius: ", rp;
print *, "toroidal magnetic field at the center: ", B0;

!Backround field and plasma settings:
print *, "flag for background: ", flag_back;
print *, "velocity of the moving frame: ", V_gal_sys;
print *, "scale factor for the Vz velocity profile: ", V_scale;
print *, "collision coefficient for electrons: ", zele;
print *, "collision coefficient for ions: ", zion;
print *, "flag for debugging mode: ", flag_debug;

!Particles settings:
print *, "particle masses: ", mass;
print *, "particle charges: ", charge;

print *, "huge_factor: ", huge_factor;

end subroutine

!--------------------------------------------------------------------
