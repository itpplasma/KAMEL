module antenna_data

!Contains antenna data and current density coefficients Jth(m,n), Jz(m,n) for the current mode

use constants, only: dp, dpc

real (dp) :: ra;       !small radius (cm) of antenna location
real (dp) :: wa;       !current density layer width
real (dp) :: I0;       !current in antenna coils (statamp)
complex (dpc) :: flab; !frequency (Hz) in the laboratory frame
integer :: dma;        !dimension of modes array
integer, dimension(:), allocatable :: modes; !array of modes (m,n)
logical :: flag_debug; !debug flag

complex(dpc), dimension(2) :: Ja_cyl

end module

!------------------------------------------------------------------------------

subroutine copy_antenna_data_to_antenna_module (as)

use constants, only: pp, dp, im;
use antenna_data;

integer(pp), intent(in) :: as
real(dp) :: flab_re, flab_im;

call set_antenna_settings_c (as, ra, wa, I0, flab_re, flab_im, dma, flag_debug);

flab = flab_re + im*flab_im;

!array 'modes' is not set!

if (flag_debug) call print_antenna_module_data ();

end subroutine

!------------------------------------------------------------------------------

subroutine set_current_density_in_antenna_module ()

use antenna_data, only: Ja_cyl;

implicit none;

call current_density (Ja_cyl); !Fourier amplutudes of the (Jth, Jz)

end subroutine

!------------------------------------------------------------------------------

subroutine calc_current_density_r_s_p (r, Ja_rsp)

use constants, only: dp, dpc, isqrt2pi
use antenna_data, only: ra, wa, Ja_cyl;

real(dp), intent(in) :: r
complex(dpc), dimension(2), intent(out) :: Ja_rsp

real(dp) :: delta

!(r,s,p) antenna current density Fourier amplutudes: optimize!
call cyl2rsp (r, Ja_cyl(1), Ja_cyl(2), Ja_rsp(1), Ja_rsp(2));

!delta-function approximation:
delta = (isqrt2pi/wa)*exp(-((r-ra)**2)/(2.0*wa**2));

Ja_rsp(1) = Ja_rsp(1)*delta;
Ja_rsp(2) = Ja_rsp(2)*delta;

end subroutine

!------------------------------------------------------------------------------

subroutine print_antenna_module_data ()

use antenna_data;

implicit none;

print *

!Antenna settings:
print *, "antenna radius: ", ra;
print *, "antenna layer width: ", wa;
print *, "antenna coils current: ", I0;
print *, "antenna lab frequency: ", flab;
print *, "dimension of modes array: ", dma;
print *, "flag for debugging mode: ", flag_debug;

end subroutine

!--------------------------------------------------------------------
