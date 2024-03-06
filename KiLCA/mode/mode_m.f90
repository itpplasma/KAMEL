!<Describes data related to one mode of perturbation.

module mode_data

use constants, only: pp;

!Contains data for the current mode

!Wave settings:
integer :: m;         !m harmonic number
integer :: n;         !n harmonic number
complex(8) :: olab;   !wave frequency in a lab frame
complex(8) :: omov;   !wave frequency in a moving frame
real(8) :: r_res;     !resonance location for the given mode
complex(8) :: det;    !determinant

!pointers to common structures:
integer(pp) :: sd_ptr;
integer(pp) :: bp_ptr;

!Directories for a given mode:
character (1024) :: path2linear;       !path to linear-data
character (1024) :: path2dispersion;   !path to dispersion-data
character (1024) :: path2poincare;     !path to poincare-data

!pointers to structures for the given mode
integer(pp) :: wd_ptr;
integer(pp) :: cp_ptr;
integer(pp) :: sp_ptr;
integer(pp) :: ep_ptr;
integer(pp) :: fp_ptr;
integer(pp) :: qp_ptr;

end module

!------------------------------------------------------------------------------

subroutine set_settings_in_mode_data_module (sd)

use constants, only: pp;
use mode_data, only: sd_ptr;

implicit none;

integer(pp), intent(in) :: sd

sd_ptr = sd;

end subroutine

!------------------------------------------------------------------------------

subroutine set_back_profiles_in_mode_data_module (bp)

use constants, only: pp;
use mode_data, only: bp_ptr;

implicit none;

integer(pp), intent(in) :: bp

bp_ptr = bp;

end subroutine

!------------------------------------------------------------------------------

subroutine set_wave_data_in_mode_data_module (wd)

use constants, only: pp;
use mode_data, only: wd_ptr;

implicit none;

integer(pp), intent(in) :: wd

wd_ptr = wd

end subroutine

!------------------------------------------------------------------------------

subroutine set_cond_profiles_in_mode_data_module (cp)

use constants, only: pp;
use mode_data, only: cp_ptr;

implicit none;

integer(pp), intent(in) :: cp

cp_ptr = cp;

end subroutine

!------------------------------------------------------------------------------

subroutine set_sysmat_profiles_in_mode_data_module (sp)

use constants, only: pp;
use mode_data, only: sp_ptr;

implicit none;

integer(pp), intent(in) :: sp

sp_ptr = sp;

end subroutine

!------------------------------------------------------------------------------

subroutine set_eig_profiles_in_mode_data_module (ep)

use constants, only: pp;
use mode_data, only: ep_ptr;

implicit none;

integer(pp), intent(in) :: ep

ep_ptr = ep;

end subroutine

!------------------------------------------------------------------------------

subroutine set_field_profiles_in_mode_data_module (fp)

use constants, only: pp;
use mode_data, only: fp_ptr;

implicit none;

integer(pp), intent(in) :: fp

fp_ptr = fp;

end subroutine

!------------------------------------------------------------------------------

subroutine set_quants_profiles_in_mode_data_module (qp)

use constants, only: pp;
use mode_data, only: qp_ptr;

implicit none;

integer(pp), intent(in) :: qp

qp_ptr = qp;

end subroutine

!------------------------------------------------------------------------------

subroutine set_wave_parameters_in_mode_data_module (m_p, n_p, olab_re, olab_im, omov_re, omov_im)

use constants, only: dp, im;
use mode_data, only: m, n, olab, omov;

implicit none;

integer,  intent(in) :: m_p, n_p;
real(dp), intent(in) :: olab_re, olab_im, omov_re, omov_im;

m = m_p;
n = n_p;

olab = olab_re + im*olab_im;
omov = omov_re + im*omov_im;

end subroutine

!--------------------------------------------------------------------

subroutine set_resonance_location_in_mode_data_module (r_res_p)

use constants, only: dp;
use mode_data, only: r_res;

implicit none;

real (dp) :: r_res_p;

r_res = r_res_p

end subroutine

!--------------------------------------------------------------------

subroutine copy_mode_paths_to_mode_data_module (md)

use constants, only: pp;
use mode_data, only: path2linear, path2dispersion, path2poincare;

implicit none;

integer(pp), intent(in) :: md;

call copy_mode_paths_from_mode_data_struct (md, path2linear, path2dispersion, path2poincare);

end subroutine

!--------------------------------------------------------------------

subroutine clear_all_data_in_mode_data_module ()

use mode_data;

implicit none;

m = 0;
n = 0;
olab = 0.0d0;
omov = 0.0d0;
r_res = 0.0d0;

path2linear = "NULL";
path2dispersion = "NULL";
path2poincare = "NULL";

sd_ptr = 0;
bp_ptr = 0;
wd_ptr = 0;
cp_ptr = 0;
sp_ptr = 0;
ep_ptr = 0;
fp_ptr = 0;
qp_ptr = 0

end subroutine

!--------------------------------------------------------------------

subroutine save_mode_det_data ()

use mode_data;

implicit none;

integer :: stat

character(1024) :: fname

write (fname,'(A)') 'det.dat'

print *, trim(fname)

!opens or creates a file mode_det.dat
open(unit=10,iostat=stat,file=trim(fname),form='formatted')
if (stat /= 0) then
    print *, 'save_mode_det_data: file open error: stat=', stat
end if

endfile (10);

write (10,*) real(olab), aimag(olab), abs(det), real(det), aimag(det)
write (10,*)

close (10);

end subroutine

!--------------------------------------------------------------------
