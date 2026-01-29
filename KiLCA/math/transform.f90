!==============================================================================
! SUBROUTINE cyl2rsp
!==============================================================================
! Transforms from (theta,z) to (perpendicular, parallel).
!------------------------------------------------------------------------------
!
! input:
!
!    r  ... radial position
!    th ... theta component
!    z  ... z component
!
! output:
!
!    per ... perpendicular component (s)
!    par ... parallel component (p)
!
!==============================================================================

subroutine cyl2rsp (r, th, z, per, par)

use constants,  only: dpc
use conduct_parameters, only: ht_, hz_
use background, only: flag_back

implicit none;

complex(dpc), intent(in)  :: th, z
complex(dpc), intent(out) :: per, par
real(dpc),    intent(in)  :: r

external :: eval_and_set_background_parameters_spec_independent

call eval_and_set_background_parameters_spec_independent (r, flag_back)

per = hz_*th - ht_*z
par = ht_*th + hz_*z

end subroutine

!==============================================================================
! subroutine rsp2cyl
!==============================================================================
! transforms from (perpendicular, parallel) to (theta,z).
!------------------------------------------------------------------------------
!
! input:
!
! input:
!
!    r   ... radial position
!    per ... perpendicular component
!    par ... parallel component
!
! output:
!
!    th ... theta component
!    z  ... z component
!==============================================================================

subroutine rsp2cyl (r, per, par, th, z)

use constants, only: dpc
use conduct_parameters, only: ht_, hz_
use background, only: flag_back

implicit none;

real(dpc),    intent(in)  :: r
complex(dpc), intent(in)  :: per, par
complex(dpc), intent(out) :: th, z

external :: eval_and_set_background_parameters_spec_independent

call eval_and_set_background_parameters_spec_independent (r, flag_back)

th = hz_*per + ht_*par
z  = hz_*par - ht_*per

end subroutine
