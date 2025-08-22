!<The module is used for Fortran - C++ interface

module core

use constants, only: pp;

integer(pp) :: cd_ptr; !>pointer to core_data object
integer(pp) :: sd_ptr; !>pointer to settings object
integer(pp) :: bp_ptr; !>pointer to background object

end module

!------------------------------------------------------------------------------

subroutine set_core_data_in_core_module (cd)

    use constants, only: pp
    use core, only: cd_ptr

    implicit none

    integer(pp), intent(in) :: cd

    cd_ptr = cd

end subroutine

!------------------------------------------------------------------------------

subroutine set_settings_in_core_module (sd)

use constants, only: pp;
use core, only: sd_ptr;

implicit none;

integer(pp), intent(in) :: sd

sd_ptr = sd;

end subroutine

!------------------------------------------------------------------------------

subroutine set_background_in_core_module (bp)

use constants, only: pp;
use core, only: bp_ptr;

implicit none;

integer(pp), intent(in) :: bp

bp_ptr = bp;

end subroutine

!------------------------------------------------------------------------------
