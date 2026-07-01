!> The basic driver program to run KiLCA code, ported from progs/main_linear.cpp.
!> Builds the project path (current directory or the first command-line
!> argument), then drives the now-Fortran core-data pipeline: create the core
!> data, register its handle in the core module, compute the mode-independent
!> data, then dispatch to the antenna or eigenmode mode-dependent path.
program main_linear
    use, intrinsic :: iso_c_binding, only: c_int, c_intptr_t, c_char
    use kilca_core_data_m, only: core_data_create_, core_data_destroy_, &
        core_data_calc_and_set_mode_independent_, &
        core_data_calc_and_set_mode_dependent_antenna_, &
        core_data_calc_and_set_mode_dependent_eigmode_
    use kilca_progs_common_m, only: get_project_path, to_cstr
    implicit none

    interface
        subroutine set_core_data_in_core_module(cd) &
            bind(C, name="set_core_data_in_core_module_")
            import :: c_intptr_t
            integer(c_intptr_t), intent(in) :: cd
        end subroutine set_core_data_in_core_module

        integer(c_int) function get_antenna_flag_eigmode() &
            bind(C, name="get_antenna_flag_eigmode_")
            import :: c_int
        end function get_antenna_flag_eigmode
    end interface

    integer(c_intptr_t) :: cd
    character(kind=c_char), allocatable :: cpath(:)

    cpath = to_cstr(get_project_path())

    cd = core_data_create_(cpath)
    call set_core_data_in_core_module(cd)

    call core_data_calc_and_set_mode_independent_(cd)

    if (get_antenna_flag_eigmode() == 0) then
        call core_data_calc_and_set_mode_dependent_antenna_(cd)
    else
        call core_data_calc_and_set_mode_dependent_eigmode_(cd)
    end if

    call core_data_destroy_(cd)
end program main_linear
