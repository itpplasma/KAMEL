!> The basic driver program to run KiLCA code, ported from progs/main_linear.cpp.
!> Builds the project path (current directory or the first command-line
!> argument), then drives the now-Fortran core-data pipeline: create the core
!> data, register its handle in the core module, compute the mode-independent
!> data, then dispatch to the antenna or eigenmode mode-dependent path.
program main_linear
    use kilca_legacy_interfaces_m, only: set_core_data_in_core_module
    use kilca_antenna_settings_m, only: get_antenna_flag_eigmode
    use, intrinsic :: iso_c_binding, only: c_intptr_t
    use kilca_core_data_m, only: core_data_create_, core_data_destroy_, &
        core_data_calc_and_set_mode_independent_, &
        core_data_calc_and_set_mode_dependent_antenna_, &
        core_data_calc_and_set_mode_dependent_eigmode_
    use kilca_progs_common_m, only: get_project_path, to_cstr
    implicit none
    integer(c_intptr_t) :: cd

    cd = core_data_create_(get_project_path())
    call set_core_data_in_core_module(cd)

    call core_data_calc_and_set_mode_independent_(cd)

    if (get_antenna_flag_eigmode() == 0) then
        call core_data_calc_and_set_mode_dependent_antenna_(cd)
    else
        call core_data_calc_and_set_mode_dependent_eigmode_(cd)
    end if

    call core_data_destroy_(cd)
end program main_linear
