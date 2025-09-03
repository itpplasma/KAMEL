subroutine calculate_total_toroidal_torque(time_index)
    use iso_fortran_env, only: dp => real64
    use baseparam_mod, only: pi
    use grid_mod, only: T_EM_phi_e, T_EM_phi_i, T_tot_phi_e, T_tot_phi_i
    use integration, only: simpson_nonequi
    use wave_code_data, only: r

    implicit none

    integer, intent(in) :: time_index
    integer :: n

    n = size(r)
    if (mod(n - 1, 2) /= 0) then
        n = n - 1 ! Ensure n is even for Simpson's rule
    end if

    call simpson_nonequi(T_tot_phi_e(time_index), r(1:n), r(1:n) * T_EM_phi_e(1:n))
    call simpson_nonequi(T_tot_phi_i(time_index), r(1:n), r(1:n) * T_EM_phi_i(1:n))

    T_tot_phi_e(time_index) = T_tot_phi_e(time_index) * 2.0_dp * pi
    T_tot_phi_i(time_index) = T_tot_phi_i(time_index) * 2.0_dp * pi

end subroutine calculate_total_toroidal_torque

subroutine write_total_toroidal_torque_to_file(time_index)
    use grid_mod, only: T_tot_phi_e, T_tot_phi_i
    use h5mod

    implicit none

    integer, intent(in) :: time_index

    call h5_init()
    call h5_open_rw(path2out, h5_id)

    h5overwrite = .true.

    h5_currentgrp = "/"//trim(h5_mode_groupname)//"/T_tot_phi_e"
    call h5_add_double_1(h5_id, trim(h5_currentgrp), T_tot_phi_e(1:time_index), &
                         lbound(T_tot_phi_e(1:time_index)), ubound(T_tot_phi_e(1:time_index)), &
                         comment="Total toroidal torque on electrons")
    h5_currentgrp = "/"//trim(h5_mode_groupname)//"/T_tot_phi_i"
    call h5_add_double_1(h5_id, trim(h5_currentgrp), T_tot_phi_i(1:time_index), &
                         lbound(T_tot_phi_i(1:time_index)), ubound(T_tot_phi_i(1:time_index)), &
                         comment="Total toroidal torque on ions")

    h5overwrite = .false.

    call h5_close(h5_id)
    call h5_deinit()

end subroutine write_total_toroidal_torque_to_file
