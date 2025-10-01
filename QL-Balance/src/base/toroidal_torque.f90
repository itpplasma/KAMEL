subroutine calculate_total_toroidal_torque(time_index)
    use iso_fortran_env, only: dp => real64
    use baseparam_mod, only: pi, rtor ! rtor is the major radius
    use grid_mod, only: T_EM_phi_e, T_EM_phi_i, T_tot_phi_e, T_tot_phi_i
    use simpson_integration, only: simpson_nonequi
    use wave_code_data, only: r

    implicit none

    integer, intent(in) :: time_index

    call simpson_nonequi(T_tot_phi_e(time_index), r, r * T_EM_phi_e)
    call simpson_nonequi(T_tot_phi_i(time_index), r, r * T_EM_phi_i)

    T_tot_phi_e(time_index) = T_tot_phi_e(time_index) * 2.0_dp * pi**2 * rtor
    T_tot_phi_i(time_index) = T_tot_phi_i(time_index) * 2.0_dp * pi**2 * rtor

end subroutine calculate_total_toroidal_torque

subroutine write_total_toroidal_torque_to_file(time_index)
    use grid_mod, only: T_tot_phi_e, T_tot_phi_i
    use KAMEL_hdf5_tools, only: HID_T, H5T_NATIVE_DOUBLE, h5_init, h5_deinit, h5_open_rw, &
                                h5_close, h5_define_unlimited_array, h5_append_double_0
    use h5mod, only: h5_id, path2out, h5_mode_groupname, h5_currentgrp

    implicit none

    integer, intent(in) :: time_index
    integer(HID_T), save :: dsetid_e, dsetid_i ! save == static

    call h5_init()
    call h5_open_rw(path2out, h5_id)

    if (time_index == 1) then
        h5_currentgrp = "/"//trim(h5_mode_groupname)//"/T_tot_phi_e"
        call h5_define_unlimited_array(h5_id, trim(h5_currentgrp), H5T_NATIVE_DOUBLE, dsetid_e)
        h5_currentgrp = "/"//trim(h5_mode_groupname)//"/T_tot_phi_i"
        call h5_define_unlimited_array(h5_id, trim(h5_currentgrp), H5T_NATIVE_DOUBLE, dsetid_i)
    end if

    call h5_append_double_0(dsetid_e, T_tot_phi_e(time_index), time_index)
    call h5_append_double_0(dsetid_i, T_tot_phi_i(time_index), time_index)

    call h5_close(h5_id)
    call h5_deinit()

end subroutine write_total_toroidal_torque_to_file
