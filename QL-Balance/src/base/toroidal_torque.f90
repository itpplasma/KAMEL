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

    T_tot_phi_e(time_index) = T_tot_phi_e(time_index) * 4.0_dp * pi**2 * rtor
    T_tot_phi_i(time_index) = T_tot_phi_i(time_index) * 4.0_dp * pi**2 * rtor

end subroutine calculate_total_toroidal_torque

subroutine write_total_toroidal_torque_to_file(time_index)
    use iso_fortran_env, only: dp => real64
    use grid_mod, only: T_tot_phi_e, T_tot_phi_i
    use KAMEL_hdf5_tools, only: h5_init, h5_deinit, h5_open_rw, h5_close, &
                                h5_obj_exists, h5_get_bounds_1, h5_get, h5_delete, h5_add
    use h5mod, only: h5_id, path2out, h5_mode_groupname

    implicit none

    integer, intent(in) :: time_index

    if (time_index < 1) error stop 'torque output time index must be positive'
    if (time_index > size(T_tot_phi_e) .or. time_index > size(T_tot_phi_i)) &
        error stop 'torque output time index exceeds storage'
    call h5_init()
    call h5_open_rw(path2out, h5_id)
    call write_history('/'//trim(h5_mode_groupname)//'/T_tot_phi_e', T_tot_phi_e(time_index))
    call write_history('/'//trim(h5_mode_groupname)//'/T_tot_phi_i', T_tot_phi_i(time_index))
    call h5_close(h5_id)
    call h5_deinit()

contains

    subroutine write_history(path, value)
        character(*), intent(in) :: path
        real(dp), intent(in) :: value
        real(dp), allocatable :: previous(:), history(:)
        integer :: lo, hi
        logical :: exists

        ! Unlimited-array handles are process-local buffers in Fortio. Resolve
        ! history by path on every write so reopening or restarting preserves it.
        call h5_obj_exists(h5_id, path, exists)
        hi = 0
        if (exists) then
            call h5_get_bounds_1(h5_id, path, lo, hi)
            if (lo /= 1 .or. hi < lo) error stop 'invalid torque history bounds'
            allocate (previous(hi))
            call h5_get(h5_id, path, previous)
        end if
        allocate (history(max(hi, time_index)), source=0.0_dp)
        if (exists) history(:hi) = previous
        history(time_index) = value
        if (exists) call h5_delete(h5_id, path)
        call h5_add(h5_id, path, history, [1], [size(history)])
    end subroutine write_history

end subroutine write_total_toroidal_torque_to_file
