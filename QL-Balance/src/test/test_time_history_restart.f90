program test_time_history_restart
    use QLBalance_kinds, only: dp
    use time_evolution, only: write_time_info_to_h5, time_ind, timstep, time, &
                              timscal, firstiterationdone
    use control_mod, only: data_verbosity
    use QLbalance_diag, only: rate_dql, timscal_dql
    use h5mod, only: path2out, h5_mode_groupname
    use KAMEL_hdf5_tools, only: HID_T, h5_create, h5_open, h5_close, h5_deinit, &
                                h5_add, h5_get, h5_get_bounds
    implicit none
    integer(HID_T) :: file_id
    real(dp) :: expected(3, 3), detailed(6, 1)

    allocate (timscal(1))
    h5_mode_groupname = 'multi_mode'
    path2out = 'test_time_history_restart.h5'
    expected(:, 1) = [1.0_dp, 0.1_dp, 0.1_dp]
    expected(:, 2) = [2.0_dp, 0.2_dp, 0.3_dp]
    expected(:, 3) = [3.0_dp, 0.4_dp, 0.7_dp]
    call h5_create(path2out, file_id)
    call h5_add(file_id, '/multi_mode/timstep_evol.dat', expected(:, :2), [1, 1], [3, 2])
    call h5_close(file_id)
    call h5_deinit()

    data_verbosity = 1
    firstiterationdone = .false.
    time_ind = 3
    timstep = 0.4_dp
    time = 0.7_dp
    call write_time_info_to_h5()
    call verify(expected)

    ! The process may subsequently write a different file with a different
    ! verbosity; neither path may inherit a stale matrix handle.
    firstiterationdone = .true.
    path2out = 'test_time_history_detailed.h5'
    call h5_create(path2out, file_id)
    call h5_close(file_id)
    call h5_deinit()
    data_verbosity = 2
    time_ind = 1
    timscal_dql = 0.15_dp
    timscal(1) = 0.25_dp
    rate_dql = 0.35_dp
    detailed(:, 1) = [1.0_dp, timstep, timscal_dql, timscal(1), rate_dql, time]
    call write_time_info_to_h5()
    call verify(detailed)

    path2out = 'test_time_history_restart.h5'
    data_verbosity = 1
    time_ind = 2
    expected(:, 2) = [2.0_dp, timstep, time]
    call write_time_info_to_h5()
    call verify(expected)
    print *, 'time history restart tests passed'
contains
    subroutine verify(wanted)
        real(dp), intent(in) :: wanted(:, :)
        real(dp), allocatable :: actual(:, :)
        integer :: lo1, lo2, hi1, hi2
        call h5_open(path2out, file_id)
        call h5_get_bounds(file_id, '/multi_mode/timstep_evol.dat', lo1, lo2, hi1, hi2)
        if (any([lo1, lo2] /= 1) .or. any([hi1, hi2] /= shape(wanted))) &
            error stop 'time history bounds changed'
        allocate (actual(hi1, hi2))
        call h5_get(file_id, '/multi_mode/timstep_evol.dat', actual)
        if (any(abs(actual - wanted) > 1.0e-12_dp)) error stop 'time history lost on reopen'
        call h5_close(file_id)
        call h5_deinit()
    end subroutine verify
end program test_time_history_restart
