program test_toroidal_torque
    use iso_fortran_env, only: dp => real64
    use baseparam_mod, only: pi, rtor
    use grid_mod, only: T_EM_phi_e, T_EM_phi_i, T_tot_phi_e, T_tot_phi_i
    use wave_code_data, only: r
    use h5mod, only: path2out, h5_mode_groupname
    use KAMEL_hdf5_tools, only: HID_T, h5_add, h5_close, h5_create, h5_deinit, &
                                h5_get, h5_open, h5_get_bounds_1

    implicit none

    real(dp), parameter :: tolerance = 1.0e-12_dp
    real(dp) :: expected_e, expected_i
    external :: calculate_total_toroidal_torque, write_total_toroidal_torque_to_file
    integer(HID_T) :: file_id
    integer :: lo, hi

    allocate (r(3))
    allocate (T_EM_phi_e(3), T_EM_phi_i(3))
    allocate (T_tot_phi_e(3), T_tot_phi_i(3))

    r = [0.0_dp, 0.5_dp, 1.0_dp]
    T_EM_phi_e = 1.0_dp
    T_EM_phi_i = -2.0_dp
    rtor = 3.0_dp

    call calculate_total_toroidal_torque(1)

    expected_e = 4.0_dp * pi**2 * rtor * 0.5_dp
    expected_i = -2.0_dp * expected_e

    call assert_close(T_tot_phi_e(1), expected_e, "electron total torque")
    call assert_close(T_tot_phi_i(1), expected_i, "ion total torque")

    ! Seed prior output independently: the first writer call models a fresh
    ! process continuing at step three, with no live unlimited-dataset handles.
    path2out = 'test_toroidal_torque_restart.h5'
    h5_mode_groupname = 'multi_mode'
    call h5_create(path2out, file_id)
    call h5_add(file_id, '/multi_mode/T_tot_phi_e', [11.0_dp, 12.0_dp], [1], [2])
    call h5_add(file_id, '/multi_mode/T_tot_phi_i', [-21.0_dp, -22.0_dp], [1], [2])
    call h5_close(file_id)
    call h5_deinit()
    T_tot_phi_e(3) = 13.0_dp
    T_tot_phi_i(3) = -23.0_dp
    call write_total_toroidal_torque_to_file(3)
    call verify_history(path2out, 'multi_mode', [11.0_dp, 12.0_dp, 13.0_dp], &
                        [-21.0_dp, -22.0_dp, -23.0_dp])

    ! A new file and another group must never inherit old handles or values.
    path2out = 'test_toroidal_torque_other.h5'
    call h5_create(path2out, file_id)
    call h5_close(file_id)
    call h5_deinit()
    T_tot_phi_e(1:2) = [31.0_dp, 32.0_dp]
    T_tot_phi_i(1:2) = [-41.0_dp, -42.0_dp]
    call write_total_toroidal_torque_to_file(1)
    call write_total_toroidal_torque_to_file(2)
    call verify_history(path2out, 'multi_mode', [31.0_dp, 32.0_dp], [-41.0_dp, -42.0_dp])
    h5_mode_groupname = 'single_mode'
    call write_total_toroidal_torque_to_file(1)
    call verify_history(path2out, 'single_mode', [31.0_dp], [-41.0_dp])
    call verify_history(path2out, 'multi_mode', [31.0_dp, 32.0_dp], [-41.0_dp, -42.0_dp])

    ! Reopen the first output and replace one entry while retaining both sides.
    path2out = 'test_toroidal_torque_restart.h5'
    h5_mode_groupname = 'multi_mode'
    call write_total_toroidal_torque_to_file(2)
    call verify_history(path2out, 'multi_mode', [11.0_dp, 32.0_dp, 13.0_dp], &
                        [-21.0_dp, -42.0_dp, -23.0_dp])
    print *, 'toroidal torque integration and restart output tests passed'

contains

    subroutine verify_history(path, group, expected_e, expected_i)
        character(*), intent(in) :: path, group
        real(dp), intent(in) :: expected_e(:), expected_i(:)
        real(dp), allocatable :: actual(:)
        integer :: index

        call h5_open(path, file_id)
        call h5_get_bounds_1(file_id, '/'//group//'/T_tot_phi_e', lo, hi)
        if (lo /= 1 .or. hi /= size(expected_e)) error stop 'electron history extent mismatch'
        allocate (actual(hi))
        call h5_get(file_id, '/'//group//'/T_tot_phi_e', actual)
        do index = 1, hi
            call assert_close(actual(index), expected_e(index), 'electron torque history')
        end do
        call h5_get_bounds_1(file_id, '/'//group//'/T_tot_phi_i', lo, hi)
        if (lo /= 1 .or. hi /= size(expected_i)) error stop 'ion history extent mismatch'
        call h5_get(file_id, '/'//group//'/T_tot_phi_i', actual)
        do index = 1, hi
            call assert_close(actual(index), expected_i(index), 'ion torque history')
        end do
        call h5_close(file_id)
        call h5_deinit()
    end subroutine verify_history

    subroutine assert_close(actual, expected, label)
        real(dp), intent(in) :: actual, expected
        character(len=*), intent(in) :: label

        if (abs(actual - expected) > tolerance * max(1.0_dp, abs(expected))) then
            print '(A,A,A,ES24.16,A,ES24.16)', "FAIL: ", label, &
                " got ", actual, " expected ", expected
            stop 1
        end if
    end subroutine assert_close

end program test_toroidal_torque
