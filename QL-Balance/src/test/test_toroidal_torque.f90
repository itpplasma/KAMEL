program test_toroidal_torque
    use iso_fortran_env, only: dp => real64
    use baseparam_mod, only: pi, rtor
    use grid_mod, only: T_EM_phi_e, T_EM_phi_i, T_tot_phi_e, T_tot_phi_i
    use wave_code_data, only: r

    implicit none

    real(dp), parameter :: tolerance = 1.0e-12_dp
    real(dp) :: expected_e, expected_i
    external :: calculate_total_toroidal_torque

    allocate(r(3))
    allocate(T_EM_phi_e(3), T_EM_phi_i(3))
    allocate(T_tot_phi_e(1), T_tot_phi_i(1))

    r = [0.0_dp, 0.5_dp, 1.0_dp]
    T_EM_phi_e = 1.0_dp
    T_EM_phi_i = -2.0_dp
    rtor = 3.0_dp

    call calculate_total_toroidal_torque(1)

    expected_e = 4.0_dp * pi**2 * rtor * 0.5_dp
    expected_i = -2.0_dp * expected_e

    call assert_close(T_tot_phi_e(1), expected_e, "electron total torque")
    call assert_close(T_tot_phi_i(1), expected_i, "ion total torque")

contains

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
