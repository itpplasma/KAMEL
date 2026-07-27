program test_paramscan
    use iso_fortran_env, only: dp => real64
    use grid_mod, only: npoic
    use paramscan_mod, only: fac_n, fac_Te, fac_Ti, fac_vz, &
        ifac_n, ifac_Te, ifac_Ti, ifac_vz, rescale_profiles
    use plasma_parameters, only: params, hold_n, hold_vz, hold_Te, hold_Ti, &
        alloc_hold_parameters
    implicit none

    real(dp), parameter :: tolerance = 1.0e-13_dp
    real(dp) :: rotation_term(3, 3)
    real(dp) :: rotation_coefficient(3)
    integer :: scan_index

    npoic = 3

    allocate(params(4, npoic))
    params(1, :) = [1.0_dp, 2.0_dp, 3.0_dp]
    params(2, :) = [10.0_dp, 20.0_dp, 30.0_dp]
    params(3, :) = [100.0_dp, 200.0_dp, 300.0_dp]
    params(4, :) = [400.0_dp, 500.0_dp, 600.0_dp]

    call alloc_hold_parameters

    allocate(fac_n(1), fac_Te(1), fac_Ti(1), fac_vz(3))
    fac_n = 1.0_dp
    fac_Te = 1.0_dp
    fac_Ti = 1.0_dp
    fac_vz = [1.0_dp, 2.0_dp, 3.0_dp]
    ifac_n = 1
    ifac_Te = 1
    ifac_Ti = 1
    rotation_coefficient = [0.5_dp, 1.0_dp, 1.5_dp]

    do scan_index = 1, size(fac_vz)
        ifac_vz = scan_index
        call rescale_profiles
        call assert_close(params(2, :), hold_vz * fac_vz(scan_index), &
            "rotation profile scaling")
        rotation_term(:, scan_index) = rotation_coefficient * params(2, :)
    end do

    call assert_close(rotation_term(:, 3) - 2.0_dp * rotation_term(:, 2) &
        + rotation_term(:, 1), [0.0_dp, 0.0_dp, 0.0_dp], &
        "linear force-balance rotation contribution")
    call assert_close(params(1, :), hold_n, "density unchanged")
    call assert_close(params(3, :), hold_Te, "electron temperature unchanged")
    call assert_close(params(4, :), hold_Ti, "ion temperature unchanged")

contains

    subroutine assert_close(actual, expected, label)
        real(dp), intent(in) :: actual(:), expected(:)
        character(len=*), intent(in) :: label

        if (any(abs(actual - expected) > tolerance * max(1.0_dp, abs(expected)))) then
            print '(A,A)', "FAIL: ", label
            print '(A,*(ES15.7,1X))', "  actual:   ", actual
            print '(A,*(ES15.7,1X))', "  expected: ", expected
            error stop 1
        end if
    end subroutine assert_close

end program test_paramscan
