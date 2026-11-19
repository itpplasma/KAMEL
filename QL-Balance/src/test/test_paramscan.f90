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

    ! Interior stencil, and a stencil clamped to the last rb point.
    call check_interpolation_at_resonance(5.3_dp)
    call check_interpolation_at_resonance(9.7_dp**1.2_dp)

contains

    !> Linear profiles on rb must be reproduced exactly at r_resonant. Using
    !> the rc average of Ercov with the rb coefficients shifts Er_res by half
    !> a cell, far outside the tolerance.
    subroutine check_interpolation_at_resonance(r_res)
        use grid_mod, only: npoib, rb, r_resonant, dqle22, Ercov
        use paramscan_mod, only: Er_res, br_abs_res_parscan, dqle22_res, &
            interpolate_Br_Dql_at_res_parscan
        use wave_code_data, only: Br, antenna_factor

        real(dp), intent(in) :: r_res
        integer :: i

        npoib = 10
        if (allocated(rb)) deallocate(rb, Ercov, dqle22, Br, r_resonant, &
            Er_res, br_abs_res_parscan, dqle22_res)
        ! Nonuniform, so rb and its midpoints rc differ cell by cell.
        allocate(rb(npoib))
        rb = [(real(i, dp)**1.2_dp, i = 1, npoib)]
        allocate(Ercov(npoib), dqle22(npoib), Br(npoib), r_resonant(1))
        Ercov = 2.0_dp + 0.5_dp * rb
        dqle22 = 1.0e3_dp * (1.0_dp + rb)
        Br = (3.0_dp, 4.0_dp) * rb
        antenna_factor = 4.0_dp
        r_resonant(1) = r_res

        ! fac_vz has more than one entry, so Er_res is computed.
        allocate(Er_res(1, 1, 1, size(fac_vz)), br_abs_res_parscan(1, 1, 1, size(fac_vz)), &
            dqle22_res(1, 1, 1, size(fac_vz)))
        ifac_vz = 2

        call interpolate_Br_Dql_at_res_parscan

        call assert_close([Er_res(1, 1, 1, 2)], [2.0_dp + 0.5_dp * r_res], &
            "Er_res of a linear Ercov")
        call assert_close([dqle22_res(1, 1, 1, 2)], [1.0e3_dp * (1.0_dp + r_res)], &
            "dqle22_res of a linear dqle22")
        call assert_close([br_abs_res_parscan(1, 1, 1, 2)], [2.0_dp * 5.0_dp * r_res], &
            "|Br| at resonance times sqrt(antenna_factor)")
    end subroutine check_interpolation_at_resonance

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
