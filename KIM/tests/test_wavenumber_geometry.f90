program test_wavenumber_geometry

    use KIM_kinds_m, only: dp
    use constants_m, only: sol
    use wavenumber_geometry_m, only: parallel_wavenumber, perpendicular_wavenumber, &
                                    exb_rotation_frequency

    implicit none

    integer, parameter :: m_mode = -7
    integer, parameter :: n_mode = 3
    real(dp), parameter :: r = 40.0_dp
    real(dp), parameter :: R0 = 165.0_dp
    real(dp), parameter :: hth = 0.6_dp
    real(dp), parameter :: hz = -0.8_dp
    real(dp), parameter :: scale = 2.5_dp
    real(dp), parameter :: Er = -0.12_dp
    real(dp), parameter :: B0 = 2.1e4_dp
    real(dp), parameter :: tol = 100.0_dp * epsilon(1.0_dp)

    real(dp) :: kp, ks, expected_kp, expected_ks
    real(dp) :: scaled_ks, old_ks, omega_E

    expected_kp = hth * m_mode / r + hz * n_mode / R0
    expected_ks = hz * m_mode / r - hth * n_mode / R0

    kp = parallel_wavenumber(m_mode, n_mode, r, R0, hth, hz)
    ks = perpendicular_wavenumber(m_mode, n_mode, r, R0, hth, hz)

    call assert_close('parallel wavenumber', kp, expected_kp, tol)
    call assert_close('perpendicular wavenumber', ks, expected_ks, tol)
    call assert_close('orthogonal rotation identity', kp**2 + ks**2, &
                      (real(m_mode, dp) / r)**2 + (real(n_mode, dp) / R0)**2, tol)

    scaled_ks = perpendicular_wavenumber(m_mode, n_mode, scale * r, scale * R0, hth, hz)
    call assert_close('length scaling', scaled_ks, ks / scale, tol)

    omega_E = exb_rotation_frequency(ks, Er, B0)
    call assert_close('ExB rotation frequency', omega_E, -sol * Er * ks / B0, tol)

    ! Both hth and n_mode are nonzero, so the historical extra division by r
    ! on the toroidal term must be distinguishable from the correct expression.
    old_ks = (real(m_mode, dp) * hz - real(n_mode, dp) * hth / R0) / r
    if (abs(old_ks - expected_ks) <= 1.0e-6_dp) then
        print *, 'negative control did not distinguish the old ks expression'
        error stop 1
    end if

    print *, 'KIM wavenumber geometry OK'

contains

    subroutine assert_close(label, actual, expected, tolerance)
        character(*), intent(in) :: label
        real(dp), intent(in) :: actual, expected, tolerance

        if (abs(actual - expected) > tolerance * max(1.0_dp, abs(expected))) then
            print *, trim(label), ': expected ', expected, ', got ', actual
            error stop 2
        end if
    end subroutine assert_close

end program test_wavenumber_geometry
