program test_local_dispersion_geometry

    use KIM_kinds_m, only: dp
    use local_dispersion_geometry_m, only: local_bessel_argument, &
                                           local_laplacian_wavenumber_squared

    implicit none

    complex(dp), parameter :: kr_squared = cmplx(2.0_dp, 3.0_dp, dp)
    real(dp), parameter :: ks = 4.0_dp
    real(dp), parameter :: kparallel = 5.0_dp
    real(dp), parameter :: rho_t = 0.2_dp
    real(dp), parameter :: tolerance = 100.0_dp * epsilon(1.0_dp)

    call assert_close('Bessel argument', &
                      local_bessel_argument(kr_squared, ks, rho_t), &
                      cmplx(0.72_dp, 0.12_dp, dp))
    call assert_close('Laplacian wavenumber', &
                      local_laplacian_wavenumber_squared(kr_squared, ks, kparallel), &
                      cmplx(43.0_dp, 3.0_dp, dp))

    print *, 'KIM local dispersion geometry OK'

contains

    subroutine assert_close(label, actual, expected)
        character(*), intent(in) :: label
        complex(dp), intent(in) :: actual, expected

        if (abs(actual - expected) > tolerance * max(1.0_dp, abs(expected))) then
            print *, trim(label), ': expected ', expected, ', got ', actual
            error stop 1
        end if
    end subroutine assert_close

end program test_local_dispersion_geometry
