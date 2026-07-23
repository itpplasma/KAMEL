program test_susceptibility_collisionless_limit

    use getIfunc_config_m, only: boole_energy_conservation
    use KIM_kinds_m, only: dp
    use species_m, only: evaluate_susceptibility
    use use_libcerf_m, only: w_of_z_F

    implicit none

    real(dp), parameter :: pi = 4.0_dp * atan(1.0_dp)
    real(dp), parameter :: z_values(3) = [0.7_dp, 5.0_dp, 10.0_dp]
    real(dp), parameter :: scales(5) = &
        [1.0_dp, 10.0_dp, 100.0_dp, 1000.0_dp, 10000.0_dp]
    complex(dp) :: numerical(0:3, 0:3), exact(0:3, 0:3)
    real(dp) :: x1, x2, z, large_argument_error
    integer :: i, conservation, iz

    do conservation = 1, 2
        boole_energy_conservation = conservation == 1
        large_argument_error = 0.0_dp
        print '(A,L1)', 'energy conservation = ', boole_energy_conservation
        do iz = 1, size(z_values)
            z = z_values(iz)
            print '(A,F5.1)', 'z = ', z
            print '(A)', 'scale  err(I00)  err(I10)  err(I21)  err(I11)  err(I31)'
            do i = 1, size(scales)
                x1 = scales(i)
                x2 = sqrt(2.0_dp) * z * x1
                call evaluate_susceptibility(x1, x2, numerical)
                call collisionless_limit(x1, x2, exact)
                print '(F7.0,5(1X,ES10.3))', scales(i), &
                    relative_error(numerical(0, 0), exact(0, 0)), &
                    relative_error(numerical(1, 0), exact(1, 0)), &
                    relative_error(numerical(2, 1), exact(2, 1)), &
                    relative_error(numerical(1, 1), exact(1, 1)), &
                    relative_error(numerical(3, 1), exact(3, 1))
                if (scales(i) == 10000.0_dp) then
                    large_argument_error = max(large_argument_error, &
                        relative_error(numerical(0, 0), exact(0, 0)), &
                        relative_error(numerical(1, 0), exact(1, 0)), &
                        relative_error(numerical(2, 1), exact(2, 1)), &
                        relative_error(numerical(1, 1), exact(1, 1)), &
                        relative_error(numerical(3, 1), exact(3, 1)))
                end if
            end do
        end do
        if (large_argument_error > 1.0e-10_dp) then
            print *, 'FAIL: large-argument susceptibility error = ', &
                large_argument_error
            error stop
        end if
    end do

    x1 = -10000.0_dp
    do iz = 1, size(z_values)
        x2 = sqrt(2.0_dp) * z_values(iz) * abs(x1)
        call evaluate_susceptibility(x1, x2, numerical)
        call collisionless_limit(x1, x2, exact)
        large_argument_error = max( &
            relative_error(numerical(0, 0), exact(0, 0)), &
            relative_error(numerical(1, 0), exact(1, 0)), &
            relative_error(numerical(2, 1), exact(2, 1)), &
            relative_error(numerical(1, 1), exact(1, 1)), &
            relative_error(numerical(3, 1), exact(3, 1)))
        if (large_argument_error > 1.0e-10_dp) then
            print *, 'FAIL: negative-x1 susceptibility error = ', &
                large_argument_error
            error stop
        end if
    end do

contains

    subroutine collisionless_limit(x1, x2, values)
        real(dp), intent(in) :: x1, x2
        complex(dp), intent(out) :: values(0:3, 0:3)
        real(dp) :: sigma, zeta
        complex(dp) :: plasma_z, common

        sigma = sign(1.0_dp, x1)
        zeta = x2 / (sqrt(2.0_dp) * x1)
        plasma_z = cmplx(0.0_dp, 1.0_dp, dp) * sqrt(pi) &
            * w_of_z_F(cmplx(sigma * zeta, 0.0_dp, dp))
        common = 1.0_dp + sigma * zeta * plasma_z

        values = cmplx(0.0_dp, 0.0_dp, dp)
        values(0, 0) = -cmplx(0.0_dp, 1.0_dp, dp) * plasma_z &
            / (sqrt(2.0_dp) * abs(x1))
        values(1, 0) = -cmplx(0.0_dp, 1.0_dp, dp) * common / x1
        values(2, 1) = -cmplx(0.0_dp, 1.0_dp, dp) &
            * (1.0_dp + 2.0_dp * zeta**2 * common) / x1
        values(1, 1) = sqrt(2.0_dp) * zeta * values(1, 0)
        values(3, 1) = sqrt(2.0_dp) * zeta * values(2, 1)
    end subroutine collisionless_limit

    real(dp) function relative_error(actual, expected) result(error)
        complex(dp), intent(in) :: actual, expected

        error = abs(actual - expected) / max(abs(expected), tiny(1.0_dp))
    end function relative_error

end program test_susceptibility_collisionless_limit
