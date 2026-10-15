program test_susceptibility_collisionless_limit

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
    integer :: i, conservation_model, iz

    do conservation_model = 0, 3
        large_argument_error = 0.0_dp
        print '(A,I0)', 'conservation model = ', conservation_model
        do iz = 1, size(z_values)
            z = z_values(iz)
            print '(A,F5.1)', 'z = ', z
            print '(A)', 'scale  err(I00)  err(I10)  err(I21)  err(I11)  err(I31)'
            do i = 1, size(scales)
                x1 = scales(i)
                x2 = sqrt(2.0_dp) * z * x1
                call evaluate_susceptibility(x1, x2, conservation_model, numerical)
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
        call evaluate_susceptibility(x1, x2, 1, numerical)
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

    call check_finite_collision_ray()

contains

    subroutine check_finite_collision_ray()
        ! Independent 60-digit quadrature of the generating integral, inserted
        ! in the energy-conserving susceptibility (Markl 2023 appendix A).
        ! These points stay below evaluate_susceptibility's Landau substitution
        ! and require the actual collisional special-function evaluation.
        real(dp), parameter :: ray_x1(3) = [3000.0_dp, 6000.0_dp, 9000.0_dp]
        complex(dp), parameter :: expected_j(3) = [ &
            cmplx(-0.09900148378348024263_dp, 0.01247667045059763900_dp, dp), &
            cmplx(-0.09900240536971374563_dp, 0.01247365120557467912_dp, dp), &
            cmplx(-0.09900271250295240704_dp, 0.01247264481483745869_dp, dp)]
        complex(dp) :: susceptibility(0:3, 0:3), response, expected
        real(dp) :: signed_x1, error
        integer :: point, harmonic

        do point = 1, size(ray_x1)
            do harmonic = -1, 1, 2
                signed_x1 = real(harmonic, dp) * ray_x1(point)
                call evaluate_susceptibility(signed_x1, -0.1_dp * signed_x1, 1, susceptibility)
                response = cmplx(0.0_dp, 1.0_dp, dp) * signed_x1 * susceptibility(1, 1)
                expected = expected_j(point)
                if (harmonic < 0) expected = conjg(expected)
                error = abs(response - expected) / abs(expected)
                print '(a,f8.0,a,es12.4)', 'finite-collision x1=', signed_x1, ' error=', error
                if (error > 1.0e-8_dp) error stop 'finite-collision response missed quadrature oracle'
            end do
        end do
    end subroutine check_finite_collision_ray

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
