program test_adaptive_normalization
    use KIM_kinds_m, only: dp
    use constants_m, only: pi
    use grid_m, only: gauss_int_nodes_Nx, gauss_int_nodes_Nxp, &
        quadpack_algorithm, quadpack_epsabs, quadpack_epsrel, &
        quadpack_key, quadpack_limit, quadpack_use_u_substitution, &
        rg_grid, rkf45_atol, rkf45_rtol, theta_integration_method
    use integrals_rkf45_m, only: init_rkf45_int, integrate_F1, integrate_F2, &
        integrate_F3, rkf45_config_t
    use integrands_rkf45_m, only: rkf45_integrand_context_t
    use quadpack_integration_m, only: init_quadpack_module

    implicit none

    integer, parameter :: n_radial = 600
    integer, parameter :: n_theta = 600
    real(dp), parameter :: relative_tolerance = 2.0e-5_dp
    real(dp), parameter :: theta_lower = 0.2_dp
    real(dp), parameter :: theta_upper = pi - theta_lower

    type(rkf45_config_t) :: integration_config
    type(rkf45_integrand_context_t) :: context
    real(dp) :: actual(3), oracle(3)
    integer :: failures

    call configure_case(integration_config, context)
    call direct_segment_oracle(context, oracle)

    failures = 0
    theta_integration_method = "RKF45"
    call run_adaptive(integration_config, context, actual)
    call check_backend("RKF45", actual, oracle, failures)

    theta_integration_method = "QUADPACK"
    quadpack_algorithm = "QAG"
    call run_adaptive(integration_config, context, actual)
    call check_backend("QUADPACK QAG", actual, oracle, failures)

    quadpack_algorithm = "QAGS"
    call run_adaptive(integration_config, context, actual)
    call check_backend("QUADPACK QAGS", actual, oracle, failures)

    if (failures /= 0) then
        write (*, '(A,I0)') "adaptive normalization failures: ", failures
        error stop 1
    end if

    write (*, '(A)') "adaptive F1-F3 normalization matches direct segment integrals"

contains

    subroutine configure_case(config, test_context)
        type(rkf45_config_t), intent(out) :: config
        type(rkf45_integrand_context_t), intent(out) :: test_context

        gauss_int_nodes_Nx = 1
        gauss_int_nodes_Nxp = 1
        rkf45_atol = 1.0e-8_dp
        rkf45_rtol = 1.0e-8_dp
        quadpack_epsabs = 1.0e-11_dp
        quadpack_epsrel = 1.0e-11_dp
        quadpack_key = 6
        quadpack_limit = 500
        quadpack_use_u_substitution = .true.

        rg_grid%npts_b = 2
        rg_grid%npts_c = 1
        allocate (rg_grid%xb(2))
        rg_grid%xb = [-0.45_dp, 0.8_dp]

        config%Nx = gauss_int_nodes_Nx
        config%Nxp = gauss_int_nodes_Nxp
        call init_rkf45_int(config)
        call init_quadpack_module()

        test_context%j = 1
        test_context%rhoT = 0.7_dp
        test_context%ks = 0.4_dp
        test_context%x = -0.3_dp
        test_context%xp = 0.8_dp
        test_context%xlm1 = -0.5_dp
        test_context%xl = test_context%x
        test_context%xlp1 = -0.1_dp
        test_context%xlpm1 = 0.6_dp
        test_context%xlp = test_context%xp
        test_context%xlpp1 = 1.0_dp
    end subroutine configure_case

    subroutine run_adaptive(config, test_context, result)
        type(rkf45_config_t), intent(in) :: config
        type(rkf45_integrand_context_t), intent(inout) :: test_context
        real(dp), intent(out) :: result(3)

        call integrate_F1(result(1), config, test_context)
        call integrate_F2(result(2), config, test_context)
        call integrate_F3(result(3), config, test_context)
    end subroutine run_adaptive

    subroutine direct_segment_oracle(test_context, result)
        type(rkf45_integrand_context_t), intent(in) :: test_context
        real(dp), intent(out) :: result(3)

        real(dp) :: c, d, d_radial, d_theta, gaussian, outside
        real(dp) :: pointwise(3), radial_point, radial_sum(3), s, theta
        real(dp) :: theta_sum(3), xbar
        integer :: radial_index, theta_index

        d = test_context%x - test_context%xp
        xbar = 0.5_dp * (test_context%x + test_context%xp)
        d_radial = (rg_grid%xb(2) - rg_grid%xb(1)) / real(n_radial, dp)
        d_theta = (theta_upper - theta_lower) / real(n_theta, dp)
        theta_sum = 0.0_dp

        do theta_index = 0, n_theta
            theta = theta_lower + real(theta_index, dp) * d_theta
            c = cos(theta)
            s = sin(theta)
            outside = exp(-test_context%ks**2 * test_context%rhoT**2 &
                - d**2 / (4.0_dp * test_context%rhoT**2 * (1.0_dp - c)))
            radial_sum = 0.0_dp

            do radial_index = 0, n_radial
                radial_point = rg_grid%xb(1) + real(radial_index, dp) * d_radial
                gaussian = exp(-(radial_point - xbar)**2 &
                    / (test_context%rhoT**2 * (1.0_dp + c)))
                pointwise(1) = outside * gaussian / (test_context%rhoT**2 * s)
                pointwise(2) = pointwise(1) * ( &
                    test_context%rhoT**2 * test_context%ks**2 + 1.0_dp / s**2 &
                    - (radial_point - xbar)**2 &
                    / (test_context%rhoT**2 * (1.0_dp + c)**2) &
                    - d**2 / (4.0_dp * test_context%rhoT**2 * (1.0_dp - c)**2))
                pointwise(3) = pointwise(1) * c * ( &
                    (radial_point - xbar)**2 &
                    / (test_context%rhoT**2 * (1.0_dp + c)**2) &
                    - d**2 / (4.0_dp * test_context%rhoT**2 * (1.0_dp - c)**2) &
                    + c / s**2)
                radial_sum = radial_sum &
                    + simpson_weight(radial_index, n_radial) * pointwise
            end do

            radial_sum = radial_sum * d_radial / 3.0_dp
            theta_sum = theta_sum &
                + simpson_weight(theta_index, n_theta) * radial_sum
        end do

        result = theta_sum * d_theta / 3.0_dp
        result = result * (test_context%xlp1 - test_context%xlm1) &
            * (test_context%xlpp1 - test_context%xlpm1)
    end subroutine direct_segment_oracle

    subroutine check_backend(backend, result, expected, failure_count)
        character(*), intent(in) :: backend
        real(dp), intent(in) :: result(3), expected(3)
        integer, intent(inout) :: failure_count

        character(len=2), parameter :: label(3) = ["F1", "F2", "F3"]
        real(dp) :: relative_error, scale, stale_error
        integer :: index

        do index = 1, 3
            scale = max(abs(expected(index)), 1.0e-12_dp)
            relative_error = abs(result(index) - expected(index)) / scale
            stale_error = abs(result(index) - 2.0_dp * pi * expected(index)) / scale
            write (*, '(A,1X,A,3(1X,A,ES13.5))') trim(backend), label(index), &
                "actual=", result(index), "oracle=", expected(index), &
                "ratio=", result(index) / expected(index)

            if (relative_error > relative_tolerance) then
                write (*, '(A,ES13.5)') "  direct-oracle relative error: ", &
                    relative_error
                failure_count = failure_count + 1
            end if
            if (stale_error < 1.0_dp) then
                write (*, '(A,ES13.5)') "  stale 2*pi mutation error: ", stale_error
                failure_count = failure_count + 1
            end if
        end do
    end subroutine check_backend

    pure function simpson_weight(index, interval_count) result(weight)
        integer, intent(in) :: index, interval_count
        real(dp) :: weight

        if (index == 0 .or. index == interval_count) then
            weight = 1.0_dp
        else if (mod(index, 2) == 0) then
            weight = 2.0_dp
        else
            weight = 4.0_dp
        end if
    end function simpson_weight

end program test_adaptive_normalization
