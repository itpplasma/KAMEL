program test_adaptive_local_cutoff
    use KIM_kinds_m, only: dp
    use constants_m, only: pi
    use grid_m, only: quadpack_algorithm, quadpack_epsabs, quadpack_epsrel, &
        quadpack_key, quadpack_limit, quadpack_use_u_substitution, rg_grid, &
        rkf45_atol, rkf45_rtol, theta_integration_method
    use integrals_rkf45_m, only: init_rkf45_int, integrate_F1, integrate_F2, &
        integrate_F3, rkf45_config_t
    use integrands_rkf45_m, only: rkf45_integrand_context_t
    use omp_lib, only: omp_set_dynamic, omp_set_num_threads
    use quadpack_integration_m, only: init_quadpack_module, &
        quadpack_integrate_F1, quadpack_integrate_F2, quadpack_integrate_F3

    implicit none

    real(dp), parameter :: comparison_tolerance = 2.0e-5_dp
    real(dp), parameter :: mutation_floor = 1.0e-4_dp
    real(dp), parameter :: sentinel_a(2) = [7.0_dp, 7.18_dp]
    real(dp), parameter :: sentinel_b(2) = [-9.0_dp, 11.0_dp]

    type(rkf45_config_t) :: config
    type(rkf45_integrand_context_t) :: context
    real(dp) :: local_oracle(3), shared_lower(3), shared_upper(3)
    integer :: failures

    call configure_case(config, context)
    call omp_set_dynamic(.false.)
    call omp_set_num_threads(2)

    quadpack_algorithm = "QAG"
    call explicit_cutoff_sum(config, context, .true., 0.0_dp, local_oracle)
    call explicit_cutoff_sum(config, context, .false., 0.05_dp, shared_lower)
    call explicit_cutoff_sum(config, context, .false., 0.20_dp, shared_upper)

    failures = 0
    call check_cutoff_sensitivity(local_oracle, shared_lower, shared_upper, failures)
    call check_backend("RKF45", "RKF45", config, context, local_oracle, failures)
    call check_backend("QAG", "QUADPACK", config, context, local_oracle, failures)
    call check_backend("QAGS", "QUADPACK", config, context, local_oracle, failures)

    if (failures /= 0) then
        write (*, '(A,I0)') "adaptive local-cutoff failures: ", failures
        error stop 1
    end if

    write (*, '(A)') "adaptive endpoints use each quadrature node's local separation"

contains

    subroutine configure_case(integration_config, test_context)
        type(rkf45_config_t), intent(out) :: integration_config
        type(rkf45_integrand_context_t), intent(out) :: test_context

        rkf45_atol = 1.0e-7_dp
        rkf45_rtol = 1.0e-7_dp
        quadpack_epsabs = 1.0e-10_dp
        quadpack_epsrel = 1.0e-10_dp
        quadpack_key = 6
        quadpack_limit = 500
        quadpack_use_u_substitution = .true.

        rg_grid%npts_b = 2
        rg_grid%npts_c = 1
        allocate (rg_grid%xb(2))
        rg_grid%xb = [-0.6_dp, 0.8_dp]

        integration_config%Nx = 2
        integration_config%Nxp = 2
        call init_rkf45_int(integration_config)
        call init_quadpack_module()

        test_context%j = 1
        test_context%rhoT = 0.25_dp
        test_context%ks = 0.4_dp
        test_context%x = sentinel_a(1)
        test_context%xp = sentinel_a(2)
        test_context%xlm1 = -0.4_dp
        test_context%xl = 0.0_dp
        test_context%xlp1 = 0.4_dp
        test_context%xlpm1 = -0.25_dp
        test_context%xlp = 0.15_dp
        test_context%xlpp1 = 0.55_dp
    end subroutine configure_case

    subroutine check_backend(algorithm, method, integration_config, test_context, &
            expected, failure_count)
        character(*), intent(in) :: algorithm, method
        type(rkf45_config_t), intent(in) :: integration_config
        type(rkf45_integrand_context_t), intent(inout) :: test_context
        real(dp), intent(in) :: expected(3)
        integer, intent(inout) :: failure_count

        real(dp) :: first(3), second(3)
        integer :: function_index

        theta_integration_method = method
        quadpack_algorithm = algorithm
        call evaluate(integration_config, test_context, sentinel_a, first)
        call evaluate(integration_config, test_context, sentinel_b, second)

        do function_index = 1, 3
            call assert_close(algorithm // " sentinel invariance", function_index, &
                first(function_index), second(function_index), &
                comparison_tolerance, failure_count)
            call assert_close(algorithm // " local-cutoff oracle", function_index, &
                first(function_index), expected(function_index), &
                comparison_tolerance, failure_count)
        end do
    end subroutine check_backend

    subroutine evaluate(integration_config, test_context, sentinel, result)
        type(rkf45_config_t), intent(in) :: integration_config
        type(rkf45_integrand_context_t), intent(inout) :: test_context
        real(dp), intent(in) :: sentinel(2)
        real(dp), intent(out) :: result(3)

        test_context%x = sentinel(1)
        test_context%xp = sentinel(2)
        call integrate_F1(result(1), integration_config, test_context)
        test_context%x = sentinel(1)
        test_context%xp = sentinel(2)
        call integrate_F2(result(2), integration_config, test_context)
        test_context%x = sentinel(1)
        test_context%xp = sentinel(2)
        call integrate_F3(result(3), integration_config, test_context)
    end subroutine evaluate

    subroutine explicit_cutoff_sum(integration_config, test_context, use_local, &
            fixed_cutoff, result)
        type(rkf45_config_t), intent(in) :: integration_config
        type(rkf45_integrand_context_t), intent(in) :: test_context
        logical, intent(in) :: use_local
        real(dp), intent(in) :: fixed_cutoff
        real(dp), intent(out) :: result(3)

        type(rkf45_integrand_context_t) :: node_context
        real(dp) :: cutoff, node_result(3), normalization, weight
        integer :: j, k

        result = 0.0_dp
        do j = 1, integration_config%Nxp
            do k = 1, integration_config%Nx
                node_context = test_context
                node_context%xp = map_node(integration_config%x_xp(j), &
                    test_context%xlpm1, test_context%xlpp1)
                node_context%x = map_node(integration_config%x_x(k), &
                    test_context%xlm1, test_context%xlp1)
                if (use_local) then
                    cutoff = local_cutoff(node_context%x, node_context%xp, &
                        node_context%rhoT)
                else
                    cutoff = fixed_cutoff
                end if

                call integrate_node(node_context, cutoff, node_result)
                weight = integration_config%w_xp(j) * integration_config%w_x(k)
                result = result + weight * node_result
            end do
        end do

        normalization = (test_context%xlp1 - test_context%xlm1) &
            * (test_context%xlpp1 - test_context%xlpm1) / 4.0_dp
        result = result * exp(-test_context%ks**2 * test_context%rhoT**2) &
            * normalization
        result(2) = result(2) * (-1.0_dp) / (8.0_dp * test_context%rhoT**4)
        result(3) = result(3) * (-1.0_dp) / (4.0_dp * test_context%rhoT**4)
    end subroutine explicit_cutoff_sum

    subroutine integrate_node(node_context, cutoff, result)
        type(rkf45_integrand_context_t), intent(in) :: node_context
        real(dp), intent(in) :: cutoff
        real(dp), intent(out) :: result(3)

        call quadpack_integrate_F1(result(1), quadpack_epsabs, quadpack_epsrel, &
            node_context, cutoff, pi - cutoff, quadpack_use_u_substitution)
        call quadpack_integrate_F2(result(2), quadpack_epsabs, quadpack_epsrel, &
            node_context, cutoff, pi - cutoff, quadpack_use_u_substitution)
        call quadpack_integrate_F3(result(3), quadpack_epsabs, quadpack_epsrel, &
            node_context, cutoff, pi - cutoff, quadpack_use_u_substitution)
    end subroutine integrate_node

    subroutine check_cutoff_sensitivity(local, lower, upper, failure_count)
        real(dp), intent(in) :: local(3), lower(3), upper(3)
        integer, intent(inout) :: failure_count

        real(dp) :: lower_error, scale, upper_error
        integer :: function_index

        do function_index = 1, 3
            scale = max(abs(local(function_index)), tiny(1.0_dp))
            lower_error = abs(local(function_index) - lower(function_index)) / scale
            upper_error = abs(local(function_index) - upper(function_index)) / scale
            write (*, '(A,I0,2(1X,A,ES12.4))') "cutoff sensitivity F", &
                function_index, "lower=", lower_error, "upper=", upper_error
            if (max(lower_error, upper_error) < mutation_floor) then
                write (*, '(A,I0)') "FAIL: shared-cutoff mutation is invisible for F", &
                    function_index
                failure_count = failure_count + 1
            end if
        end do
    end subroutine check_cutoff_sensitivity

    subroutine assert_close(label, function_index, actual, expected, tolerance, &
            failure_count)
        character(*), intent(in) :: label
        integer, intent(in) :: function_index
        real(dp), intent(in) :: actual, expected, tolerance
        integer, intent(inout) :: failure_count

        real(dp) :: error, scale

        scale = max(abs(expected), tiny(1.0_dp))
        error = abs(actual - expected) / scale
        write (*, '(A,1X,A,I0,3(1X,A,ES13.5))') trim(label), "F", &
            function_index, "actual=", actual, "expected=", expected, &
            "relative_error=", error
        if (error > tolerance) then
            failure_count = failure_count + 1
        end if
    end subroutine assert_close

    pure function local_cutoff(x, xp, rho_thermal) result(cutoff)
        real(dp), intent(in) :: x, xp, rho_thermal
        real(dp) :: cutoff

        cutoff = max(0.05_dp, min(0.20_dp, &
            abs(x - xp) / (sqrt(50.0_dp) * max(abs(rho_thermal), tiny(1.0_dp)))))
    end function local_cutoff

    pure function map_node(node, lower, upper) result(mapped)
        real(dp), intent(in) :: node, lower, upper
        real(dp) :: mapped

        mapped = 0.5_dp * ((upper - lower) * node + upper + lower)
    end function map_node

end program test_adaptive_local_cutoff
