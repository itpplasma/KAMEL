program test_integration_methods
    ! Test program to compare RKF45 and QUADPACK integration methods
    ! Verifies consistency between methods for various test integrands
    
    use KIM_kinds_m, only: dp
    use integrals_rkf45_m, only: rkf45_config_t, init_rkf45_int, &
                                 integrate_F1, integrate_F2, integrate_F3
    use integrands_rkf45_m, only: rkf45_integrand_context_t
    use grid_m, only: theta_integration_method, rkf45_atol, rkf45_rtol, &
                      quadpack_epsabs, quadpack_epsrel, quadpack_use_u_substitution, quadpack_algorithm, &
                      gauss_int_nodes_Nx, gauss_int_nodes_Nxp, rg_grid
    use quadpack_integration_m, only: init_quadpack_module
    use omp_lib
    
    implicit none
    
    ! Test parameters
    integer, parameter :: n_tests = 5
    real(dp), parameter :: test_tolerance = 1.0d-8
    
    ! Integration variables
    type(rkf45_config_t) :: rkf45_conf
    type(rkf45_integrand_context_t) :: context
    real(dp) :: result_rkf45, result_quadpack
    real(dp) :: rel_error, abs_error
    integer :: test_case, i
    logical :: all_tests_passed
    real(dp) :: start_time, end_time, time_rkf45, time_quadpack
    
    ! Test case parameters
    real(dp), dimension(n_tests) :: rhoT_values
    real(dp), dimension(n_tests) :: ks_values
    real(dp), dimension(n_tests) :: x_values
    real(dp), dimension(n_tests) :: xp_values
    
    ! Initialize test output
    write(*,*) "==============================================="
    write(*,*) "Integration Method Comparison Test"
    write(*,*) "==============================================="
    write(*,*)
    
    ! Set integration parameters
    gauss_int_nodes_Nx = 20
    gauss_int_nodes_Nxp = 20
    rkf45_atol = 1.0d-10
    rkf45_rtol = 1.0d-10
    quadpack_epsabs = 1.0d-10
    quadpack_epsrel = 1.0d-10
    quadpack_use_u_substitution = .true.
    
    ! Initialize minimal rg_grid used by integrands (two-cell domain [0,1])
    rg_grid%npts_b = 2
    rg_grid%npts_c = 1
    allocate(rg_grid%xb(rg_grid%npts_b))
    rg_grid%xb(1) = 0.0d0
    rg_grid%xb(2) = 1.0d0

    ! Initialize integration modules
    rkf45_conf%Nx = gauss_int_nodes_Nx
    rkf45_conf%Nxp = gauss_int_nodes_Nxp
    call init_rkf45_int(rkf45_conf)
    call init_quadpack_module()
    
    ! Set up test cases
    rhoT_values = [0.1d0, 0.5d0, 1.0d0, 2.0d0, 5.0d0]
    ks_values = [0.1d0, 0.5d0, 1.0d0, 2.0d0, 3.0d0]
    x_values = [0.1d0, 0.3d0, 0.5d0, 0.7d0, 0.9d0]
    xp_values = [0.2d0, 0.4d0, 0.6d0, 0.8d0, 0.95d0]
    
    all_tests_passed = .true.
    
    ! Test F1 integrand (QAG)
    write(*,*) "Testing F1 integrand..."
    write(*,*) "----------------------------------------"
    
    do test_case = 1, n_tests
        ! Set up context for test case
        context%j = 1
        context%rhoT = rhoT_values(test_case)
        context%ks = ks_values(test_case)
        context%x = x_values(test_case)
        context%xp = xp_values(test_case)
        context%xlm1 = 0.0d0
        context%xl = 0.5d0
        context%xlp1 = 1.0d0
        context%xlpm1 = 0.0d0
        context%xlp = 0.5d0
        context%xlpp1 = 1.0d0
        
        ! Test with RKF45
        theta_integration_method = "RKF45"
        start_time = omp_get_wtime()
        call integrate_F1(result_rkf45, rkf45_conf, context)
        end_time = omp_get_wtime()
        time_rkf45 = end_time - start_time
        
        ! Test with QUADPACK (QAG)
        theta_integration_method = "QUADPACK"
        quadpack_algorithm = "QAG"
        start_time = omp_get_wtime()
        call integrate_F1(result_quadpack, rkf45_conf, context)
        end_time = omp_get_wtime()
        time_quadpack = end_time - start_time
        
        ! Calculate errors
        abs_error = abs(result_quadpack - result_rkf45)
        if (abs(result_rkf45) > 1.0d-15) then
            rel_error = abs_error / abs(result_rkf45)
        else
            rel_error = abs_error
        end if
        
        ! Print results
        write(*,'(A,I2,A)') "Test case ", test_case, ":"
        write(*,'(A,ES12.5,A,ES12.5,A,ES12.5,A,ES12.5)') &
            "  Parameters: rhoT=", rhoT_values(test_case), &
            ", ks=", ks_values(test_case), &
            ", x=", x_values(test_case), &
            ", xp=", xp_values(test_case)
        write(*,'(A,ES16.9)') "  RKF45 result:    ", result_rkf45
        write(*,'(A,ES16.9)') "  QUADPACK result: ", result_quadpack
        write(*,'(A,ES12.5)') "  Absolute error:  ", abs_error
        write(*,'(A,ES12.5)') "  Relative error:  ", rel_error
        write(*,'(A,F8.5,A,F8.5,A)') "  Time: RKF45=", time_rkf45*1000, &
                                     "ms, QUADPACK=", time_quadpack*1000, "ms"
        
        if (rel_error > test_tolerance) then
            write(*,*) "  WARNING: Error exceeds tolerance!"
            all_tests_passed = .false.
        else
            write(*,*) "  PASSED"
        end if
        write(*,*)
    end do
    
    ! Test F2 integrand (QAG)
    write(*,*) "Testing F2 integrand..."
    write(*,*) "----------------------------------------"
    
    do test_case = 1, n_tests
        ! Set up context
        context%j = 1
        context%rhoT = rhoT_values(test_case)
        context%ks = ks_values(test_case)
        context%x = x_values(test_case)
        context%xp = xp_values(test_case)
        context%xlm1 = 0.0d0
        context%xl = 0.5d0
        context%xlp1 = 1.0d0
        context%xlpm1 = 0.0d0
        context%xlp = 0.5d0
        context%xlpp1 = 1.0d0
        
        ! Test with RKF45
        theta_integration_method = "RKF45"
        start_time = omp_get_wtime()
        call integrate_F2(result_rkf45, rkf45_conf, context)
        end_time = omp_get_wtime()
        time_rkf45 = end_time - start_time
        
        ! Test with QUADPACK (QAG)
        theta_integration_method = "QUADPACK"
        quadpack_algorithm = "QAG"
        start_time = omp_get_wtime()
        call integrate_F2(result_quadpack, rkf45_conf, context)
        end_time = omp_get_wtime()
        time_quadpack = end_time - start_time
        
        ! Calculate errors
        abs_error = abs(result_quadpack - result_rkf45)
        if (abs(result_rkf45) > 1.0d-15) then
            rel_error = abs_error / abs(result_rkf45)
        else
            rel_error = abs_error
        end if
        
        ! Print results
        write(*,'(A,I2,A)') "Test case ", test_case, ":"
        write(*,'(A,ES12.5,A,ES12.5,A,ES12.5,A,ES12.5)') &
            "  Parameters: rhoT=", rhoT_values(test_case), &
            ", ks=", ks_values(test_case), &
            ", x=", x_values(test_case), &
            ", xp=", xp_values(test_case)
        write(*,'(A,ES16.9)') "  RKF45 result:    ", result_rkf45
        write(*,'(A,ES16.9)') "  QUADPACK result: ", result_quadpack
        write(*,'(A,ES12.5)') "  Absolute error:  ", abs_error
        write(*,'(A,ES12.5)') "  Relative error:  ", rel_error
        write(*,'(A,F8.5,A,F8.5,A)') "  Time: RKF45=", time_rkf45*1000, &
                                     "ms, QUADPACK=", time_quadpack*1000, "ms"
        
        if (rel_error > test_tolerance) then
            write(*,*) "  WARNING: Error exceeds tolerance!"
            all_tests_passed = .false.
        else
            write(*,*) "  PASSED"
        end if
        write(*,*)
    end do
    
    ! Test F3 integrand (QAG)
    write(*,*) "Testing F3 integrand..."
    write(*,*) "----------------------------------------"
    
    do test_case = 1, n_tests
        ! Set up context
        context%j = 1
        context%rhoT = rhoT_values(test_case)
        context%ks = ks_values(test_case)
        context%x = x_values(test_case)
        context%xp = xp_values(test_case)
        context%xlm1 = 0.0d0
        context%xl = 0.5d0
        context%xlp1 = 1.0d0
        context%xlpm1 = 0.0d0
        context%xlp = 0.5d0
        context%xlpp1 = 1.0d0
        
        ! Test with RKF45
        theta_integration_method = "RKF45"
        start_time = omp_get_wtime()
        call integrate_F3(result_rkf45, rkf45_conf, context)
        end_time = omp_get_wtime()
        time_rkf45 = end_time - start_time
        
        ! Test with QUADPACK (QAG)
        theta_integration_method = "QUADPACK"
        quadpack_algorithm = "QAG"
        start_time = omp_get_wtime()
        call integrate_F3(result_quadpack, rkf45_conf, context)
        end_time = omp_get_wtime()
        time_quadpack = end_time - start_time
        
        ! Calculate errors
        abs_error = abs(result_quadpack - result_rkf45)
        if (abs(result_rkf45) > 1.0d-15) then
            rel_error = abs_error / abs(result_rkf45)
        else
            rel_error = abs_error
        end if
        
        ! Print results
        write(*,'(A,I2,A)') "Test case ", test_case, ":"
        write(*,'(A,ES12.5,A,ES12.5,A,ES12.5,A,ES12.5)') &
            "  Parameters: rhoT=", rhoT_values(test_case), &
            ", ks=", ks_values(test_case), &
            ", x=", x_values(test_case), &
            ", xp=", xp_values(test_case)
        write(*,'(A,ES16.9)') "  RKF45 result:    ", result_rkf45
        write(*,'(A,ES16.9)') "  QUADPACK result: ", result_quadpack
        write(*,'(A,ES12.5)') "  Absolute error:  ", abs_error
        write(*,'(A,ES12.5)') "  Relative error:  ", rel_error
        write(*,'(A,F8.5,A,F8.5,A)') "  Time: RKF45=", time_rkf45*1000, &
                                     "ms, QUADPACK=", time_quadpack*1000, "ms"
        
        if (rel_error > test_tolerance) then
            write(*,*) "  WARNING: Error exceeds tolerance!"
            all_tests_passed = .false.
        else
            write(*,*) "  PASSED"
        end if
        write(*,*)
    end do
    
    ! Repeat tests with QUADPACK QAGS algorithm
    write(*,*)
    write(*,*) "Testing QUADPACK QAGS algorithm..."
    write(*,*) "----------------------------------------"

    ! F1 with QAGS (single case suffices to exercise branch)
    context%j = 1
    context%rhoT = 0.5d0
    context%ks = 1.0d0
    context%x = 0.4d0
    context%xp = 0.6d0
    context%xlm1 = 0.0d0
    context%xl = 0.5d0
    context%xlp1 = 1.0d0
    context%xlpm1 = 0.0d0
    context%xlp = 0.5d0
    context%xlpp1 = 1.0d0

    theta_integration_method = "RKF45"
    call integrate_F1(result_rkf45, rkf45_conf, context)
    theta_integration_method = "QUADPACK"
    quadpack_algorithm = "QAGS"
    call integrate_F1(result_quadpack, rkf45_conf, context)
    abs_error = abs(result_quadpack - result_rkf45)
    if (abs(result_rkf45) > 1.0d-15) then
        rel_error = abs_error / abs(result_rkf45)
    else
        rel_error = abs_error
    end if
    write(*,'(A,ES16.9)') "  RKF45 result (F1):    ", result_rkf45
    write(*,'(A,ES16.9)') "  QUADPACK QAGS (F1):   ", result_quadpack
    write(*,'(A,ES12.5)') "  Relative error:       ", rel_error
    if (rel_error > test_tolerance) then
        write(*,*) "  WARNING: QAGS error exceeds tolerance!"
        all_tests_passed = .false.
    else
        write(*,*) "  QAGS F1 PASSED"
    end if

    ! Summary
    write(*,*) "==============================================="
    if (all_tests_passed) then
        write(*,*) "ALL TESTS PASSED"
    else
        write(*,*) "SOME TESTS FAILED"
    end if
    write(*,*) "==============================================="
    
    if (.not. all_tests_passed) then
        stop 1
    end if
    
end program test_integration_methods
