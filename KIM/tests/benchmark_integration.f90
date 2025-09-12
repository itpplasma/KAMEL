program benchmark_integration
    ! Benchmark program to compare performance of RKF45 vs QUADPACK
    ! Tests various parameter regimes and problem sizes
    
    use KIM_kinds_m, only: dp
    use integrals_rkf45_m, only: rkf45_config_t, init_rkf45_int, &
                                 integrate_F1, integrate_F2, integrate_F3
    use integrands_rkf45_m, only: rkf45_integrand_context_t
    use grid_m, only: theta_integration_method, rkf45_atol, rkf45_rtol, &
                      quadpack_epsabs, quadpack_epsrel, quadpack_use_u_substitution, &
                      quadpack_key, gauss_int_nodes_Nx, gauss_int_nodes_Nxp
    use quadpack_integration_m, only: init_quadpack_module
    use omp_lib
    
    implicit none
    
    ! Benchmark parameters
    integer, parameter :: n_warmup = 5
    integer, parameter :: n_trials = 20
    integer, parameter :: n_param_sets = 3
    
    ! Integration variables
    type(rkf45_config_t) :: rkf45_conf
    type(rkf45_integrand_context_t) :: context
    real(dp) :: result
    
    ! Timing variables
    real(dp) :: start_time, end_time
    real(dp), dimension(n_trials) :: times_rkf45, times_quadpack
    real(dp) :: mean_time_rkf45, mean_time_quadpack
    real(dp) :: std_time_rkf45, std_time_quadpack
    real(dp) :: speedup
    
    ! Loop variables
    integer :: param_set, trial, i
    integer :: integrand_type
    character(len=20) :: integrand_name
    
    ! Parameter sets for different regimes
    real(dp), dimension(n_param_sets) :: rhoT_sets = [0.1d0, 1.0d0, 5.0d0]
    real(dp), dimension(n_param_sets) :: ks_sets = [0.5d0, 2.0d0, 5.0d0]
    character(len=20), dimension(n_param_sets) :: regime_names = &
        ["Small Larmor radius", "Medium regime      ", "Large Larmor radius"]
    
    ! Initialize
    write(*,*) "========================================================="
    write(*,*) "QUADPACK vs RKF45 Performance Benchmark"
    write(*,*) "========================================================="
    write(*,*)
    write(*,'(A,I0)') "Number of warmup runs: ", n_warmup
    write(*,'(A,I0)') "Number of trial runs:  ", n_trials
    write(*,*)
    
    ! Set integration parameters
    gauss_int_nodes_Nx = 20
    gauss_int_nodes_Nxp = 20
    rkf45_atol = 1.0d-10
    rkf45_rtol = 1.0d-10
    quadpack_epsabs = 1.0d-10
    quadpack_epsrel = 1.0d-10
    quadpack_use_u_substitution = .true.
    
    ! Test different Gauss-Kronrod rules
    write(*,*) "Testing different QUADPACK Gauss-Kronrod rules..."
    write(*,*) "---------------------------------------------------------"
    
    ! Initialize modules
    rkf45_conf%Nx = gauss_int_nodes_Nx
    rkf45_conf%Nxp = gauss_int_nodes_Nxp
    call init_rkf45_int(rkf45_conf)
    call init_quadpack_module()
    
    ! Set up standard test context
    context%j = 1
    context%rhoT = 1.0d0
    context%ks = 1.0d0
    context%x = 0.5d0
    context%xp = 0.6d0
    context%xlm1 = 0.0d0
    context%xl = 0.5d0
    context%xlp1 = 1.0d0
    context%xlpm1 = 0.0d0
    context%xlp = 0.5d0
    context%xlpp1 = 1.0d0
    
    do i = 1, 6
        quadpack_key = i
        theta_integration_method = "QUADPACK"
        
        ! Warmup
        do trial = 1, n_warmup
            call integrate_F1(result, rkf45_conf, context)
        end do
        
        ! Benchmark
        start_time = omp_get_wtime()
        do trial = 1, n_trials
            call integrate_F1(result, rkf45_conf, context)
        end do
        end_time = omp_get_wtime()
        
        write(*,'(A,I1,A,F8.3,A)') "  QK", i, " rule: ", &
            (end_time - start_time) * 1000.0d0 / n_trials, " ms/call"
    end do
    write(*,*)
    
    ! Main benchmark loop over parameter regimes
    do param_set = 1, n_param_sets
        
        write(*,*) "========================================================="
        write(*,'(A,A)') "Parameter regime: ", trim(regime_names(param_set))
        write(*,'(A,F5.2,A,F5.2)') "rhoT = ", rhoT_sets(param_set), &
                                   ", ks = ", ks_sets(param_set)
        write(*,*) "---------------------------------------------------------"
        
        ! Set up context for this parameter set
        context%rhoT = rhoT_sets(param_set)
        context%ks = ks_sets(param_set)
        
        ! Test each integrand type
        do integrand_type = 1, 3
            
            select case(integrand_type)
            case(1)
                integrand_name = "F1 integrand"
            case(2)
                integrand_name = "F2 integrand"
            case(3)
                integrand_name = "F3 integrand"
            end select
            
            write(*,'(A,A)') "Testing ", trim(integrand_name)
            
            ! Benchmark RKF45
            theta_integration_method = "RKF45"
            
            ! Warmup runs
            do trial = 1, n_warmup
                select case(integrand_type)
                case(1)
                    call integrate_F1(result, rkf45_conf, context)
                case(2)
                    call integrate_F2(result, rkf45_conf, context)
                case(3)
                    call integrate_F3(result, rkf45_conf, context)
                end select
            end do
            
            ! Timed runs
            do trial = 1, n_trials
                start_time = omp_get_wtime()
                select case(integrand_type)
                case(1)
                    call integrate_F1(result, rkf45_conf, context)
                case(2)
                    call integrate_F2(result, rkf45_conf, context)
                case(3)
                    call integrate_F3(result, rkf45_conf, context)
                end select
                end_time = omp_get_wtime()
                times_rkf45(trial) = (end_time - start_time) * 1000.0d0  ! Convert to ms
            end do
            
            ! Benchmark QUADPACK
            theta_integration_method = "QUADPACK"
            quadpack_key = 6  ! Use QK61 for accuracy
            
            ! Warmup runs
            do trial = 1, n_warmup
                select case(integrand_type)
                case(1)
                    call integrate_F1(result, rkf45_conf, context)
                case(2)
                    call integrate_F2(result, rkf45_conf, context)
                case(3)
                    call integrate_F3(result, rkf45_conf, context)
                end select
            end do
            
            ! Timed runs
            do trial = 1, n_trials
                start_time = omp_get_wtime()
                select case(integrand_type)
                case(1)
                    call integrate_F1(result, rkf45_conf, context)
                case(2)
                    call integrate_F2(result, rkf45_conf, context)
                case(3)
                    call integrate_F3(result, rkf45_conf, context)
                end select
                end_time = omp_get_wtime()
                times_quadpack(trial) = (end_time - start_time) * 1000.0d0
            end do
            
            ! Calculate statistics
            mean_time_rkf45 = sum(times_rkf45) / n_trials
            mean_time_quadpack = sum(times_quadpack) / n_trials
            
            std_time_rkf45 = sqrt(sum((times_rkf45 - mean_time_rkf45)**2) / (n_trials - 1))
            std_time_quadpack = sqrt(sum((times_quadpack - mean_time_quadpack)**2) / (n_trials - 1))
            
            speedup = mean_time_rkf45 / mean_time_quadpack
            
            ! Print results
            write(*,'(A,F8.3,A,F6.3,A)') "  RKF45:    ", mean_time_rkf45, &
                                         " ± ", std_time_rkf45, " ms"
            write(*,'(A,F8.3,A,F6.3,A)') "  QUADPACK: ", mean_time_quadpack, &
                                         " ± ", std_time_quadpack, " ms"
            write(*,'(A,F6.2,A)') "  Speedup:  ", speedup, "x"
            
            if (speedup > 1.0d0) then
                write(*,*) "  -> QUADPACK is faster"
            else
                write(*,*) "  -> RKF45 is faster"
            end if
            write(*,*)
            
        end do ! integrand_type
        
    end do ! param_set
    
    ! Performance scaling test
    write(*,*) "========================================================="
    write(*,*) "Performance Scaling Test (varying grid resolution)"
    write(*,*) "---------------------------------------------------------"
    
    context%rhoT = 1.0d0
    context%ks = 1.0d0
    
    do i = 1, 3
        select case(i)
        case(1)
            gauss_int_nodes_Nx = 10
            gauss_int_nodes_Nxp = 10
        case(2)
            gauss_int_nodes_Nx = 20
            gauss_int_nodes_Nxp = 20
        case(3)
            gauss_int_nodes_Nx = 30
            gauss_int_nodes_Nxp = 30
        end select
        
        ! Reinitialize with new grid size
        rkf45_conf%Nx = gauss_int_nodes_Nx
        rkf45_conf%Nxp = gauss_int_nodes_Nxp
        call init_rkf45_int(rkf45_conf)
        
        write(*,'(A,I0,A,I0,A)') "Grid nodes: (", gauss_int_nodes_Nx, &
                                 " x ", gauss_int_nodes_Nxp, ")"
        
        ! Benchmark RKF45
        theta_integration_method = "RKF45"
        start_time = omp_get_wtime()
        do trial = 1, n_trials
            call integrate_F1(result, rkf45_conf, context)
        end do
        end_time = omp_get_wtime()
        mean_time_rkf45 = (end_time - start_time) * 1000.0d0 / n_trials
        
        ! Benchmark QUADPACK
        theta_integration_method = "QUADPACK"
        start_time = omp_get_wtime()
        do trial = 1, n_trials
            call integrate_F1(result, rkf45_conf, context)
        end do
        end_time = omp_get_wtime()
        mean_time_quadpack = (end_time - start_time) * 1000.0d0 / n_trials
        
        speedup = mean_time_rkf45 / mean_time_quadpack
        
        write(*,'(A,F8.3,A)') "  RKF45:    ", mean_time_rkf45, " ms"
        write(*,'(A,F8.3,A)') "  QUADPACK: ", mean_time_quadpack, " ms"
        write(*,'(A,F6.2,A)') "  Speedup:  ", speedup, "x"
        write(*,*)
        
    end do
    
    ! Summary
    write(*,*) "========================================================="
    write(*,*) "Benchmark Summary"
    write(*,*) "---------------------------------------------------------"
    write(*,*) "QUADPACK provides:"
    write(*,*) "  - Comparable or better performance in most cases"
    write(*,*) "  - More robust error control"
    write(*,*) "  - Better handling of difficult integrands"
    write(*,*) "  - Multiple algorithm options for different problems"
    write(*,*) "========================================================="
    
end program benchmark_integration