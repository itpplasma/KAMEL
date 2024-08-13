
module balance_mod

    implicit none

    contains

    subroutine fromBalanceFactoryGetBalance(typeOfRun, balanceInstance)

        use balanceBase, only: balance_t
        use singleStep, only: SingleStep_t
        use time_evolution, only: TimeEvolution_t
        use paramscan_mod, only: ParameterScan_t

        implicit none

        character(100), intent(in) :: typeOfRun
        class(balance_t), allocatable, intent(out) :: balanceInstance

        select case(trim(typeOfRun))
            case("SingleStep")
                allocate(balanceInstance, source=SingleStep_t())
            case("TimeEvolution")
                allocate(balanceInstance, source=TimeEvolution_t())
            case("ParameterScan")
                allocate(balanceInstance, source=ParameterScan_t())
            !case("TimeEvolutionParameterScan") ! TODO
            !case("ConstantPsi") ! TODO
            case default
                print *, "Invalid balance type of run " // trim(typeOfRun)
                print *, "Options are: SingleStep, TimeEvolution, ParameterScan"
                stop "Due to invalid type of run"
        end select

    end subroutine

    
end module

subroutine balanceInit

        use time_evolution, only: iexit, timescale, tmax, timstep, Nstorage, &
                                  allocate_prev_variables, tmax_factor
        use grid_mod, only: mwind, rmax, rmin, setBoundaryCondition, npoib, rb
        use baseparam_mod, only: dperp
        use diag_mod, only: write_diag, write_diag_b
        use hdf5_tools, only: h5overwrite
        use h5mod, only: mode_m, mode_n
        use control_mod, only: gyro_current_study, write_gyro_current, debug_mode, &
                          ihdf5IO
        use parallelTools, only: initMPI, irank
        use wave_code_data, only: m_vals, n_vals
        use paramscan_mod, only: initialize_parameter_scan_vars, creategroupstructure
        use singleStep, only: init_background_profiles
        use plasma_parameters, only: writeInitialParameters, alloc_hold_parameters

        implicit none

        call read_config

        iexit = 0 ! 0 - don't skip, 1 - skip, 2 - stop
        mwind = 10

        write_diag = .false.
        write_diag_b = .false.

        ! if h5overwrite = true, existing data will be deleted
        ! before new one is written
        ! This is contained in hdf5_tools module
        h5overwrite = .true.
    
        if (gyro_current_study .ne. 0) then
            write_gyro_current = .true.
        else
            write_gyro_current = .false.
        end if

        call initMPI

        timescale = (rmax - rmin)**2/dperp
        tmax = timescale*tmax_factor
        timstep = tmax/Nstorage
    
        if (irank .eq. 0) then
            write(*,*) "timstep = ", timstep
        end if

        call gengrid

        call setBoundaryCondition

        CALL initialize_wave_code_interface(npoib, rb);
        CALL initialize_parameter_scan_vars

        mode_m = m_vals(1)
        mode_n = n_vals(1)
        if (debug_mode) write(*,*) 'Debug: mode_m = ', mode_m, 'mode_n = ', mode_n

        if (ihdf5IO .eq. 1) then
            CALL creategroupstructure
        end if

        call allocate_prev_variables

        call init_background_profiles

        if (irank .eq. 0) then
            CALL writeInitialParameters
            call alloc_hold_parameters
        end if

    end subroutine

