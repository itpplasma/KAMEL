module time_evolution_stellarator

    use control_mod
    use h5mod
    use balance_base, only: balance_t
    use time_evolution
    use QLBalance_kinds, only: dp

    implicit none

    logical :: set_momentum_source_to_zero, set_Q_neo_to_zero
    integer :: update_transport_coefficients
    logical :: reduce_time_step, turn_off_heat_sources

    integer :: hyst_mod_stage = 0
    integer :: ant_fac_flattop
    real(dp) :: t_flattop_begin, t_flattop_end

    real(dp) :: hyst_mod_amp_fac, hyst_mod_freq, hyst_mod_phase, bif_crit_flattop
    real(dp) :: delta_t_flattop

    type, extends(balance_t) :: time_evolution_stellarator_t
        contains
            procedure :: init_balance => initTimeEvolution
            procedure :: run_balance => runTimeEvolution
    end type

    contains

    subroutine initTimeEvolution(this)

        use recstep_mod, only: tol
        use transp_coeffs_mod, only: rescale_transp_coeffs_by_ant_fac
        use grid_mod, only: mwind, rmax, rmin, set_boundary_condition, npoib, rb
        use baseparam_mod, only: dperp, tol_max
        use QLbalance_diag, only: write_diag, write_diag_b
        use KAMEL_hdf5_tools, only: h5overwrite
        use h5mod, only: mode_m, mode_n
        use control_mod, only: gyro_current_study, write_gyro_current, &
                        ihdf5IO
        use logger_m, only: log_debug
        use wave_code_data, only: m_vals, n_vals
        use plasma_parameters, only: write_initial_parameters, alloc_hold_parameters, &
                                params, params_begbeg, init_background_profiles
        implicit none

        class(time_evolution_stellarator_t), intent(inout) :: this
        this%runType = "TimeEvolution"

        call read_stell_config
        call print_stell_config

        iexit = 0 ! 0 - dont skip, 1 - skip, 2 - stop
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

        !call read_config
        timescale = (rmax - rmin)**2 / dperp
        tmax = timescale * tmax_factor
        timstep = tmax / Nstorage
        time = 0.0d0
        tol = tol_max

        call gengrid
        call set_boundary_condition

        call initialize_wave_code_interface(npoib, rb)

        mode_m = m_vals(1)
        mode_n = n_vals(1)
        if (ihdf5IO .eq. 1) then
            call create_group_structure_timeevol
        end if
        call log_debug('mode_m/mode_n set')

        call allocate_prev_variables
        call init_background_profiles
        call write_initial_parameters

        call alloc_Br_Dqle_for_timeevol
        call calc_geometric_parameter_profiles
        call initialize_get_dql
        call initialize_D_one_over_nu
        call initialize_antenna_factor
        call det_balance_eqs_source_terms_stell

        if (.not. suppression_mode) call write_kin_prof_data_to_disk

        call allocate_timscal_and_params
        timstep = timstep*tol
        scratch = .true.

        call reset_timstep_arr_w_timstep

        call get_dql
        call rescale_transp_coeffs_by_ant_fac
        call hold_prev_transp_coeffs

        params_begbeg = params

    end subroutine

    subroutine runTimeEvolution(this)

        use baseparam_mod, only: factolmax, factolred
        use recstep_mod, only: tol
        use plasma_parameters, only: params, params_beg, params_begbeg, limit_temps_from_below
        use restart_mod, only: redostep
        use transp_coeffs_mod, only: rescale_transp_coeffs_by_ant_fac
        use recstep_mod, only: timstep_arr
        use grid_mod, only: Ipar
        use writeData_m, only: writefort9999_stellarator

        implicit none

        class(time_evolution_stellarator_t), intent(inout) :: this
        integer :: iredo = 0

        do time_ind = 1, Nstorage

            call copy_kin_profs_to_yprev
            redostep = .false.

            if (update_transport_coefficients .eq. 1) then
                call get_D_one_over_nu
            else if (update_transport_coefficients .eq. 2) then
                call rescale_D_one_over_nu_to_new_n_and_T
            end if

            call get_dql
            Ipar_time(time_ind) = Ipar
            call rescale_transp_coeffs_by_ant_fac
            call interp_Br_Dql_at_resonance_timeevol
            call determine_Dql_diagnostic

            call write_br_dqle22_time_data
            call message_Br_Dqle_values

            if (data_verbosity >= 2) then
                call writefort9999_stellarator(dqle11_prev, dqli11_prev)
            end if

            if (.true.) then
                call hold_prev_transp_coeffs
                params_begbeg = params
            else
                call redoTimeStep
            end if

            do ! evolve in time, reduce time step if too large
                iredo = iredo + 1
                params_beg = params

                print *, ""
                call log_debug("Timstep before evolvestep")
                call evolvestep_stell(timstep, eps)

                call limit_temps_from_below

                call calc_params_num_and_denom
                call smooth_params_num_and_denom
                call determine_timscal

                if (maxval(timscal) .lt. tol * factolmax) then
                    exit
                end if

                timstep_arr = timstep_arr * factolred
                params = params_beg

                call log_debug("Redoing step: Maxval(timscal) > tol * factolmax")
                if (iredo > 300) then
                    stop "Redoing step: Maxval(timscal) is not lesser than tol * factolmax after 100 redos"
                end if
            end do

            call rescale_time_step_array
            call set_time_step
            call stop_if_time_step_too_small

            call reset_timstep_arr_w_timstep
            call write_time_info

            call relax_plasma_parameters

            timstep_arr = 0.0d0
            call evolvestep_stell(timstep, eps)
            timstep_arr = timstep
            time = time + timstep

            call log_debug("msg_time_info")
            if (.not. suppression_mode) call write_kin_profile_at_time_index
            call set_first_iteration_true
            call check_linear_discr_pen_ratio
            call stop_if_antenna_fac_max_reached

            call ramp_coil
        end do

    end subroutine

    subroutine read_stell_config

        implicit none

        NAMELIST /BALANCE_STELL/ set_momentum_source_to_zero, set_Q_neo_to_zero, &
            update_transport_coefficients, reduce_time_step, turn_off_heat_sources, &
            hyst_mod_amp_fac, hyst_mod_freq, hyst_mod_phase, bif_crit_flattop, delta_t_flattop

        ! read the parameters from namelist file
        open (22, file='stell_conf.nml');
        read (22, NML=BALANCE_STELL)
        close (22);

    end subroutine

    subroutine print_stell_config

        implicit none

        print *, "============================================================================="
        print *, "    Stellarator mode configuration"
        print *, "        set_momentum_source_to_zero           = ", set_momentum_source_to_zero
        print *, "        set_Q_neo_to_zero                     = ", set_Q_neo_to_zero
        print *, "        update_transport_coefficients         = ", update_transport_coefficients
        print *, "        reduce_time_step                      = ", reduce_time_step
        print *, "        turn_off_heat_sources                 = ", turn_off_heat_sources
        print *, "        hyst_mod_amp_fac                      = ", hyst_mod_amp_fac
        print *, "        hyst_mod_freq                         = ", hyst_mod_freq
        print *, "        hyst_mod_phase                        = ", hyst_mod_phase
        print *, "        bif_crit_flattop                      = ", bif_crit_flattop
        print *, "        delta_t_flattop                       = ", delta_t_flattop
        print *, "============================================================================="

    end subroutine


end module
