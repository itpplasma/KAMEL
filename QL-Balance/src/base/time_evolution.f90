module time_evolution

    use control_mod
    use h5mod
    use balance_base, only: balance_t
    use QLBalance_kinds, only: dp
    use logger_m, only: log_debug, log_info
    use periodic_amplitude_state_m, only: periodic_amplitudes

    implicit none

    logical :: br_stopping ! trigger Br stopping criterion
    logical :: discr_reached = .false. ! variable to say if discrepancy to linear regression
    logical :: scratch

    integer :: Nstorage
    integer :: ramp_up_mode !> control ramp up mode of the RMP coil current amplitude
    integer :: save_prof_time_step
    integer :: iexit ! used for ramp-up skipping of saving
    integer :: ramp_up_down = 0 !> used in hysteresis mode, tells if ramp-up (0) or ramp-down (1)
    logical :: periodic_response_trial = .false.
    integer :: last_periodic_rejections = 0
    real(dp) :: accepted_current_floor = 1.0e-30_dp, accepted_max_scale = 1.0e12_dp
    real(dp) :: accepted_antenna_factor = 1.0_dp
    real(dp), allocatable :: periodic_transport_unscaled(:,:)
    integer :: periodic_completed_steps = 0
    integer :: time_ind

    real(dp) :: tmax, timescale
    real(dp) :: br_beta = 0
    real(dp) :: br_predicted = 0.0_dp

    real(dp) :: tmax_factor!, antenna_factor
    real(dp) :: stop_time_step
    real(dp) :: timstep_min
    real(dp) :: t_max_ramp_up = 1e-2 !> 10ms ramp up until antenna_factor_max is reached
    real(dp) :: timstep
    real(dp) :: time
    real(dp) :: t_hysteresis_turn = 0
    real(dp) :: constant_time_step = 0.0_dp
    logical :: set_constant_time_step = .false.

    real(dp), dimension(:), allocatable :: timscal

    real(dp), DIMENSION(:), ALLOCATABLE :: yprev
    real(dp), DIMENSION(:), ALLOCATABLE :: dqle11_prev
    real(dp), DIMENSION(:), ALLOCATABLE :: dqle12_prev
    real(dp), DIMENSION(:), ALLOCATABLE :: dqle21_prev
    real(dp), DIMENSION(:), ALLOCATABLE :: dqle22_prev
    real(dp), DIMENSION(:), ALLOCATABLE :: dqli11_prev
    real(dp), DIMENSION(:), ALLOCATABLE :: dqli12_prev
    real(dp), DIMENSION(:), ALLOCATABLE :: dqli21_prev
    real(dp), DIMENSION(:), ALLOCATABLE :: dqli22_prev

    real(dp) :: antenna_factor_max
    real(dp) :: antenna_max_stopping

    !needed for interpolation of br abs and stopping criterion
    real(dp), DIMENSION(:), ALLOCATABLE :: br_abs
    complex(dp), DIMENSION(:), ALLOCATABLE :: br_formfactor
    complex(dp), DIMENSION(:), ALLOCATABLE :: br_vac_res
    real(dp), DIMENSION(:), ALLOCATABLE :: br_abs_time
    real(dp), DIMENSION(:), ALLOCATABLE :: br_abs_antenna_factor
    real(dp), DIMENSION(:), ALLOCATABLE :: dqle22_res_time
    real(dp), DIMENSION(:), ALLOCATABLE :: dae22_res_time
    real(dp), DIMENSION(:), ALLOCATABLE :: bif_criterion
    complex(dp), dimension(:), allocatable :: Ipar_time

    logical :: firstiterationdone = .false. !Some steps in saving the data to hdf5 file
    !need to be done only the first time iteration

    integer(HID_T) :: time_dataset_id !> variable to save the time dataset id

    !integer :: time_ind

    type, extends(balance_t) :: TimeEvolution_t
        contains
            procedure :: init_balance => initTimeEvolution
            procedure :: run_balance => runTimeEvolution
    end type

    private :: initTimeEvolution
    private :: runTimeEvolution

    contains

    subroutine initTimeEvolution(this)

        use recstep_mod, only: tol
        use periodic_checkpoint_m, only: pending_periodic_restart, reset_periodic_restart
        use kim_wave_code_adapter_m, only: kim_periodic_mode_selected
        use transp_coeffs_mod, only: rescale_transp_coeffs_by_ant_fac, &
            compute_antenna_factor_from_Ipar
        use grid_mod, only: mwind, rmax, rmin, set_boundary_condition, npoib, rb
        use baseparam_mod, only: dperp, tol_max
        use QLbalance_diag, only: write_diag, write_diag_b
        use KAMEL_hdf5_tools, only: h5overwrite
        use h5mod, only: mode_m, mode_n
        use control_mod, only: gyro_current_study, write_gyro_current, &
                        ihdf5IO
        use wave_code_data, only: m_vals, n_vals
        use plasma_parameters, only: write_initial_parameters, alloc_hold_parameters, &
                                params, params_begbeg, init_background_profiles
        use resonances_mod, only: write_resonant_radii_to_hdf5
        use logger_m, only: log_debug

        implicit none

        class(TimeEvolution_t), intent(inout) :: this
        this%runType = "TimeEvolution"

        periodic_completed_steps = 0
        time_ind = 0
        if (readfromtimestep == 0) then
            call periodic_amplitudes%reset()
            call reset_periodic_restart()
        end if
        iexit = 0 ! 0 - dont skip, 1 - skip, 2 - stop
        mwind = 10
        write_diag = .false.
        write_diag_b = .false.

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
        call write_periodic_workflow_provenance
        call log_debug('mode_m/mode_n set')

        call write_resonant_radii_to_hdf5

        call allocate_prev_variables
        call init_background_profiles
        call write_initial_parameters
        !call alloc_hold_parameters

        call calc_geometric_parameter_profiles
        call initialize_get_dql
        call initialize_antenna_factor
        call det_balance_eqs_source_terms

        ! Preserve the established step-zero profile output. There is no
        ! accepted evolution step to checkpoint until the first commit.
        if (.not. suppression_mode) call write_kin_prof_data_to_disk

        call allocate_timscal_and_params
        timstep = timstep*tol
        scratch = .true.

        call reset_timstep_arr_w_timstep
        !timstep_arr = timstep
        !tim_stack = timstep_arr

        call alloc_Br_Dqle_for_timeevol
        if (trim(wave_code) == 'KIM' .and. kim_periodic_mode_selected()) then
            if (pending_periodic_restart%active) then
                call apply_periodic_restart(pending_periodic_restart)
                call restore_periodic_accepted_response(continuation=.true.)
            else
                call initialize_periodic_response()
            end if
        else
            call get_dql
            call compute_antenna_factor_from_Ipar
            call rescale_transp_coeffs_by_ant_fac
        end if
        call hold_prev_transp_coeffs

        params_begbeg = params

    end subroutine

    subroutine runTimeEvolution(this)
        class(TimeEvolution_t), intent(inout) :: this

        block
            integer :: first_step, last_step
            first_step = periodic_completed_steps+1
            last_step = periodic_completed_steps+Nstorage
            do time_ind = first_step, last_step
                call doStep(this)
            end do
        end block
    end subroutine runTimeEvolution

    subroutine doStep(this, stat)
        use baseparam_mod, only: factolmax, factolred
        use kim_wave_code_adapter_m, only: kim_periodic_mode_selected
        use plasma_parameters, only: params, params_beg, params_begbeg, limit_temps_from_below
        use logger_m, only: log_debug
        use recstep_mod, only: timstep_arr, tol
        use restart_mod, only: redostep
        use transp_coeffs_mod, only: rescale_transp_coeffs_by_ant_fac, &
            compute_antenna_factor_from_Ipar
        use writeData_m, only: writefort9999

        implicit none

        class(TimeEvolution_t), intent(inout) :: this

        integer :: iredo
        integer, intent(out), optional :: stat

        if (present(stat)) stat = 0
        if (trim(wave_code) == 'KIM' .and. kim_periodic_mode_selected()) then
            call do_periodic_step(this, stat)
            return
        end if

        call copy_kin_profs_to_yprev
        redostep = .false.

        call get_dql
        call compute_antenna_factor_from_Ipar
        call rescale_transp_coeffs_by_ant_fac
        call interp_Br_Dql_at_resonance_timeevol
        call determine_Dql_diagnostic

        call write_br_dqle22_time_data
        call message_Br_Dqle_values

        if (data_verbosity >= 2) then
            call writefort9999(dqle11_prev, dqli11_prev)
        end if

        if (.true.) then
            call hold_prev_transp_coeffs
            params_begbeg = params
        else
            call redoTimeStep
        end if

        iredo = 0
        do ! redo step loop
            iredo = iredo + 1
            params_beg = params

            print *, ""
            call log_debug("Timstep before evolvestep")

            call evolvestep(timstep, eps)

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
            if (iredo > 100) then
                stop "Redoing step: Maxval(timscal) is not lesser than tol * factolmax " // &
                        "after 100 redos"
            end if
        end do

        call rescale_time_step_array
        call set_time_step
        call stop_if_time_step_too_small
        call reset_timstep_arr_w_timstep
        call write_time_info
        call relax_plasma_parameters

        timstep_arr = 0.0d0
        call evolvestep(timstep, eps)
        timstep_arr = timstep
        time = time + timstep

        call log_debug("msg_time_info")
        if (.not. suppression_mode) call write_kin_profile_at_time_index
        call set_first_iteration_true
        call calculate_total_toroidal_torque(time_ind)
        call write_total_toroidal_torque_to_file(time_ind)
        call check_linear_discr_pen_ratio
        call stop_if_antenna_fac_max_reached

        call ramp_coil
    end subroutine doStep

    subroutine sync_periodic_amplitude_trial()
        use kim_wave_code_adapter_m, only: kim_current_records, kim_periodic_scale_modes, &
            kim_periodic_current_unit, kim_periodic_scale_status
        use wave_code_data, only: I_par_toroidal
        complex(dp), allocatable :: residual(:)
        integer :: i

        if (.not. allocated(kim_periodic_scale_modes)) &
            error stop 'periodic response is missing during synchronization'
        allocate(residual(size(kim_periodic_scale_modes)))
        do i = 1, size(residual)
            residual(i) = kim_current_records(i)%residual
        end do
        call periodic_amplitudes%begin_trial(kim_periodic_scale_modes, &
            kim_periodic_current_unit, residual, kim_periodic_scale_status, &
            I_par_toroidal, kim_current_relaxation)
    end subroutine sync_periodic_amplitude_trial

    subroutine cache_periodic_transport()
        use grid_mod, only: dqle11, dqle12, dqle21, dqle22, &
            dqli11, dqli12, dqli21, dqli22
        if (allocated(periodic_transport_unscaled)) deallocate(periodic_transport_unscaled)
        allocate(periodic_transport_unscaled(8,size(dqle11)))
        periodic_transport_unscaled(1,:) = dqle11
        periodic_transport_unscaled(2,:) = dqle12
        periodic_transport_unscaled(3,:) = dqle21
        periodic_transport_unscaled(4,:) = dqle22
        periodic_transport_unscaled(5,:) = dqli11
        periodic_transport_unscaled(6,:) = dqli12
        periodic_transport_unscaled(7,:) = dqli21
        periodic_transport_unscaled(8,:) = dqli22
    end subroutine cache_periodic_transport

    subroutine apply_periodic_transport()
        use grid_mod, only: dqle11, dqle12, dqle21, dqle22, &
            dqli11, dqli12, dqli21, dqli22
        use wave_code_data, only: antenna_factor
        use transp_coeffs_mod, only: compute_antenna_factor_from_Ipar
        if (.not. allocated(periodic_transport_unscaled)) &
            error stop 'missing accepted periodic transport'
        call compute_antenna_factor_from_Ipar
        dqle11 = periodic_transport_unscaled(1,:)*antenna_factor
        dqle12 = periodic_transport_unscaled(2,:)*antenna_factor
        dqle21 = periodic_transport_unscaled(3,:)*antenna_factor
        dqle22 = periodic_transport_unscaled(4,:)*antenna_factor
        dqli11 = periodic_transport_unscaled(5,:)*antenna_factor
        dqli12 = periodic_transport_unscaled(6,:)*antenna_factor
        dqli21 = periodic_transport_unscaled(7,:)*antenna_factor
        dqli22 = periodic_transport_unscaled(8,:)*antenna_factor
    end subroutine apply_periodic_transport

    subroutine initialize_periodic_response(stat)
        use QLbalance_diag, only: timscal_dql, timscal_dqli, rate_dql
        use wave_code_data, only: antenna_factor
        use transp_coeffs_mod, only: compute_antenna_factor_from_Ipar, &
            rescale_transp_coeffs_by_ant_fac
        integer, intent(out), optional :: stat
        integer :: ierr
        periodic_response_trial = .true.
        call get_dql
        periodic_response_trial = .false.
        call sync_periodic_amplitude_trial()
        call periodic_amplitudes%accept(stat=ierr)
        if (ierr == 0) then
            call publish_periodic_response()
            call cache_periodic_transport()
            call apply_periodic_transport()
            timscal_dql = 0.0_dp
            timscal_dqli = 0.0_dp
            rate_dql = 0.0_dp
            accepted_current_floor = kim_current_floor
            accepted_max_scale = kim_current_max_scale
            accepted_antenna_factor = antenna_factor
        end if
        if (present(stat)) then
            stat = ierr
        elseif (ierr /= 0) then
            error stop 'initial periodic response failed normalization guards'
        end if
    end subroutine initialize_periodic_response

    subroutine restore_periodic_accepted_response(stat, continuation)
        use kim_wave_code_adapter_m, only: kim_use_accepted_amplitudes, &
            kim_periodic_scale_status
        use wave_code_data, only: I_par_toroidal, antenna_factor
        use transp_coeffs_mod, only: compute_antenna_factor_from_Ipar, &
            rescale_transp_coeffs_by_ant_fac
        integer, intent(out), optional :: stat
        logical, intent(in), optional :: continuation
        real(dp) :: next_target, next_relaxation, next_floor, next_max, next_antenna
        integer :: ierr
        next_target = I_par_toroidal
        next_relaxation = kim_current_relaxation
        next_floor = kim_current_floor
        next_max = kim_current_max_scale
        next_antenna = antenna_factor
        call periodic_amplitudes%reject()
        I_par_toroidal = periodic_amplitudes%accepted_target_current
        kim_current_relaxation = periodic_amplitudes%accepted_relaxation
        kim_current_floor = accepted_current_floor
        kim_current_max_scale = accepted_max_scale
        antenna_factor = accepted_antenna_factor
        call kim_use_accepted_amplitudes()
        ! Reconstruct all fields, currents, background caches and transport
        ! from the restored accepted profiles and prescribed accepted scales.
        periodic_response_trial = .true.
        call get_dql
        periodic_response_trial = .false.
        ierr = 0
        if (any(kim_periodic_scale_status > 0)) ierr = 1
        if (present(continuation)) then
            if (continuation) then
                I_par_toroidal = next_target
                kim_current_relaxation = next_relaxation
                kim_current_floor = next_floor
                kim_current_max_scale = next_max
                antenna_factor = next_antenna
            end if
        end if
        call cache_periodic_transport()
        call apply_periodic_transport()
        if (present(stat)) then
            stat = ierr
        elseif (ierr /= 0) then
            error stop 'could not reconstruct accepted periodic response'
        end if
    end subroutine restore_periodic_accepted_response

    subroutine do_periodic_step(this, stat)
        use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
        use baseparam_mod, only: factolmax, factolred
        use plasma_parameters, only: params, params_beg, params_begbeg, &
            limit_temps_from_below
        use recstep_mod, only: timstep_arr, tim_stack, tol
        use restart_mod, only: redostep
        use wave_code_data, only: antenna_factor
        use transp_coeffs_mod, only: compute_antenna_factor_from_Ipar, &
            rescale_transp_coeffs_by_ant_fac
        class(TimeEvolution_t), intent(inout) :: this
        integer, intent(out), optional :: stat
        real(dp), allocatable :: accepted_profiles(:,:), old_steps(:), old_stack(:)
        real(dp) :: old_dt, old_time, trial_dt
        logical :: old_scratch
        integer :: attempt, ierr, restore_status

        if (present(stat)) stat = 0
        last_periodic_rejections = 0
        if (.not. periodic_amplitudes%initialized) then
            call initialize_periodic_response(ierr)
            if (ierr /= 0) then
                call failed_step(ierr)
                return
            end if
        end if
        call apply_periodic_transport()
        accepted_antenna_factor = antenna_factor
        accepted_profiles = params
        old_steps = timstep_arr
        old_stack = tim_stack
        old_dt = timstep
        old_time = time
        old_scratch = scratch
        trial_dt = timstep
        call copy_kin_profs_to_yprev
        params_begbeg = params
        call hold_prev_transp_coeffs
        redostep = .false.

        ! The accepted response already belongs to these profiles. In
        ! particular, do not replace the initial unit drive before advancing.
        do attempt = 1, 100
            params = accepted_profiles
            params_beg = params
            timstep = trial_dt
            timstep_arr = trial_dt
            call evolvestep(timstep, eps)
            call limit_temps_from_below
            call calc_params_num_and_denom
            call smooth_params_num_and_denom
            call determine_timscal
            if (.not. all(ieee_is_finite(params)) .or. any(params(1,:) <= 0.0_dp) .or. &
                .not. all(ieee_is_finite(timscal)) .or. &
                maxval(timscal) >= tol*factolmax) then
                last_periodic_rejections = last_periodic_rejections + 1
                params = accepted_profiles
                call periodic_amplitudes%reject()
                trial_dt = trial_dt*factolred
                if (trial_dt < timstep_min) exit
                cycle
            end if

            call rescale_time_step_array
            call set_time_step
            call reset_timstep_arr_w_timstep
            call relax_plasma_parameters
            timstep_arr = 0.0_dp
            call evolvestep(timstep, eps)
            timstep_arr = timstep

            ! Stage the response for the candidate profiles before either is
            ! committed. A guarded response cannot replace the accepted one.
            periodic_response_trial = .true.
            call get_dql
            periodic_response_trial = .false.
            call sync_periodic_amplitude_trial()
            call periodic_amplitudes%accept(stat=ierr)
            if (ierr /= 0) then
                last_periodic_rejections = last_periodic_rejections + 1
                exit
            end if
            call publish_periodic_response()
            call cache_periodic_transport()
            call apply_periodic_transport()
            accepted_current_floor = kim_current_floor
            accepted_max_scale = kim_current_max_scale
            accepted_antenna_factor = antenna_factor
            periodic_completed_steps = time_ind
            time = old_time + trial_dt
            call write_time_info
            call interp_Br_Dql_at_resonance_timeevol
            call determine_Dql_diagnostic
            call write_br_dqle22_time_data
            call message_Br_Dqle_values
            call set_first_iteration_true
            call calculate_total_toroidal_torque(time_ind)
            call write_total_toroidal_torque_to_file(time_ind)
            call check_linear_discr_pen_ratio
            call stop_if_antenna_fac_max_reached
            call ramp_coil
            if (.not. suppression_mode) call write_kin_profile_at_time_index
            call stop_if_time_step_too_small
            return
        end do

        params = accepted_profiles
        params_beg = accepted_profiles
        params_begbeg = accepted_profiles
        time = old_time
        timstep = old_dt
        timstep_arr = old_steps
        tim_stack = old_stack
        scratch = old_scratch
        call restore_periodic_accepted_response(restore_status)
        if (restore_status /= 0) error stop 'failed to restore accepted periodic response'
        call hold_prev_transp_coeffs
        call failed_step(1)
    contains
        subroutine failed_step(code)
            integer, intent(in) :: code
            if (present(stat)) then
                stat = code
            else
                error stop 'periodic profile step rejected; accepted state restored'
            end if
        end subroutine failed_step
    end subroutine do_periodic_step

    subroutine publish_periodic_response()
        use writeData_m, only: write_fields_currs_transp_coefs_to_h5, &
            write_D_one_over_nu_to_h5
        use periodic_transport_benchmark_m, only: write_transport_benchmark
        use periodic_current_diagnostics_m, only: write_periodic_current_diagnostics
        if (suppression_mode) return
        if (modulo(time_ind, save_prof_time_step) /= 0) return
        call write_fields_currs_transp_coefs_to_h5(time_ind)
        call write_transport_benchmark(time_ind)
        call write_periodic_current_diagnostics(time_ind)
        call write_D_one_over_nu_to_h5(time_ind)
    end subroutine publish_periodic_response

    subroutine capture_periodic_restart(payload)
        use periodic_checkpoint_m, only: periodic_restart_t
        use grid_mod, only: rc, rb, Ercov, dery_equisource
        use QLbalance_diag, only: timscal_dql, timscal_dqli, rate_dql
        use plasma_parameters, only: params
        use recstep_mod, only: timstep_arr, tim_stack, y_stack, nstack, tol
        use wave_code_data, only: q, dPhi0, Vth, I_par_toroidal, antenna_factor
        use kim_wave_code_adapter_m, only: kim_get_boundary_drive
        type(periodic_restart_t), intent(out) :: payload
        interface
            subroutine periodic_ramp_checkpoint(payload, restore)
                use periodic_checkpoint_m, only: periodic_restart_t
                type(periodic_restart_t), intent(inout) :: payload
                logical, intent(in) :: restore
            end subroutine
        end interface
        payload%active = .true.
        payload%accepted_step = periodic_completed_steps
        payload%timscal_dql = timscal_dql
        payload%timscal_dqli = timscal_dqli
        payload%rate_dql = rate_dql
        payload%dery_equisource = dery_equisource
        payload%accepted_floor = accepted_current_floor
        payload%accepted_max_scale = accepted_max_scale
        payload%time = time
        payload%next_dt = timstep
        payload%tolerance = tol
        payload%ramp_mode = ramp_up_mode
        payload%ramp_up_down = ramp_up_down
        payload%iexit = iexit
        payload%br_beta = br_beta
        payload%br_predicted = br_predicted
        payload%t_hysteresis_turn = t_hysteresis_turn
        payload%antenna_factor = antenna_factor
        payload%antenna_factor_max = antenna_factor_max
        payload%antenna_max_stopping = antenna_max_stopping
        payload%t_max_ramp_up = t_max_ramp_up
        payload%target_current = I_par_toroidal
        payload%current_floor = kim_current_floor
        payload%max_scale = kim_current_max_scale
        payload%relaxation = kim_current_relaxation
        payload%constant_dt = constant_time_step
        payload%constant_dt_enabled = set_constant_time_step
        payload%br_stopping = br_stopping
        payload%discr_reached = discr_reached
        payload%scratch = scratch
        payload%tmax = tmax
        payload%timescale = timescale
        payload%tmax_factor = tmax_factor
        payload%stop_time_step = stop_time_step
        payload%timstep_min = timstep_min
        payload%params = params
        payload%rc = rc
        payload%rb = rb
        payload%q = q
        payload%er = -dPhi0
        payload%vth = Vth
        payload%timstep_arr = timstep_arr
        payload%tim_stack = tim_stack
        if (allocated(y_stack)) then
            payload%y_stack = y_stack
            payload%nstack = nstack
        end if
        payload%boundary_br = kim_get_boundary_drive()
        payload%br_abs = br_abs(1:periodic_completed_steps)
        payload%br_abs_time = br_abs_time(1:periodic_completed_steps)
        payload%br_abs_antenna_factor = br_abs_antenna_factor(1:periodic_completed_steps)
        payload%dqle22_res_time = dqle22_res_time(1:periodic_completed_steps)
        payload%dae22_res_time = dae22_res_time(1:periodic_completed_steps)
        payload%bif_criterion = bif_criterion(1:periodic_completed_steps)
        payload%Ipar_time = Ipar_time(1:periodic_completed_steps)
        payload%br_formfactor = br_formfactor(1:periodic_completed_steps)
        payload%br_vac_res = br_vac_res(1:periodic_completed_steps)
        call periodic_ramp_checkpoint(payload, .false.)
    end subroutine capture_periodic_restart

    subroutine apply_periodic_restart(payload)
        use periodic_checkpoint_m, only: periodic_restart_t
        use grid_mod, only: rc, rb, Ercov, dery_equisource
        use QLbalance_diag, only: timscal_dql, timscal_dqli, rate_dql
        use plasma_parameters, only: params
        use recstep_mod, only: timstep_arr, tim_stack, y_stack, nstack, tol
        use wave_code_data, only: q, dPhi0, Vth, I_par_toroidal, antenna_factor
        use kim_wave_code_adapter_m, only: kim_set_boundary_drive
        type(periodic_restart_t), intent(in) :: payload
        type(periodic_restart_t) :: ramp_payload
        integer :: count
        interface
            subroutine periodic_ramp_checkpoint(payload, restore)
                use periodic_checkpoint_m, only: periodic_restart_t
                type(periodic_restart_t), intent(inout) :: payload
                logical, intent(in) :: restore
            end subroutine
        end interface
        if (.not. payload%active) error stop 'cannot apply an inactive periodic checkpoint'
        if (size(rc) /= size(payload%rc) .or. size(rb) /= size(payload%rb)) &
            error stop 'periodic restart grid size changed'
        if (any(rc /= payload%rc) .or. any(rb /= payload%rb)) &
            error stop 'periodic restart grid changed'
        if (any(shape(params) /= shape(payload%params))) &
            error stop 'periodic restart profile shape changed'
        params = payload%params
        q = payload%q
        dPhi0 = -payload%er
        Ercov = payload%er
        Vth = payload%vth
        timstep_arr = payload%timstep_arr
        if (allocated(payload%tim_stack)) tim_stack = payload%tim_stack
        if (allocated(y_stack)) deallocate(y_stack)
        if (allocated(payload%y_stack)) y_stack = payload%y_stack
        nstack = payload%nstack
        periodic_completed_steps = payload%accepted_step
        timscal_dql = payload%timscal_dql
        timscal_dqli = payload%timscal_dqli
        rate_dql = payload%rate_dql
        dery_equisource = payload%dery_equisource
        time = payload%time
        timstep = payload%next_dt
        tol = payload%tolerance
        ramp_up_mode = payload%ramp_mode
        ramp_up_down = payload%ramp_up_down
        iexit = payload%iexit
        br_beta = payload%br_beta
        br_predicted = payload%br_predicted
        t_hysteresis_turn = payload%t_hysteresis_turn
        antenna_factor = payload%antenna_factor
        antenna_factor_max = payload%antenna_factor_max
        antenna_max_stopping = payload%antenna_max_stopping
        t_max_ramp_up = payload%t_max_ramp_up
        I_par_toroidal = payload%target_current
        kim_current_floor = payload%current_floor
        kim_current_max_scale = payload%max_scale
        kim_current_relaxation = payload%relaxation
        constant_time_step = payload%constant_dt
        set_constant_time_step = payload%constant_dt_enabled
        br_stopping = payload%br_stopping
        discr_reached = payload%discr_reached
        scratch = payload%scratch
        tmax = payload%tmax
        timescale = payload%timescale
        tmax_factor = payload%tmax_factor
        stop_time_step = payload%stop_time_step
        timstep_min = payload%timstep_min
        time_ind = periodic_completed_steps
        call kim_set_boundary_drive(payload%boundary_br)
        accepted_current_floor = payload%accepted_floor
        accepted_max_scale = payload%accepted_max_scale
        accepted_antenna_factor = antenna_factor
        count = periodic_completed_steps
        if (size(br_abs) < count) error stop 'periodic restart diagnostic history is too short'
        br_abs(1:count) = payload%br_abs
        br_abs_time(1:count) = payload%br_abs_time
        br_abs_antenna_factor(1:count) = payload%br_abs_antenna_factor
        dqle22_res_time(1:count) = payload%dqle22_res_time
        dae22_res_time(1:count) = payload%dae22_res_time
        bif_criterion(1:count) = payload%bif_criterion
        Ipar_time(1:count) = payload%Ipar_time
        br_formfactor(1:count) = payload%br_formfactor
        br_vac_res(1:count) = payload%br_vac_res
        ramp_payload = payload
        call periodic_ramp_checkpoint(ramp_payload, .true.)
    end subroutine apply_periodic_restart

    subroutine write_periodic_workflow_provenance
        use control_mod, only: wave_code, kim_run_type, kim_profiles_from_balance, kim_config_path, &
            kim_n_modes, kim_m_list, kim_n_list, kim_electron_transport_model, kim_ion_transport_model, &
            kim_bparallel_source, kim_benchmark_mode
        use wave_code_data, only: I_par_toroidal
        use control_mod, only: ihdf5IO
        use h5mod, only: h5_id, h5_mode_groupname, path2out
        use KAMEL_hdf5_tools, only: h5_init, h5_open_rw, h5_close, h5_deinit, h5_add_string, &
            h5_add_double_1, h5_create_parent_groups
        use periodic_amplitude_state_m, only: periodic_normalization_version, periodic_phase_policy
        real(dp), allocatable :: modes_m(:), modes_n(:)
        character(len=1024) :: group

        if (trim(wave_code) /= 'KIM' .or. trim(kim_run_type) /= 'electrostatic_periodic') return
        if (ihdf5IO /= 1) return
        allocate(modes_m(kim_n_modes), modes_n(kim_n_modes))
        modes_m = real(kim_m_list(1:kim_n_modes), dp)
        modes_n = real(kim_n_list(1:kim_n_modes), dp)
        group = "/"//trim(h5_mode_groupname)//"/periodic_workflow"
        call h5_init()
        call h5_open_rw(path2out, h5_id)
        call h5_create_parent_groups(h5_id, trim(group)//"/")
        call h5_add_string(h5_id, trim(group)//"/wave_code", trim(wave_code))
        call h5_add_string(h5_id, trim(group)//"/kim_run_type", trim(kim_run_type))
        call h5_add_string(h5_id, trim(group)//"/electron_transport_model", trim(kim_electron_transport_model))
        call h5_add_string(h5_id, trim(group)//"/ion_transport_model", trim(kim_ion_transport_model))
        call h5_add_string(h5_id, trim(group)//"/bparallel_source", trim(kim_bparallel_source))
        call h5_add_string(h5_id, trim(group)//"/benchmark_mode", trim(kim_benchmark_mode))
        call h5_add_string(h5_id, trim(group)//"/phase_policy", periodic_phase_policy)
        call h5_add_string(h5_id, trim(group)//"/fourier_convention", "exp(i*(kr*r + ell*theta - omega*t))")
        call h5_add_string(h5_id, trim(group)//"/field_order", "Phi,Br,Bparallel")
        call h5_add_string(h5_id, trim(group)//"/transition_contract", "compact-C1-common-window")
        call h5_add_string(h5_id, trim(group)//"/kim_config_path", trim(kim_config_path))
        call h5_add_string(h5_id, trim(group)//"/algebra_generator_sha256", &
            "a7591175092dd15b54ddd0eaf294f990f3441a90f3bcd1f0459092e4bf36891e")
        call h5_add_double_1(h5_id, trim(group)//"/mode_m", modes_m, (/1/), (/kim_n_modes/))
        call h5_add_double_1(h5_id, trim(group)//"/mode_n", modes_n, (/1/), (/kim_n_modes/))
        call h5_add_double_1(h5_id, trim(group)//"/target_current", [I_par_toroidal], (/1/), (/1/))
        call h5_add_double_1(h5_id, trim(group)//"/normalization_version", &
            [real(periodic_normalization_version, dp)], (/1/), (/1/))
        call h5_add_double_1(h5_id, trim(group)//"/profiles_from_balance", &
            [merge(1.0_dp, 0.0_dp, kim_profiles_from_balance)], (/1/), (/1/))
        call h5_close(h5_id)
        call h5_deinit()
        deallocate(modes_m, modes_n)
    end subroutine write_periodic_workflow_provenance

    subroutine allocate_prev_variables

        use recstep_mod, only: timstep_arr, tim_stack
        use grid_mod, only: neqset, npoib

        implicit none

        allocate (yprev(neqset))
        allocate (dqle11_prev(npoib))
        allocate (dqle12_prev(npoib))
        allocate (dqle21_prev(npoib))
        allocate (dqle22_prev(npoib))
        allocate (dqli11_prev(npoib))
        allocate (dqli12_prev(npoib))
        allocate (dqli21_prev(npoib))
        allocate (dqli22_prev(npoib))
        allocate (timstep_arr(neqset), tim_stack(neqset))

    end subroutine

    subroutine copy_kin_profs_to_yprev

        use grid_mod, only: npoi, nbaleqs
        use plasma_parameters, only: params

        implicit none

        integer :: ipoi, ieq, k

        call log_debug('yprev loop')
        do ipoi = 1, npoi
            do ieq = 1, nbaleqs
                k = nbaleqs*(ipoi - 1) + ieq
                yprev(k) = params(ieq, ipoi)
            end do
        end do

    end subroutine

    subroutine initialize_antenna_factor

        use wave_code_data, only: antenna_factor
        implicit none

        antenna_factor_max = antenna_factor
        if (ramp_up_mode .eq. 4) then
            antenna_factor = 0d0
        else
            antenna_factor = 1.d-4
        end if

    end subroutine

    subroutine alloc_Br_Dqle_for_timeevol

        use grid_mod, only: T_tot_phi_e, T_tot_phi_i

        use periodic_checkpoint_m, only: pending_periodic_restart
        implicit none
        integer :: history_size

        history_size = Nstorage
        if (pending_periodic_restart%active) &
            history_size = history_size+pending_periodic_restart%accepted_step
        allocate(br_abs(history_size))
        allocate(br_formfactor(history_size))
        allocate(br_vac_res(history_size))
        allocate(br_abs_antenna_factor(history_size))
        allocate(br_abs_time(history_size))
        allocate(dqle22_res_time(history_size))
        allocate(dae22_res_time(history_size))
        allocate(bif_criterion(history_size))
        allocate(Ipar_time(history_size))
        allocate(T_tot_phi_e(history_size))
        allocate(T_tot_phi_i(history_size))

    end subroutine

    subroutine hold_prev_transp_coeffs

        use grid_mod, only: dqle11, dqle12, dqle21, dqle22, &
                            dqli11, dqli12, dqli21, dqli22
        implicit none

        dqle11_prev = dqle11
        dqle12_prev = dqle12
        dqle21_prev = dqle21
        dqle22_prev = dqle22
        dqli11_prev = dqli11
        dqli12_prev = dqli12
        dqli21_prev = dqli21
        dqli22_prev = dqli22

    end subroutine

    subroutine allocate_timscal_and_params

        use grid_mod, only: npoi, npoic, nbaleqs, dummy
        use plasma_parameters, only: params_beg, params_num, params_denom, params_begbeg

        implicit none

        allocate(timscal(npoi), dummy(npoic))
        allocate(params_beg(nbaleqs, npoic), params_num(nbaleqs, npoic))
        allocate(params_denom(nbaleqs, npoic))
        allocate(params_begbeg(nbaleqs, npoic))

    end subroutine


    !> @brief subroutine write_br_dqle22_time_data. Writes radial magnetic field perturbation evaluated at the resonant
    !> surface, the antenna factor, the time and Dqle22 evaluated at the resonant surface for a given
    !> time step to the hdf5 file
    !> @param[in] i Integer of time step to which the data will be saved. Goes from 1:i.
    !> @param[in] br_abs_time Time value of the time evolution.
    !> @param[in] br_abs_antenna_factor Value of the antenna factor, i.e. the RMP coil current.
    !> @param[in] br_abs Absolute value of the radial magnetic field evaluated at the resonant surface in question.
    !> @param[in] dqle22_res_time Value of Dqle22 evaluated at the resonant surface during the time evolution.
    subroutine write_br_dqle22_time_data

        use control_mod
        use baseparam_mod
        use h5mod
        use KAMEL_hdf5_tools
        use wave_code_data, only: antenna_factor
        use resonances_mod, only: r_res

        implicit none

        call log_debug("writing out br time evolution data")

        if (ihdf5IO .eq. 1) then
        !if (.false.) then
            CALL h5_init()
            CALL h5_open_rw(path2out, h5_id)

            h5overwrite = .true.

            h5_currentgrp = "/"//trim(h5_mode_groupname) //"/br_abs_time"
            CALL h5_add_double_1(h5_id, trim(h5_currentgrp), br_abs_time(1:time_ind), &
                lbound(br_abs_time(1:time_ind)), ubound(br_abs_time(1:time_ind)))

            h5_currentgrp = "/"//trim(h5_mode_groupname) //"/br_vac_res"
            CALL h5_add_double_1(h5_id, trim(h5_currentgrp), abs(br_vac_res(1:time_ind)), &
                lbound(br_vac_res(1:time_ind)), ubound(br_vac_res(1:time_ind)))

            h5_currentgrp = "/"//trim(h5_mode_groupname) //"/br_abs_antenna_factor"
            CALL h5_add_double_1(h5_id, trim(h5_currentgrp), br_abs_antenna_factor(1:time_ind), &
                lbound(br_abs_antenna_factor(1:time_ind)), ubound(br_abs_antenna_factor(1:time_ind)))

            h5_currentgrp = "/"//trim(h5_mode_groupname) //"/br_abs_res"
            CALL h5_add_double_1(h5_id, trim(h5_currentgrp), br_abs(1:time_ind), &
                lbound(br_abs(1:time_ind)), ubound(br_abs(1:time_ind)))

            h5_currentgrp = "/"//trim(h5_mode_groupname) //"/dqle22_res_time"
            CALL h5_add_double_1(h5_id, trim(h5_currentgrp), dqle22_res_time(1:time_ind), &
                lbound(dqle22_res_time(1:time_ind)), ubound(dqle22_res_time(1:time_ind)))

            h5_currentgrp = "/"//trim(h5_mode_groupname) //"/dae22_res_time"
            CALL h5_add_double_1(h5_id, trim(h5_currentgrp), dae22_res_time(1:time_ind), &
                lbound(dae22_res_time(1:time_ind)), ubound(dae22_res_time(1:time_ind)))

            h5_currentgrp = "/"//trim(h5_mode_groupname) //"/bifurcation_criterion"
            CALL h5_add_double_1(h5_id, trim(h5_currentgrp), bif_criterion(1:time_ind), &
                lbound(bif_criterion(1:time_ind)), ubound(bif_criterion(1:time_ind)))

            h5_currentgrp = "/"//trim(h5_mode_groupname) //"/br_formfactor"
            CALL h5_add_complex_1(h5_id, trim(h5_currentgrp), br_formfactor(1:time_ind), &
                lbound(real(br_formfactor(1:time_ind))), ubound(real(br_formfactor(1:time_ind))))

            h5_currentgrp = "/"//trim(h5_mode_groupname) //"/Ipar"
            CALL h5_add_complex_1(h5_id, trim(h5_currentgrp), Ipar_time(1:time_ind), &
                lbound(real(Ipar_time(1:time_ind))), ubound(real(Ipar_time(1:time_ind))))

            h5overwrite = .false.

            CALL h5_close(h5_id)
            CALL h5_deinit()

        else
            open (777, file='br_abs_res.dat', position='append')
            write (777, *) time_ind, time, antenna_factor, br_abs(time_ind)
            close (777)
        end if
    end subroutine ! write_br_dqle22_time_data


    subroutine check_linear_discr_pen_ratio
        use writeData_m, only: write_fields_currs_transp_coefs_to_h5

        implicit none

        if (time_ind .gt. 50 .and. .not. discr_reached) then
            ! calculate beta only once
            if (br_beta .eq. 0) then
                ! Calculate slope from the data until this time step. (Simple linear regression algorithm)
                br_beta = sum(br_abs(3:time_ind)*br_abs_time(3:time_ind))/sum(br_abs_time(1:time_ind)**2)
                write(*,*) 'br_beta = ', br_beta
            end if
            ! calculate the from the linear regression predicted value of Br_abs
            br_predicted = br_beta*br_abs_time(time_ind)
            write(*,*) "Delta = ", abs(br_abs(time_ind) - br_predicted)
            if (abs(br_abs(time_ind) - br_predicted) .gt. 0.1) then
                write(*,*) 'discrepancy to linearly predicted value of Br_abs_res > delta'
                if (modulo(time_ind, save_prof_time_step) .ne. 0) then
                    if (suppression_mode .eqv. .false.) then
                        call write_fields_currs_transp_coefs_to_h5(time_ind)
                    end if
                end if
                if (suppression_mode .eqv. .false.) then
                    CALL write_kin_prof_data_to_disk
                end if
                if (br_stopping) then

                    call write_reason_for_stop_to_h5("discrepancy to " //&
                        "linearly predicted value of Br_abs_res > delta")
                    CALL write_br_dqle22_time_data!, br_abs_time, br_abs_antenna_factor, br_abs, dqle22_res_time)
                    stop "Finished time evolution: br_stopping"

                else
                    call write_br_discrepancy_reached_info
                    discr_reached = .true.
                end if
            end if
        end if

    end subroutine

    subroutine write_br_discrepancy_reached_info

        use h5mod

        implicit none

        CALL h5_init()
        CALL h5_open_rw(path2out, h5_id)
        CALL h5_add_string(h5_id, trim(h5_mode_groupname)// &
            '/info', 'discrepancy to linearly predicted value of Br_abs_res > delta')
        CALL h5_add_double_1(h5_id, trim(h5_mode_groupname)// &
            '/discrep_time', (/time_ind*1.d0, time/), (/1/), (/2/))
        CALL h5_close(h5_id)
        CALL h5_deinit()

    end subroutine


    !> @brief subroutine write_kin_prof_data_to_disk(time_ind). Writes the profile data to hdf5 files.
    !> Formerly, this data was written to fort.1xxx ascii files.
    !> This routine was added because of the change that only every
    !>  "save_prof_time_step"th timestep is written. If the program is to be stopped
    !> because a stopping criterion was met, the profiles should be written for that
    !> last time step. Because this occurs more than once, it is more convenient to
    !> summarize this in a subroutine.
    !> @author Markus Markl
    !> @date 12.03.2021
    !> @param[in] time_ind Current step of the time evolution. Used to name the fort.1000 group in which
    !> the data is written.
    subroutine write_kin_prof_data_to_disk

        use grid_mod
        use plasma_parameters
        use control_mod
        use baseparam_mod
        use h5mod
        use KAMEL_hdf5_tools
        use wave_code_data, only: Vth, m_vals, n_vals
        use kim_wave_code_adapter_m, only: kim_periodic_mode_selected
        use periodic_checkpoint_m, only: periodic_restart_t, write_periodic_checkpoint

        implicit none
        type(periodic_restart_t) :: checkpoint
        integer :: ipoi

        if (ihdf5IO .eq. 1) then
            call log_debug("Write kinetic profiles")
            do ipoi = 1, npoic
                sqg_bthet_overcavg(ipoi) = 0.5d0*(sqrt_g_times_B_theta_over_c(ipoi) &
                                                + sqrt_g_times_B_theta_over_c(ipoi + 1))
                Ercovavg(ipoi) = 0.5d0*(Ercov(ipoi) + Ercov(ipoi + 1))
            end do
            ! h5_mode_groupname
            !h5_currentgrp = "/"//trim(h5_mode_groupname) &
            h5_currentgrp = trim(h5_mode_groupname) &
                            //"/KinProfiles"

            CALL h5_init()
            CALL h5_open_rw(path2out, h5_id)

            write (h5_currentgrp, "(A,A,I0,A)") trim(h5_currentgrp), &
                "/", 1000 + time_ind, "/"

            call log_debug("h5_currentgrp " // trim(h5_currentgrp))
            call log_debug("defining KinProfiles/1000 group")

            CALL h5_obj_exists(h5_id, trim(h5_currentgrp), h5_exists_log)
            if (.not. h5_exists_log) then
                CALL h5_create_parent_groups(h5_id, trim(h5_currentgrp))
            end if

            CALL h5_add_float_1(h5_id, trim(h5_currentgrp)//"rc", &
                                real(rc), lbound(rc), ubound(rc))

            CALL h5_add_float_1(h5_id, trim(h5_currentgrp)//"n", &
                                real(params(1, :)), lbound(params(1, :)), ubound(params(1, :)))

            CALL h5_add_float_1(h5_id, trim(h5_currentgrp)//"Vz", &
                                real(params(2, :)), lbound(params(2, :)), ubound(params(2, :)))

            CALL h5_add_float_1(h5_id, trim(h5_currentgrp)//"Te", &
                                real(params(3, :)/ev), lbound(params(3, :)), ubound(params(3, :)))

            CALL h5_add_float_1(h5_id, trim(h5_currentgrp)//"Ti", &
                                real(params(4, :)/ev), lbound(params(4, :)), ubound(params(4, :)))

            CALL h5_add_float_1(h5_id, trim(h5_currentgrp)//"Er", &
                                real(Ercovavg), lbound(Ercovavg), ubound(Ercovavg))

            CALL h5_add_float_1(h5_id, trim(h5_currentgrp)//"sqg_btheta_overc", &
                                real(sqg_bthet_overcavg), lbound(sqg_bthet_overcavg), &
                                ubound(sqg_bthet_overcavg))

            CALL h5_add_float_1(h5_id, trim(h5_currentgrp)//"Vth", &
                                real(Vth), lbound(Vth), ubound(Vth))

            if (trim(wave_code) == 'KIM' .and. kim_periodic_mode_selected() .and. &
                periodic_amplitudes%initialized .and. time_ind > 0) then
                call capture_periodic_restart(checkpoint)
                call write_periodic_checkpoint(h5_id, trim(h5_currentgrp), &
                    m_vals, n_vals, periodic_amplitudes, checkpoint)
            end if

            CALL h5_close(h5_id)
            CALL h5_deinit()
            call log_debug("finished writing KinProfiles")

        else
            do ipoi = 1, npoic
                write (1000 + time_ind, *) rc(ipoi), params(1:2, ipoi) &
                    , params(3, ipoi)/ev &
                    , params(4, ipoi)/ev &
                    , 0.5d0*(Ercov(ipoi) + Ercov(ipoi + 1)) &
                    , 0.5d0*(sqrt_g_times_B_theta_over_c(ipoi) + &
                            sqrt_g_times_B_theta_over_c(ipoi + 1))
            end do
            close (1000 + time_ind)
        end if

    end subroutine write_kin_prof_data_to_disk


    subroutine write_kin_profile_at_time_index

        implicit none

        call log_debug("Write kinetic profiles at time index")
        if (modulo(time_ind, save_prof_time_step) .eq. 0) then
            call write_kin_prof_data_to_disk
        end if

    end subroutine


    subroutine interp_Br_Dql_at_resonance_timeevol

        use PolyLagrangeInterpolation
        use grid_mod, only: npoib, r_resonant, rb, dqle22, dae22
        use wave_code_data, only: antenna_factor, Br

        implicit none

        integer :: indResRadius, ind_begin_interp, ind_end_interp

        call binsrc(rb, 1, npoib, r_resonant(1), indResRadius)
        call get_ind_Lagr_interp(indResRadius, ind_begin_interp, ind_end_interp)
        call plag_coeff(nlagr, nder, r_resonant(1), rb(ind_begin_interp:ind_end_interp), coef)

        br_abs(time_ind) = sum(coef(0, :)*abs(Br(ind_begin_interp:ind_end_interp)))!*sqrt(antenna_factor)
        dqle22_res_time(time_ind) = sum(coef(0, :)*dqle22(ind_begin_interp:ind_end_interp))
        dae22_res_time(time_ind) = sum(coef(0, :)*dae22(ind_begin_interp:ind_end_interp))
        bif_criterion(time_ind) = dqle22_res_time(time_ind)/dae22_res_time(time_ind)

        if (bif_criterion(time_ind) .ge. 1.0d0) then
            write(*,*) "!!! Bifurcation criterion met !!!"
            write(*,*) "bif_criterion = ", bif_criterion(time_ind)
        end if

        ! save the time for the improved stopping criterion
        br_abs_time(time_ind) = time
        br_abs_antenna_factor(time_ind) = antenna_factor

    end subroutine


    subroutine message_Br_Dqle_values

        use wave_code_data, only: antenna_factor

        implicit none

        write(*,*) " "
        write(*,*) "+ + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + +"
        WRITE(*,'(A9,F10.4,A16,I4)') '  time = ', br_abs_time(time_ind), " s, time_ind = ", time_ind
        WRITE(*,'(A23,E10.5,A11,F6.2,A12)') '    Antenna_factor   = ', antenna_factor, " which are ", &
            antenna_factor/antenna_factor_max*100, "% of the max"
        WRITE(*,'(A23,F10.5,A2)') '    Br abs res * C_mn= ', br_abs(time_ind) * SQRT(antenna_factor), " G"
        WRITE(*,'(A23,F10.5,A2)') '    Br abs res       = ', br_abs(time_ind), " G"
        WRITE(*,'(A23,F10.3,A7)') '    Dqle22 res       = ', dqle22_res_time(time_ind), " cm^2/s"
        WRITE(*,'(A23,F10.5)')    '    bif crit         = ', bif_criterion(time_ind)
        write(*,*) '   Form factor      = ', abs(br_formfactor(time_ind))
        write(*,*) '   Ipar             = ', abs(Ipar_time(time_ind))
        write(*,*) "+ + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + +"
        write(*,*) " "

    end subroutine

    subroutine stop_if_time_step_too_small

        use h5mod

        implicit none

        character(*), parameter :: reason = 'timestep < stop time step'

        if (timstep .lt. stop_time_step .and. time .gt. 1.0d-3) then
            write(*,*) 'stop: timestep smaller than stop limit'
            if (suppression_mode .eqv. .false.) then
                call write_kin_prof_data_to_disk
            end if

            call write_reason_for_stop_to_h5(reason)
            stop
        end if

    end subroutine


    subroutine calc_params_num_and_denom

        use plasma_parameters, only: params_num, params_denom, params, params_beg

        implicit none

        params_num = (params - params_beg)**2
        params_denom = params**2 + params_beg**2

    end subroutine

    subroutine smooth_params_num_and_denom

        use grid_mod, only: npoi, nbaleqs, mwind, dummy
        use plasma_parameters, only: params_num, params_denom
        use QLBalance_kinds, only: dp

        implicit none

        integer :: ieq
        real(dp), allocatable :: row_buffer(:)

        allocate(row_buffer(npoi))
        do ieq = 1, nbaleqs
            row_buffer = params_num(ieq, :)
            call smooth_array_gauss(npoi, mwind, row_buffer, dummy)
            params_num(ieq, :) = dummy
            row_buffer = params_denom(ieq, :)
            call smooth_array_gauss(npoi, mwind, row_buffer, dummy)
            params_denom(ieq, :) = dummy
        end do
        deallocate(row_buffer)

    end subroutine

    subroutine determine_timscal

        use grid_mod, only: npoi, rc
        use plasma_parameters, only: params_num, params_denom
        use baseparam_mod, only: factolmax
        use recstep_mod, only: tol

        implicit none

        integer :: ipoi

        do ipoi = 1, npoi
            if (rc(ipoi) .lt. 0.95d0*rc(npoi)) then
                timscal(ipoi) = sum(sqrt(params_num(3:4, ipoi)/params_denom(3:4, ipoi)))
            else
                timscal(ipoi) = 1d-30 !0.d0
            end if
        end do

        call log_debug("determine_timscal complete")

    end subroutine

    subroutine rescale_time_step_array

        use grid_mod, only: npoi, nbaleqs
        use recstep_mod, only: tim_stack, timstep_arr, tol
        use QLbalance_diag, only: timscal_dql

        implicit none

        integer :: ipoi, ieq, k

        timscal = timscal + timscal_dql

        do ipoi = 1, npoi
            do ieq = 1, nbaleqs
                k = nbaleqs*(ipoi - 1) + ieq
                !timstep_arr(k)=timstep_arr(k)/timscal(ipoi)*tol
                timstep_arr(k) = timstep_arr(k)/max(timscal(ipoi), epsilon(1.d0))*tol
                ! steady state solution:
                !if (ieq .gt. 1 .and. r(ipoi) .gt. rsepar-0.5d0) then
                !    timstep_arr(k) = 0d0
                !end if
            end do
        end do

        timstep_arr = timstep_arr*timescale/(timstep_arr + timescale)
        if (scratch) then
            scratch = .false.
            tim_stack = timstep_arr
        end if
        timstep_arr = 2.d0*timstep_arr*tim_stack/(timstep_arr + tim_stack)

    end subroutine

    subroutine set_time_step

        use recstep_mod, only: timstep_arr, tol

        implicit none

        timstep = minval(timstep_arr)

        if (.not. set_constant_time_step) then
            ! limit time step from below:
            timstep = max(timstep, timstep_min)
            ! limit timestep from above:
            !if (ramp_up_mode .ne. 0) timstep = min(timstep,0.1)
            !timstep = min(timstep,0.005)
        else
            ! use for constant time step:
            timstep = constant_time_step
            write(*,*) "constant time step = ", timstep
        end if

        call log_debug('timstep set')

    end subroutine

    subroutine reset_timstep_arr_w_timstep

        use recstep_mod, only: tim_stack, timstep_arr

        implicit none

        timstep_arr = 0.0d0
        timstep_arr = timstep
        tim_stack = timstep_arr

    end subroutine

    subroutine write_time_info

        implicit none

        if (ihdf5IO .eq. 1) then
            call write_time_info_to_h5
        else
            call write_time_info_to_txt
        end if

    end subroutine

    subroutine write_time_info_to_txt

        use QLbalance_diag, only: rate_dql, timscal_dql

        implicit none

        open (4321, file='timstep_evol.dat', position='append')
        write (4321, *) time_ind, timstep, timscal_dql, timscal(1), rate_dql, time
        close (4321)

    end subroutine

    subroutine write_time_info_to_h5

        use QLbalance_diag, only: rate_dql, timscal_dql

        implicit none
        real(dp), allocatable :: previous(:,:), history(:,:)
        integer :: rows, lo1, lo2, hi1, hi2
        logical :: exists

        if (time_ind < 1) error stop 'time history index must be positive'
        rows = 3
        if (data_verbosity >= 2) rows = 6
        h5_currentgrp = '/'//trim(h5_mode_groupname)//'/timstep_evol.dat'
        call h5_init()
        call h5_open_rw(path2out, h5_id)
        call h5_obj_exists(h5_id, trim(h5_currentgrp), exists)
        hi2 = 0
        if (exists) then
            call h5_get_bounds(h5_id, trim(h5_currentgrp), lo1, lo2, hi1, hi2)
            if (lo1 /= 1 .or. lo2 /= 1 .or. hi1 /= rows .or. hi2 < 1) &
                error stop 'time history shape or verbosity changed'
            allocate(previous(rows,hi2))
            call h5_get(h5_id, trim(h5_currentgrp), previous)
        end if
        ! Reopening a file cannot reuse Fortio's process-local unlimited buffer.
        ! Resolve its stored columns by path before replacing this one dataset.
        allocate(history(rows,max(hi2,time_ind)), source=0.0_dp)
        if (exists) history(:,:hi2) = previous
        if (data_verbosity >= 2) then
            history(:,time_ind) = [real(time_ind,dp), timstep, &
                timscal_dql, timscal(1), rate_dql, time]
        else
            history(:,time_ind) = [real(time_ind,dp), timstep, time]
        end if
        if (exists) call h5_delete(h5_id, trim(h5_currentgrp))
        call h5_add(h5_id, trim(h5_currentgrp), history, [1,1], shape(history))
        call h5_close(h5_id)
        call h5_deinit()

    end subroutine

    subroutine relax_plasma_parameters

        use grid_mod, only: npoi, nbaleqs
        use plasma_parameters, only: params
        use baseparam_mod, only: urelax

        implicit none

        integer :: ipoi, ieq, k

        do ipoi = 1, npoi
            do ieq = 1, nbaleqs
                k = nbaleqs*(ipoi - 1) + ieq
                params(ieq, ipoi) = yprev(k)*urelax + params(ieq, ipoi)*(1.d0 - urelax)
            end do
        end do

    end subroutine

    subroutine msg_time_info

        implicit none

        write(*,*) ' '
        write(*,*) 'Debug: i = ', int2(time_ind), 'time = ', real(time)
        write(*,*) ' '

    end subroutine

    subroutine set_first_iteration_true

        implicit none

        if (firstiterationdone .eqv. .false.) firstiterationdone = .true.

    end subroutine


    subroutine determine_Dql_diagnostic

        use grid_mod, only: dqle11, dqli11
        use QLbalance_diag

        implicit none

        timscal_dql = maxval(abs(dqle11_prev - dqle11))/ &
            max(maxval(abs(dqle11_prev) + abs(dqle11)), tiny(1.0_dp))
        ind_dqle = maxloc(abs(dqle11_prev - dqle11))
        timscal_dqli = maxval(abs(dqli11_prev - dqli11))/ &
            max(maxval(abs(dqli11_prev) + abs(dqli11)), tiny(1.0_dp))
        ind_dqli = maxloc(abs(dqli11_prev - dqli11))
        rate_dql = timscal_dql/timstep

    end subroutine

    subroutine create_group_structure_timeevol

        use control_mod
        use wave_code_data, only: m_vals, n_vals
        use periodic_checkpoint_m, only: periodic_mode_group
        use h5mod

        implicit none

        call log_debug("Creating group structure for TimeEvol")

        h5_mode_groupname = periodic_mode_group(m_vals, n_vals)

        CALL h5_init()
        CALL h5_open_rw(path2out, h5_id)

        if (.not. suppression_mode) then
            call log_debug("h5_mode_groupname " // trim(h5_mode_groupname))
            CALL h5_create_parent_groups(h5_id, trim(h5_mode_groupname) //'/')
            CALL h5_create_parent_groups(h5_id, trim(h5_mode_groupname)//"/KinProfiles/")

            call create_group_if_not_existent(trim(h5_mode_groupname)//"/LinearProfiles/")
            call create_group_if_not_existent("/init_params")
        else
            call log_debug("h5_mode_groupname: " // trim(h5_mode_groupname))
            call create_group_if_not_existent(trim(h5_mode_groupname))
            call create_group_if_not_existent("/init_params")
        end if

        CALL h5_close(h5_id)
        CALL h5_deinit()

        call log_debug("finished creating group structure for TimeEvol")
    end subroutine



    subroutine redoTimeStep

        use recstep_mod, only: timstep_arr
        use grid_mod, only: npoic, rc, Ercov
        use plasma_parameters, only: params, params_begbeg
        use baseparam_mod, only: eV, factolmax
        use restart_mod

        implicit none

        integer :: ipoi

        print *, 'redo step with old DQL'

        call hold_prev_transp_coeffs
        iunit_redo = 137

        open (iunit_redo, file='params_redostep.after')
        do ipoi = 1, npoic
            write (iunit_redo, *) rc(ipoi), params(1:2, ipoi) &
                , params(3, ipoi)/ev &
                , params(4, ipoi)/ev &
                , 0.5d0*(Ercov(ipoi) + Ercov(ipoi + 1))
        end do
        close (iunit_redo)
        params = params_begbeg
        open (iunit_redo, file='params_redostep.before')
        do ipoi = 1, npoic
            write (iunit_redo, *) rc(ipoi), params(1:2, ipoi) &
                , params(3, ipoi)/ev &
                , params(4, ipoi)/ev &
                , 0.5d0*(Ercov(ipoi) + Ercov(ipoi + 1))
        end do
        close (iunit_redo)

        timstep = timstep/factolmax
        timstep_arr = timstep

    end subroutine

end module
