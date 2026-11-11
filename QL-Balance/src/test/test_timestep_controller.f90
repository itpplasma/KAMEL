program test_timestep_controller
    !! Adaptive time-step controller of the QL-Balance time evolution on a
    !! synthetic grid. No wave code is involved: with irf = 0 the production
    !! doStep skips get_dql and advances the profiles with the production
    !! implicit solver using zero quasilinear transport and the constant
    !! anomalous diffusion read from the fixture input file.
    !!
    !! The expected values are produced by an oracle that calls the step
    !! pieces in the order doStep does (evolvestep, determine_timscal,
    !! redo with factolred) and by the documented update rule evaluated with
    !! plain arithmetic.
    use control_mod, only: wave_code, kim_run_type, type_of_run, irf, suppression_mode, &
                           misalign_diffusion, ihdf5IO, data_verbosity, eps, temperature_limit, &
                           write_gyro_current, gyro_current_study, jpar_method
    use wave_code_data, only: dim_mn, m_vals, n_vals, r, n, Te, Ti, q, Vth, Vz, &
                              dPhi0, I_par_toroidal, antenna_factor
    use grid_mod, only: rb, r_resonant, gg_width, gg_factor, npoib, npoi, rmin, rmax, npoimin, &
                        iboutype, set_boundary_condition, mwind, dery_equisource, source_term, &
                        dae22, rc, rb_cut_in, rb_cut_out, re_cut_out
    use plasma_parameters, only: params, params_b, params_beg, params_num, params_denom, &
                                 limit_temps_from_below, init_background_profiles
    use baseparam_mod, only: ev, am, rtor, Z_i, btor, rsepar, factolmax, factolred, urelax
    use QLBalance_diag, only: timscal_dql
    use paramscan_mod, only: viscosity_factor
    use resonances_mod, only: prop, numres, r_res, width_res, ampl_res
    use recstep_mod, only: tol, timstep_arr, tim_stack
    use time_evolution, only: TimeEvolution_t, doStep, allocate_prev_variables, &
                              allocate_timscal_and_params, alloc_Br_Dqle_for_timeevol, &
                              reset_timstep_arr_w_timstep, set_time_step, determine_timscal, &
                              calc_params_num_and_denom, smooth_params_num_and_denom, &
                              stop_if_time_step_too_small, timscal, timstep, timstep_min, &
                              stop_time_step, timescale, time, time_ind, Nstorage, scratch, &
                              set_constant_time_step, constant_time_step, ramp_up_mode, &
                              antenna_factor_max, antenna_max_stopping, br_stopping, &
                              firstiterationdone, br_abs, br_formfactor, br_vac_res, &
                              br_abs_antenna_factor, br_abs_time, dqle22_res_time, &
                              dae22_res_time, bif_criterion, save_prof_time_step, tmax, t_max_ramp_up
    use h5mod, only: path2out, path2inp, h5_mode_groupname, h5_id, group_id_1, h5_create, &
                     h5_define_group, h5_close_group, h5_close, h5_deinit, h5_add, h5_open, &
                     h5_obj_exists
    implicit none

    type(TimeEvolution_t) :: evolution
    character(1024) :: argument, executable, command
    integer :: stat, child_status, command_status
    real(8), allocatable :: profile_before_rejection(:, :)
    external :: gengrid, allocate_wave_code_data, calc_geometric_parameter_profiles
    external :: evolvestep, det_balance_eqs_source_terms, initialize_get_dql

    call get_command_argument(1, argument)
    call get_command_argument(0, executable)

    call build_fixture()

    if (trim(argument) == 'stop') then
        call run_stop_child()
        stop
    end if

    call check_clock_over_accepted_steps()
    call check_constant_step_retry()
    call check_step_floor()
    call check_stop_guard()
    call check_controller_watches_temperatures_only()

    print *, 'PASS: adaptive time-step controller'

contains

    subroutine require(condition, message)
        logical, intent(in) :: condition
        character(*), intent(in) :: message
        if (.not. condition) then
            print *, 'FAIL: ', message
            error stop 1
        end if
    end subroutine require

    logical function close_to(actual, expected, rel)
        real(8), intent(in) :: actual, expected, rel
        close_to = abs(actual - expected) <= rel*abs(expected)
    end function close_to

    subroutine build_fixture()
        real(8), allocatable :: input_r(:), diffusion(:)

        rmin = 2.0d0
        rmax = 60.0d0
        npoimin = 99
        iboutype = 2
        gg_factor = 0.0d0
        gg_width = 58.0d0
        prop = .false.
        numres = 0
        call gengrid()
        call set_boundary_condition()
        numres = 1
        allocate (r_res(1), width_res(1), ampl_res(1), r_resonant(1))
        r_res = 31.0d0
        r_resonant = r_res
        width_res = 1.0d0
        ampl_res = 0.0d0
        dim_mn = 1
        allocate (m_vals(1), n_vals(1), input_r(npoib))
        m_vals = -6
        n_vals = 2
        input_r = rb
        call allocate_wave_code_data(npoib, input_r)
        n = 2.0d13*(1.1d0 - r/100.0d0)
        Te = 1.0d3*(1.2d0 - r/100.0d0)
        Ti = Te
        q = 1.5d0 + 3.0d0*(r - 2.0d0)/58.0d0
        dPhi0 = 0.5d0
        Vth = 0.0d0
        Vz = 0.0d0
        rtor = 165.0d0
        am = 2.0d0
        Z_i = 1.0d0
        btor = -18000.0d0
        rsepar = 60.0d0
        viscosity_factor = 1.0d0
        mwind = 1
        dery_equisource = 0.0d0
        source_term = 0.0d0
        call init_background_profiles()
        params_b(1, :) = n
        params_b(2, :) = Vz/rtor
        params_b(3, :) = Te*ev
        params_b(4, :) = Ti*ev
        call calc_geometric_parameter_profiles()

        wave_code = 'KiLCA'
        kim_run_type = ''
        type_of_run = 'TimeEvolution'
        irf = 0
        suppression_mode = .true.
        misalign_diffusion = .false.
        jpar_method = 'conductivity'
        ihdf5IO = 1
        data_verbosity = 0
        write_gyro_current = .false.
        gyro_current_study = 0
        I_par_toroidal = 0.0d0
        antenna_factor = 1.0d0
        temperature_limit = 20.0d0
        eps = 1.0d-8
        factolmax = 2.0d0
        factolred = 0.5d0
        urelax = 0.0d0
        tol = 0.1d0
        timstep = 1.0d-4
        timstep_min = 1.0d-14
        stop_time_step = 0.0d0
        timescale = 1.0d-2
        time = 0.0d0
        time_ind = 1
        tmax = 1.0d0
        t_max_ramp_up = 1.0d12
        Nstorage = 16
        save_prof_time_step = 100
        scratch = .true.
        set_constant_time_step = .false.
        constant_time_step = timstep
        ramp_up_mode = 4
        antenna_factor_max = 1.0d0
        antenna_max_stopping = 1.0d6
        br_stopping = .false.
        firstiterationdone = .false.
        call allocate_prev_variables()
        call allocate_timscal_and_params()
        call reset_timstep_arr_w_timstep()

        rb_cut_in = 0.0d0
        rb_cut_out = maxval(rb) + 1.0d0
        re_cut_out = rb_cut_out + 1.0d0
        diffusion = spread(1500.0d0, 1, npoib)
        path2inp = 'timestep_controller_input.h5'
        call h5_create(trim(path2inp), h5_id)
        call h5_define_group(h5_id, 'da_estimation', group_id_1)
        call h5_close_group(group_id_1)
        call h5_add(h5_id, '/da_estimation/r', rb, [1], [npoib])
        call h5_add(h5_id, '/da_estimation/Da', diffusion, [1], [npoib])
        call h5_close(h5_id)
        call h5_deinit()
        path2out = 'timestep_controller.h5'
        if (trim(argument) == 'stop') path2out = 'timestep_controller_stop.h5'
        h5_mode_groupname = 'timestep_controller'
        call h5_create(trim(path2out), h5_id)
        call h5_define_group(h5_id, trim(h5_mode_groupname), group_id_1)
        call h5_close_group(group_id_1)
        call h5_close(h5_id)
        call h5_deinit()

        irf = 2
        call initialize_get_dql()
        irf = 0
        call det_balance_eqs_source_terms()
        call alloc_Br_Dqle_for_timeevol()
        br_abs = 0.0d0
        br_formfactor = 0.0d0
        br_vac_res = 0.0d0
        br_abs_antenna_factor = 0.0d0
        br_abs_time = 0.0d0
        dqle22_res_time = 0.0d0
        dae22_res_time = 0.0d0
        bif_criterion = 0.0d0
        dae22 = 1.0d0
    end subroutine build_fixture

    !> Replays the redo loop of doStep on the current profiles and returns the
    !> step size of the accepted solve, the number of rejections and the
    !> accepted profiles. The caller's profiles are restored.
    subroutine oracle_accepted_solve(step_used, rejections, accepted_params, accepted_timscal)
        real(8), intent(out) :: step_used, accepted_params(:, :), accepted_timscal(:)
        integer, intent(out) :: rejections
        real(8) :: entry_arr(size(timstep_arr))
        real(8), allocatable :: entry_params(:, :), entry_beg(:, :)
        integer :: iredo

        entry_arr = timstep_arr
        entry_params = params
        entry_beg = params_beg
        iredo = 0
        do
            iredo = iredo + 1
            params_beg = params
            call evolvestep(timstep, eps)
            call limit_temps_from_below()
            call calc_params_num_and_denom()
            call smooth_params_num_and_denom()
            call determine_timscal()
            if (maxval(timscal) < tol*factolmax) exit
            timstep_arr = timstep_arr*factolred
            params = params_beg
            call require(iredo <= 100, 'oracle did not converge')
        end do
        step_used = timstep_arr(1)
        rejections = iredo - 1
        accepted_params = params
        accepted_timscal = timscal
        timstep_arr = entry_arr
        params = entry_params
        params_beg = entry_beg
    end subroutine oracle_accepted_solve

    !> Next step size by the documented rule: timstep_arr/timscal*tol, limited
    !> by timescale, harmonically averaged with the previous step (tim_stack),
    !> floored at timstep_min.
    real(8) function documented_next_step(step_used, accepted_timscal, previous_stack, first)
        real(8), intent(in) :: step_used, accepted_timscal(:), previous_stack
        logical, intent(in) :: first
        real(8) :: candidate, s
        s = maxval(accepted_timscal + timscal_dql)
        candidate = step_used/max(s, epsilon(1.0d0))*tol
        candidate = candidate*timescale/(candidate + timescale)
        if (.not. first) candidate = 2.0d0*candidate*previous_stack/(candidate + previous_stack)
        documented_next_step = max(candidate, timstep_min)
    end function documented_next_step

    subroutine check_clock_over_accepted_steps()
        real(8) :: entry_steps(3), step_used, expected_clock, expected_next
        real(8), allocatable :: oracle_params(:, :), oracle_timscal(:)
        integer :: i, rejections, total_rejections
        logical :: first

        allocate (oracle_params(size(params, 1), size(params, 2)), oracle_timscal(size(timscal)))
        ! Entry steps: a small step (accepted at once), then steps large enough
        ! to violate tol*factolmax and force the redo loop.
        entry_steps = [1.0d-4, 1.0d3, 1.0d4]
        expected_clock = 0.0d0
        total_rejections = 0
        first = .true.
        do i = 1, size(entry_steps)
            if (i == 3) profile_before_rejection = params
            timstep = entry_steps(i)
            call reset_timstep_arr_w_timstep()
            call oracle_accepted_solve(step_used, rejections, oracle_params, oracle_timscal)
            call require(close_to(step_used, entry_steps(i)*factolred**rejections, 1.0d-14), &
                         'oracle step is not the entry step reduced by factolred per rejection')
            if (i == 1) call require(rejections == 0, 'small first step was rejected')
            if (i == 3) call require(rejections > 0, 'large step was not rejected')
            total_rejections = total_rejections + rejections

            time_ind = i
            expected_clock = expected_clock + step_used
            expected_next = documented_next_step(step_used, oracle_timscal, tim_stack(1), first)
            call doStep(evolution, stat)
            call require(stat == 0, 'doStep failed')

            call require(maxval(abs(params - oracle_params)/abs(oracle_params)) < 1.0d-10, &
                         'doStep did not advance the profiles with the accepted retry step')
            call require(close_to(timstep, expected_next, 1.0d-10), &
                         'next step does not follow the documented update rule')
            print '(a,i0,a,i0,3(a,es12.5))', 'step ', i, ' rejections ', rejections, &
                ' used ', step_used, ' next ', timstep, ' clock-expected ', time - expected_clock
            call require(close_to(time, expected_clock, 1.0d-10), &
                         'clock differs from the sum of the steps used by the accepted solves')
            first = .false.
        end do
        call require(total_rejections >= 1, 'fixture did not exercise the redo loop')
    end subroutine check_clock_over_accepted_steps

    subroutine check_constant_step_retry()
        real(8) :: step_used, clock_before
        real(8), allocatable :: oracle_params(:, :), oracle_timscal(:)
        integer :: rejections

        allocate (oracle_params(size(params, 1), size(params, 2)), oracle_timscal(size(timscal)))
        params = profile_before_rejection
        time_ind = 1
        set_constant_time_step = .true.
        constant_time_step = 1.0d4
        timstep = constant_time_step
        call reset_timstep_arr_w_timstep()
        call oracle_accepted_solve(step_used, rejections, oracle_params, oracle_timscal)
        call require(rejections > 0, 'constant step was not rejected')
        clock_before = time
        time_ind = time_ind + 1
        call doStep(evolution, stat)
        call require(stat == 0, 'rejected constant step did not terminate')
        call require(maxval(abs(params - oracle_params)/abs(oracle_params)) < 1.0d-10, &
                     'rejected constant step was not retried with a reduced step')
        call require(timstep == constant_time_step, 'constant step mode changed the step size')
        call require(close_to(time - clock_before, step_used, 1.0d-10), &
                     'clock after a rejected constant step differs from the step used')
        set_constant_time_step = .false.
    end subroutine check_constant_step_retry

    subroutine check_step_floor()
        timstep_arr = 1.0d-30
        timstep_min = 1.0d-9
        call set_time_step()
        call require(timstep == timstep_min, 'step below the floor was not raised to timstep_min')
        timstep_arr = 3.0d-7
        timstep_arr(5) = 2.0d-7
        call set_time_step()
        call require(timstep == 2.0d-7, 'step above the floor is not the minimum of timstep_arr')
        set_constant_time_step = .true.
        constant_time_step = 4.0d-7
        call set_time_step()
        call require(timstep == 4.0d-7, 'constant step mode ignored constant_time_step')
        set_constant_time_step = .false.
        timstep_min = 1.0d-14
    end subroutine check_step_floor

    subroutine check_stop_guard()
        logical :: marker_exists, group_has_reason

        ! No stop: step below the limit during the first millisecond, or above it later.
        stop_time_step = 1.0d-9
        timstep = 1.0d-12
        time = 5.0d-4
        call stop_if_time_step_too_small()
        timstep = 1.0d-6
        time = 1.0d0
        call stop_if_time_step_too_small()
        stop_time_step = 0.0d0

        ! Stop: a separate process reaching the guard must terminate before
        ! it can write the marker file.
        open (unit=77, file='stop_guard_marker', status='replace')
        close (77, status='delete')
        call require(index(trim(executable), "'") == 0, 'test executable path contains an apostrophe')
        command = "'"//trim(executable)//"' stop"
        call execute_command_line(trim(command), exitstat=child_status, cmdstat=command_status)
        call require(command_status == 0 .and. child_status == 0, 'stop child failed to run')
        inquire (file='stop_guard_marker', exist=marker_exists)
        call require(.not. marker_exists, 'step below stop_time_step did not stop the run')
        call h5_open('timestep_controller_stop.h5', h5_id)
        call h5_obj_exists(h5_id, 'timestep_controller/stopping_criterion', group_has_reason)
        call h5_close(h5_id)
        call h5_deinit()
        call require(group_has_reason, 'stop reason was not recorded in the output file')
    end subroutine check_stop_guard

    subroutine run_stop_child()
        stop_time_step = 1.0d-9
        timstep = 1.0d-12
        time = 1.0d0
        call stop_if_time_step_too_small()
        open (unit=77, file='stop_guard_marker', status='replace')
        close (77)
    end subroutine run_stop_child

    subroutine check_controller_watches_temperatures_only()
        integer :: ipoi
        real(8) :: expected

        params_num = 0.0d0
        params_denom = 1.0d0
        params_num(1, :) = 0.5d0
        params_num(2, :) = 0.5d0
        call determine_timscal()
        call require(all(timscal <= 1.0d-29), &
                     'density/rotation change influenced the step controller')
        params_num(3, :) = 0.04d0
        params_num(4, :) = 0.09d0
        call determine_timscal()
        expected = sqrt(0.04d0) + sqrt(0.09d0)
        do ipoi = 1, npoi
            if (rc(ipoi) < 0.95d0*rc(npoi)) then
                call require(close_to(timscal(ipoi), expected, 1.0d-12), &
                             'temperature change is not the sum of the Te and Ti measures')
            else
                call require(timscal(ipoi) <= 1.0d-29, 'edge region contributes to the controller')
            end if
        end do
    end subroutine check_controller_watches_temperatures_only

end program test_timestep_controller
