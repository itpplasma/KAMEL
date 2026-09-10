program test_periodic_time_evolution
    !! Actual periodic KIM transport coupled to QL's implicit sparse advance.
    use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
    use control_mod, only: wave_code, kim_config_path, kim_profiles_from_balance, &
                      type_of_run, kim_run_type, kim_ion_transport_model, kim_transport_benchmark, &
                           irf, suppression_mode, misalign_diffusion, ihdf5IO, data_verbosity, &
                      write_gyro_current, gyro_current_study, jpar_method, eps, temperature_limit, &
                           kim_current_floor, kim_current_max_scale, kim_current_relaxation
    use kim_wave_code_adapter_m, only: kim_initialize, kim_vac_Br, kim_periodic_scale_modes, &
                                       kim_periodic_current_unit, kim_periodic_scale_status, &
                                       kim_embedding_metadata, kim_transition_weights
    use periodic_amplitude_state_m, only: periodic_amplitudes, periodic_amplitude_state_t
    use periodic_checkpoint_m, only: pending_periodic_restart, periodic_mode_group
    use wave_code_data, only: dim_mn, m_vals, n_vals, r, n, Te, Ti, q, Vth, Vz, &
                         dPhi0, I_par_toroidal, antenna_factor, Es, Ep, Er, Et, Ez, Br, Bp, Jpe, Jpi
    use grid_mod, only: rb, rc, r_resonant, gg_width, gg_factor, npoib, npoic, npoi, &
                        rmin, rmax, npoimin, iboutype, set_boundary_condition, mwind, neqset, &
                        dery_equisource, source_term, dqle11, dqle12, dqle21, dqle22, &
                        dqli11, dqli12, dqli21, dqli22
    use plasma_parameters, only: params, params_b, init_background_profiles
    use baseparam_mod, only: ev, am, rtor, Z_i, btor, rsepar, factolmax, factolred, urelax
    use QLBalance_diag, only: timscal_dql, timscal_dqli, rate_dql
    use paramscan_mod, only: viscosity_factor
    use resonances_mod, only: prop, numres, r_res, width_res, ampl_res
    use recstep_mod, only: tol, timstep_arr, tim_stack
    use time_evolution, only: TimeEvolution_t, doStep, initialize_periodic_response, &
                              restore_periodic_accepted_response, last_periodic_rejections, &
                              apply_periodic_restart, determine_Dql_diagnostic, &
                 allocate_prev_variables, allocate_timscal_and_params, alloc_Br_Dqle_for_timeevol, &
                   hold_prev_transp_coeffs, reset_timstep_arr_w_timstep, Nstorage, time_ind, time, &
                              timstep, timstep_min, timescale, tmax, save_prof_time_step, scratch, &
                     set_constant_time_step, constant_time_step, ramp_up_mode, antenna_factor_max, &
                              antenna_max_stopping, br_stopping, firstiterationdone, &
                            br_abs, br_formfactor, br_vac_res, br_abs_antenna_factor, br_abs_time, &
                              dqle22_res_time, dae22_res_time, bif_criterion, Ipar_time
    use h5mod, only: path2out, path2inp, path2time, h5_mode_groupname, h5_id, group_id_1, &
                     h5_create, h5_define_group, h5_close_group, h5_close, h5_deinit, h5_add, &
                     h5_open, h5_get, h5_obj_exists
    implicit none

    type(TimeEvolution_t) :: evolution
    type(periodic_amplitude_state_t) :: accepted_state
    real(8), allocatable :: input_r(:), initial_profiles(:, :), one_advance(:, :)
    real(8), allocatable :: accepted_profiles(:, :), accepted_transport(:, :)
    real(8) :: accepted_time, accepted_timestep, accepted_floor, expected_advance
    integer :: stat, child_status, command_status
    character(1024) :: argument, executable, command
    logical :: restarted
    external :: gengrid, allocate_wave_code_data, calc_geometric_parameter_profiles
    external :: evolvestep
    external :: initialize_get_dql, det_balance_eqs_source_terms
    external :: read_background_profiles_h5_timeevol, interp_background_profiles

    call get_command_argument(2, argument)
    restarted = trim(argument) == 'restart'
    call get_command_argument(0, executable)

    ! Use the production finite-volume grid and interpolation stencils. The
    ! fixture supplies its known resonance directly, so no external q file is
    ! needed by grid generation.
    rmin = 2.0d0
    rmax = 60.0d0
    npoimin = 199
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
    n = 2.0d13 * (1.1d0 - r / 100.0d0)
    Te = 1.0d3 * (1.2d0 - r / 100.0d0)
    Ti = Te
    q = 1.5d0 + 3.0d0 * (r - 2.0d0) / 58.0d0
    dPhi0 = 0.5d0
    Vth = 0.0d0
    Vz = -8.0d4 * 165.0d0
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
    params_b(2, :) = Vz / rtor
    params_b(3, :) = Te * ev
    params_b(4, :) = Ti * ev
    call calc_geometric_parameter_profiles()

    wave_code = 'KIM'
    kim_run_type = 'electrostatic_periodic'
    kim_ion_transport_model = 'finite_larmor_radius'
    kim_transport_benchmark = .false.
    kim_profiles_from_balance = .true.
    type_of_run = 'TimeEvolution'
    irf = 1
    suppression_mode = .true.
    misalign_diffusion = .false.
    jpar_method = 'conductivity'
    ihdf5IO = 1
    data_verbosity = 0
    write_gyro_current = .false.
    gyro_current_study = 0
    I_par_toroidal = 0.25d0
    antenna_factor = 1.0d0
    kim_current_floor = 1.0d-30
    kim_current_max_scale = 1.0d12
    kim_current_relaxation = 0.4d0
    temperature_limit = 20.0d0
    eps = 1.0d-8
    factolmax = 2.0d0
    factolred = 0.5d0
    urelax = 0.0d0
    tol = 0.1d0
    timstep = 1.0d-8
    timstep_min = 1.0d-14
    timescale = 1.0d-2
    time = 0.0d0
    time_ind = 1
    tmax = 1.0d0
    Nstorage = 16
    save_prof_time_step = 100
    scratch = .true.
    set_constant_time_step = .true.
    constant_time_step = timstep
    ramp_up_mode = 0
    antenna_factor_max = 1.0d0
    antenna_max_stopping = 1.0d6
    br_stopping = .false.
    firstiterationdone = .false.
    call allocate_prev_variables()
    call allocate_timscal_and_params()
    call reset_timstep_arr_w_timstep()
    call create_inputs_and_output()

    kim_config_path = 'KIM_config_periodic_time_evolution.nml'
    call write_config(trim(kim_config_path))
    if (restarted) then
        path2time = 'periodic_time_evolution.h5'
        call read_background_profiles_h5_timeevol(2)
        call require(pending_periodic_restart%active, 'real profile loader missed the checkpoint')
        call interp_background_profiles()
        call init_background_profiles()
        call calc_geometric_parameter_profiles()
    end if
    call kim_initialize(npoib, input_r)
    kim_vac_Br = (1.0d0, 0.0d0)
    call initialize_get_dql()
    call det_balance_eqs_source_terms()
    call require(maxval(abs(dery_equisource)) > 0.0d0, 'production fixed source vanished')
    call alloc_Br_Dqle_for_timeevol()
    br_abs = 0.0d0
    br_formfactor = 0.0d0
    br_vac_res = 0.0d0
    br_abs_antenna_factor = 0.0d0
    br_abs_time = 0.0d0
    dqle22_res_time = 0.0d0
    dae22_res_time = 0.0d0
    bif_criterion = 0.0d0
    Ipar_time = 0.0d0
    if (restarted) call require(size(br_abs) == Nstorage + &
                                pending_periodic_restart%accepted_step, &
                                'restart history capacity omitted the accepted prefix')
    if (restarted) then
        call require(any(dery_equisource /= pending_periodic_restart%dery_equisource), &
                     'restart fixture did not regenerate a different equilibrium source')
        timscal_dql = -1000.0d0
        timscal_dqli = -2000.0d0
        rate_dql = -3000.0d0
        call apply_periodic_restart(pending_periodic_restart)
        call require(timscal_dql == pending_periodic_restart%timscal_dql .and. &
                     timscal_dqli == pending_periodic_restart%timscal_dqli .and. &
                     rate_dql == pending_periodic_restart%rate_dql, &
                     'restart did not restore adaptive transport controller state')
        call require(all(dery_equisource == pending_periodic_restart%dery_equisource), &
                     'restart did not restore the original fixed equilibrium source')
        call restore_periodic_accepted_response(stat, continuation=.true.)
        call require(stat == 0, 'could not rebuild the checkpoint response in a fresh process')
        call check_accepted_response()
        time_ind = 3
        suppression_mode = .false.
        save_prof_time_step = 1
        call doStep(evolution, stat)
        call require(stat == 0, 'actual restarted doStep failed')
        call check_accepted_response()
        call compare_saved_continuation()
        print *, 'PASS: cold restart reproduces the actual periodic doStep trajectory'
        stop
    end if
    call initialize_periodic_response(stat)
    call require(stat == 0, 'initial periodic response failed')
    call require(all(periodic_amplitudes%accepted == (1.0d0, 0.0d0)), &
                 'first accepted periodic response is not the constant unit drive')
    call hold_prev_transp_coeffs()
    initial_profiles = params

    ! A direct call to the same real sparse solver establishes the expected
    ! first advance using the accepted unit-drive transport. doStep must use
    ! exactly that response, rather than normalize again before advancing.
    call evolvestep(timstep, eps)
    one_advance = params
    params = initial_profiles
    call restore_periodic_accepted_response(stat)
    call require(stat == 0, 'could not restore the unit-drive initial response')
    call reset_timstep_arr_w_timstep()
    suppression_mode = .false.
    save_prof_time_step = 1
    call doStep(evolution, stat)
    call require(stat == 0, 'first actual periodic doStep failed')
    call require(last_periodic_rejections == 0, 'small first step unexpectedly rejected')
    call compare_profiles(params, one_advance, 'first doStep changed its frozen transport advance')
    call require(maxval(abs((params - initial_profiles) / initial_profiles)) > 1.0d-12, &
                 'the production implicit step did not evolve the profiles')
    call check_accepted_response()

    ! Force a genuine adaptive rejection using a large proposed profile step.
    ! The accepted retry must still leave a matched physical response/state.
    time_ind = time_ind + 1
    set_constant_time_step = .false.
    tol = 1.0d-6
    timstep = 1.0d-3
    accepted_time = time
    accepted_timestep = timstep
    call reset_timstep_arr_w_timstep()
    call doStep(evolution, stat)
    call require(stat == 0, 'adaptive periodic step did not recover after rejection')
    call require(last_periodic_rejections > 0, 'large profile step did not exercise rejection')
    expected_advance = accepted_timestep * factolred**last_periodic_rejections
    call require(abs((time - accepted_time) - expected_advance) <= 1.0d-12 * expected_advance, &
                 'accepted clock increment differs from the successfully retried profile advance')
    call check_accepted_response()

    ! Write through the actual accepted-profile checkpoint path, then advance
    ! uninterrupted once. A fresh executable below must reproduce this step
    ! after going through the actual restart profile reader and response rebuild.
    call check_saved_step(2, .true.)
    time_ind = 3
    call doStep(evolution, stat)
    call require(stat == 0, 'uninterrupted continuation step failed')
    call check_accepted_response()
    call save_continuation()

    ! Invalid candidate normalization must not commit profiles, fields,
    ! transport, controls, or amplitude history from the failed candidate.
    accepted_profiles = params
    call capture_transport(accepted_transport)
    accepted_state = periodic_amplitudes
    accepted_time = time
    accepted_timestep = timstep
    accepted_floor = kim_current_floor
    kim_current_floor = 1.0d100
    time_ind = time_ind + 1
    call doStep(evolution, stat)
    call require(stat /= 0, 'invalid candidate normalization was silently accepted')
    call check_saved_step(time_ind, .false.)
    call compare_profiles(params, accepted_profiles, 'failed candidate changed accepted profiles')
    call compare_transport(accepted_transport)
    call require(time == accepted_time, 'failed candidate advanced accepted time')
    call require(kim_current_floor == accepted_floor, 'failed candidate did not restore controls')
    call require(all(periodic_amplitudes%accepted == accepted_state%accepted), &
                 'failed candidate changed accepted complex amplitude')
    call require(all(periodic_amplitudes%trial == periodic_amplitudes%accepted), &
                 'failed candidate left a live amplitude trial')
    call require(all(kim_periodic_scale_modes == periodic_amplitudes%accepted), &
                 'failed candidate left transport at a different complex amplitude')
    call require(index(trim(executable), "'") == 0, 'test executable path contains an apostrophe')
    command = "'"//trim(executable)//"' '"//trim(kim_config_path)//"' restart"
    call execute_command_line(trim(command), exitstat=child_status, cmdstat=command_status)
    call require(command_status == 0 .and. child_status == 0, 'cold restart subprocess failed')
    dqle11 = 0.0d0
    dqli11 = 0.0d0
    call hold_prev_transp_coeffs()
    call determine_Dql_diagnostic()
    call require(ieee_is_finite(timscal_dql) .and. ieee_is_finite(timscal_dqli) .and. &
                 ieee_is_finite(rate_dql), 'zero transport made the adaptive diagnostic non-finite')
    print *, 'PASS: actual periodic doStep advances, rejects and restores coupled plasma/response'

contains

    subroutine create_inputs_and_output()
        use grid_mod, only: rb_cut_in, rb_cut_out, re_cut_out
        real(8), allocatable :: diffusion(:)
        rb_cut_in = 0.0d0
        rb_cut_out = maxval(rb) + 1.0d0
        re_cut_out = rb_cut_out + 1.0d0
        diffusion = spread(1500.0d0, 1, npoib)
        path2inp = 'periodic_time_input.h5'
        call h5_create(trim(path2inp), h5_id)
        call h5_define_group(h5_id, 'da_estimation', group_id_1)
        call h5_close_group(group_id_1)
        call h5_add(h5_id, '/da_estimation/r', rb, [1], [npoib])
        call h5_add(h5_id, '/da_estimation/Da', diffusion, [1], [npoib])
        call h5_close(h5_id)
        call h5_deinit()
        path2out = 'periodic_time_evolution.h5'
        if (restarted) path2out = 'periodic_time_restart.h5'
        h5_mode_groupname = periodic_mode_group(m_vals, n_vals)
        call h5_create(trim(path2out), h5_id)
        call h5_define_group(h5_id, trim(h5_mode_groupname), group_id_1)
        call h5_close_group(group_id_1)
        call h5_close(h5_id)
        call h5_deinit()
    end subroutine create_inputs_and_output

    subroutine check_saved_step(index, expected)
        integer, intent(in) :: index
        logical, intent(in) :: expected
        character(256) :: group
        logical :: found
        call h5_open(trim(path2out), h5_id)
        write (group, '(A,A,A,I0)') '/', trim(h5_mode_groupname), '/LinearProfiles/', index
        call h5_obj_exists(h5_id, trim(group), found)
        call require(found .eqv. expected, 'linear response output violates step acceptance')
        write (group, '(A,A,A,I0)') '/', trim(h5_mode_groupname), '/CurrentNormalization/', index
        call h5_obj_exists(h5_id, trim(group), found)
        call require(found .eqv. expected, 'normalization output violates step acceptance')
        write (group, '(A,A,A,I0)') '/', trim(h5_mode_groupname), '/KinProfiles/', 1000 + index
        call h5_obj_exists(h5_id, trim(group), found)
        call require(found .eqv. expected, 'profile checkpoint violates step acceptance')
        call h5_close(h5_id)
        call h5_deinit()
    end subroutine check_saved_step

    subroutine save_continuation()
        real(8), allocatable :: transport(:, :), response(:, :, :)
        call capture_transport(transport)
        call capture_response(response)
        call h5_create('periodic_time_expected.h5', h5_id)
        call h5_add(h5_id, 'profiles', params, [1, 1], shape(params))
        call h5_add(h5_id, 'transport', transport, [1, 1], shape(transport))
        call h5_add(h5_id, 'response', response, [1, 1, 1], shape(response))
        call h5_add(h5_id, 'unit_current', [real(kim_periodic_current_unit(1)), &
                                            aimag(kim_periodic_current_unit(1))], [1], [2])
        call h5_add(h5_id, 'embedding', kim_embedding_metadata(:, 1), [1], [4])
        call h5_add(h5_id, 'weights', kim_transition_weights(:, 1), [1], [npoib])
        call h5_add(h5_id, 'normalization_status', kim_periodic_scale_status(1))
        call h5_add(h5_id, 'fixed_source', dery_equisource, [1], [neqset])
        call h5_add(h5_id, 'amplitude', [real(periodic_amplitudes%accepted(1)), &
                                         aimag(periodic_amplitudes%accepted(1))], [1], [2])
        call h5_add(h5_id, 'time', time)
        call h5_add(h5_id, 'next_dt', timstep)
        call h5_close(h5_id)
        call h5_deinit()
    end subroutine save_continuation

    subroutine compare_saved_continuation()
        real(8), allocatable :: profiles(:, :), transport(:, :), response(:, :, :)
        real(8), allocatable :: actual(:, :, :), weights(:), fixed_source(:)
        real(8) :: amplitude(2), unit_current(2), embedding(4), saved_time, next_dt
        integer :: status, component
        allocate (profiles(4, npoic), transport(8, npoib), response(9, npoib, 2))
        allocate (weights(npoib), fixed_source(neqset))
        call h5_open('periodic_time_expected.h5', h5_id)
        call h5_get(h5_id, 'profiles', profiles)
        call h5_get(h5_id, 'transport', transport)
        call h5_get(h5_id, 'response', response)
        call h5_get(h5_id, 'unit_current', unit_current)
        call h5_get(h5_id, 'embedding', embedding)
        call h5_get(h5_id, 'weights', weights)
        call h5_get(h5_id, 'normalization_status', status)
        call h5_get(h5_id, 'fixed_source', fixed_source)
        call h5_get(h5_id, 'amplitude', amplitude)
        call h5_get(h5_id, 'time', saved_time)
        call h5_get(h5_id, 'next_dt', next_dt)
        call h5_close(h5_id)
        call h5_deinit()
        call compare_profiles(params, profiles, 'cold restart changed the evolved profiles')
        call compare_transport(transport)
        call capture_response(actual)
        do component = 1, 9
            call require(maxval(abs(actual(component, :, :) - response(component, :, :))) &
                         <= 1.0d-8 * max(maxval(abs(response(component, :, :))), 1.0d-250), &
                         'cold restart changed a physical field or species current')
        end do
        call require(all(dery_equisource == fixed_source), 'cold restart changed fixed sources')
        call require(status == kim_periodic_scale_status(1), 'cold restart changed guard status')
        call require(abs(kim_periodic_current_unit(1) &
                         - cmplx(unit_current(1), unit_current(2), 8)) &
                <= 1.0d-10 * abs(kim_periodic_current_unit(1)), 'cold restart changed unit current')
        call require(maxval(abs(kim_embedding_metadata(:, 1) - embedding)) <= 1.0d-10, &
                     'cold restart moved the periodic embedding')
        call require(maxval(abs(kim_transition_weights(:, 1) - weights)) <= 1.0d-10, &
                     'cold restart changed periodic transition weights')
        call require(abs(time - saved_time) <= 1.0d-12 * abs(saved_time), &
                     'cold restart changed accepted time')
        call require(abs(timstep - next_dt) <= 1.0d-12 * abs(next_dt), &
                     'cold restart changed the next adaptive timestep')
        call require(abs(periodic_amplitudes%accepted(1) - cmplx(amplitude(1), amplitude(2), 8)) &
                     <= 1.0d-10 * abs(periodic_amplitudes%accepted(1)), &
                     'cold restart changed accepted complex amplitude')
    end subroutine compare_saved_continuation

    subroutine capture_response(values)
        real(8), allocatable, intent(out) :: values(:, :, :)
        complex(8) :: response(9, npoib)
        response = transpose(reshape([Es, Ep, Er, Et, Ez, Br, Bp, Jpe, Jpi], [npoib, 9]))
        allocate (values(9, npoib, 2))
        values(:, :, 1) = real(response)
        values(:, :, 2) = aimag(response)
    end subroutine capture_response

    subroutine check_accepted_response()
        call require(all(ieee_is_finite(params)), 'non-finite accepted profiles')
        call require(all(params(1, :) > 0.0d0) .and. all(params(3:4, :) > 0.0d0), &
                     'accepted profiles lost physical density/temperature')
        call require(all(periodic_amplitudes%accepted_status == 0), &
                     'accepted response contains a normalization guard failure')
        call require(all(periodic_amplitudes%accepted == periodic_amplitudes%trial), &
                     'successful doStep did not commit its complex amplitude')
        call require(all(kim_periodic_scale_modes == periodic_amplitudes%accepted), &
                     'accepted transport and amplitude state disagree')
        call require(all(ieee_is_finite(dqle11)) .and. all(ieee_is_finite(dqli11)), &
                     'non-finite accepted transport')
    end subroutine check_accepted_response

    subroutine capture_transport(values)
        real(8), allocatable, intent(out) :: values(:, :)
        allocate (values(8, npoib))
        values(1, :) = dqle11
        values(2, :) = dqle12
        values(3, :) = dqle21
        values(4, :) = dqle22
        values(5, :) = dqli11
        values(6, :) = dqli12
        values(7, :) = dqli21
        values(8, :) = dqli22
    end subroutine capture_transport

    subroutine compare_transport(expected)
        real(8), intent(in) :: expected(:, :)
        real(8), allocatable :: actual(:, :)
        integer :: component
        call capture_transport(actual)
        do component = 1, 8
            call require(maxval(abs(actual(component, :) - expected(component, :))) <= &
                         1.0d-8 * max(maxval(abs(expected(component, :))), 1.0d-250), &
                         'failed candidate left stale transport')
        end do
    end subroutine compare_transport

    subroutine compare_profiles(actual, expected, message)
        real(8), intent(in) :: actual(:, :), expected(:, :)
        character(*), intent(in) :: message
        real(8) :: relative
        relative = maxval(abs(actual - expected) / max(abs(expected), 1.0d-250))
        if (relative > 1.0d-10) then
            print *, 'profile relative error:', relative
            call require(.false., message)
        end if
    end subroutine compare_profiles

    subroutine require(condition, message)
        logical, intent(in) :: condition
        character(*), intent(in) :: message
        if (.not. condition) then
            print *, 'FAIL: ', message
            error stop 1
        end if
    end subroutine require
    subroutine write_config(path)
        character(*), intent(in) :: path
        integer :: unit
        open (newunit=unit, file=path, status='replace', action='write')
        write (unit, '(A)') &
            '&KIM_CONFIG', &
            ' number_of_ion_species = 1', &
            ' artificial_debye_case = 0', &
            " type_of_run = 'electrostatic_periodic'", &
            " collision_model = 'FokkerPlanck'", &
            " ion_collision_model = 'FokkerPlanck'", &
            ' read_species_from_namelist = .false.', &
            " plasma_type = 'D'", &
            ' turn_off_electrons = .false.', &
            ' turn_off_ions = .false.', &
            ' rescale_density = .false.', &
            ' ion_flr_scale_factor = 1.0', &
            '/', &
            '&WKB_DISPERSION', '/', &
            '&KIM_IO', &
            " profile_location = './'", &
            " output_path = './periodic_time_evolution_output/'", &
            ' hdf5_input = .false.', &
            ' hdf5_output = .false.', &
            ' log_level = 3', &
            ' data_verbosity = 0', &
            ' calculate_asymptotics = .false.', &
            ' write_diagnostics_dat = .false.', &
            '/', &
            '&KIM_SETUP', &
            ' btor = -18000.0', &
            ' R0 = 165.0', &
            ' m_mode = -6', &
            ' n_mode = 2', &
            ' omega = 0.0', &
            ' spline_base = 1', &
            ' type_br_field = 12', &
            ' collisions_off = .false.', &
            ' set_profiles_constant = 0', &
            ' bc_type = 3', &
            ' mphi_max = 0', &
            ' Br_boundary_re = 1.0', &
            ' Br_boundary_im = 0.0', &
            '/', &
            '&KIM_GRID', &
            " grid_spacing_rg = 'equidistant'", &
            " grid_spacing_xl = 'equidistant'", &
            ' l_space_dim = 32', &
            ' rg_space_dim = 32', &
            " theta_integration = 'GaussLegendre'", &
            ' Larmor_skip_factor = 5', &
            ' gauss_int_nodes_Ntheta = 17', &
            ' gauss_int_nodes_Nx = 31', &
            ' gauss_int_nodes_Nxp = 30', &
            ' r_plas = 50.0', &
            ' r_min = 10.0', &
            ' width_res = 0.5', &
            ' ampl_res = 15.0', &
            '/', &
            '&KIM_PROFILES', '/', &
            '&KIM_PERIODIC', &
            ' periodic_dr_asis_scale = 5.0', &
            ' periodic_dr_tr_scale = 10.0', &
            ' periodic_kmax_scale = 0.4', &
            ' periodic_n_rg = 64', &
            ' periodic_Bparallel_ratio = (0.0, 0.0)', &
            '/'
        close (unit)
    end subroutine write_config
end program test_periodic_time_evolution
