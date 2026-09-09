program test_periodic_target_current
    !! Actual two-mode KIM responses, normalized before QL transport assembly.
    use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
    use kim_wave_code_adapter_m, only: kim_initialize, kim_run_for_all_modes, &
        kim_get_wave_fields, kim_get_wave_vectors, kim_get_current_densities, &
        kim_D_ion_modes, kim_embedding_metadata, &
        kim_periodic_scale_modes, kim_periodic_current_unit, kim_periodic_scale_status
    use control_mod, only: wave_code, kim_config_path, kim_profiles_from_balance, &
        type_of_run, kim_run_type, ihdf5IO, &
        kim_current_relaxation, kim_current_floor, kim_current_max_scale
    use wave_code_data, only: dim_mn, m_vals, n_vals, r, n, Te, Ti, q, &
        Vth, Vz, dPhi0, Es, Ep, Er, Et, Ez, Br, Bp, B0, nue, om_E, ks, &
        Jpe, Jpi, I_par_toroidal, antenna_factor
    use periodic_current_diagnostics_m, only: write_periodic_current_diagnostics
    use h5mod, only: h5_create, h5_close, h5_open, h5_deinit, h5_get, &
        h5_id, h5_mode_groupname, path2out, h5_obj_exists, h5_get_bounds, &
        h5_define_group, h5_close_group, group_id_1
    use setup_m, only: Br_boundary_re, Br_boundary_im
    use baseparam_mod, only: ev, c, e_mass
    use grid_mod, only: rb, r_resonant, gg_width
    use QLBalance_diag, only: i_mn_loop
    use transp_coeffs_mod, only: compute_antenna_factor_from_Ipar
    implicit none

    integer, parameter :: npts = 401, nmodes = 2
    real(8), parameter :: target = 0.25d0, speed_of_light = 2.99792458d10
    type :: response_t
        complex(8) :: fields(7, npts, nmodes), currents(2, npts, nmodes)
        complex(8) :: unit_current(nmodes), scale(nmodes)
        real(8) :: ion(2, 2, npts, nmodes), electron(2, 2, npts, nmodes)
    end type response_t
    type(response_t) :: baseline, doubled, zero_drive, changed_drive
    real(8) :: input_r(npts)
    integer :: point, mode, active
    logical :: exists
    external :: allocate_wave_code_data

    dim_mn = nmodes
    allocate (m_vals(nmodes), n_vals(nmodes))
    m_vals = [-6, -7]
    n_vals = [2, 2]
    do point = 1, npts
        input_r(point) = 2.0d0 + 58.0d0*real(point - 1, 8)/real(npts - 1, 8)
    end do
    call allocate_wave_code_data(npts, input_r)
    n = 2.0d13*(1.1d0 - r/100.0d0)
    Te = 1.0d3*(1.2d0 - r/100.0d0)
    Ti = Te
    ! Both resonant windows remain inside the configured [10,50] plasma grid.
    q = 1.5d0 + 3.0d0*(r - 2.0d0)/58.0d0
    dPhi0 = 0.5d0
    Vth = 0.0d0
    Vz = 0.0d0
    allocate (rb(npts), r_resonant(nmodes))
    rb = r
    gg_width = maxval(r) - minval(r)

    wave_code = 'KIM'
    kim_run_type = 'electrostatic_periodic'
    kim_profiles_from_balance = .true.
    type_of_run = 'SingleStep'
    kim_config_path = 'KIM_config_periodic_target.nml'
    call write_config(trim(kim_config_path))
    call kim_initialize(npts, input_r)

    ihdf5IO = 1
    kim_current_relaxation = 1.0d0
    path2out = 'periodic_target_current.h5'
    h5_mode_groupname = 'target_current'
    call h5_create(trim(path2out), h5_id)
    call h5_define_group(h5_id, trim(h5_mode_groupname), group_id_1)
    call h5_close_group(group_id_1)
    call h5_close(h5_id)
    call h5_deinit()

    call solve_and_capture(target, (1.0d0, 0.0d0), baseline)
    call check_saved_current(1, target, baseline)
    call require(abs(baseline%unit_current(1) - baseline%unit_current(2)) > &
        1.0d-3*maxval(abs(baseline%unit_current)), 'two-mode fixture has identical unit currents')
    call solve_and_capture(2.0d0*target, (1.0d0, 0.0d0), doubled)
    call compare_scaled_response(doubled, baseline, 2.0d0)
    call check_saved_current(2, 2.0d0*target, doubled)

    ! A positive target defines its own unit-drive solve. The arbitrary input
    ! drive, including exactly zero, must neither suppress it nor survive as
    ! an extra factor. Its configured value must be restored after both modes.
    call solve_and_capture(target, (0.0d0, 0.0d0), zero_drive)
    call compare_scaled_response(zero_drive, baseline, 1.0d0)
    call check_saved_current(3, target, zero_drive)
    call solve_and_capture(target, (2.0d0, -0.7d0), changed_drive)
    call compare_scaled_response(changed_drive, baseline, 1.0d0)
    call check_saved_current(4, target, changed_drive)

    call check_guarded_response()

    I_par_toroidal = 0.0d0
    antenna_factor = 9.0d0
    call compute_antenna_factor_from_Ipar()
    call require(antenna_factor == 9.0d0, 'zero target erased the manual antenna factor')
    call kim_run_for_all_modes()
    do mode = 1, nmodes
        call kim_get_wave_fields(mode)
        call require(count(r >= kim_embedding_metadata(1, mode) .and. &
            r <= kim_embedding_metadata(2, mode)) > 0, 'manual fixture has no core samples')
        call require(maxval(abs(Br - cmplx(2.0d0, -0.7d0, 8)), &
            mask=r >= kim_embedding_metadata(1, mode) .and. &
            r <= kim_embedding_metadata(2, mode)) < 1.0d-12, &
            'manual periodic solve ignored its configured magnetic drive')
    end do
    ! Reusing a saved index after disabling normalization must remove its
    ! previous guard record, in addition to avoiding new manual records.
    call write_periodic_current_diagnostics(5)
    call write_periodic_current_diagnostics(6)
    call h5_open(trim(path2out), h5_id)
    call h5_get(h5_id, '/target_current/CurrentNormalization/5/active', active)
    call require(active == 0, 'manual drive retained an active normalization record')
    call h5_obj_exists(h5_id, &
        '/target_current/CurrentNormalization/5/mode_1/target_current', exists)
    call require(.not. exists, 'manual drive retained a stale same-index target')
    call h5_obj_exists(h5_id, &
        '/target_current/CurrentNormalization/5/mode_2/normalized_jpar', exists)
    call require(.not. exists, 'manual drive retained a stale same-index current profile')
    call h5_obj_exists(h5_id, '/target_current/CurrentNormalization/6', exists)
    call require(.not. exists, 'manual drive wrote a stale target-normalization record')
    call h5_close(h5_id)
    call h5_deinit()
    call require(Br_boundary_re == 2.0d0 .and. Br_boundary_im == -0.7d0, &
        'manual solve changed configured Br')
    I_par_toroidal = -1.0d0
    antenna_factor = 4.0d0
    call compute_antenna_factor_from_Ipar()
    call require(antenna_factor == 4.0d0, 'negative target erased the manual antenna factor')
    print *, 'PASS: two periodic modes reach target current with linear/quadratic scaling'

contains

    subroutine solve_and_capture(requested, configured_drive, result)
        real(8), intent(in) :: requested
        complex(8), intent(in) :: configured_drive
        type(response_t), intent(out) :: result
        real(8) :: vt(npts), core_lo, core_hi
        complex(8) :: measured_current
        integer :: mode

        I_par_toroidal = requested
        Br_boundary_re = real(configured_drive)
        Br_boundary_im = aimag(configured_drive)
        call kim_run_for_all_modes()
        call require(Br_boundary_re == real(configured_drive) .and. &
            Br_boundary_im == aimag(configured_drive), 'normalization changed configured Br')
        call require(all(kim_periodic_scale_status == 0), 'positive target was not normalized')
        result%unit_current = kim_periodic_current_unit
        result%scale = kim_periodic_scale_modes
        do mode = 1, nmodes
            i_mn_loop = mode
            core_lo = kim_embedding_metadata(1, mode)
            core_hi = kim_embedding_metadata(2, mode)
            r_resonant(mode) = 0.5d0*(core_lo + core_hi)
            call kim_get_wave_fields(mode)
            call kim_get_wave_vectors(mode)
            call kim_get_current_densities(mode)
            result%fields(:, :, mode) = transpose(reshape([Es, Ep, Er, Et, Ez, Br, Bp], &
                [npts, 7]))
            result%currents(1, :, mode) = Jpe
            result%currents(2, :, mode) = Jpi
            result%ion(:, :, :, mode) = kim_D_ion_modes(:, :, :, mode)
            ! Use the established electron transport routine with the scaled
            ! physical E and Br, just as QL-Balance does after selecting a mode.
            om_E = ks*c*dPhi0/B0
            vt = sqrt(Te*ev/e_mass)
            call calc_transport_coeffs_ornuhl(npts, vt, nue, &
                result%electron(1, 1, :, mode), result%electron(1, 2, :, mode), &
                result%electron(2, 1, :, mode), result%electron(2, 2, :, mode))
            call require(all(ieee_is_finite(result%ion(:, :, :, mode))), 'non-finite ion tensor')
            call require(all(ieee_is_finite(result%electron(:, :, :, mode))), &
                'non-finite electron tensor')
            call require(maxval(abs(result%ion(:, :, :, mode))) > tiny(1.0d0), &
                'target-current ion tensor vanished')
            call require(maxval(abs(result%electron(:, :, :, mode))) > tiny(1.0d0), &
                'target-current electron tensor vanished')

            ! Independently integrate the currents actually supplied to
            ! QL-Balance. This never uses the stored unit current or scale.
            measured_current = clipped_output_integral(Jpe + Jpi, core_lo, core_hi) &
                *(2.0d0*acos(-1.0d0)/speed_of_light)
            print *, 'mode, target, independently integrated current:', mode, requested, &
                measured_current
            call require(abs(measured_current - cmplx(requested, 0.0d0, 8)) < &
                0.01d0*requested, 'embedded current misses the target beyond grid accuracy')
            antenna_factor = 9.0d0
            call compute_antenna_factor_from_Ipar()
            call require(antenna_factor == 1.0d0, 'positive target would be antenna-scaled twice')
        end do
    end subroutine solve_and_capture

    subroutine check_guarded_response()
        real(8) :: unit_current(2), achieved(2), relative, saved_floor
        integer :: mode, status
        character(160) :: group

        I_par_toroidal = target
        saved_floor = kim_current_floor
        kim_current_floor = 1.0d100
        call kim_run_for_all_modes()
        call require(all(kim_periodic_scale_status == 2), 'actual floor guard did not activate')
        do mode = 1, nmodes
            call kim_get_wave_fields(mode)
            call kim_get_current_densities(mode)
            call require(all([Es, Ep, Er, Et, Ez, Br, Bp, Jpe, Jpi] == (0.0d0, 0.0d0)), &
                'guarded response supplied nonzero or non-finite fields/current')
            call require(all(kim_D_ion_modes(:, :, :, mode) == 0.0d0), &
                'guarded response supplied a nonzero or non-finite ion tensor')
        end do
        call write_periodic_current_diagnostics(5)
        call h5_open(trim(path2out), h5_id)
        do mode = 1, nmodes
            write(group, '(A,I0,A)') '/target_current/CurrentNormalization/5/mode_', mode, '/'
            call h5_get(h5_id, trim(group)//'status', status)
            call h5_get(h5_id, trim(group)//'achieved_current', achieved)
            call h5_get(h5_id, trim(group)//'unit_current', unit_current)
            call h5_get(h5_id, trim(group)//'relative_residual', relative)
            call require(status == 2 .and. all(achieved == 0.0d0) .and. relative == 1.0d0, &
                'saved guard output misrepresents its suppressed response')
            call require(all(ieee_is_finite(unit_current)) .and. &
                maxval(abs(unit_current)) > 0.0d0, 'guard lost the finite raw unit current')
        end do
        call h5_close(h5_id)
        call h5_deinit()
        kim_current_floor = saved_floor
    end subroutine check_guarded_response

    subroutine check_saved_current(time_index, requested, response)
        integer, intent(in) :: time_index
        real(8), intent(in) :: requested
        type(response_t), intent(in) :: response
        real(8), allocatable :: radius(:), profile(:, :), unit_profile(:, :)
        real(8) :: bounds(2), achieved(2), residual(2), unit_current(2), scale(2)
        real(8) :: saved_target, relative_residual, relaxation, current_floor, max_scale_ratio
        complex(8) :: measured, measured_unit, saved_achieved
        character(160) :: group
        integer :: mode, status, lower_bound, upper_bound

        call write_periodic_current_diagnostics(time_index)
        ! Repeating a saved output index must remain safe and deterministic.
        call write_periodic_current_diagnostics(time_index)
        call h5_open(trim(path2out), h5_id)
        do mode = 1, nmodes
            write (group, '(A,I0,A,I0,A)') '/target_current/CurrentNormalization/', &
                time_index, '/mode_', mode, '/'
            call h5_get_bounds(h5_id, trim(group)//'r', lower_bound, upper_bound)
            allocate (radius(upper_bound - lower_bound + 1))
            allocate (profile(size(radius), 2), unit_profile(size(radius), 2))
            call h5_get(h5_id, trim(group)//'r', radius)
            call h5_get(h5_id, trim(group)//'core_bounds', bounds)
            call h5_get(h5_id, trim(group)//'normalized_jpar', profile)
            call h5_get(h5_id, trim(group)//'unit_jpar', unit_profile)
            call h5_get(h5_id, trim(group)//'unit_current', unit_current)
            call h5_get(h5_id, trim(group)//'achieved_current', achieved)
            call h5_get(h5_id, trim(group)//'residual', residual)
            call h5_get(h5_id, trim(group)//'relative_residual', relative_residual)
            call h5_get(h5_id, trim(group)//'target_current', saved_target)
            call h5_get(h5_id, trim(group)//'scale', scale)
            call h5_get(h5_id, trim(group)//'status', status)
            call h5_get(h5_id, trim(group)//'relaxation', relaxation)
            call h5_get(h5_id, trim(group)//'current_floor', current_floor)
            call h5_get(h5_id, trim(group)//'max_scale_ratio', max_scale_ratio)
            call require(relaxation == kim_current_relaxation .and. &
                current_floor == kim_current_floor .and. max_scale_ratio == kim_current_max_scale, &
                'saved guard and relaxation settings differ from the solve configuration')
            call require(saved_target == requested .and. status == 0, &
                'saved target/status does not describe the actual solve')
            call require(size(profile, 1) == size(radius) .and. size(profile, 2) == 2, &
                'saved complex current profile has incorrect shape')
            call require(all(ieee_is_finite(profile)), 'non-finite saved normalized current')
            ! Integrate the saved native current profile independently. Using
            ! its piecewise-linear r*J interpolant avoids a second embedding
            ! error and checks the production quadrature's clipped core.
            measured = native_integral(radius, cmplx(profile(:, 1), profile(:, 2), 8), &
                bounds(1), bounds(2))*(2.0d0*acos(-1.0d0)/speed_of_light)
            measured_unit = native_integral(radius, &
                cmplx(unit_profile(:, 1), unit_profile(:, 2), 8), bounds(1), bounds(2))
            saved_achieved = cmplx(achieved(1), achieved(2), 8)
            call require(abs(measured - cmplx(requested, 0.0d0, 8)) < 1.0d-10*requested, &
                'saved physical current profile does not integrate to target')
            call require(abs(saved_achieved - measured) < 1.0d-10*requested, &
                'saved achieved current disagrees with native profile integration')
            call require(abs(cmplx(residual(1), residual(2), 8) &
                - (measured - requested)) < 1.0d-10*requested, &
                'saved complex residual disagrees with native profile integration')
            call require(abs(relative_residual - abs(measured - requested)/requested) &
                < 1.0d-10, 'saved relative residual is inconsistent')
            call require(abs(measured_unit - cmplx(unit_current(1), unit_current(2), 8)) &
                < 1.0d-10*abs(measured_unit), 'saved unit-current quadrature is inconsistent')
            call require(abs(measured_unit - response%unit_current(mode)) &
                < 1.0d-10*abs(measured_unit), 'saved record belongs to a different mode')
            call require(abs(cmplx(scale(1), scale(2), 8) - response%scale(mode)) &
                < 1.0d-10*abs(response%scale(mode)), 'saved complex scale is inconsistent')
            deallocate (radius, profile, unit_profile)
        end do
        call h5_close(h5_id)
        call h5_deinit()
    end subroutine check_saved_current

    function native_integral(radius, current, core_lo, core_hi) result(integral)
        real(8), intent(in) :: radius(:), core_lo, core_hi
        complex(8), intent(in) :: current(:)
        complex(8) :: integral
        integral = native_primitive(radius, current, core_hi) &
            - native_primitive(radius, current, core_lo)
    end function native_integral

    function native_primitive(radius, current, endpoint) result(value)
        real(8), intent(in) :: radius(:), endpoint
        complex(8), intent(in) :: current(:)
        complex(8) :: value, left_value, slope
        real(8) :: dx, h
        integer :: i
        value = (0.0d0, 0.0d0)
        do i = 1, size(radius) - 1
            dx = min(endpoint, radius(i + 1)) - radius(i)
            if (dx <= 0.0d0) exit
            h = radius(i + 1) - radius(i)
            left_value = radius(i)*current(i)
            slope = (radius(i + 1)*current(i + 1) - left_value)/h
            value = value + dx*left_value + 0.5d0*dx**2*slope
        end do
    end function native_primitive

    function clipped_output_integral(current, core_lo, core_hi) result(integral)
        complex(8), intent(in) :: current(npts)
        real(8), intent(in) :: core_lo, core_hi
        complex(8) :: integral, left_current, right_current
        real(8) :: left, right, midpoint, h
        integer :: i

        integral = (0.0d0, 0.0d0)
        do i = 1, npts - 1
            left = max(r(i), core_lo)
            right = min(r(i + 1), core_hi)
            if (right <= left) cycle
            h = r(i + 1) - r(i)
            left_current = current(i) + (current(i + 1) - current(i))*(left - r(i))/h
            right_current = current(i) + (current(i + 1) - current(i))*(right - r(i))/h
            midpoint = 0.5d0*(left + right)
            ! Simpson's rule is exact for r times the linear current on each
            ! clipped interval; neither endpoint can borrow outside-core area.
            integral = integral + (right - left)/6.0d0*(left*left_current &
                + 2.0d0*midpoint*(left_current + right_current) + right*right_current)
        end do
    end function clipped_output_integral

    subroutine compare_scaled_response(actual, reference, amplitude)
        type(response_t), intent(in) :: actual, reference
        real(8), intent(in) :: amplitude
        real(8), parameter :: tolerance = 2.0d-8
        integer :: mode, component, row, column

        do mode = 1, nmodes
            do component = 1, size(actual%fields, 1)
                call require(maxval(abs(actual%fields(component, :, mode) &
                    - amplitude*reference%fields(component, :, mode))) <= &
                    tolerance*amplitude*max(1.0d-250, &
                    maxval(abs(reference%fields(component, :, mode)))), &
                    'physical field component does not scale linearly with target current')
            end do
            do component = 1, size(actual%currents, 1)
                call require(maxval(abs(actual%currents(component, :, mode) &
                    - amplitude*reference%currents(component, :, mode))) <= &
                    tolerance*amplitude*max(1.0d-250, &
                    maxval(abs(reference%currents(component, :, mode)))), &
                    'species current does not scale linearly with target current')
            end do
            do column = 1, 2
                do row = 1, 2
                    call require(maxval(abs(actual%ion(row, column, :, mode) &
                        - amplitude**2*reference%ion(row, column, :, mode))) <= &
                        tolerance*amplitude**2*max(1.0d-250, &
                        maxval(abs(reference%ion(row, column, :, mode)))), &
                        'ion tensor entry does not scale quadratically with target current')
                    call require(maxval(abs(actual%electron(row, column, :, mode) &
                        - amplitude**2*reference%electron(row, column, :, mode))) <= &
                        tolerance*amplitude**2*max(1.0d-250, &
                        maxval(abs(reference%electron(row, column, :, mode)))), &
                        'electron tensor entry does not scale quadratically with target current')
                end do
            end do
            call require(abs(actual%unit_current(mode) - reference%unit_current(mode)) < &
                tolerance*abs(reference%unit_current(mode)), 'unit current depends on input drive')
            call require(abs(actual%scale(mode) - amplitude*reference%scale(mode)) < &
                tolerance*amplitude*abs(reference%scale(mode)), 'incorrect complex drive scale')
        end do
    end subroutine compare_scaled_response

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
            " output_path = './periodic_target_output/'", &
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
end program test_periodic_target_current
