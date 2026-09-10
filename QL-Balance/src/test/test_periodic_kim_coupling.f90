program test_periodic_kim_coupling
    !! Exercise the actual periodic KIM solver through the QL-Balance adapter.
    !! The half-open Fourier grid must embed without changing input profiles,
    !! and repeating a stationary solve must reproduce the same response.
    use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
    use kim_wave_code_adapter_m, only: kim_initialize, kim_run_for_all_modes, &
        kim_get_wave_fields, kim_D_ion_modes, kim_transition_weights, &
        kim_embedding_metadata, kim_update_profiles, kim_get_wave_vectors
    use control_mod, only: wave_code, kim_config_path, kim_profiles_from_balance, &
        type_of_run, kim_run_type, kim_transport_benchmark, kim_ion_transport_model
    use wave_code_data, only: dim_mn, m_vals, n_vals, r, n, Te, Ti, q, &
        Vth, Vz, dPhi0, Es, Ep, Er, Et, Ez, Br, Bp, B0, nue, nui
    use plasma_parameters, only: params_b
    use baseparam_mod, only: ev, rtor
    use periodic_transport_benchmark_m, only: select_periodic_ion_transport, &
        reset_transport_benchmark, write_transport_benchmark
    implicit none

    integer, parameter :: npts = 401
    real(8) :: input_r(npts), input_profiles(npts, 5)
    real(8) :: initial_tensor(2, 2, npts), initial_weights(npts)
    real(8) :: initial_nue(npts), initial_nui(npts), initial_b0(npts)
    complex(8) :: initial_es(npts), initial_br(npts)
    integer :: i
    external :: allocate_wave_code_data

    dim_mn = 1
    allocate(m_vals(1), n_vals(1))
    m_vals = -6
    n_vals = 2
    do i = 1, npts
        input_r(i) = 2.0d0 + 58.0d0*real(i - 1, 8)/real(npts - 1, 8)
    end do
    call allocate_wave_code_data(npts, input_r)
    n = 2.0d13*(1.1d0 - r/100.0d0)
    Te = 1.0d3*(1.2d0 - r/100.0d0)
    Ti = Te
    q = 1.5d0 + 2.5d0*(r - 2.0d0)/58.0d0
    dPhi0 = 0.5d0
    Vth = 0.0d0
    Vz = 0.0d0
    input_profiles(:, 1) = n
    input_profiles(:, 2) = Te
    input_profiles(:, 3) = Ti
    input_profiles(:, 4) = q
    input_profiles(:, 5) = dPhi0

    wave_code = 'KIM'
    kim_run_type = 'electrostatic_periodic'
    kim_profiles_from_balance = .true.
    type_of_run = 'SingleStep'
    kim_config_path = 'KIM_config_periodic_coupling.nml'
    call write_config(trim(kim_config_path))
    call kim_initialize(npts, input_r)
    call kim_run_for_all_modes()
    call kim_get_wave_fields(1)
    call check_response()
    call check_input_profiles()
    call check_transport_benchmark()
    initial_tensor = kim_D_ion_modes(:, :, :, 1)
    initial_weights = kim_transition_weights(:, 1)
    initial_es = Es
    initial_br = Br
    initial_nue = nue
    initial_nui = nui
    initial_b0 = B0

    call kim_run_for_all_modes()
    call kim_get_wave_fields(1)
    call check_response()
    call check_input_profiles()
    call require(maxval(abs(kim_transition_weights(:, 1) - initial_weights)) < 1.0d-10, &
        'repeated stationary solve moved the periodic support')
    call require(maxval(abs(Br - initial_br)) < 1.0d-10, &
        'repeated stationary solve changed the magnetic response')
    call require(maxval(abs(Es - initial_es)) < &
        1.0d-8*max(maxval(abs(initial_es)), tiny(1.0d0)), &
        'repeated stationary solve changed the electric response')
    call require(maxval(abs(kim_D_ion_modes(:, :, :, 1) - initial_tensor)) < &
        1.0d-8*max(maxval(abs(initial_tensor)), tiny(1.0d0)), &
        'repeated stationary solve changed the ion tensor')
    call require(maxval(abs(B0 - initial_b0)) < 1.0d-10*maxval(abs(initial_b0)), &
        'repeated solve changed the global magnetic background')
    call require(maxval(abs(nue - initial_nue)) < 1.0d-10*maxval(abs(initial_nue)), &
        'repeated solve changed the global electron collision profile')
    call require(maxval(abs(nui - initial_nui)) < 1.0d-10*maxval(abs(initial_nui)), &
        'repeated solve changed the global ion collision profile')

    ! A fresh transport profile must replace the saved background, rather than
    ! be silently overwritten by the stationary-solve restoration mechanism.
    rtor = 165.0d0
    allocate(params_b(4, npts))
    params_b(1, :) = 1.5d0*input_profiles(:, 1)
    params_b(2, :) = 0.0d0
    params_b(3, :) = input_profiles(:, 2)*ev
    params_b(4, :) = input_profiles(:, 3)*ev
    input_profiles(:, 1) = params_b(1, :)
    input_profiles(:, 2) = params_b(3, :)/ev
    input_profiles(:, 3) = params_b(4, :)/ev
    call kim_update_profiles()
    call kim_run_for_all_modes()
    call kim_get_wave_fields(1)
    call check_response()
    call check_input_profiles()
    call require(all(nue > 1.2d0*initial_nue) .and. all(nue < 1.8d0*initial_nue), &
        'electron collision profile did not follow the density increase')
    call require(all(nui > 1.2d0*initial_nui) .and. all(nui < 1.8d0*initial_nui), &
        'ion collision profile did not follow the density increase')
    call require(maxval(abs(kim_D_ion_modes(:, :, :, 1) - initial_tensor)) > &
        1.0d-5*maxval(abs(initial_tensor)), 'density feedback did not change the ion tensor')
    print *, 'PASS: periodic KIM fields and tensor reach QL-Balance reproducibly'

contains

    subroutine check_transport_benchmark()
        use grid_mod, only: rb, r_resonant, gg_width
        use QLBalance_diag, only: i_mn_loop
        use wave_code_data, only: om_E, ks
        use baseparam_mod, only: c, p_mass, am
        use h5mod, only: h5_create, h5_close, h5_open, h5_deinit, h5_get, &
            h5_id, h5_mode_groupname, path2out, h5_obj_exists
        real(8) :: selected(2, 2, npts), disabled(2, 2, npts)
        real(8) :: old(2, 2, npts), new(2, 2, npts), residual(2, 2, npts)
        real(8) :: relative(2, 2, npts), vt(npts), expected(2, 2, npts)
        complex(8) :: saved_es(npts), saved_br(npts)
        character(32) :: models(2)
        character(256) :: group
        logical :: exists
        integer :: model

        allocate(rb(npts), r_resonant(1))
        rb = r
        r_resonant = sum(kim_embedding_metadata(1:2, 1))/2.0d0
        gg_width = maxval(r) - minval(r)
        i_mn_loop = 1
        call kim_get_wave_vectors(1)
        om_E = ks*c*dPhi0/B0
        am = 2.0d0
        vt = sqrt(Ti*ev/(p_mass*am))
        saved_es = Es
        saved_br = Br
        models = [character(32) :: 'finite_larmor_radius', 'drift_kinetic']
        path2out = 'periodic_transport_benchmark.h5'
        h5_mode_groupname = 'coupling'
        call h5_create(trim(path2out), h5_id)
        call h5_close(h5_id)
        call h5_deinit()

        do model = 1, 2
            kim_ion_transport_model = models(model)
            kim_transport_benchmark = .false.
            call reset_transport_benchmark()
            call select_periodic_ion_transport(1, vt, nui, disabled)
            kim_transport_benchmark = .true.
            call select_periodic_ion_transport(1, vt, nui, selected)
            call require(all(selected == disabled), 'benchmark changed selected ion transport')
            call require(all(Es == saved_es) .and. all(Br == saved_br), &
                'benchmark changed physical fields')
            call write_transport_benchmark(model)
            ! Repeated output at the same time index must replace the existing group.
            call write_transport_benchmark(model)
            write(group, '(A,I0,A)') '/coupling/TransportBenchmark/', model, '/mode_1/'
            call h5_open(trim(path2out), h5_id)
            call h5_get(h5_id, trim(group)//'drift_kinetic', old)
            call h5_get(h5_id, trim(group)//'finite_larmor_radius', new)
            call h5_get(h5_id, trim(group)//'absolute_residual', residual)
            call h5_get(h5_id, trim(group)//'relative_residual', relative)
            call h5_close(h5_id)
            call h5_deinit()
            call require(all(new == kim_D_ion_modes(:, :, :, 1)), &
                'saved benchmark is not the actual spectral ion tensor')
            expected = abs(new - old)
            call require(all(residual == expected), 'incorrect saved absolute residual')
            where (max(abs(old), abs(new)) > 0.0d0)
                expected = expected/max(abs(old), abs(new))
            elsewhere
                expected = 0.0d0
            end where
            call require(maxval(abs(relative - expected)) < 1.0d-14, &
                'incorrect saved relative residual')
            if (model == 1) call require(all(selected == new), 'FLR selection changed')
            if (model == 2) call require(all(selected == old), 'drift selection changed')
        end do
        kim_transport_benchmark = .false.
        call reset_transport_benchmark()
        call write_transport_benchmark(3)
        call h5_open(trim(path2out), h5_id)
        call h5_obj_exists(h5_id, '/coupling/TransportBenchmark/3', exists)
        call h5_close(h5_id)
        call h5_deinit()
        call require(.not. exists, 'disabled benchmark wrote stale data')
        kim_ion_transport_model = 'finite_larmor_radius'
        deallocate(rb, r_resonant)
    end subroutine check_transport_benchmark


    subroutine check_response()
        real(8) :: tensor(2, 2), tensor_scale, symmetric_cross, determinant, weight
        real(8) :: core_lo, core_hi, requested_width, sampled_width
        complex(8) :: fields(7)
        integer :: j

        call require(allocated(kim_D_ion_modes), 'adapter did not return ion tensor')
        call require(all(ieee_is_finite(B0)) .and. all(B0 > 1.7d4) .and. &
            all(B0 < 1.9d4), 'global magnetic field corrupted by periodic background')
        call require(all(ieee_is_finite(nue)) .and. all(nue > 0.0d0), &
            'global electron collision profile is not finite and positive')
        call require(all(ieee_is_finite(nui)) .and. all(nui > 0.0d0), &
            'global ion collision profile is not finite and positive')
        ! These smooth density/temperature profiles vary by less than a factor
        ! three, so their collision profiles cannot develop orders of magnitude
        ! of artificial radial variation from extrapolating a narrow window.
        call require(maxval(nue) < 10.0d0*minval(nue), 'unphysical electron collision range')
        call require(maxval(nui) < 10.0d0*minval(nui), 'unphysical ion collision range')
        call require(allocated(kim_embedding_metadata), 'embedding metadata missing')
        core_lo = kim_embedding_metadata(1, 1)
        core_hi = kim_embedding_metadata(2, 1)
        requested_width = kim_embedding_metadata(3, 1)
        sampled_width = kim_embedding_metadata(4, 1)
        call require(sampled_width > 0.0d0 .and. sampled_width < requested_width, &
            'half-open Fourier grid did not limit the available transition width')
        call require(all(ieee_is_finite(kim_D_ion_modes)), 'non-finite ion tensor')
        call require(all(ieee_is_finite(kim_transition_weights)), 'non-finite transition')
        call require(all(kim_transition_weights >= 0.0d0) .and. &
            all(kim_transition_weights <= 1.0d0), 'transition weights outside [0,1]')
        call require(any(kim_transition_weights(:, 1) == 1.0d0), &
            'test grid does not sample the trusted core')
        call require(any(kim_transition_weights(:, 1) > 0.0d0 .and. &
            kim_transition_weights(:, 1) < 1.0d0), 'test grid does not sample transition')
        call require(any(kim_transition_weights(:, 1) == 0.0d0), &
            'test grid does not sample outside support')
        call require(maxval(abs(Es)) > tiny(1.0d0), 'driven electric response vanished')
        call require(maxval(abs(kim_D_ion_modes)) > tiny(1.0d0), 'ion response vanished')
        call require(all(ieee_is_finite(real(Es))) .and. all(ieee_is_finite(aimag(Es))), &
            'non-finite perpendicular electric response')

        do j = 1, npts
            weight = kim_transition_weights(j, 1)
            if (r(j) <= core_lo - sampled_width .or. r(j) >= core_hi + sampled_width) &
                call require(weight == 0.0d0, 'transition exceeds its reported sampled support')
            fields = [Es(j), Ep(j), Er(j), Et(j), Ez(j), Br(j), Bp(j)]
            call require(all(ieee_is_finite(real(fields))) .and. &
                all(ieee_is_finite(aimag(fields))), 'non-finite embedded physical field')
            ! Br is prescribed as one Gauss in this actual KIM solve. Its
            ! compact embedding is therefore independently known at every point.
            call require(abs(Br(j) - cmplx(weight, 0.0d0, 8)) < 1.0d-12, &
                'prescribed radial field was not embedded with the field weight')
            if (weight == 0.0d0) then
                call require(abs(Es(j)) + abs(Ep(j)) + abs(Er(j)) + abs(Et(j)) + &
                    abs(Ez(j)) + abs(Br(j)) + abs(Bp(j)) == 0.0d0, &
                    'physical perturbation leaked outside compact support')
                call require(all(kim_D_ion_modes(:, :, j, 1) == 0.0d0), &
                    'ion transport leaked outside compact support')
            end if
            ! A real transport tensor may have an antisymmetric cross part.
            ! Dissipation depends on its symmetric part, whose principal
            ! minors must be nonnegative up to relative numerical roundoff.
            tensor = kim_D_ion_modes(:, :, j, 1)
            tensor_scale = maxval(abs(tensor))
            if (tensor_scale <= tiny(1.0d0)) cycle
            tensor = tensor/tensor_scale
            symmetric_cross = 0.5d0*(tensor(1, 2) + tensor(2, 1))
            determinant = tensor(1, 1)*tensor(2, 2) - symmetric_cross**2
            if (min(tensor(1, 1), tensor(2, 2), determinant) < -1.0d-8) then
                print *, 'Non-positive transport at r = ', r(j), ' tensor = ', tensor
                error stop 'periodic ion tensor has a negative dissipative direction'
            end if
        end do
    end subroutine check_response

    subroutine check_input_profiles()
        call require(all(r == input_r), 'periodic solve changed QL-Balance radial grid')
        call require(all(n == input_profiles(:, 1)), 'periodic solve changed input density')
        call require(all(Te == input_profiles(:, 2)), 'periodic solve changed electron temperature')
        call require(all(Ti == input_profiles(:, 3)), 'periodic solve changed ion temperature')
        call require(all(q == input_profiles(:, 4)), 'periodic solve changed safety factor')
        call require(all(dPhi0 == input_profiles(:, 5)), 'periodic solve changed equilibrium Er')
    end subroutine check_input_profiles

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
        open(newunit=unit, file=path, status='replace', action='write')
        write(unit, '(A)') &
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
            " output_path = './periodic_coupling_output/'", &
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
            ' periodic_n_rg = 32', &
            ' periodic_Bparallel_ratio = (0.0, 0.0)', &
            '/'
        close(unit)
    end subroutine write_config
end program test_periodic_kim_coupling
