program test_periodic_multimode_assembly
    !! Drive the actual get_dql profile-refresh, KIM solve and accumulation path.
    use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
    use kim_wave_code_adapter_m, only: kim_initialize, kim_vac_Br, kim_transition_weights
    use control_mod, only: wave_code, kim_config_path, kim_profiles_from_balance, &
        type_of_run, kim_run_type, kim_ion_transport_model, kim_transport_benchmark, &
        irf, suppression_mode, misalign_diffusion, ihdf5IO, data_verbosity, &
        write_gyro_current, gyro_current_study, jpar_method
    use wave_code_data, only: dim_mn, m_vals, n_vals, r, n, Te, Ti, q, &
        Vth, Vz, dPhi0, I_par_toroidal, antenna_factor
    use grid_mod, only: rb, r_resonant, gg_width, nbaleqs, npoib, npoic, npoi_der, &
        ipbeg, ipend, deriv_coef, reint_coef, sqrt_g_times_B_theta_over_c, Ercov, mwind, &
        dqle11, dqle12, dqle21, dqle22, dqli11, dqli12, dqli21, dqli22, &
        de11, de12, de21, de22, di11, di12, di21, di22, rb_cut_in, rb_cut_out, re_cut_out
    use plasma_parameters, only: params, params_b, ddr_params_nl
    use setup_m, only: Br_boundary_re, Br_boundary_im
    use baseparam_mod, only: ev, am, rtor, Z_i
    use time_evolution, only: time_ind, save_prof_time_step
    use h5mod, only: path2out, h5_mode_groupname, h5_id, group_id_1, &
        h5_create, h5_define_group, h5_close_group, h5_close, h5_deinit
    implicit none

    integer, parameter :: npts = 401
    real(8) :: input_r(npts), single_a(8, npts), single_b(8, npts), combined(8, npts)
    real(8) :: phased(8, npts), reverse(8, npts)
    logical :: union_support(npts), overlap(npts)
    integer :: point
    external :: allocate_wave_code_data, get_dql

    dim_mn = 2
    allocate(m_vals(2), n_vals(2))
    m_vals = [-6, -9]
    n_vals = [2, 3]
    do point = 1, npts
        input_r(point) = 2.0d0 + 58.0d0*real(point - 1, 8)/real(npts - 1, 8)
    end do
    call allocate_wave_code_data(npts, input_r)
    n = 2.0d13
    Te = 1.0d3
    Ti = 1.0d3
    q = 1.5d0 + 3.0d0*(r - 2.0d0)/58.0d0
    dPhi0 = 0.5d0
    Vth = 0.0d0
    Vz = 0.0d0
    rtor = 165.0d0
    am = 2.0d0
    Z_i = 1.0d0

    ! Constant thermodynamic profiles permit exact collocated interpolation
    ! and zero derivatives. get_dql still executes both interpolation and
    ! input smoothing, then refreshes KIM from these physical QL profiles.
    nbaleqs = 4
    npoib = npts
    npoic = npts
    npoi_der = 1
    mwind = 1
    allocate(rb(npts), ipbeg(npts), ipend(npts), deriv_coef(1, npts), reint_coef(1, npts))
    allocate(sqrt_g_times_B_theta_over_c(npts), Ercov(npts))
    allocate(params(4, npts), params_b(4, npts), ddr_params_nl(4, npts))
    allocate(dqle11(npts), dqle12(npts), dqle21(npts), dqle22(npts))
    allocate(dqli11(npts), dqli12(npts), dqli21(npts), dqli22(npts))
    allocate(de11(npts), de12(npts), de21(npts), de22(npts))
    allocate(di11(npts), di12(npts), di21(npts), di22(npts))
    rb = r
    ipbeg = [(point, point=1, npts)]
    ipend = ipbeg
    deriv_coef = 0.0d0
    reint_coef = 1.0d0
    sqrt_g_times_B_theta_over_c = 0.0d0
    Ercov = 0.0d0
    params(1, :) = n
    params(2, :) = 0.0d0
    params(3, :) = Te*ev
    params(4, :) = Ti*ev
    params_b = params
    ddr_params_nl = 0.0d0
    gg_width = 58.0d0
    rb_cut_in = 0.0d0
    rb_cut_out = 61.0d0
    re_cut_out = 62.0d0

    wave_code = 'KIM'
    kim_run_type = 'electrostatic_periodic'
    kim_ion_transport_model = 'finite_larmor_radius'
    kim_transport_benchmark = .false.
    kim_profiles_from_balance = .true.
    type_of_run = 'SingleStep'
    irf = 1
    suppression_mode = .true.
    misalign_diffusion = .false.
    jpar_method = 'conductivity'
    ihdf5IO = 1
    data_verbosity = 0
    write_gyro_current = .false.
    gyro_current_study = 0
    I_par_toroidal = 0.0d0
    antenna_factor = 1.0d0
    time_ind = 1
    save_prof_time_step = 2
    kim_config_path = 'KIM_config_periodic_multimode_assembly.nml'
    call write_config(trim(kim_config_path))
    call kim_initialize(npts, input_r)

    ! The direct-current diagnostic invoked by get_dql expects an existing
    ! run file even when detailed current output is disabled.
    path2out = 'periodic_multimode_assembly.h5'
    h5_mode_groupname = 'assembly'
    call h5_create(trim(path2out), h5_id)
    call h5_define_group(h5_id, trim(h5_mode_groupname), group_id_1)
    call h5_close_group(group_id_1)
    call h5_close(h5_id)
    call h5_deinit()

    call run_case([-6], [2], (1.0d0, 0.0d0), single_a)
    call run_case([-9], [3], (1.0d0, 0.0d0), single_b)
    call run_case([-6, -9], [2, 3], (1.0d0, 0.0d0), combined)
    overlap = kim_transition_weights(:, 1) > 0.0d0 .and. &
        kim_transition_weights(:, 2) > 0.0d0
    union_support = any(kim_transition_weights > 0.0d0, dim=2)
    call require(count(overlap) > 4, 'distinct modes do not have resolved overlapping support')
    call require(maxval(abs(single_a), mask=spread(overlap, 1, 8)) > 0.0d0 .and. &
        maxval(abs(single_b), mask=spread(overlap, 1, 8)) > 0.0d0, &
        'overlap fixture has a vanishing mode response')
    call compare(combined, single_a + single_b, 'incoherent sum of independently solved modes')
    call require(all(pack(combined, spread(.not. union_support, 1, 8)) == 0.0d0), &
        'get_dql spread transport outside the union of periodic supports')

    call run_case([-9, -6], [3, 2], (1.0d0, 0.0d0), reverse)
    call compare(reverse, combined, 'production accumulation depends on mode order')
    call run_case([-6, -9], [2, 3], (0.6d0, 0.8d0), phased)
    call compare(phased, combined, 'production transport depends on an overall drive phase')
    print *, 'PASS: actual get_dql adds overlapping periodic mode transport incoherently'

contains

    subroutine run_case(modes_m, modes_n, drive, tensor)
        integer, intent(in) :: modes_m(:), modes_n(:)
        complex(8), intent(in) :: drive
        real(8), intent(out) :: tensor(8, npts)
        integer :: mode

        dim_mn = size(modes_m)
        deallocate(m_vals, n_vals)
        allocate(m_vals(dim_mn), n_vals(dim_mn))
        m_vals = modes_m
        n_vals = modes_n
        if (allocated(r_resonant)) deallocate(r_resonant)
        allocate(r_resonant(dim_mn))
        do mode = 1, dim_mn
            r_resonant(mode) = 2.0d0 + 58.0d0*(-real(m_vals(mode), 8)/n_vals(mode) &
                - 1.5d0)/3.0d0
        end do
        Br_boundary_re = real(drive)
        Br_boundary_im = aimag(drive)
        if (allocated(kim_vac_Br)) deallocate(kim_vac_Br)
        allocate(kim_vac_Br(npts, dim_mn))
        ! The periodic source prescribes a constant vacuum radial field.
        ! This only initializes get_dql's diagnostic vacuum form factor;
        ! every plasma response and transport coefficient is solver-produced.
        kim_vac_Br = drive
        call get_dql()
        tensor(1, :) = dqle11
        tensor(2, :) = dqle12
        tensor(3, :) = dqle21
        tensor(4, :) = dqle22
        tensor(5, :) = dqli11
        tensor(6, :) = dqli12
        tensor(7, :) = dqli21
        tensor(8, :) = dqli22
        call require(all(ieee_is_finite(tensor)), 'actual get_dql produced non-finite transport')
        call require(maxval(abs(tensor(1:4, :))) > 0.0d0, 'electron transport vanished')
        call require(maxval(abs(tensor(5:8, :))) > 0.0d0, 'ion transport vanished')
    end subroutine run_case

    subroutine compare(actual, expected, message)
        real(8), intent(in) :: actual(8, npts), expected(8, npts)
        character(*), intent(in) :: message
        real(8) :: error, reference
        integer :: component
        do component = 1, 8
            error = maxval(abs(actual(component, :) - expected(component, :)))
            reference = maxval(abs(expected(component, :)))
            if (error > 2.0d-8*max(reference, 1.0d-250)) then
                print *, 'component, absolute error, reference:', component, error, reference
                call require(.false., message)
            end if
        end do
    end subroutine compare

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
            " output_path = './periodic_multimode_assembly_output/'", &
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
end program test_periodic_multimode_assembly
