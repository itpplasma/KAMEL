program test_state_diagnostics_refresh
    !! After a time step, QL-Balance needs the state-dependent diagnostics (RMP
    !! torque densities, E_r, boundary profiles, ...) at the accepted profiles.
    !! They used to come from an implicit step with dt = 0, i.e. a full Jacobian
    !! probe and a sparse solve of the identity. refresh_state_diagnostics
    !! evaluates only the state part of the right-hand side.
    !!
    !! Asserted: the refresh leaves the profiles untouched and gives exactly the
    !! diagnostics of the dt = 0 step evaluated at the same profiles.
    use, intrinsic :: ieee_arithmetic, only: ieee_value, ieee_quiet_nan
    use QLBalance_kinds, only: dp
    use grid_mod, only: rb, npoib, npoic, rmin, rmax, npoimin, iboutype, gg_factor, &
                        gg_width, set_boundary_condition, mwind, neqset, &
                        dery_equisource, source_term, rb_cut_in, rb_cut_out, re_cut_out, &
                        dqle11, dqle12, dqle21, dqle22, dqli11, dqli12, dqli21, dqli22, &
                        torque_density_e, torque_density_i, Ercov, dni22, polforce, &
                        qlheat_e, qlheat_i
    use plasma_parameters, only: params, params_b, ddr_params_nl, init_background_profiles
    use baseparam_mod, only: am, rtor, Z_i, btor, rsepar
    use paramscan_mod, only: viscosity_factor
    use resonances_mod, only: prop, numres
    use control_mod, only: ihdf5IO
    use wave_code_data, only: r, n, Te, Ti, q, Vth, Vz
    use recstep_mod, only: timstep_arr
    use rhs_balance_m, only: refresh_state_diagnostics
    use h5mod, only: path2inp, h5_id, group_id_1, h5_create, h5_define_group, &
                     h5_close_group, h5_close, h5_deinit, h5_add

    implicit none

    type :: diagnostics_t
        real(dp), allocatable :: torque_e(:), torque_i(:), er(:), dni(:), polforce(:)
        real(dp), allocatable :: qlheat_e(:), qlheat_i(:), params_b(:, :), ddr_params(:, :)
    end type diagnostics_t

    type(diagnostics_t) :: after_zero_step, after_refresh
    real(dp), allocatable :: input_r(:), state(:, :)
    real(dp) :: nan
    integer :: num_failed
    external :: gengrid, allocate_wave_code_data, calc_geometric_parameter_profiles
    external :: evolvestep

    num_failed = 0
    nan = ieee_value(1.0_dp, ieee_quiet_nan)
    call setup_plasma()
    allocate (timstep_arr(neqset))
    state = params

    timstep_arr = 0.0_dp
    call evolvestep(0.0_dp, 1.0e-8_dp)
    after_zero_step = capture()
    print '(A,4ES10.2)', '  dt = 0 step moved the profiles by (max rel., n Vphi Te Ti):', &
        maxval(abs(params - state), dim=2) / maxval(abs(state), dim=2)

    params = state
    call spoil_diagnostics()
    call refresh_state_diagnostics()
    after_refresh = capture()

    call require(all(params == state), 'refresh changed the profiles')
    call require_same(after_refresh%torque_e, after_zero_step%torque_e, 'torque_density_e')
    call require_same(after_refresh%torque_i, after_zero_step%torque_i, 'torque_density_i')
    call require_same(after_refresh%er, after_zero_step%er, 'Ercov')
    call require_same(after_refresh%dni, after_zero_step%dni, 'dni22')
    call require_same(after_refresh%polforce, after_zero_step%polforce, 'polforce')
    call require_same(after_refresh%qlheat_e, after_zero_step%qlheat_e, 'qlheat_e')
    call require_same(after_refresh%qlheat_i, after_zero_step%qlheat_i, 'qlheat_i')
    call require_same(reshape(after_refresh%params_b, [size(params_b)]), &
                      reshape(after_zero_step%params_b, [size(params_b)]), 'params_b')
    call require_same(reshape(after_refresh%ddr_params, [size(ddr_params_nl)]), &
                      reshape(after_zero_step%ddr_params, [size(ddr_params_nl)]), &
                      'ddr_params_nl')
    call require(maxval(abs(after_zero_step%torque_e)) > 0.0_dp, &
                 'fixture has no RMP torque, so the comparison is empty')

    if (num_failed > 0) then
        print *, 'TESTS FAILED: ', num_failed
        error stop 1
    end if
    print *, 'PASS: refresh_state_diagnostics reproduces the dt = 0 step diagnostics'

contains

    function capture() result(d)
        type(diagnostics_t) :: d
        d%torque_e = torque_density_e
        d%torque_i = torque_density_i
        d%er = Ercov
        d%dni = dni22
        d%polforce = polforce
        d%qlheat_e = qlheat_e
        d%qlheat_i = qlheat_i
        d%params_b = params_b
        d%ddr_params = ddr_params_nl
    end function capture

    subroutine spoil_diagnostics()
        torque_density_e = nan
        torque_density_i = nan
        Ercov = nan
        dni22 = nan
        polforce = nan
        qlheat_e = nan
        qlheat_i = nan
        params_b = nan
        ddr_params_nl = nan
    end subroutine spoil_diagnostics

    subroutine require_same(actual, expected, name)
        real(dp), intent(in) :: actual(:), expected(:)
        character(*), intent(in) :: name
        integer :: i
        if (all(actual == expected)) return
        i = maxloc(abs(actual - expected), dim=1, mask=actual == actual)
        if (i == 0) i = 1
        print '(A,A,A,I0,A,ES24.16,A,ES24.16)', '  FAIL: ', name, ' differs, e.g. at ', i, &
            ': ', actual(i), ' vs ', expected(i)
        num_failed = num_failed + 1
    end subroutine require_same

    subroutine setup_plasma()
        real(dp), allocatable :: diffusion(:), bump(:)

        rmin = 2.0_dp
        rmax = 60.0_dp
        npoimin = 300
        iboutype = 2
        gg_factor = 0.0_dp
        gg_width = 58.0_dp
        prop = .false.
        numres = 0
        call gengrid()
        call set_boundary_condition()

        allocate (input_r(npoib))
        input_r = rb
        call allocate_wave_code_data(npoib, input_r)
        n = 2.0e13_dp * (1.1_dp - r / 100.0_dp)
        Te = 1.0e3_dp * (1.2_dp - r / 100.0_dp)
        Ti = 0.8_dp * Te
        q = 1.5_dp + 3.0_dp * (r - 2.0_dp) / 58.0_dp
        Vth = 3.0e5_dp * r / 60.0_dp
        Vz = -8.0e4_dp * 165.0_dp * (1.0_dp - r / 120.0_dp)
        rtor = 165.0_dp
        am = 2.0_dp
        Z_i = 1.0_dp
        btor = -18000.0_dp
        rsepar = 60.0_dp
        mwind = 1
        dery_equisource = 0.0_dp
        source_term = 0.0_dp
        call init_background_profiles()
        call calc_geometric_parameter_profiles()

        viscosity_factor = 1.0_dp

        ! Localised quasilinear transport around a resonance at r = 31 cm.
        bump = 5.0e4_dp * exp(-((rb - 31.0_dp) / 0.3_dp)**2)
        dqle11 = bump
        dqle12 = 1.5_dp * bump
        dqle21 = 1.5_dp * bump
        dqle22 = 3.75_dp * bump
        dqli11 = 0.5_dp * bump
        dqli12 = 0.75_dp * bump
        dqli21 = 0.75_dp * bump
        dqli22 = 1.875_dp * bump

        ihdf5IO = 1
        rb_cut_in = 0.0_dp
        rb_cut_out = maxval(rb) + 1.0_dp
        re_cut_out = rb_cut_out + 1.0_dp
        diffusion = spread(1.0e4_dp, 1, npoib)
        path2inp = 'state_diagnostics_refresh_input.h5'
        call h5_create(trim(path2inp), h5_id)
        call h5_define_group(h5_id, 'da_estimation', group_id_1)
        call h5_close_group(group_id_1)
        call h5_add(h5_id, '/da_estimation/r', rb, [1], [npoib])
        call h5_add(h5_id, '/da_estimation/Da', diffusion, [1], [npoib])
        call h5_close(h5_id)
        call h5_deinit()
    end subroutine setup_plasma

    subroutine require(condition, message)
        logical, intent(in) :: condition
        character(*), intent(in) :: message
        if (.not. condition) then
            print *, 'FAIL: ', message
            num_failed = num_failed + 1
        end if
    end subroutine require

end program test_state_diagnostics_refresh
