program test_evolvestep_stationary
    !! Accuracy of the sparse solve in the implicit step (I - dt A) y' = y + dt q.
    !!
    !! The fixture is made stationary by the equilibrium source
    !! (det_balance_eqs_source_terms sets dery_equisource = -f(y)), so the exact
    !! solution is y' = y. Quasilinear transport is strong and localised on a fine
    !! grid, as near a resonant surface, and the unknowns span ~20 orders of
    !! magnitude in CGS units (n ~ 1e13 cm^-3, T ~ 1e-9 erg).
    !!
    !! Asserted: the row-relative residual of the computed y' in the system built
    !! at y is at round-off level. Printed only: the departure of y' from y, which
    !! is the residual amplified by the condition number of (I - dt A) and grows
    !! with dt for this stiff fixture.
    use QLBalance_kinds, only: dp
    use grid_mod, only: rb, npoib, npoic, rmin, rmax, npoimin, iboutype, gg_factor, &
                        gg_width, set_boundary_condition, mwind, neqset, nbaleqs, &
                        dery_equisource, source_term, rb_cut_in, rb_cut_out, re_cut_out, &
                        dqle11, dqle12, dqle21, dqle22, dqli11, dqli12, dqli21, dqli22
    use plasma_parameters, only: params, init_background_profiles
    use baseparam_mod, only: am, rtor, Z_i, btor, rsepar, ev
    use paramscan_mod, only: viscosity_factor
    use resonances_mod, only: prop, numres
    use control_mod, only: ihdf5IO
    use wave_code_data, only: r, n, Te, Ti, q, Vth, Vz
    use recstep_mod, only: timstep_arr
    use matrix_mod, only: nz, irow, icol, amat, rhsvec
    use rhs_balance_m, only: initialize_rhs, rhs_balance
    use h5mod, only: path2inp, h5_id, group_id_1, h5_create, h5_define_group, &
                     h5_close_group, h5_close, h5_deinit, h5_add

    implicit none

    real(dp), parameter :: residual_tol = 1.0e-14_dp
    real(dp), parameter :: dt_values(3) = [1.0e-2_dp, 1.0e-1_dp, 5.0e-1_dp]
    real(dp), allocatable :: input_r(:), state(:, :), y_new(:), y_old(:), dyv(:)
    real(dp), allocatable :: lhs(:), mag(:)
    real(dp) :: max_rel(4), res_rel
    integer :: num_failed, it, ieq, npoi, k
    external :: gengrid, allocate_wave_code_data, calc_geometric_parameter_profiles
    external :: det_balance_eqs_source_terms, evolvestep

    num_failed = 0
    call setup_plasma()
    viscosity_factor = 1.0_dp
    Vth = 3.0e5_dp * r / 60.0_dp
    npoi = npoic
    allocate (timstep_arr(neqset), dyv(neqset))

    call det_balance_eqs_source_terms()
    state = params(:, 1:npoi)
    y_old = reshape(state, [neqset])

    do it = 1, size(dt_values)
        params(:, 1:npoi) = state
        timstep_arr = dt_values(it)
        call evolvestep(dt_values(it), 1.0e-8_dp)
        y_new = reshape(params(:, 1:npoi), [neqset])
        do ieq = 1, nbaleqs
            max_rel(ieq) = maxval(abs(params(ieq, 1:npoi) - state(ieq, :))) / &
                           maxval(abs(state(ieq, :)))
        end do

        params(:, 1:npoi) = state
        call initialize_rhs(y_old, dyv)
        call rhs_balance(0.0_dp, y_old, dyv)
        lhs = y_new
        mag = abs(y_new) + abs(y_old) + dt_values(it) * abs(rhsvec)
        do k = 1, nz
            lhs(irow(k)) = lhs(irow(k)) - dt_values(it) * amat(k) * y_new(icol(k))
            mag(irow(k)) = mag(irow(k)) + dt_values(it) * abs(amat(k) * y_new(icol(k)))
        end do
        res_rel = maxval(abs(lhs - y_old - dt_values(it) * rhsvec) / mag)

        print '(A,ES9.2,A,ES10.3,A,4ES10.2)', '  dt =', dt_values(it), &
            '  residual:', res_rel, '  departure (n, Vphi, Te, Ti):', max_rel
        call require(res_rel <= residual_tol, 'sparse solve is not accurate to round-off')
    end do

    if (num_failed > 0) then
        print *, 'TESTS FAILED: ', num_failed
        error stop 1
    end if
    print *, 'PASS: the implicit step is solved to round-off accuracy'

contains

    subroutine setup_plasma()
        real(dp), allocatable :: diffusion(:), bump(:)

        rmin = 2.0_dp
        rmax = 60.0_dp
        npoimin = 1500
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
        Vth = 0.0_dp
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

        ! Switch off momentum diffusion so that the torque is the only momentum term.
        viscosity_factor = 0.0_dp

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
        path2inp = 'evolvestep_stationary_input.h5'
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

end program test_evolvestep_stationary
