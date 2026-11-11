program test_rhs_source_consistency
    !! The implicit balance step advances dy/dt = A(y)*y + q(y), where A comes
    !! from rhs_balance's frozen-coefficient column probing and q is returned by
    !! rhs_balance. Whatever the split between A and q, their sum at the current
    !! state must reproduce the physical right-hand side f(y), with every term
    !! counted once [B1, #199].
    !!
    !! The momentum equation is used because, with the viscosity switched off,
    !! its only right-hand side is the RMP torque density
    !!   f_2 = Z/(n <g_phiphi>) * T_EM / (m_i),
    !!   T_EM = sqrt(g) B^theta/c * e * (Gamma^ql_e - Z Gamma^ql_i)
    !!   [Markl2023 (23), Heyn2014 (64), (67)],
    !! so the reference value is computed here directly from the transport
    !! coefficients and the boundary-interpolated profiles, without the
    !! production assembly routines. The second case adds poloidal rotation,
    !! whose constant contribution -sqrt(g) B^theta/c q V_pol/r to E0r
    !! [Markl2023 (26)] does not respond to the Jacobian probe.
    use QLBalance_kinds, only: dp
    use grid_mod, only: rb, npoib, npoic, rmin, rmax, npoimin, iboutype, gg_factor, &
                        gg_width, set_boundary_condition, mwind, neqset, nbaleqs, &
                        dery_equisource, source_term, ipbeg, ipend, deriv_coef, reint_coef, &
                        sqrt_g_times_B_theta_over_c, gpp_av, rb_cut_in, rb_cut_out, &
                        re_cut_out, dqle11, dqle12, dqle21, dqle22, dqli11, dqli12, dqli21, &
                        dqli22
    use plasma_parameters, only: params, init_background_profiles
    use baseparam_mod, only: am, rtor, Z_i, btor, rsepar, e_charge, p_mass
    use paramscan_mod, only: viscosity_factor
    use resonances_mod, only: prop, numres
    use control_mod, only: ihdf5IO
    use wave_code_data, only: r, n, Te, Ti, q, Vth, Vz
    use matrix_mod, only: nz, irow, icol, amat, rhsvec
    use rhs_balance_m, only: initialize_rhs, rhs_balance
    use h5mod, only: path2inp, h5_id, group_id_1, h5_create, h5_define_group, &
                     h5_close_group, h5_close, h5_deinit, h5_add

    implicit none

    real(dp), parameter :: rel_tol = 1.0e-8_dp
    real(dp), allocatable :: input_r(:)
    integer :: num_failed
    external :: gengrid, allocate_wave_code_data, calc_geometric_parameter_profiles

    num_failed = 0
    call setup_plasma()

    Vth = 0.0_dp
    call check_torque_counted_once('no poloidal rotation')

    Vth = 3.0e5_dp * r / 60.0_dp
    call check_torque_counted_once('with poloidal rotation')

    ! Electron energy row: without anomalous transport and without electron QL heat
    ! flux coefficients, the electron heat flux vanishes and the row reduces to
    !   dTe/dt = (-Te dn/dt + H_e / (3/2)) / n,  H_e = e E0r Gamma^ql_e,
    ! so the QL heating can be checked against the particle row.
    dqle21 = 0.0_dp
    dqle22 = 0.0_dp
    call write_anomalous_diffusion(0.0_dp)
    call check_electron_heating_counted_once('electron QL heating')

    if (num_failed > 0) then
        print *, 'TESTS FAILED: ', num_failed
        error stop 1
    end if
    print *, 'PASS: A(y)*y + q(y) reproduces the RMP torque once'

contains

    subroutine check_torque_counted_once(label)
        character(*), intent(in) :: label
        real(dp), allocatable :: y(:), dy(:), a_times_y(:), rate(:), torque_rate(:)
        integer :: ipoi, k, npoi, i_mom
        real(dp) :: scale, ratio, max_err, max_err_q

        npoi = npoic
        allocate (y(neqset), dy(neqset), a_times_y(neqset), rate(neqset))
        y = reshape(params(:, 1:npoi), [neqset])

        call initialize_rhs(y, dy)
        call rhs_balance(0.0_dp, y, dy)

        a_times_y = 0.0_dp
        do k = 1, nz
            a_times_y(irow(k)) = a_times_y(irow(k)) + amat(k) * y(icol(k))
        end do
        rate = a_times_y + rhsvec

        call reference_torque_rate(torque_rate)
        scale = maxval(abs(torque_rate))
        call require(scale > 0.0_dp, label//': fixture produces no RMP torque')

        max_err = 0.0_dp
        max_err_q = 0.0_dp
        ratio = 0.0_dp
        do ipoi = 1, npoi
            i_mom = nbaleqs * (ipoi - 1) + 2
            max_err = max(max_err, abs(rate(i_mom) - torque_rate(ipoi)))
            max_err_q = max(max_err_q, abs(rhsvec(i_mom) - torque_rate(ipoi)))
            ratio = ratio + rate(i_mom) * torque_rate(ipoi)
        end do
        ratio = ratio / sum(torque_rate**2)
        print '(2A)', '  case: ', label
        ! Diagnostic only: the split between A and q is a design choice.
        print '(A,ES12.4)', '    q(y) momentum rows vs reference torque, max rel. error: ', &
            max_err_q / scale
        print '(A,ES12.4)', '    A*y+q momentum rows vs reference torque, max rel. error: ', &
            max_err / scale
        print '(A,F10.6)', '    least-squares ratio (A*y+q)/torque: ', ratio
        call require(max_err <= rel_tol * scale, &
                     label//': A(y)*y + q(y) does not reproduce the RMP torque once')
    end subroutine check_torque_counted_once

    subroutine check_electron_heating_counted_once(label)
        character(*), intent(in) :: label
        real(dp), allocatable :: y(:), dy(:), a_times_y(:), rate(:), heat_rate(:), expected(:)
        integer :: ipoi, k, npoi, i_par, i_te
        real(dp) :: scale, ratio_num, ratio_den, max_err

        npoi = npoic
        allocate (y(neqset), dy(neqset), a_times_y(neqset), rate(neqset), expected(npoi))
        y = reshape(params(:, 1:npoi), [neqset])

        call initialize_rhs(y, dy)
        call rhs_balance(0.0_dp, y, dy)

        a_times_y = 0.0_dp
        do k = 1, nz
            a_times_y(irow(k)) = a_times_y(irow(k)) + amat(k) * y(icol(k))
        end do
        rate = a_times_y + rhsvec

        call reference_electron_heat_rate(heat_rate)
        max_err = 0.0_dp
        ratio_num = 0.0_dp
        ratio_den = 0.0_dp
        do ipoi = 1, npoi
            i_par = nbaleqs * (ipoi - 1) + 1
            i_te = nbaleqs * (ipoi - 1) + 3
            expected(ipoi) = (-params(3, ipoi) * rate(i_par) + heat_rate(ipoi) / 1.5_dp) / &
                             params(1, ipoi)
            ! heating part of the row as applied by the step, and as expected
            ratio_num = ratio_num + (rate(i_te) * params(1, ipoi) * 1.5_dp + &
                                     1.5_dp * params(3, ipoi) * rate(i_par)) * heat_rate(ipoi)
            ratio_den = ratio_den + heat_rate(ipoi)**2
        end do
        scale = maxval(abs(expected))
        call require(maxval(abs(heat_rate)) > 0.0_dp, label//': fixture produces no QL heating')
        do ipoi = 1, npoi
            max_err = max(max_err, abs(rate(nbaleqs * (ipoi - 1) + 3) - expected(ipoi)))
        end do
        print '(2A)', '  case: ', label
        print '(A,ES12.4)', '    A*y+q Te rows vs expected, max rel. error: ', max_err / scale
        print '(A,F10.6)', '    least-squares ratio applied/expected QL heating: ', &
            ratio_num / ratio_den
        call require(max_err <= rel_tol * scale, &
                     label//': A(y)*y + q(y) does not apply the QL heating once')
    end subroutine check_electron_heating_counted_once

    subroutine setup_plasma()
        real(dp), allocatable :: bump(:)

        rmin = 2.0_dp
        rmax = 60.0_dp
        npoimin = 120
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
        bump = 1.0e3_dp * exp(-((rb - 31.0_dp) / 2.0_dp)**2)
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
        path2inp = 'rhs_source_consistency_input.h5'
        call write_anomalous_diffusion(1.0e4_dp)
    end subroutine setup_plasma

    subroutine reference_torque_rate(rate_out)
        !! Torque contribution to d(Vphi/R)/dt on cell centres, from the textbook
        !! flux-force relations Gamma = -n (D11 A1 + D12 A2) [Heyn2014 (35)].
        real(dp), allocatable, intent(out) :: rate_out(:)
        real(dp), allocatable :: torque_b(:)
        real(dp) :: nb, vphib, teb, tib, dn, dte, dti, er, a1e, a2e, a1i, a2i
        real(dp) :: gamma_e, gamma_i
        integer :: ib, ic, lo, hi

        allocate (torque_b(npoib), rate_out(npoic))
        do ib = 1, npoib
            lo = ipbeg(ib)
            hi = ipend(ib)
            nb = sum(params(1, lo:hi) * reint_coef(:, ib))
            vphib = sum(params(2, lo:hi) * reint_coef(:, ib))
            teb = sum(params(3, lo:hi) * reint_coef(:, ib))
            tib = sum(params(4, lo:hi) * reint_coef(:, ib))
            dn = sum(params(1, lo:hi) * deriv_coef(:, ib))
            dte = sum(params(3, lo:hi) * deriv_coef(:, ib))
            dti = sum(params(4, lo:hi) * deriv_coef(:, ib))
            er = (tib / nb * dn + dti) / (Z_i * e_charge) + &
                 sqrt_g_times_B_theta_over_c(ib) * (vphib - q(ib) * Vth(ib) / rb(ib))
            a1e = dn / nb - 1.5_dp * dte / teb + e_charge * er / teb
            a2e = dte / teb
            a1i = dn / nb - 1.5_dp * dti / tib - Z_i * e_charge * er / tib
            a2i = dti / tib
            gamma_e = -(dqle11(ib) * a1e + dqle12(ib) * a2e) * nb
            gamma_i = -(dqli11(ib) * a1i + dqli12(ib) * a2i) * nb / Z_i
            torque_b(ib) = sqrt_g_times_B_theta_over_c(ib) * e_charge * &
                           (gamma_e - Z_i * gamma_i)
        end do
        do ic = 1, npoic
            rate_out(ic) = 0.5_dp * (torque_b(ic) + torque_b(ic + 1)) / (am * p_mass) * &
                           Z_i / params(1, ic) * 2.0_dp / (gpp_av(ic) + gpp_av(ic + 1))
        end do
    end subroutine reference_torque_rate

    subroutine reference_electron_heat_rate(rate_out)
        !! Electron QL heating density H_e = e E0r Gamma^ql_e on cell centres, with
        !! the sign convention of the production code (+e E0r Gamma_e; the sign is
        !! an open physics question and not what this test checks).
        real(dp), allocatable, intent(out) :: rate_out(:)
        real(dp), allocatable :: heat_b(:)
        real(dp) :: nb, vphib, teb, tib, dn, dte, dti, er, a1e, a2e, gamma_e
        integer :: ib, ic, lo, hi

        allocate (heat_b(npoib), rate_out(npoic))
        do ib = 1, npoib
            lo = ipbeg(ib)
            hi = ipend(ib)
            nb = sum(params(1, lo:hi) * reint_coef(:, ib))
            vphib = sum(params(2, lo:hi) * reint_coef(:, ib))
            teb = sum(params(3, lo:hi) * reint_coef(:, ib))
            tib = sum(params(4, lo:hi) * reint_coef(:, ib))
            dn = sum(params(1, lo:hi) * deriv_coef(:, ib))
            dte = sum(params(3, lo:hi) * deriv_coef(:, ib))
            dti = sum(params(4, lo:hi) * deriv_coef(:, ib))
            er = (tib / nb * dn + dti) / (Z_i * e_charge) + &
                 sqrt_g_times_B_theta_over_c(ib) * (vphib - q(ib) * Vth(ib) / rb(ib))
            a1e = dn / nb - 1.5_dp * dte / teb + e_charge * er / teb
            a2e = dte / teb
            gamma_e = -(dqle11(ib) * a1e + dqle12(ib) * a2e) * nb
            heat_b(ib) = e_charge * er * gamma_e
        end do
        do ic = 1, npoic
            rate_out(ic) = 0.5_dp * (heat_b(ic) + heat_b(ic + 1))
        end do
    end subroutine reference_electron_heat_rate

    subroutine write_anomalous_diffusion(value)
        real(dp), intent(in) :: value
        real(dp), allocatable :: diffusion(:)
        diffusion = spread(value, 1, npoib)
        call h5_create(trim(path2inp), h5_id)
        call h5_define_group(h5_id, 'da_estimation', group_id_1)
        call h5_close_group(group_id_1)
        call h5_add(h5_id, '/da_estimation/r', rb, [1], [npoib])
        call h5_add(h5_id, '/da_estimation/Da', diffusion, [1], [npoib])
        call h5_close(h5_id)
        call h5_deinit()
    end subroutine write_anomalous_diffusion

    subroutine require(condition, message)
        logical, intent(in) :: condition
        character(*), intent(in) :: message
        if (.not. condition) then
            print *, 'FAIL: ', message
            num_failed = num_failed + 1
        end if
    end subroutine require

end program test_rhs_source_consistency
