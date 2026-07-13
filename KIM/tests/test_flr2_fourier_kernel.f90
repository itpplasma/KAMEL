program test_flr2_fourier_kernel
    ! Plumbing-validation harness for the fused FLR2 Fourier kernel.
    !
    ! Proves two things before any kernel math is written:
    !   (a) the new flr2_fourier_kernel_m module compiles and links into a
    !       standalone CTest program, and
    !   (b) an electrostatic FokkerPlanck run populates a global `plasma`
    !       whose per-species susceptibility functions I00(j,0) and
    !       backgrounds rho_L(j), lambda_D(j) are available on the
    !       guiding-centre grid rg_grid%xb.
    !
    ! Setup mirrors the in-memory-profile path of test_kim_diagnostics.f90:
    !   write_test_namelist -> nml_config_path -> profiles_in_memory ->
    !   kim_init -> set_profiles_from_arrays -> electrostatic %init().
    ! The electrostatic init calls set_plasma_quantities(plasma), which runs
    ! calculate_plasma_backs (rho_L, lambda_D) and
    ! calculate_thermodynamic_forces_and_susc (I00) on rg_grid%xb. No %run()
    ! is needed to read the boundary-grid susceptibility functions.

    use KIM_kinds_m, only: dp
    use flr2_fourier_kernel_m, only: hatG_rho_phi, hatG_rho_B
    use flr2_fourier_kernel_m, only: hatG_rho_phi_diag_sp, hatG_rho_B_diag_sp
    use flr2_fourier_kernel_m, only: fp_ks_model, FP_KPERP_FULL, FP_KR_ONLY_APPROX

    implicit none

    call test_independent_rho_phi_oracle()
    call test_populated_plasma_and_kernel_stub()
    call test_diagonal_matches_inline()
    call test_collapse_rho_phi()
    call test_phase_rho_phi()
    call test_collapse_rho_phi_ks2()
    call test_zero_wavenumber_limit()
    call test_collapse_rho_B()
    call test_collapse_rho_B_ks2()
    call test_zero_wavenumber_limit_rho_B()
    call test_scaled_bessel_large_arg()
    call test_large_bcross_kernel_finite()

    print *, 'All tests PASSED'
    stop 0

contains

    subroutine test_independent_rho_phi_oracle()
        use, intrinsic :: ieee_arithmetic, only: ieee_value, ieee_quiet_nan
        use flr2_fourier_kernel_m, only: FP_KPERP_FULL, FP_KR_ONLY_APPROX
        use flr2_fourier_kernel_m, only: flr_arg_pair, fp_rho_phi_harmonic, fp_rho_B_harmonic
        use flr2_fourier_kernel_m, only: fp_rho_phi_debye

        integer :: unit, ios, case_id, charge_sign, mphi, include_debye
        integer :: status, nrows
        real(dp) :: lambda_D, vT, omega_c, nu, ks, kr, krp, rg, A1, A2
        real(dp) :: kpar, omega_E, Vpar, omega, x1, x2, rho_L
        real(dp) :: bplus_ref, bcross_ref, phase_re, phase_im
        real(dp) :: I00_re, I00_im, I20_re, I20_im, I01_re, I01_im, I21_re, I21_im
        real(dp) :: F0_ref, F2_ref
        real(dp) :: dynamic_re, dynamic_im, debye_re, debye_im
        real(dp) :: total_re, total_im, bplus, bcross, scaled_error
        real(dp) :: rho_B_re, rho_B_im
        real(dp) :: nan_value, max_scaled_error
        complex(dp) :: I00, I20, I01, I21, actual_dynamic, actual_debye, actual_rho_B, expected
        complex(dp) :: representative(-2:2), representative_total
        character(len=2048) :: line
        logical :: seen_charge(-1:1), seen_large, seen_zero_flr, seen_debye

        nrows = 0
        max_scaled_error = 0.0_dp
        seen_charge = .false.
        seen_large = .false.
        seen_zero_flr = .false.
        seen_debye = .false.
        representative = (0.0_dp, 0.0_dp)
        open(newunit=unit, file='rho_phi_kernels.dat', status='old', action='read', iostat=ios)
        if (ios /= 0) error stop 'cannot open rho_phi_kernels.dat'

        do
            read(unit, '(A)', iostat=ios) line
            if (ios < 0) exit
            if (ios > 0) error stop 'cannot read rho-Phi oracle line'
            if (len_trim(line) == 0 .or. line(1:1) == '#') cycle
            read(line, *, iostat=ios) case_id, charge_sign, lambda_D, vT, omega_c, nu, &
                ks, kr, krp, rg, mphi, A1, A2, kpar, omega_E, Vpar, omega, &
                include_debye, x1, x2, rho_L, bplus_ref, bcross_ref, phase_re, &
                phase_im, I00_re, I00_im, I20_re, I20_im, I01_re, I01_im, &
                I21_re, I21_im, F0_ref, F2_ref, &
                dynamic_re, dynamic_im, debye_re, debye_im, total_re, total_im, &
                rho_B_re, rho_B_im
            if (ios /= 0) error stop 'malformed rho-Phi oracle row'

            call flr_arg_pair(ks, kr, krp, rho_L, FP_KPERP_FULL, bplus, bcross, status)
            if (status /= 0) error stop 'full finite-radius FLR arguments rejected oracle row'
            call assert_oracle_real('bplus', bplus, bplus_ref, case_id, max_scaled_error)
            call assert_oracle_real('bcross', bcross, bcross_ref, case_id, max_scaled_error)

            I00 = cmplx(I00_re, I00_im, kind=dp)
            I20 = cmplx(I20_re, I20_im, kind=dp)
            call fp_rho_phi_harmonic(lambda_D, vT, omega_c, nu, ks, kr, krp, rg, &
                mphi, A1, A2, I00, I20, FP_KPERP_FULL, actual_dynamic, status)
            if (status /= 0) error stop 'rho-Phi harmonic rejected oracle row'
            expected = cmplx(dynamic_re, dynamic_im, kind=dp)
            call assert_oracle_complex('dynamic', actual_dynamic, expected, case_id, max_scaled_error)

            call fp_rho_phi_debye(lambda_D, kr, krp, rg, include_debye == 1 .and. mphi == 0, actual_debye)
            expected = cmplx(debye_re, debye_im, kind=dp)
            call assert_oracle_complex('Debye', actual_debye, expected, case_id, max_scaled_error)
            expected = cmplx(total_re, total_im, kind=dp)
            call assert_oracle_complex('total', actual_dynamic + actual_debye, expected, &
                case_id, max_scaled_error)

            I01 = cmplx(I01_re, I01_im, kind=dp)
            I21 = cmplx(I21_re, I21_im, kind=dp)
            call fp_rho_B_harmonic(lambda_D, vT, omega_c, nu, ks, kr, krp, rg, &
                mphi, A1, A2, I01, I21, FP_KPERP_FULL, actual_rho_B, status)
            if (status /= 0) error stop 'rho-B harmonic rejected oracle row'
            expected = cmplx(rho_B_re, rho_B_im, kind=dp)
            call assert_oracle_complex('rho-B', actual_rho_B, expected, case_id, max_scaled_error)

            nrows = nrows + 1
            if (charge_sign >= -1 .and. charge_sign <= 1) seen_charge(charge_sign) = .true.
            if (bcross_ref >= 500.0_dp) seen_large = .true.
            if (rho_L <= 1.0e-8_dp) seen_zero_flr = .true.
            if (include_debye == 1) seen_debye = .true.
            if (case_id >= 9 .and. case_id <= 13) representative(mphi) = actual_dynamic
        end do
        close(unit)

        if (nrows /= 13) error stop 'rho-Phi oracle row coverage changed'
        if (.not. seen_charge(-1) .or. .not. seen_charge(1)) error stop 'oracle lacks both charge signs'
        if (.not. seen_large) error stop 'oracle lacks large-Bessel case'
        if (.not. seen_zero_flr) error stop 'oracle lacks zero-FLR case'
        if (.not. seen_debye) error stop 'oracle lacks Debye case'

        representative_total = sum(representative)
        scaled_error = abs(representative_total - representative(0))/abs(representative_total)
        if (scaled_error >= 5.0e-3_dp) error stop 'mphi=0 truncation exceeds verified bound'
        scaled_error = (abs(representative(-2)) + abs(representative(2)))/abs(representative_total)
        if (scaled_error >= 5.0e-5_dp) error stop '|mphi|=2 shell exceeds verified bound'

        nan_value = ieee_value(0.0_dp, ieee_quiet_nan)
        call flr_arg_pair(nan_value, 1.0_dp, 2.0_dp, 0.5_dp, FP_KPERP_FULL, &
            bplus, bcross, status)
        if (status == 0) error stop 'full finite-radius model accepted undefined ks'
        call flr_arg_pair(nan_value, 1.0_dp, 2.0_dp, 0.5_dp, FP_KR_ONLY_APPROX, &
            bplus, bcross, status)
        if (status /= 0) error stop 'named radial-only approximation used undefined ks'

        print *, 'PASS: ', nrows, ' independent finite-radius rho-Phi oracle rows'
        print *, 'Maximum scaled kernel-oracle error: ', max_scaled_error
    end subroutine test_independent_rho_phi_oracle

    subroutine assert_oracle_real(label, actual, expected, case_id, max_error)
        character(len=*), intent(in) :: label
        real(dp), intent(in) :: actual, expected
        integer, intent(in) :: case_id
        real(dp), intent(inout) :: max_error
        real(dp) :: error

        error = abs(actual - expected)/max(1.0_dp, abs(expected))
        max_error = max(max_error, error)
        if (error > 1.0e-12_dp) then
            print *, 'FAIL: ', trim(label), ' oracle case ', case_id
            print *, ' actual/expected/error = ', actual, expected, error
            error stop
        end if
    end subroutine assert_oracle_real

    subroutine assert_oracle_complex(label, actual, expected, case_id, max_error)
        character(len=*), intent(in) :: label
        complex(dp), intent(in) :: actual, expected
        integer, intent(in) :: case_id
        real(dp), intent(inout) :: max_error
        real(dp) :: error

        error = abs(actual - expected)/max(1.0_dp, abs(expected))
        max_error = max(max_error, error)
        if (error > 1.0e-12_dp) then
            print *, 'FAIL: ', trim(label), ' oracle case ', case_id
            print *, ' actual   = ', actual
            print *, ' expected = ', expected
            print *, ' scaled error = ', error
            error stop
        end if
    end subroutine assert_oracle_complex

    subroutine test_diagonal_matches_inline()
        ! Characterization test: the shared per-species diagonal integrand
        ! functions must reproduce, bit-for-bit within tolerance, the original
        ! inline expressions in calc_hatK_Phi_in_Fourier. The INLINE reference
        ! below is the permanent anchor (verbatim copy of the source-of-truth
        ! per-species rho-Phi Debye+FP term and signed rho-B FP term). It must
        ! not be refactored to call the functions it validates.
        use config_m, only: profiles_in_memory, nml_config_path
        use config_m, only: artificial_debye_case, turn_off_ions, turn_off_electrons
        use constants_m, only: com_unit, sol
        use fortnum_special, only: bessel_in
        use species_m, only: plasma, set_profiles_from_arrays, plasma_t
        use grid_m, only: rg_grid
        use kim_base_m, only: kim_t
        use kim_mod_m, only: from_kim_factory_get_kim

        integer, parameter :: npts = 101

        real(dp) :: r_prof(npts), n_prof(npts), Te_prof(npts)
        real(dp) :: Ti_prof(npts), q_prof(npts), Er_prof(npts)
        class(kim_t), allocatable :: kim_instance
        real(dp) :: kr_arr(3), kr, b, ks
        integer :: i, j, sp, jj_arr(3), k
        complex(dp) :: phi_inline, phi_func, B_inline, B_func
        real(dp) :: tol_phi, tol_B

        call make_test_profiles(npts, r_prof, n_prof, Te_prof, Ti_prof, &
                                q_prof, Er_prof)

        call write_test_namelist('./KIM_config_fourier_diag_test.nml')
        nml_config_path = './KIM_config_fourier_diag_test.nml'

        profiles_in_memory = .true.
        call kim_init()
        call set_profiles_from_arrays(r_prof, n_prof, Te_prof, Ti_prof, &
                                      q_prof, Er_prof, npts)

        call from_kim_factory_get_kim('electrostatic', kim_instance)
        call kim_instance%init()

        ! This characterization block intentionally exercises the named
        ! historical radial-only approximation copied below.
        fp_ks_model = FP_KR_ONLY_APPROX

        kr_arr = [0.1_dp, 1.0_dp, 5.0_dp]
        ! A few interior guiding-centre boundary-grid indices.
        jj_arr = [rg_grid%npts_b / 4, rg_grid%npts_b / 2, (3 * rg_grid%npts_b) / 4]

        do i = 1, size(kr_arr)
            kr = kr_arr(i)
            do k = 1, size(jj_arr)
                j = jj_arr(k)

                ! (i) INLINE reference: verbatim per-species expressions,
                ! summed over the non-turned-off species.
                phi_inline = (0.0_dp, 0.0_dp)
                B_inline = (0.0_dp, 0.0_dp)
                do sp = 0, plasma%n_species - 1
                    if (turn_off_ions .and. sp >= 1) cycle
                    if (turn_off_electrons .and. sp == 0) cycle

                    ks = plasma%ks(j)
                    b = (kr**2.0d0) * plasma%spec(sp)%rho_L(j)**2.0d0

                    if (artificial_debye_case <= 1) then
                        phi_inline = phi_inline - 1.0d0 / plasma%spec(sp)%lambda_D(j)**2.0d0
                    end if

                    if (artificial_debye_case == 0 .or. artificial_debye_case == 2) then
                        phi_inline = phi_inline + 1.0d0 / plasma%spec(sp)%lambda_D(j)**2.0d0 &
                            * com_unit * plasma%spec(sp)%vT(j)**2.0d0 * plasma%ks(j) &
                            / (plasma%spec(sp)%omega_c(j) * plasma%spec(sp)%nu(j)) * exp(-b) * &
                            (&
                                plasma%spec(sp)%I00(j, 0) * (&
                                    bessel_in(0, b) * (plasma%spec(sp)%A1(j) + plasma%spec(sp)%A2(j) * (1-b)) &
                                    + plasma%spec(sp)%A2(j) * b * bessel_in(-1, b) &
                                )&
                                + 0.5d0 * plasma%spec(sp)%I20(j, 0) * plasma%spec(sp)%A2(j) * bessel_in(0, b) &
                            )
                        B_inline = B_inline - 1.0d0 / plasma%spec(sp)%lambda_D(j)**2.0d0 * plasma%spec(sp)%vT(j)**3.0d0 &
                            / (plasma%spec(sp)%omega_c(j) * plasma%spec(sp)%nu(j) * sol) * exp(-b) * &
                            (&
                                plasma%spec(sp)%I01(j, 0) * (&
                                    bessel_in(0, b) * (plasma%spec(sp)%A1(j) + plasma%spec(sp)%A2(j) * (1-b)) &
                                    + plasma%spec(sp)%A2(j) * b * bessel_in(-1, b) &
                                )&
                                + 0.5d0 * plasma%spec(sp)%I21(j, 0) * plasma%spec(sp)%A2(j) * bessel_in(0, b) &
                            )
                    end if
                end do

                ! (ii) Shared functions: sum over the same non-turned-off species.
                phi_func = (0.0_dp, 0.0_dp)
                B_func = (0.0_dp, 0.0_dp)
                do sp = 0, plasma%n_species - 1
                    if (turn_off_ions .and. sp >= 1) cycle
                    if (turn_off_electrons .and. sp == 0) cycle
                    phi_func = phi_func + hatG_rho_phi_diag_sp(plasma, sp, kr, j)
                    B_func = B_func + hatG_rho_B_diag_sp(plasma, sp, kr, j)
                end do

                tol_phi = 1.0e-12_dp * (1.0_dp + abs(phi_inline))
                tol_B = 1.0e-12_dp * (1.0_dp + abs(B_inline))

                if (abs(phi_inline - phi_func) >= tol_phi) then
                    print *, 'FAIL: rho-Phi diagonal mismatch at kr=', kr, ' j=', j
                    print *, '  inline   = ', phi_inline
                    print *, '  fromfunc = ', phi_func
                    print *, '  |diff|   = ', abs(phi_inline - phi_func), ' tol = ', tol_phi
                    error stop
                end if
                if (abs(B_inline - B_func) >= tol_B) then
                    print *, 'FAIL: rho-B diagonal mismatch at kr=', kr, ' j=', j
                    print *, '  inline   = ', B_inline
                    print *, '  fromfunc = ', B_func
                    print *, '  |diff|   = ', abs(B_inline - B_func), ' tol = ', tol_B
                    error stop
                end if
            end do
        end do

        print *, 'PASS: shared diagonal integrand matches inline reference'
    end subroutine

    subroutine test_collapse_rho_phi()
        ! Collapse-by-construction: the fused two-wavenumber kernel evaluated on
        ! the diagonal k_r = k'_r must equal the per-species diagonal integrand
        ! summed over the non-turned-off species and divided by 4*pi (the
        ! source-of-truth normalization). Exact to machine precision because the
        ! diagonal and the fused kernel share the same per-species core.
        use config_m, only: profiles_in_memory, nml_config_path
        use config_m, only: turn_off_ions, turn_off_electrons
        use flr2_fourier_kernel_m, only: hatG_rho_phi, hatG_rho_phi_diag_sp
        use constants_m, only: pi
        use species_m, only: plasma, set_profiles_from_arrays
        use grid_m, only: rg_grid
        use kim_base_m, only: kim_t
        use kim_mod_m, only: from_kim_factory_get_kim

        integer, parameter :: npts = 101

        real(dp) :: r_prof(npts), n_prof(npts), Te_prof(npts)
        real(dp) :: Ti_prof(npts), q_prof(npts), Er_prof(npts)
        class(kim_t), allocatable :: kim_instance
        real(dp) :: kr_arr(3), kr, tol
        integer :: i, j, sp, jj_arr(3), k
        complex(dp) :: dref, gfused

        call make_test_profiles(npts, r_prof, n_prof, Te_prof, Ti_prof, &
                                q_prof, Er_prof)

        call write_test_namelist('./KIM_config_fourier_collapse_test.nml')
        nml_config_path = './KIM_config_fourier_collapse_test.nml'

        profiles_in_memory = .true.
        call kim_init()
        call set_profiles_from_arrays(r_prof, n_prof, Te_prof, Ti_prof, &
                                      q_prof, Er_prof, npts)

        call from_kim_factory_get_kim('electrostatic', kim_instance)
        call kim_instance%init()

        fp_ks_model = FP_KR_ONLY_APPROX

        kr_arr = [0.1_dp, 1.0_dp, 5.0_dp]
        jj_arr = [rg_grid%npts_b / 4, rg_grid%npts_b / 2, (3 * rg_grid%npts_b) / 4]

        do i = 1, size(kr_arr)
            kr = kr_arr(i)
            do k = 1, size(jj_arr)
                j = jj_arr(k)

                dref = (0.0_dp, 0.0_dp)
                do sp = 0, plasma%n_species - 1
                    if (turn_off_ions .and. sp >= 1) cycle
                    if (turn_off_electrons .and. sp == 0) cycle
                    dref = dref + hatG_rho_phi_diag_sp(plasma, sp, kr, j)
                end do
                dref = dref / (4.0d0 * pi)

                gfused = hatG_rho_phi(plasma, kr, kr, j)

                tol = 5.0e-12_dp * (1.0_dp + abs(dref))
                if (abs(gfused - dref) >= tol) then
                    print *, 'FAIL: rho-Phi collapse mismatch at kr=', kr, ' j=', j
                    print *, '  diag/4pi = ', dref
                    print *, '  fused    = ', gfused
                    print *, '  |diff|   = ', abs(gfused - dref), ' tol = ', tol
                    error stop
                end if
            end do
        end do

        print *, 'PASS: fused rho-Phi kernel collapses to diagonal at kr=krp'
    end subroutine

    subroutine test_phase_rho_phi()
        ! Off-diagonal phase structure: for k_r /= k'_r the fused kernel carries
        ! the Fourier phase exp(i*(kr-krp)*rg). Because b_+ and b_x are symmetric
        ! under (kr,krp) -> (krp,kr), the per-species cores are equal, so the only
        ! difference between hatG(kr,krp) and hatG(krp,kr) is that phase factor.
        ! Stripping each of its own phase must leave identical (phase-free,
        ! still complex) core values. Tolerance is looser than the collapse
        ! test: this compares two independent complex-exp/Bessel evaluations,
        ! not a construction identity.
        use config_m, only: profiles_in_memory, nml_config_path
        use flr2_fourier_kernel_m, only: hatG_rho_phi
        use constants_m, only: com_unit
        use species_m, only: plasma, set_profiles_from_arrays
        use grid_m, only: rg_grid
        use kim_base_m, only: kim_t
        use kim_mod_m, only: from_kim_factory_get_kim

        integer, parameter :: npts = 101

        real(dp) :: r_prof(npts), n_prof(npts), Te_prof(npts)
        real(dp) :: Ti_prof(npts), q_prof(npts), Er_prof(npts)
        class(kim_t), allocatable :: kim_instance
        real(dp) :: kr_arr(3), krp_arr(3), kr, krp, rg, tol
        integer :: i, j, jj_arr(3), k
        complex(dp) :: ci, g_fwd, g_rev, stripped_fwd, stripped_rev

        call make_test_profiles(npts, r_prof, n_prof, Te_prof, Ti_prof, &
                                q_prof, Er_prof)

        call write_test_namelist('./KIM_config_fourier_phase_test.nml')
        nml_config_path = './KIM_config_fourier_phase_test.nml'

        profiles_in_memory = .true.
        call kim_init()
        call set_profiles_from_arrays(r_prof, n_prof, Te_prof, Ti_prof, &
                                      q_prof, Er_prof, npts)

        call from_kim_factory_get_kim('electrostatic', kim_instance)
        call kim_instance%init()

        fp_ks_model = FP_KR_ONLY_APPROX
        ci = com_unit

        kr_arr = [0.1_dp, 1.0_dp, 2.0_dp]
        krp_arr = [1.0_dp, 5.0_dp, 0.5_dp]
        jj_arr = [rg_grid%npts_b / 4, rg_grid%npts_b / 2, (3 * rg_grid%npts_b) / 4]

        do i = 1, size(kr_arr)
            kr = kr_arr(i)
            krp = krp_arr(i)
            do k = 1, size(jj_arr)
                j = jj_arr(k)
                rg = rg_grid%xb(j)

                g_fwd = hatG_rho_phi(plasma, kr, krp, j)
                g_rev = hatG_rho_phi(plasma, krp, kr, j)

                stripped_fwd = g_fwd * exp(-ci * (kr - krp) * rg)
                stripped_rev = g_rev * exp(-ci * (krp - kr) * rg)

                tol = 1.0e-10_dp * (1.0_dp + abs(g_fwd))
                if (abs(stripped_fwd - stripped_rev) >= tol) then
                    print *, 'FAIL: rho-Phi phase mismatch at kr=', kr, ' krp=', krp, ' j=', j
                    print *, '  stripped_fwd = ', stripped_fwd
                    print *, '  stripped_rev = ', stripped_rev
                    print *, '  |diff|       = ', abs(stripped_fwd - stripped_rev), ' tol = ', tol
                    error stop
                end if
            end do
        end do

        print *, 'PASS: fused rho-Phi kernel carries the Fourier phase only'
    end subroutine

    subroutine test_collapse_rho_phi_ks2()
        ! Structural collapse for the full finite-radius model. This model
        ! keeps k_s in the FLR arguments, so it does NOT
        ! collapse to the k_s-dropped diagonal source-of-truth. Instead, on the
        ! diagonal k_r = k'_r the two FLR arguments both reduce to the full
        ! perpendicular argument bf = (k_s^2 + k_r^2) * rho_L^2, so the fused
        ! kernel must equal the per-species core evaluated at (bf, bf), summed
        ! over the non-turned-off species and divided by 4*pi. This validates the
        ! sqrt/2 algebra of the .true. branch without needing a k_s = 0 setup.
        use config_m, only: profiles_in_memory, nml_config_path
        use config_m, only: turn_off_ions, turn_off_electrons
        use flr2_fourier_kernel_m, only: hatG_rho_phi, core_rho_phi_sp
        use constants_m, only: pi
        use species_m, only: plasma, set_profiles_from_arrays
        use grid_m, only: rg_grid
        use kim_base_m, only: kim_t
        use kim_mod_m, only: from_kim_factory_get_kim

        integer, parameter :: npts = 101

        real(dp) :: r_prof(npts), n_prof(npts), Te_prof(npts)
        real(dp) :: Ti_prof(npts), q_prof(npts), Er_prof(npts)
        class(kim_t), allocatable :: kim_instance
        real(dp) :: kr_arr(3), kr, bf, tol
        integer :: i, j, sp, jj_arr(3), k
        complex(dp) :: ref, gfused

        call make_test_profiles(npts, r_prof, n_prof, Te_prof, Ti_prof, &
                                q_prof, Er_prof)

        call write_test_namelist('./KIM_config_fourier_collapse_ks2_test.nml')
        nml_config_path = './KIM_config_fourier_collapse_ks2_test.nml'

        profiles_in_memory = .true.
        call kim_init()
        call set_profiles_from_arrays(r_prof, n_prof, Te_prof, Ti_prof, &
                                      q_prof, Er_prof, npts)

        call from_kim_factory_get_kim('electrostatic', kim_instance)
        call kim_instance%init()

        fp_ks_model = FP_KPERP_FULL

        kr_arr = [0.1_dp, 1.0_dp, 5.0_dp]
        jj_arr = [rg_grid%npts_b / 4, rg_grid%npts_b / 2, (3 * rg_grid%npts_b) / 4]

        do i = 1, size(kr_arr)
            kr = kr_arr(i)
            do k = 1, size(jj_arr)
                j = jj_arr(k)

                ref = (0.0_dp, 0.0_dp)
                do sp = 0, plasma%n_species - 1
                    if (turn_off_ions .and. sp >= 1) cycle
                    if (turn_off_electrons .and. sp == 0) cycle
                    bf = (plasma%ks(j)**2 + kr**2) * plasma%spec(sp)%rho_L(j)**2
                    ref = ref + core_rho_phi_sp(plasma, sp, bf, bf, j)
                end do
                ref = ref / (4.0d0 * pi)

                gfused = hatG_rho_phi(plasma, kr, kr, j)

                tol = 5.0e-12_dp * (1.0_dp + abs(ref))
                if (abs(gfused - ref) >= tol) then
                    print *, 'FAIL: rho-Phi ks2 collapse mismatch at kr=', kr, ' j=', j
                    print *, '  core(bf,bf)/4pi = ', ref
                    print *, '  fused           = ', gfused
                    print *, '  |diff|          = ', abs(gfused - ref), ' tol = ', tol
                    error stop
                end if
            end do
        end do

        fp_ks_model = FP_KR_ONLY_APPROX

        print *, 'PASS: full rho-Phi model collapses to core(bf,bf) at kr=krp'
    end subroutine

    subroutine test_zero_wavenumber_limit()
        ! Analytic b -> 0 anchor in the named radial-only approximation: at
        ! kr = krp = 0,
        ! both FLR arguments vanish (b_+ = b_x = 0). The Bessel/exp assembly then
        ! reduces to a closed form INDEPENDENT of the Bessel routine, because
        ! exp(-b_+) = 1, (1 - b_+) = 1, I_0(0) = 1, I_{-1}(0) = 0, and the
        ! A2 * b_x * I_{-1}(b_x) term vanishes with b_x = 0. The Fourier phase
        ! exp(i*0*rg) = 1. So the per-species FP core collapses to
        !   1/lambda_D^2 * i * vT^2 * ks / (omega_c*nu) * ( I00*(A1+A2) + 0.5*I20*A2 )
        ! and the Debye term stays -1/lambda_D^2. We hard-code I_0(0)=1 and
        ! I_{-1}(0)=0 here rather than calling bessel_in, so this test validates
        ! the kernel's argument assembly against math, not against the Bessel
        ! routine it shares with the collapse/phase tests. Only valid for the
        ! radial-only model (the full model has arguments rho_L^2*ks^2 /= 0).
        use config_m, only: profiles_in_memory, nml_config_path
        use config_m, only: artificial_debye_case, turn_off_ions, turn_off_electrons
        use flr2_fourier_kernel_m, only: hatG_rho_phi
        use constants_m, only: com_unit, pi
        use species_m, only: plasma, set_profiles_from_arrays
        use grid_m, only: rg_grid
        use kim_base_m, only: kim_t
        use kim_mod_m, only: from_kim_factory_get_kim

        integer, parameter :: npts = 101

        real(dp) :: r_prof(npts), n_prof(npts), Te_prof(npts)
        real(dp) :: Ti_prof(npts), q_prof(npts), Er_prof(npts)
        class(kim_t), allocatable :: kim_instance
        real(dp) :: tol, ks
        integer :: j, sp, jj_arr(3), k
        complex(dp) :: closed, gfused

        call make_test_profiles(npts, r_prof, n_prof, Te_prof, Ti_prof, &
                                q_prof, Er_prof)

        call write_test_namelist('./KIM_config_fourier_zero_wavenumber_test.nml')
        nml_config_path = './KIM_config_fourier_zero_wavenumber_test.nml'

        profiles_in_memory = .true.
        call kim_init()
        call set_profiles_from_arrays(r_prof, n_prof, Te_prof, Ti_prof, &
                                      q_prof, Er_prof, npts)

        call from_kim_factory_get_kim('electrostatic', kim_instance)
        call kim_instance%init()

        ! This closed form only holds for the named radial-only approximation.
        fp_ks_model = FP_KR_ONLY_APPROX

        jj_arr = [rg_grid%npts_b / 4, rg_grid%npts_b / 2, (3 * rg_grid%npts_b) / 4]

        do k = 1, size(jj_arr)
            j = jj_arr(k)

            ! Per-species b -> 0 closed form, summed over the non-turned-off
            ! species, with I_0(0)=1 and I_{-1}(0)=0 folded in analytically.
            closed = (0.0_dp, 0.0_dp)
            do sp = 0, plasma%n_species - 1
                if (turn_off_ions .and. sp >= 1) cycle
                if (turn_off_electrons .and. sp == 0) cycle

                ks = plasma%ks(j)

                if (artificial_debye_case <= 1) then
                    closed = closed - 1.0d0 / plasma%spec(sp)%lambda_D(j)**2.0d0
                end if

                if (artificial_debye_case == 0 .or. artificial_debye_case == 2) then
                    closed = closed + 1.0d0 / plasma%spec(sp)%lambda_D(j)**2.0d0 &
                        * com_unit * plasma%spec(sp)%vT(j)**2.0d0 * ks &
                        / (plasma%spec(sp)%omega_c(j) * plasma%spec(sp)%nu(j)) * &
                        (&
                            plasma%spec(sp)%I00(j, 0) &
                                * (plasma%spec(sp)%A1(j) + plasma%spec(sp)%A2(j)) &
                            + 0.5d0 * plasma%spec(sp)%I20(j, 0) * plasma%spec(sp)%A2(j) &
                        )
                end if
            end do
            closed = closed / (4.0d0 * pi)

            gfused = hatG_rho_phi(plasma, 0.0_dp, 0.0_dp, j)

            tol = 1.0e-12_dp * (1.0_dp + abs(closed))
            if (abs(gfused - closed) >= tol) then
                print *, 'FAIL: rho-Phi b->0 limit mismatch at j=', j
                print *, '  closed = ', closed
                print *, '  fused  = ', gfused
                print *, '  |diff| = ', abs(gfused - closed), ' tol = ', tol
                error stop
            end if
        end do

        print *, 'PASS: fused rho-Phi kernel matches b->0 closed form'
    end subroutine

    subroutine test_collapse_rho_B()
        ! Collapse-by-construction for the radial-only rho-B approximation: the
        ! fused two-wavenumber hatG_rho_B evaluated on the diagonal k_r = k'_r
        ! must equal the per-species signed diagonal integrand summed over the
        ! non-turned-off species and divided by 4*pi. Exact to machine precision
        ! because both share the same per-species core.
        use config_m, only: profiles_in_memory, nml_config_path
        use config_m, only: turn_off_ions, turn_off_electrons
        use flr2_fourier_kernel_m, only: hatG_rho_B, hatG_rho_B_diag_sp
        use constants_m, only: pi
        use species_m, only: plasma, set_profiles_from_arrays
        use grid_m, only: rg_grid
        use kim_base_m, only: kim_t
        use kim_mod_m, only: from_kim_factory_get_kim

        integer, parameter :: npts = 101

        real(dp) :: r_prof(npts), n_prof(npts), Te_prof(npts)
        real(dp) :: Ti_prof(npts), q_prof(npts), Er_prof(npts)
        class(kim_t), allocatable :: kim_instance
        real(dp) :: kr_arr(3), kr, tol
        integer :: i, j, sp, jj_arr(3), k
        complex(dp) :: dref, gfused

        call make_test_profiles(npts, r_prof, n_prof, Te_prof, Ti_prof, &
                                q_prof, Er_prof)

        call write_test_namelist('./KIM_config_fourier_collapse_B_test.nml')
        nml_config_path = './KIM_config_fourier_collapse_B_test.nml'

        profiles_in_memory = .true.
        call kim_init()
        call set_profiles_from_arrays(r_prof, n_prof, Te_prof, Ti_prof, &
                                      q_prof, Er_prof, npts)

        call from_kim_factory_get_kim('electrostatic', kim_instance)
        call kim_instance%init()

        fp_ks_model = FP_KR_ONLY_APPROX

        kr_arr = [0.1_dp, 1.0_dp, 5.0_dp]
        jj_arr = [rg_grid%npts_b / 4, rg_grid%npts_b / 2, (3 * rg_grid%npts_b) / 4]

        do i = 1, size(kr_arr)
            kr = kr_arr(i)
            do k = 1, size(jj_arr)
                j = jj_arr(k)

                dref = (0.0_dp, 0.0_dp)
                do sp = 0, plasma%n_species - 1
                    if (turn_off_ions .and. sp >= 1) cycle
                    if (turn_off_electrons .and. sp == 0) cycle
                    dref = dref + hatG_rho_B_diag_sp(plasma, sp, kr, j)
                end do
                dref = dref / (4.0d0 * pi)

                gfused = hatG_rho_B(plasma, kr, kr, j)

                tol = 5.0e-12_dp * (1.0_dp + abs(dref))
                if (abs(gfused - dref) >= tol) then
                    print *, 'FAIL: rho-B collapse mismatch at kr=', kr, ' j=', j
                    print *, '  diag/4pi = ', dref
                    print *, '  fused    = ', gfused
                    print *, '  |diff|   = ', abs(gfused - dref), ' tol = ', tol
                    error stop
                end if
            end do
        end do

        print *, 'PASS: fused rho-B kernel collapses to diagonal at kr=krp'
    end subroutine

    subroutine test_collapse_rho_B_ks2()
        ! Structural collapse for the full finite-radius rho-B model. On the
        ! diagonal k_r = k'_r both FLR
        ! arguments reduce to the full perpendicular argument
        ! bf = (k_s^2 + k_r^2) * rho_L^2, so the fused kernel must equal the
        ! per-species core evaluated at (bf, bf), summed over the non-turned-off
        ! species and divided by 4*pi.
        use config_m, only: profiles_in_memory, nml_config_path
        use config_m, only: turn_off_ions, turn_off_electrons
        use flr2_fourier_kernel_m, only: hatG_rho_B, core_rho_B_sp
        use constants_m, only: pi
        use species_m, only: plasma, set_profiles_from_arrays
        use grid_m, only: rg_grid
        use kim_base_m, only: kim_t
        use kim_mod_m, only: from_kim_factory_get_kim

        integer, parameter :: npts = 101

        real(dp) :: r_prof(npts), n_prof(npts), Te_prof(npts)
        real(dp) :: Ti_prof(npts), q_prof(npts), Er_prof(npts)
        class(kim_t), allocatable :: kim_instance
        real(dp) :: kr_arr(3), kr, bf, tol
        integer :: i, j, sp, jj_arr(3), k
        complex(dp) :: ref, gfused

        call make_test_profiles(npts, r_prof, n_prof, Te_prof, Ti_prof, &
                                q_prof, Er_prof)

        call write_test_namelist('./KIM_config_fourier_collapse_B_ks2_test.nml')
        nml_config_path = './KIM_config_fourier_collapse_B_ks2_test.nml'

        profiles_in_memory = .true.
        call kim_init()
        call set_profiles_from_arrays(r_prof, n_prof, Te_prof, Ti_prof, &
                                      q_prof, Er_prof, npts)

        call from_kim_factory_get_kim('electrostatic', kim_instance)
        call kim_instance%init()

        fp_ks_model = FP_KPERP_FULL

        kr_arr = [0.1_dp, 1.0_dp, 5.0_dp]
        jj_arr = [rg_grid%npts_b / 4, rg_grid%npts_b / 2, (3 * rg_grid%npts_b) / 4]

        do i = 1, size(kr_arr)
            kr = kr_arr(i)
            do k = 1, size(jj_arr)
                j = jj_arr(k)

                ref = (0.0_dp, 0.0_dp)
                do sp = 0, plasma%n_species - 1
                    if (turn_off_ions .and. sp >= 1) cycle
                    if (turn_off_electrons .and. sp == 0) cycle
                    bf = (plasma%ks(j)**2 + kr**2) * plasma%spec(sp)%rho_L(j)**2
                    ref = ref + core_rho_B_sp(plasma, sp, bf, bf, j)
                end do
                ref = ref / (4.0d0 * pi)

                gfused = hatG_rho_B(plasma, kr, kr, j)

                tol = 5.0e-12_dp * (1.0_dp + abs(ref))
                if (abs(gfused - ref) >= tol) then
                    print *, 'FAIL: rho-B ks2 collapse mismatch at kr=', kr, ' j=', j
                    print *, '  core(bf,bf)/4pi = ', ref
                    print *, '  fused           = ', gfused
                    print *, '  |diff|          = ', abs(gfused - ref), ' tol = ', tol
                    error stop
                end if
            end do
        end do

        fp_ks_model = FP_KR_ONLY_APPROX

        print *, 'PASS: full rho-B model collapses to core(bf,bf) at kr=krp'
    end subroutine

    subroutine test_zero_wavenumber_limit_rho_B()
        ! Analytic b -> 0 anchor for the radial-only rho-B approximation: at
        ! kr = krp = 0 both FLR arguments vanish (b_+ = b_x = 0). With
        ! I_0(0) = 1, I_{-1}(0) = 0, exp(-b_+) = 1, (1 - b_+) = 1 and the
        ! A2 * b_x * I_{-1}(b_x) term vanishing, and the Fourier phase = 1, the
        ! per-species signed FP core collapses to
        !   -1/lambda_D^2 * vT^3/(omega_c*nu*sol) * ( I01*(A1+A2) + 0.5*I21*A2 ).
        ! We fold I_0(0)=1 and I_{-1}(0)=0 in analytically (no bessel_in call),
        ! so this validates the kernel's argument assembly against math, not the
        ! Bessel routine. rho-B has NO Debye term. Only valid for this named
        ! approximation.
        use config_m, only: profiles_in_memory, nml_config_path
        use config_m, only: artificial_debye_case, turn_off_ions, turn_off_electrons
        use flr2_fourier_kernel_m, only: hatG_rho_B
        use constants_m, only: pi, sol
        use species_m, only: plasma, set_profiles_from_arrays
        use grid_m, only: rg_grid
        use kim_base_m, only: kim_t
        use kim_mod_m, only: from_kim_factory_get_kim

        integer, parameter :: npts = 101

        real(dp) :: r_prof(npts), n_prof(npts), Te_prof(npts)
        real(dp) :: Ti_prof(npts), q_prof(npts), Er_prof(npts)
        class(kim_t), allocatable :: kim_instance
        real(dp) :: tol
        integer :: j, sp, jj_arr(3), k
        complex(dp) :: closed, gfused

        call make_test_profiles(npts, r_prof, n_prof, Te_prof, Ti_prof, &
                                q_prof, Er_prof)

        call write_test_namelist('./KIM_config_fourier_zero_wavenumber_B_test.nml')
        nml_config_path = './KIM_config_fourier_zero_wavenumber_B_test.nml'

        profiles_in_memory = .true.
        call kim_init()
        call set_profiles_from_arrays(r_prof, n_prof, Te_prof, Ti_prof, &
                                      q_prof, Er_prof, npts)

        call from_kim_factory_get_kim('electrostatic', kim_instance)
        call kim_instance%init()

        ! Closed form derived above only holds for the radial-only approximation.
        fp_ks_model = FP_KR_ONLY_APPROX

        jj_arr = [rg_grid%npts_b / 4, rg_grid%npts_b / 2, (3 * rg_grid%npts_b) / 4]

        do k = 1, size(jj_arr)
            j = jj_arr(k)

            ! Per-species b -> 0 closed form, summed over the non-turned-off
            ! species, with I_0(0)=1 and I_{-1}(0)=0 folded in analytically.
            ! rho-B carries no Debye term.
            closed = (0.0_dp, 0.0_dp)
            do sp = 0, plasma%n_species - 1
                if (turn_off_ions .and. sp >= 1) cycle
                if (turn_off_electrons .and. sp == 0) cycle

                if (artificial_debye_case == 0 .or. artificial_debye_case == 2) then
                    closed = closed - 1.0d0 / plasma%spec(sp)%lambda_D(j)**2.0d0 &
                        * plasma%spec(sp)%vT(j)**3.0d0 &
                        / (plasma%spec(sp)%omega_c(j) * plasma%spec(sp)%nu(j) * sol) * &
                        (&
                            plasma%spec(sp)%I01(j, 0) &
                                * (plasma%spec(sp)%A1(j) + plasma%spec(sp)%A2(j)) &
                            + 0.5d0 * plasma%spec(sp)%I21(j, 0) * plasma%spec(sp)%A2(j) &
                        )
                end if
            end do
            closed = closed / (4.0d0 * pi)

            gfused = hatG_rho_B(plasma, 0.0_dp, 0.0_dp, j)

            tol = 1.0e-12_dp * (1.0_dp + abs(closed))
            if (abs(gfused - closed) >= tol) then
                print *, 'FAIL: rho-B b->0 limit mismatch at j=', j
                print *, '  closed = ', closed
                print *, '  fused  = ', gfused
                print *, '  |diff| = ', abs(gfused - closed), ' tol = ', tol
                error stop
            end if
        end do

        print *, 'PASS: fused rho-B kernel matches b->0 closed form'
    end subroutine

    subroutine test_scaled_bessel_large_arg()
        ! Overflow-safe scaled Bessel products in the large-argument regime.
        ! At b_+ = b_x = 1000 the unscaled I_0(1000) overflows to +Inf (near
        ! the ~710 argument limit), so the naive exp(-b_+)*I_0(b_x) is Inf*0.
        ! The asymptotic branch must instead return finite values matching the
        ! analytic limits exp(-b)*I_0(b) -> 1/sqrt(2*pi*b) and I_{-1}/I_0 -> 1.
        use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
        use flr2_fourier_kernel_m, only: scaled_bessel_pair
        use constants_m, only: pi

        real(dp) :: b, sI0, sIm1, analytic_sI0, tol

        b = 1000.0_dp
        call scaled_bessel_pair(b, b, sI0, sIm1)

        if (.not. ieee_is_finite(sI0)) then
            print *, 'FAIL: scaled_bessel_pair sI0 not finite at b=1000: ', sI0
            error stop
        end if
        if (.not. ieee_is_finite(sIm1)) then
            print *, 'FAIL: scaled_bessel_pair sIm1 not finite at b=1000: ', sIm1
            error stop
        end if

        ! exp(-b)*I_0(b) -> 1/sqrt(2*pi*b) as b -> infinity.
        analytic_sI0 = 1.0_dp / sqrt(2.0_dp * pi * b)
        tol = 1.0e-3_dp * analytic_sI0
        if (abs(sI0 - analytic_sI0) >= tol) then
            print *, 'FAIL: scaled sI0 off analytic limit at b=1000'
            print *, '  sI0      = ', sI0
            print *, '  analytic = ', analytic_sI0
            print *, '  |diff|   = ', abs(sI0 - analytic_sI0), ' tol = ', tol
            error stop
        end if

        ! I_{-1}(b)/I_0(b) -> 1, so sIm1 is positive and within ~2% of sI0.
        if (.not. (sIm1 > 0.0_dp)) then
            print *, 'FAIL: scaled sIm1 must be positive at b=1000, got ', sIm1
            error stop
        end if
        if (abs(sIm1 - sI0) >= 0.02_dp * sI0) then
            print *, 'FAIL: scaled sIm1 not within 2% of sI0 at b=1000'
            print *, '  sI0    = ', sI0
            print *, '  sIm1   = ', sIm1
            print *, '  |diff| = ', abs(sIm1 - sI0), ' tol = ', 0.02_dp * sI0
            error stop
        end if

        print *, 'PASS: scaled_bessel_pair finite and matches large-arg asymptotics'
    end subroutine

    subroutine test_large_bcross_kernel_finite()
        ! End-to-end overflow guard: the per-species cores must stay finite when
        ! driven with a forced large FLR argument (b_x = 1000) that would
        ! overflow the naive exp(-b_+)*I_0(b_x) path. Uses the public
        ! core_*_sp entry points with an interior guiding-centre index and the
        ! default artificial_debye_case = 0 (FP term active).
        use config_m, only: profiles_in_memory, nml_config_path
        use config_m, only: artificial_debye_case
        use flr2_fourier_kernel_m, only: core_rho_phi_sp, core_rho_B_sp
        use species_m, only: plasma, set_profiles_from_arrays
        use grid_m, only: rg_grid
        use kim_base_m, only: kim_t
        use kim_mod_m, only: from_kim_factory_get_kim

        integer, parameter :: npts = 101

        real(dp) :: r_prof(npts), n_prof(npts), Te_prof(npts)
        real(dp) :: Ti_prof(npts), q_prof(npts), Er_prof(npts)
        class(kim_t), allocatable :: kim_instance
        real(dp) :: big
        integer :: j
        complex(dp) :: g_phi, g_B

        call make_test_profiles(npts, r_prof, n_prof, Te_prof, Ti_prof, &
                                q_prof, Er_prof)

        call write_test_namelist('./KIM_config_fourier_large_bcross_test.nml')
        nml_config_path = './KIM_config_fourier_large_bcross_test.nml'

        profiles_in_memory = .true.
        call kim_init()
        call set_profiles_from_arrays(r_prof, n_prof, Te_prof, Ti_prof, &
                                      q_prof, Er_prof, npts)

        call from_kim_factory_get_kim('electrostatic', kim_instance)
        call kim_instance%init()

        ! The namelist sets artificial_debye_case = 0, so the FP branch (which
        ! evaluates the Bessel products) is active.
        if (.not. (artificial_debye_case == 0 .or. artificial_debye_case == 2)) then
            print *, 'FAIL: FP branch inactive; artificial_debye_case = ', artificial_debye_case
            error stop
        end if

        j = rg_grid%npts_b / 2
        big = 1000.0_dp

        g_phi = core_rho_phi_sp(plasma, 0, big, big, j)
        g_B = core_rho_B_sp(plasma, 0, big, big, j)

        if (.not. is_finite_complex(g_phi)) then
            print *, 'FAIL: core_rho_phi_sp not finite at bcross=1000: ', g_phi
            error stop
        end if
        if (.not. is_finite_complex(g_B)) then
            print *, 'FAIL: core_rho_B_sp not finite at bcross=1000: ', g_B
            error stop
        end if

        print *, 'PASS: core_*_sp stay finite at large bcross (overflow guard)'
        print *, '  core_rho_phi_sp = ', g_phi
        print *, '  core_rho_B_sp   = ', g_B
    end subroutine

    subroutine test_populated_plasma_and_kernel_stub()
        use config_m, only: profiles_in_memory, nml_config_path
        use species_m, only: plasma, set_profiles_from_arrays
        use grid_m, only: rg_grid
        use kim_base_m, only: kim_t
        use kim_mod_m, only: from_kim_factory_get_kim

        integer, parameter :: npts = 101

        real(dp) :: r_prof(npts), n_prof(npts), Te_prof(npts)
        real(dp) :: Ti_prof(npts), q_prof(npts), Er_prof(npts)
        class(kim_t), allocatable :: kim_instance
        complex(dp) :: I00_val, G_stub
        real(dp) :: rho_L_val, lambda_D_val
        integer :: j

        call make_test_profiles(npts, r_prof, n_prof, Te_prof, Ti_prof, &
                                q_prof, Er_prof)

        call write_test_namelist('./KIM_config_fourier_test.nml')
        nml_config_path = './KIM_config_fourier_test.nml'

        profiles_in_memory = .true.
        call kim_init()
        call set_profiles_from_arrays(r_prof, n_prof, Te_prof, Ti_prof, &
                                      q_prof, Er_prof, npts)

        call from_kim_factory_get_kim('electrostatic', kim_instance)
        call kim_instance%init()

        ! Interior guiding-centre boundary-grid index.
        j = rg_grid%npts_b / 2

        if (.not. allocated(plasma%spec(0)%I00)) then
            print *, 'FAIL: plasma%spec(0)%I00 not allocated after init'
            error stop
        end if
        if (.not. allocated(plasma%spec(0)%rho_L)) then
            print *, 'FAIL: plasma%spec(0)%rho_L not allocated after init'
            error stop
        end if
        if (.not. allocated(plasma%spec(0)%lambda_D)) then
            print *, 'FAIL: plasma%spec(0)%lambda_D not allocated after init'
            error stop
        end if

        I00_val = plasma%spec(0)%I00(j, 0)
        rho_L_val = plasma%spec(0)%rho_L(j)
        lambda_D_val = plasma%spec(0)%lambda_D(j)

        if (.not. is_finite_complex(I00_val)) then
            print *, 'FAIL: plasma%spec(0)%I00(j,0) is not finite: ', I00_val
            error stop
        end if
        if (.not. (rho_L_val > 0.0_dp)) then
            print *, 'FAIL: plasma%spec(0)%rho_L(j) must be positive, got ', rho_L_val
            error stop
        end if
        if (.not. (lambda_D_val > 0.0_dp)) then
            print *, 'FAIL: plasma%spec(0)%lambda_D(j) must be positive, got ', lambda_D_val
            error stop
        end if

        ! Drive the (stubbed) kernel to prove it links and is callable with a
        ! populated plasma. Diagonal argument k_r = k'_r = 1.
        G_stub = hatG_rho_phi(plasma, 1.0_dp, 1.0_dp, j)

        print *, 'PASS: electrostatic init populated plasma on rg grid'
        print *, '  index j        = ', j, ' of npts_b = ', rg_grid%npts_b
        print *, '  I00(j,0)       = ', I00_val
        print *, '  rho_L(j)       = ', rho_L_val
        print *, '  lambda_D(j)    = ', lambda_D_val
        print *, '  hatG_rho_phi   = ', G_stub
    end subroutine

    logical function is_finite_complex(z) result(ok)
        use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
        complex(dp), intent(in) :: z

        ok = ieee_is_finite(real(z, dp)) .and. ieee_is_finite(aimag(z))
    end function is_finite_complex

    subroutine make_test_profiles(npts, r_prof, n_prof, Te_prof, Ti_prof, &
                                  q_prof, Er_prof)
        ! Analytic profiles on r in [2, 60] cm with q in [1.5, 2.5];
        ! q crosses |m/n| = 2 at r ~ 31 cm, inside the solver domain
        ! [r_min, r_plas] = [10, 50].
        integer, intent(in) :: npts
        real(dp), intent(out) :: r_prof(npts), n_prof(npts), Te_prof(npts)
        real(dp), intent(out) :: Ti_prof(npts), q_prof(npts), Er_prof(npts)
        integer :: i

        do i = 1, npts
            r_prof(i) = 2.0_dp + 58.0_dp * real(i - 1, dp) / real(npts - 1, dp)
            n_prof(i) = 2.0e13_dp * (1.1_dp - r_prof(i) / 100.0_dp)
            Te_prof(i) = 1.0e3_dp * (1.2_dp - r_prof(i) / 100.0_dp)
            Ti_prof(i) = Te_prof(i)
            q_prof(i) = 1.5_dp + (r_prof(i) - 2.0_dp) / 58.0_dp
            Er_prof(i) = -0.5_dp
        end do
    end subroutine

    subroutine write_test_namelist(path)
        ! Minimal electrostatic FokkerPlanck configuration on a coarse grid,
        ! written in the exact group order kim_read_config reads them. The
        ! FokkerPlanck collision model is required so the susceptibility
        ! functions (I00, ...) are computed during set_plasma_quantities.
        character(len=*), intent(in) :: path
        integer :: iunit

        open(newunit=iunit, file=path, status='replace', action='write')
        write(iunit, '(A)') '&KIM_CONFIG'
        write(iunit, '(A)') ' number_of_ion_species = 1'
        write(iunit, '(A)') ' artificial_debye_case = 0'
        write(iunit, '(A)') " type_of_run = 'electrostatic'"
        write(iunit, '(A)') " collision_model = 'FokkerPlanck'"
        write(iunit, '(A)') ' read_species_from_namelist = .false.'
        write(iunit, '(A)') ' turn_off_ions = .false.'
        write(iunit, '(A)') ' turn_off_electrons = .false.'
        write(iunit, '(A)') " plasma_type = 'D'"
        write(iunit, '(A)') ' rescale_density = .false.'
        write(iunit, '(A)') ' number_density_rescale = 1.0'
        write(iunit, '(A)') ' ion_flr_scale_factor = 1.0'
        write(iunit, '(A)') '/'
        write(iunit, '(A)') '&WKB_DISPERSION'
        write(iunit, '(A)') '/'
        write(iunit, '(A)') '&KIM_IO'
        write(iunit, '(A)') " profile_location = './'"
        write(iunit, '(A)') " output_path = './out_fourier_test/'"
        write(iunit, '(A)') ' hdf5_input = .false.'
        write(iunit, '(A)') ' hdf5_output = .false.'
        write(iunit, '(A)') ' log_level = 3'
        write(iunit, '(A)') ' data_verbosity = 0'
        write(iunit, '(A)') ' calculate_asymptotics = .false.'
        write(iunit, '(A)') " h5_out_file = ''"
        write(iunit, '(A)') ' write_diagnostics_dat = .false.'
        write(iunit, '(A)') '/'
        write(iunit, '(A)') '&KIM_SETUP'
        write(iunit, '(A)') ' btor = -18000.0'
        write(iunit, '(A)') ' R0 = 165.0'
        write(iunit, '(A)') ' m_mode = -2'
        write(iunit, '(A)') ' n_mode = 1'
        write(iunit, '(A)') ' omega = 0.0'
        write(iunit, '(A)') ' spline_base = 1'
        write(iunit, '(A)') ' type_br_field = 12'
        write(iunit, '(A)') ' collisions_off = .false.'
        write(iunit, '(A)') ' set_profiles_constant = 0'
        write(iunit, '(A)') ' bc_type = 3'
        write(iunit, '(A)') ' mphi_max = 0'
        write(iunit, '(A)') ' Br_boundary_re = 1.0'
        write(iunit, '(A)') ' Br_boundary_im = 0.0'
        write(iunit, '(A)') '/'
        write(iunit, '(A)') '&KIM_GRID'
        write(iunit, '(A)') " grid_spacing_rg = 'equidistant'"
        write(iunit, '(A)') " grid_spacing_xl = 'equidistant'"
        write(iunit, '(A)') ' l_space_dim = 32'
        write(iunit, '(A)') ' rg_space_dim = 32'
        write(iunit, '(A)') " theta_integration = 'GaussLegendre'"
        write(iunit, '(A)') " theta_integration_method = ''"
        write(iunit, '(A)') ' Larmor_skip_factor = 5'
        write(iunit, '(A)') ' gauss_int_nodes_Ntheta = 17'
        write(iunit, '(A)') ' gauss_int_nodes_Nx = 31'
        write(iunit, '(A)') ' gauss_int_nodes_Nxp = 30'
        write(iunit, '(A)') ' r_plas = 50.0'
        write(iunit, '(A)') ' r_min = 10.0'
        write(iunit, '(A)') ' width_res = 0.5'
        write(iunit, '(A)') ' ampl_res = 15.0'
        write(iunit, '(A)') ' hrmax_scaling = 1.0'
        write(iunit, '(A)') ' rkf45_atol = 1.0e-9'
        write(iunit, '(A)') ' rkf45_rtol = 1.0e-6'
        write(iunit, '(A)') ' kernel_taper_skip_threshold = 1.0e-6'
        write(iunit, '(A)') " quadpack_algorithm = 'QAG'"
        write(iunit, '(A)') ' quadpack_key = 6'
        write(iunit, '(A)') ' quadpack_limit = 500'
        write(iunit, '(A)') ' quadpack_epsabs = 1.0e-10'
        write(iunit, '(A)') ' quadpack_epsrel = 1.0e-10'
        write(iunit, '(A)') ' quadpack_use_u_substitution = .true.'
        write(iunit, '(A)') '/'
        write(iunit, '(A)') '&KIM_PROFILES'
        write(iunit, '(A)') '/'
        close(iunit)
    end subroutine

end program
