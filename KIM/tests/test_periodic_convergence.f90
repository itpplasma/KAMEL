program test_periodic_convergence
    ! Phase-3 Task 1: exercises the reusable periodic core
    ! rt_electrostatic_periodic_m::compute_periodic_delta_phi extracted from the
    ! electrostatic_periodic run-type.
    !
    ! The convergence harness (P3.2-3.4) needs to evaluate delta-Phi on a FIXED
    ! diagnostic grid independent of the window discretization (n_rg, M, dx), so
    ! that Cauchy residuals are comparable across parameter refinements. This
    ! test pins the contract of that core: given explicit window parameters, a
    ! constant Br drive, and an explicit output grid, it builds -> assembles ->
    ! solves -> reconstructs and returns delta-Phi on the SUPPLIED grid.
    !
    ! Setup mirrors test_kim_solver_periodic.f90 / test_periodic_assembly.f90:
    ! in-memory analytic profiles, (m,n) = (-6, 2) so q is resonant at q = 3,
    ! type_br_field = 12 (constant Br over the window).
    !
    ! Assertions on compute_periodic_delta_phi with a FIXED 41-point diagnostic
    ! grid (independent of n_rg = 96 / M = 24):
    !   info == 0 (non-singular solve);
    !   size(dPhi) == size(r_diag) == 41 (evaluated on the SUPPLIED grid);
    !   all dPhi finite;
    !   maxval(abs(dPhi)) > 0 (constant Br drives a response).

    use KIM_kinds_m, only: dp

    implicit none

    call test_reference_scale_validation()
    call test_core_on_fixed_grid()
    call test_nrg_convergence()
    call test_M_convergence()
    call test_dr_deformation_scan()

    print *, 'All tests PASSED'
    stop 0

contains

    subroutine test_reference_scale_validation()
        use, intrinsic :: ieee_arithmetic, only: ieee_quiet_nan, ieee_value
        use species_m, only: plasma_t
        use rt_electrostatic_periodic_m, only: select_periodic_reference_scale, &
                                               PERIODIC_SCALE_NO_ACTIVE, &
                                               PERIODIC_SCALE_INVALID_RHO

        type(plasma_t) :: test_plasma
        real(dp) :: rho_ref
        integer :: sp, reference_species, info

        test_plasma%n_species = 3
        test_plasma%grid_size = 4
        allocate(test_plasma%r_grid(4), test_plasma%spec(0:2))
        test_plasma%r_grid = [0.0_dp, 1.0_dp, 2.0_dp, 3.0_dp]
        do sp = 0, 2
            allocate(test_plasma%spec(sp)%rho_L(4))
            test_plasma%spec(sp)%rho_L = 0.1_dp * real(sp + 1, dp)
        end do

        call select_periodic_reference_scale(test_plasma, 1.5_dp, .false., .false., &
                                             rho_ref, reference_species, info)
        if (info /= PERIODIC_SCALE_NO_ACTIVE) then
            print *, 'FAIL: selector did not reject a plasma with no active species'
            error stop
        end if

        test_plasma%spec(2)%rho_L(2) = ieee_value(0.0_dp, ieee_quiet_nan)
        call select_periodic_reference_scale(test_plasma, 1.5_dp, .true., .true., &
                                             rho_ref, reference_species, info)
        if (info /= PERIODIC_SCALE_INVALID_RHO) then
            print *, 'FAIL: selector accepted a non-finite active-species rho_L'
            print *, '  info =', info, ' reference species =', reference_species
            error stop
        end if
        print *, 'PASS: reference-scale selector rejects invalid active species'
    end subroutine test_reference_scale_validation

    subroutine test_core_on_fixed_grid()
        use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
        use config_m, only: profiles_in_memory, nml_config_path
        use species_m, only: set_profiles_from_arrays
        use kim_base_m, only: kim_t
        use kim_mod_m, only: from_kim_factory_get_kim
        use kim_resonances_m, only: r_res
        use rt_electrostatic_periodic_m, only: compute_periodic_delta_phi

        integer, parameter :: npts = 201
        integer, parameter :: M = 24, n_rg = 96, ndiag = 41

        real(dp) :: r_prof(npts), n_prof(npts), Te_prof(npts)
        real(dp) :: Ti_prof(npts), q_prof(npts), Er_prof(npts)
        class(kim_t), allocatable :: kim_instance

        complex(dp), allocatable :: dPhi(:)
        real(dp) :: r_diag(ndiag)
        real(dp) :: rm, rhoL_rm, dx_asis, dx_tr
        integer :: i, info

        call make_test_profiles(npts, r_prof, n_prof, Te_prof, Ti_prof, &
                                q_prof, Er_prof)

        call write_test_namelist('./KIM_config_periodic_conv_test.nml')
        nml_config_path = './KIM_config_periodic_conv_test.nml'

        profiles_in_memory = .true.
        call kim_init()
        call set_profiles_from_arrays(r_prof, n_prof, Te_prof, Ti_prof, &
                                      q_prof, Er_prof, npts)

        call from_kim_factory_get_kim('electrostatic', kim_instance)
        call kim_instance%init()

        ! Locate the resonance (q = |m/n| = 3) on the global plasma.
        call prepare_resonances
        if (.not. (r_res > 0.0_dp)) then
            print *, 'FAIL: resonance not found, r_res = ', r_res
            error stop
        end if
        rm = r_res
        print *, 'resonant radius rm = ', rm

        ! Reference Larmor radius at rm, selected from the same active species
        ! set used by the production run type.
        rhoL_rm = interp_rho_L(rm)
        if (.not. (rhoL_rm > 0.0_dp)) then
            print *, 'FAIL: rho_L(rm) not positive, = ', rhoL_rm
            error stop
        end if
        call assert_convergence_uses_reference_scale(rm, rhoL_rm)
        print *, 'rho_L(rm) [active-species reference] = ', rhoL_rm

        dx_asis =  5.0_dp * rhoL_rm
        dx_tr   = 10.0_dp * rhoL_rm

        ! FIXED diagnostic grid: 41 equidistant points on [rm - 3 rho_L, rm + 3 rho_L],
        ! chosen INDEPENDENTLY of n_rg / M / dx (the whole point of the core).
        do i = 1, ndiag
            r_diag(i) = (rm - 3.0_dp * rhoL_rm) &
                      + 6.0_dp * rhoL_rm * real(i - 1, dp) / real(ndiag - 1, dp)
        end do

        call compute_periodic_delta_phi(rm, dx_asis, dx_tr, M, n_rg, &
                                        (1.0_dp, 0.0_dp), r_diag, dPhi, info)

        if (info /= 0) then
            print *, 'FAIL: compute_periodic_delta_phi info /= 0:', info
            error stop
        end if
        print *, 'PASS: compute_periodic_delta_phi info == 0'

        if (size(dPhi) /= ndiag) then
            print *, 'FAIL: size(dPhi) /=', ndiag, ' got ', size(dPhi)
            error stop
        end if
        print *, 'PASS: dPhi evaluated on the supplied grid, size ', ndiag

        do i = 1, ndiag
            if (.not. ieee_is_finite(real(dPhi(i), dp)) .or. &
                .not. ieee_is_finite(aimag(dPhi(i)))) then
                print *, 'FAIL: non-finite dPhi at', i
                error stop
            end if
        end do
        print *, 'PASS: all dPhi finite'

        if (.not. (maxval(abs(dPhi)) > 0.0_dp)) then
            print *, 'FAIL: dPhi all zero (constant Br should drive response)'
            error stop
        end if
        print *, 'PASS: dPhi non-zero, max|dPhi| =', maxval(abs(dPhi))
    end subroutine test_core_on_fixed_grid

    subroutine assert_convergence_uses_reference_scale(rm, convergence_rho)
        use config_m, only: turn_off_electrons, turn_off_ions
        use species_m, only: plasma
        use rt_electrostatic_periodic_m, only: select_periodic_reference_scale, &
                                               PERIODIC_SCALE_OK

        real(dp), intent(in) :: rm, convergence_rho
        real(dp) :: selected_rho, tol
        integer :: reference_species, info

        call select_periodic_reference_scale(plasma, rm, .not. turn_off_electrons, &
                                             .not. turn_off_ions, selected_rho, &
                                             reference_species, info)
        if (info /= PERIODIC_SCALE_OK) then
            print *, 'FAIL: convergence harness could not select a physical scale'
            error stop
        end if
        tol = 1.0e-10_dp * selected_rho
        if (abs(convergence_rho - selected_rho) > tol) then
            print *, 'FAIL: convergence scans do not use the selected species scale'
            print *, '  convergence rho_L =', convergence_rho
            print *, '  selected species =', reference_species
            print *, '  selected rho_L =', selected_rho
            error stop
        end if
    end subroutine assert_convergence_uses_reference_scale

    !> Phase-3 Task 2: r_g QUADRATURE convergence of the periodic solver.
    !>
    !> With M and the window (dx_asis, dx_tr) FIXED, refine only the number of
    !> radial quadrature nodes n_rg and reconstruct delta-Phi on a FIXED
    !> diagnostic grid. The periodic-trapezoidal rule is spectrally accurate for
    !> the near-periodic window integrand, so the Cauchy residual between
    !> successive refinements must fall fast and beat 1e-2. A failure to converge
    !> is a real finding about the quadrature / seam handling, not a reason to
    !> loosen the tolerance.
    subroutine test_nrg_convergence()
        use config_m, only: profiles_in_memory, nml_config_path
        use species_m, only: set_profiles_from_arrays
        use kim_base_m, only: kim_t
        use kim_mod_m, only: from_kim_factory_get_kim
        use kim_resonances_m, only: r_res
        use rt_electrostatic_periodic_m, only: compute_periodic_delta_phi

        integer, parameter :: npts = 201
        integer, parameter :: M = 32, ndiag = 41, nseq = 4
        integer, parameter :: n_rg_seq(nseq) = [48, 96, 192, 384]

        real(dp) :: r_prof(npts), n_prof(npts), Te_prof(npts)
        real(dp) :: Ti_prof(npts), q_prof(npts), Er_prof(npts)
        class(kim_t), allocatable :: kim_instance

        complex(dp), allocatable :: dPhi(:)
        complex(dp) :: dphi_k(ndiag, nseq)
        real(dp) :: r_diag(ndiag)
        real(dp) :: res_L2(nseq), res_max(nseq)
        real(dp) :: rm, rhoL_rm, dx_asis, dx_tr
        real(dp) :: denom_L2, denom_max
        integer :: i, k, info
        logical :: decreasing

        print *, ''
        print *, '=== test_nrg_convergence: r_g quadrature convergence ==='

        call make_test_profiles(npts, r_prof, n_prof, Te_prof, Ti_prof, &
                                q_prof, Er_prof)

        call write_test_namelist('./KIM_config_periodic_conv_test.nml')
        nml_config_path = './KIM_config_periodic_conv_test.nml'

        profiles_in_memory = .true.
        call kim_init()
        call set_profiles_from_arrays(r_prof, n_prof, Te_prof, Ti_prof, &
                                      q_prof, Er_prof, npts)

        call from_kim_factory_get_kim('electrostatic', kim_instance)
        call kim_instance%init()

        call prepare_resonances
        if (.not. (r_res > 0.0_dp)) then
            print *, 'FAIL: resonance not found, r_res = ', r_res
            error stop
        end if
        rm = r_res

        rhoL_rm = interp_rho_L(rm)
        if (.not. (rhoL_rm > 0.0_dp)) then
            print *, 'FAIL: rho_L(rm) not positive, = ', rhoL_rm
            error stop
        end if
        print *, 'resonant radius rm       = ', rm
        print *, 'rho_L(rm) [active-species reference] = ', rhoL_rm

        ! FIXED window (M and dx_asis / dx_tr do NOT change with n_rg): the only
        ! knob refined here is the r_g quadrature node count.
        dx_asis =  5.0_dp * rhoL_rm
        dx_tr   = 10.0_dp * rhoL_rm
        print *, 'FIXED M                  = ', M
        print *, 'FIXED dx_asis (5 rho_L)  = ', dx_asis
        print *, 'FIXED dx_tr  (10 rho_L)  = ', dx_tr
        print *, 'window L/2 (15 rho_L)    = ', dx_asis + dx_tr

        ! FIXED diagnostic grid: 41 equidistant points on [rm - 3 rho_L,
        ! rm + 3 rho_L], well inside L/2 = 15 rho_L and INDEPENDENT of n_rg, so
        ! the Cauchy residuals are comparable across refinements.
        do i = 1, ndiag
            r_diag(i) = (rm - 3.0_dp * rhoL_rm) &
                      + 6.0_dp * rhoL_rm * real(i - 1, dp) / real(ndiag - 1, dp)
        end do

        ! Refine n_rg only; store delta-Phi on the fixed grid for each n_rg.
        do k = 1, nseq
            call compute_periodic_delta_phi(rm, dx_asis, dx_tr, M, n_rg_seq(k), &
                                            (1.0_dp, 0.0_dp), r_diag, dPhi, info)
            if (info /= 0) then
                print *, 'FAIL: compute_periodic_delta_phi info /= 0 at n_rg =', &
                    n_rg_seq(k), ' info =', info
                error stop
            end if
            if (size(dPhi) /= ndiag) then
                print *, 'FAIL: size(dPhi) /=', ndiag, ' got ', size(dPhi)
                error stop
            end if
            dphi_k(:, k) = dPhi
        end do
        print *, 'PASS: all four solves returned info == 0'

        ! Cauchy relative residual between successive refinements.
        res_L2  = 0.0_dp
        res_max = 0.0_dp
        do k = 2, nseq
            denom_L2  = sqrt(sum(abs(dphi_k(:, k))**2))
            denom_max = maxval(abs(dphi_k(:, k)))
            res_L2(k)  = sqrt(sum(abs(dphi_k(:, k) - dphi_k(:, k-1))**2)) / denom_L2
            res_max(k) = maxval(abs(dphi_k(:, k) - dphi_k(:, k-1))) / denom_max
        end do

        print *, ''
        print *, '--- n_rg convergence (fixed M, window, diagnostic grid) ---'
        print '(A)', '   n_rg pair          res_L2            res_max'
        do k = 2, nseq
            print '(3X,I5," ->",I5,4X,ES16.8,2X,ES16.8)', &
                n_rg_seq(k-1), n_rg_seq(k), res_L2(k), res_max(k)
        end do
        print *, '-----------------------------------------------------------'

        ! Assert: the residual sequence DECREASES and the finest beats 1e-2.
        decreasing = (res_L2(4) < res_L2(3)) .and. (res_L2(3) < res_L2(2))
        if (.not. decreasing) then
            print *, 'FAIL: res_L2 sequence does not decrease monotonically:'
            print *, '   res_L2 =', res_L2(2), res_L2(3), res_L2(4)
            print *, '   -> real finding about r_g quadrature / seam handling;'
            print *, '      investigate, do NOT loosen the threshold blindly.'
            error stop 'test_nrg_convergence: residual not decreasing'
        end if
        if (.not. (res_L2(4) < 1.0e-2_dp)) then
            print *, 'FAIL: finest res_L2(4) =', res_L2(4), ' not < 1.0e-2'
            print *, '   res_L2 =', res_L2(2), res_L2(3), res_L2(4)
            print *, '   -> quadrature has NOT converged; investigate, do NOT'
            print *, '      loosen the threshold blindly.'
            error stop 'test_nrg_convergence: quadrature not converged'
        end if

        print *, 'PASS: res_L2 decreasing and res_L2(finest) < 1.0e-2'
        print *, '=== test_nrg_convergence PASSED ==='
    end subroutine test_nrg_convergence

    !> Phase-3 Task 3: FOURIER BASIS truncation convergence of the periodic
    !> solver.
    !>
    !> With the r_g quadrature (n_rg = 512, converged in P3.2) and the window
    !> (dx_asis, dx_tr) FIXED, refine only the number of Fourier modes M (the
    !> basis has 2M+1 modes) and reconstruct delta-Phi on a FIXED diagnostic
    !> grid. The localized solution's Fourier coefficients decay for
    !> |k_m| >> 1/rho_L, so the Cauchy residual between successive M must fall
    !> and beat 1e-2 at the finest tested M. A failure to converge is a real finding about
    !> the basis truncation or the k-resolution, not a reason to loosen the
    !> tolerance (bump n_rg first: the exp(i(k_m - k_m')r) quadrature is exact
    !> for |m - m'| < n_rg by Nyquist, and 2M = 256 < 512, so n_rg = 512 resolves
    !> M up to 128).
    subroutine test_M_convergence()
        use config_m, only: profiles_in_memory, nml_config_path
        use species_m, only: set_profiles_from_arrays
        use kim_base_m, only: kim_t
        use kim_mod_m, only: from_kim_factory_get_kim
        use kim_resonances_m, only: r_res
        use rt_electrostatic_periodic_m, only: compute_periodic_delta_phi

        integer, parameter :: npts = 201
        integer, parameter :: n_rg = 512, ndiag = 41, nseq = 5
        integer, parameter :: M_seq(nseq) = [32, 48, 64, 96, 128]

        real(dp) :: r_prof(npts), n_prof(npts), Te_prof(npts)
        real(dp) :: Ti_prof(npts), q_prof(npts), Er_prof(npts)
        class(kim_t), allocatable :: kim_instance

        complex(dp), allocatable :: dPhi(:)
        complex(dp) :: dphi_k(ndiag, nseq)
        real(dp) :: r_diag(ndiag)
        real(dp) :: res_L2(nseq), res_max(nseq)
        real(dp) :: rm, rhoL_rm, dx_asis, dx_tr
        real(dp) :: denom_L2, denom_max
        integer :: i, k, info
        logical :: decreasing

        print *, ''
        print *, '=== test_M_convergence: Fourier basis truncation convergence ==='

        call make_test_profiles(npts, r_prof, n_prof, Te_prof, Ti_prof, &
                                q_prof, Er_prof)

        call write_test_namelist('./KIM_config_periodic_conv_test.nml')
        nml_config_path = './KIM_config_periodic_conv_test.nml'

        profiles_in_memory = .true.
        call kim_init()
        call set_profiles_from_arrays(r_prof, n_prof, Te_prof, Ti_prof, &
                                      q_prof, Er_prof, npts)

        call from_kim_factory_get_kim('electrostatic', kim_instance)
        call kim_instance%init()

        call prepare_resonances
        if (.not. (r_res > 0.0_dp)) then
            print *, 'FAIL: resonance not found, r_res = ', r_res
            error stop
        end if
        rm = r_res

        rhoL_rm = interp_rho_L(rm)
        if (.not. (rhoL_rm > 0.0_dp)) then
            print *, 'FAIL: rho_L(rm) not positive, = ', rhoL_rm
            error stop
        end if
        print *, 'resonant radius rm       = ', rm
        print *, 'rho_L(rm) [active-species reference] = ', rhoL_rm

        ! FIXED quadrature and window (only M is refined here). n_rg = 512 is the
        ! P3.2-converged node count; it resolves the k-modes for M up to 128.
        dx_asis =  5.0_dp * rhoL_rm
        dx_tr   = 10.0_dp * rhoL_rm
        print *, 'FIXED n_rg               = ', n_rg
        print *, 'FIXED dx_asis (5 rho_L)  = ', dx_asis
        print *, 'FIXED dx_tr  (10 rho_L)  = ', dx_tr
        print *, 'window L/2 (15 rho_L)    = ', dx_asis + dx_tr

        ! FIXED diagnostic grid: 41 equidistant points on [rm - 3 rho_L,
        ! rm + 3 rho_L], INDEPENDENT of M, so the Cauchy residuals are
        ! comparable across basis refinements.
        do i = 1, ndiag
            r_diag(i) = (rm - 3.0_dp * rhoL_rm) &
                      + 6.0_dp * rhoL_rm * real(i - 1, dp) / real(ndiag - 1, dp)
        end do

        ! Refine M only; store delta-Phi on the fixed grid for each M.
        do k = 1, nseq
            call compute_periodic_delta_phi(rm, dx_asis, dx_tr, M_seq(k), n_rg, &
                                            (1.0_dp, 0.0_dp), r_diag, dPhi, info)
            if (info /= 0) then
                print *, 'FAIL: compute_periodic_delta_phi info /= 0 at M =', &
                    M_seq(k), ' info =', info
                error stop
            end if
            if (size(dPhi) /= ndiag) then
                print *, 'FAIL: size(dPhi) /=', ndiag, ' got ', size(dPhi)
                error stop
            end if
            dphi_k(:, k) = dPhi
        end do
        print *, 'PASS: all Fourier-refinement solves returned info == 0'

        ! Cauchy relative residual between successive refinements.
        res_L2  = 0.0_dp
        res_max = 0.0_dp
        do k = 2, nseq
            denom_L2  = sqrt(sum(abs(dphi_k(:, k))**2))
            denom_max = maxval(abs(dphi_k(:, k)))
            res_L2(k)  = sqrt(sum(abs(dphi_k(:, k) - dphi_k(:, k-1))**2)) / denom_L2
            res_max(k) = maxval(abs(dphi_k(:, k) - dphi_k(:, k-1))) / denom_max
        end do

        print *, ''
        print *, '--- M convergence (fixed n_rg, window, diagnostic grid) ---'
        print '(A)', '   M pair             res_L2            res_max'
        do k = 2, nseq
            print '(3X,I5," ->",I5,4X,ES16.8,2X,ES16.8)', &
                M_seq(k-1), M_seq(k), res_L2(k), res_max(k)
        end do
        print *, '-----------------------------------------------------------'

        ! Assert: the residual sequence DECREASES and the finest beats 1e-2.
        decreasing = .true.
        do k = 3, nseq
            decreasing = decreasing .and. (res_L2(k) < res_L2(k - 1))
        end do
        if (.not. decreasing) then
            print *, 'FAIL: res_L2 sequence does not decrease monotonically:'
            print *, '   res_L2 =', res_L2(2:nseq)
            print *, '   -> real finding about basis truncation / k-resolution;'
            print *, '      investigate (bump n_rg first), do NOT loosen the'
            print *, '      threshold blindly.'
            error stop 'test_M_convergence: residual not decreasing'
        end if
        if (.not. (res_L2(nseq) < 1.0e-2_dp)) then
            print *, 'FAIL: finest res_L2 =', res_L2(nseq), ' not < 1.0e-2'
            print *, '   res_L2 =', res_L2(2:nseq)
            print *, '   -> Fourier basis has NOT converged; investigate (bump'
            print *, '      n_rg first), do NOT loosen the threshold blindly.'
            error stop 'test_M_convergence: basis not converged'
        end if

        print *, 'PASS: res_L2 decreasing and res_L2(finest) < 1.0e-2'
        print *, '=== test_M_convergence PASSED ==='
    end subroutine test_M_convergence

    !> Phase-3 Task 4: PERIODIZATION-DEFORMATION error bar (design section 5.1).
    !>
    !> This is the payoff of Phase 3: it quantifies the PHYSICAL error introduced
    !> by periodizing the background, distinct from the numerical convergence
    !> pinned by Tasks 2-3. As the as-is half-width dx_asis grows, LESS of the
    !> background is deformed near the resonant layer, so the resonant-layer
    !> delta-Phi approaches the true (undeformed) solution. The Cauchy residual
    !> between successive dx_asis is the periodization-deformation error bar.
    !>
    !> CRITICAL (design section 4.1): the NUMERICAL resolution is held constant
    !> across the scan so the residual measures PHYSICAL deformation, not
    !> numerical under-resolution. Hold dx_tr fixed at 10 rho_L so dx_asis is
    !> the only physical knob. As dx_asis grows, L = 2*(dx_asis + dx_tr) grows,
    !> so M and n_rg are SCALED with L:
    !>   M    = ceiling( (27/rhoL) * L / (2 pi) ) -> fixed converged k_max
    !>   n_rg = max(ceiling(16 L/rhoL), 4 M)      -> converged quadrature/Nyquist
    !> The cutoff follows the preceding M scan, where M=128 on the 15-rho_L
    !> half-window reduced the Cauchy residual below 1e-2. This scan therefore
    !> measures deformation rather than the known under-resolution at k_max=5/rho_L.
    !>
    !> SOFT GATE (this REPORTS the error bar, it does NOT tightly pass/fail): the
    !> test error-stops ONLY on a genuine defect -- non-finite / all-zero dPhi, a
    !> non-finite residual, or a residual sequence that GROWS without bound (a
    !> growing residual means the deformation is NOT converging: a real finding).
    !> A decreasing / plateauing residual is the physical expectation; the plateau
    !> value (residual at the largest dx_asis) is the reported deformation error
    !> bar. A small-but-finite plateau does NOT hard-fail.
    subroutine test_dr_deformation_scan()
        use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
        use config_m, only: profiles_in_memory, nml_config_path
        use species_m, only: set_profiles_from_arrays
        use kim_base_m, only: kim_t
        use kim_mod_m, only: from_kim_factory_get_kim
        use kim_resonances_m, only: r_res
        use rt_electrostatic_periodic_m, only: compute_periodic_delta_phi

        integer, parameter :: npts = 201
        integer, parameter :: ndiag = 41, nseq = 4
        real(dp), parameter :: dx_asis_scale(nseq) = &
            [3.0_dp, 5.0_dp, 8.0_dp, 12.0_dp]

        real(dp) :: r_prof(npts), n_prof(npts), Te_prof(npts)
        real(dp) :: Ti_prof(npts), q_prof(npts), Er_prof(npts)
        class(kim_t), allocatable :: kim_instance

        complex(dp), allocatable :: dPhi(:)
        complex(dp) :: dphi_k(ndiag, nseq)
        real(dp) :: r_diag(ndiag)
        real(dp) :: res_L2(nseq), res_max(nseq)
        real(dp) :: rm, rhoL_rm, dx_asis, dx_tr, L, k_max, pi
        real(dp) :: denom_L2, denom_max
        integer :: M_seq(nseq), n_rg_seq(nseq)
        integer :: i, k, info
        logical :: growing

        pi = acos(-1.0_dp)

        print *, ''
        print *, '=== test_dr_deformation_scan: PERIODIZATION-DEFORMATION error bar ==='

        call make_test_profiles(npts, r_prof, n_prof, Te_prof, Ti_prof, &
                                q_prof, Er_prof)

        call write_test_namelist('./KIM_config_periodic_conv_test.nml')
        nml_config_path = './KIM_config_periodic_conv_test.nml'

        profiles_in_memory = .true.
        call kim_init()
        call set_profiles_from_arrays(r_prof, n_prof, Te_prof, Ti_prof, &
                                      q_prof, Er_prof, npts)

        call from_kim_factory_get_kim('electrostatic', kim_instance)
        call kim_instance%init()

        call prepare_resonances
        if (.not. (r_res > 0.0_dp)) then
            print *, 'FAIL: resonance not found, r_res = ', r_res
            error stop
        end if
        rm = r_res

        rhoL_rm = interp_rho_L(rm)
        if (.not. (rhoL_rm > 0.0_dp)) then
            print *, 'FAIL: rho_L(rm) not positive, = ', rhoL_rm
            error stop
        end if
        print *, 'resonant radius rm       = ', rm
        print *, 'rho_L(rm) [active-species reference] = ', rhoL_rm

        ! FIXED diagnostic grid: 41 equidistant points on [rm - 2 rho_L,
        ! rm + 2 rho_L]. The SAME physical grid is used for EVERY dx_asis, so the
        ! Cauchy residuals are directly comparable. 2 rho_L is safely INSIDE the
        ! as-is region for ALL dx_asis in the scan (the smallest dx_asis =
        ! 3 rho_L > 2 rho_L), so the diagnostic layer always sees the undeformed
        ! background and the residual isolates the periodization deformation.
        do i = 1, ndiag
            r_diag(i) = (rm - 2.0_dp * rhoL_rm) &
                      + 4.0_dp * rhoL_rm * real(i - 1, dp) / real(ndiag - 1, dp)
        end do

        ! Scan ONLY the as-is half-width. Keep dx_tr fixed at the run-type's
        ! configured 10-rho_L transition scale, and hold numerical resolution
        ! constant by scaling M and n_rg with L. The fixed k_max = 27/rho_L is
        ! justified by test_M_convergence above.
        do k = 1, nseq
            dx_asis = dx_asis_scale(k) * rhoL_rm
            dx_tr   = 10.0_dp * rhoL_rm
            L       = 2.0_dp * (dx_asis + dx_tr)
            k_max   = 27.0_dp / rhoL_rm
            M_seq(k)    = ceiling(k_max * L / (2.0_dp * pi))
            n_rg_seq(k) = max(ceiling(16.0_dp * L / rhoL_rm), 4 * M_seq(k))

            print '(A,F6.1,A,ES14.6,A,I5,A,I6)', &
                '   dx_asis = ', dx_asis_scale(k), ' rho_L  (=', dx_asis, &
                ')  ->  M = ', M_seq(k), '  n_rg = ', n_rg_seq(k)

            call compute_periodic_delta_phi(rm, dx_asis, dx_tr, M_seq(k), &
                                            n_rg_seq(k), (1.0_dp, 0.0_dp), &
                                            r_diag, dPhi, info)
            if (info /= 0) then
                print *, 'FAIL: compute_periodic_delta_phi info /= 0 at dx_asis =', &
                    dx_asis_scale(k), ' rho_L, info =', info
                error stop
            end if
            if (size(dPhi) /= ndiag) then
                print *, 'FAIL: size(dPhi) /=', ndiag, ' got ', size(dPhi)
                error stop
            end if
            do i = 1, ndiag
                if (.not. ieee_is_finite(real(dPhi(i), dp)) .or. &
                    .not. ieee_is_finite(aimag(dPhi(i)))) then
                    print *, 'FAIL: non-finite dPhi at i =', i, &
                        ' dx_asis =', dx_asis_scale(k), ' rho_L'
                    error stop 'test_dr_deformation_scan: non-finite dPhi'
                end if
            end do
            if (.not. (maxval(abs(dPhi)) > 0.0_dp)) then
                print *, 'FAIL: dPhi all zero at dx_asis =', dx_asis_scale(k), ' rho_L'
                error stop 'test_dr_deformation_scan: all-zero dPhi'
            end if
            dphi_k(:, k) = dPhi
        end do
        print *, 'PASS: all solves returned info == 0 with finite, non-zero dPhi'

        ! Cauchy relative residual between successive dx_asis (k = 2, 3, 4).
        res_L2  = 0.0_dp
        res_max = 0.0_dp
        do k = 2, nseq
            denom_L2  = sqrt(sum(abs(dphi_k(:, k))**2))
            denom_max = maxval(abs(dphi_k(:, k)))
            res_L2(k)  = sqrt(sum(abs(dphi_k(:, k) - dphi_k(:, k-1))**2)) / denom_L2
            res_max(k) = maxval(abs(dphi_k(:, k) - dphi_k(:, k-1))) / denom_max
            if (.not. ieee_is_finite(res_L2(k)) .or. &
                .not. ieee_is_finite(res_max(k))) then
                print *, 'FAIL: non-finite residual at k =', k, &
                    ' res_L2 =', res_L2(k), ' res_max =', res_max(k)
                error stop 'test_dr_deformation_scan: non-finite residual'
            end if
        end do

        print *, ''
        print *, '==================================================================='
        print *, '   PERIODIZATION-DEFORMATION ERROR BAR (design section 5.1)'
        print *, '   fixed diagnostic layer [rm-2 rho_L, rm+2 rho_L]; M, n_rg scaled'
        print *, '   with L so k_max and r_g spacing are constant across the scan'
        print *, '==================================================================='
        print '(A)', '   dx_asis pair [rho_L]      res_L2            res_max'
        do k = 2, nseq
            print '(3X,F5.1," ->",F5.1,6X,ES16.8,2X,ES16.8)', &
                dx_asis_scale(k-1), dx_asis_scale(k), res_L2(k), res_max(k)
        end do
        print *, '-------------------------------------------------------------------'
        print '(A,ES16.8)', &
            '   DEFORMATION ERROR BAR (plateau, res_L2 at largest dx_asis)  = ', &
            res_L2(nseq)
        print '(A,ES16.8)', &
            '   DEFORMATION ERROR BAR (plateau, res_max at largest dx_asis) = ', &
            res_max(nseq)
        print *, '==================================================================='

        ! SOFT GATE: only error-stop on a genuinely GROWING residual sequence (the
        ! deformation is NOT converging -- a finding to investigate). Non-finite /
        ! all-zero / non-finite-residual defects are already caught above. A
        ! decreasing / plateauing sequence is the physical expectation and is
        ! reported, NOT failed, however small the finite plateau.
        growing = (res_L2(3) > res_L2(2)) .and. (res_L2(4) > res_L2(3))
        if (growing) then
            print *, 'FINDING: res_L2 sequence GROWS across the dx_asis scan:'
            print *, '   res_L2 =', res_L2(2), res_L2(3), res_L2(4)
            print *, '   -> periodization deformation is NOT converging as the'
            print *, '      as-is window grows; this is a real finding to'
            print *, '      investigate, not a numerical-resolution artifact'
            print *, '      (M and n_rg were scaled with L to hold resolution).'
            error stop 'test_dr_deformation_scan: deformation residual growing'
        end if

        if (res_L2(nseq) <= res_L2(2)) then
            print *, 'PASS: deformation residual decreases / plateaus (converging);'
            print '(A,ES14.6)', &
                '       reported periodization-deformation error bar = ', res_L2(nseq)
        else
            print *, 'NOTE: non-monotone but non-growing residual; plateau reported'
            print '(A,ES14.6)', &
                '       reported periodization-deformation error bar = ', res_L2(nseq)
        end if
        print *, '=== test_dr_deformation_scan PASSED ==='
    end subroutine test_dr_deformation_scan

    !> Active-species reference scale selected by the production run type.
    real(dp) function interp_rho_L(x) result(rhoLx)
        use config_m, only: turn_off_electrons, turn_off_ions
        use species_m, only: plasma
        use rt_electrostatic_periodic_m, only: select_periodic_reference_scale, &
                                               PERIODIC_SCALE_OK
        real(dp), intent(in) :: x
        integer :: reference_species, info

        call select_periodic_reference_scale(plasma, x, .not. turn_off_electrons, &
                                             .not. turn_off_ions, rhoLx, &
                                             reference_species, info)
        if (info /= PERIODIC_SCALE_OK) then
            print *, 'FAIL: could not select convergence reference scale, info =', info
            error stop
        end if
    end function interp_rho_L

    subroutine make_test_profiles(npts, r_prof, n_prof, Te_prof, Ti_prof, &
                                  q_prof, Er_prof)
        ! Analytic profiles on r in [2, 60] cm; q crosses 3 inside [10, 50] so
        ! the (m,n) = (-6, 2) mode (|m/n| = 3) has a resonance there.
        integer, intent(in) :: npts
        real(dp), intent(out) :: r_prof(npts), n_prof(npts), Te_prof(npts)
        real(dp), intent(out) :: Ti_prof(npts), q_prof(npts), Er_prof(npts)
        integer :: i

        do i = 1, npts
            r_prof(i) = 2.0_dp + 58.0_dp * real(i - 1, dp) / real(npts - 1, dp)
            n_prof(i) = 2.0e13_dp * (1.1_dp - r_prof(i) / 100.0_dp)
            Te_prof(i) = 1.0e3_dp * (1.2_dp - r_prof(i) / 100.0_dp)
            Ti_prof(i) = Te_prof(i)
            q_prof(i) = 1.5_dp + 2.5_dp * (r_prof(i) - 2.0_dp) / 58.0_dp
            Er_prof(i) = -0.5_dp
        end do
    end subroutine make_test_profiles

    subroutine write_test_namelist(path)
        ! Minimal electrostatic FokkerPlanck configuration; m_mode = -6,
        ! n_mode = 2 makes q resonant at q = 3, type_br_field = 12 (constant Br).
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
        write(iunit, '(A)') " output_path = './out_periodic_conv_test/'"
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
        write(iunit, '(A)') ' m_mode = -6'
        write(iunit, '(A)') ' n_mode = 2'
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
    end subroutine write_test_namelist

end program test_periodic_convergence
