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

    call test_core_on_fixed_grid()
    call test_nrg_convergence()
    call test_M_convergence()

    print *, 'All tests PASSED'
    stop 0

contains

    subroutine test_core_on_fixed_grid()
        use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
        use config_m, only: profiles_in_memory, nml_config_path
        use species_m, only: plasma, set_profiles_from_arrays
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

        ! Representative electron Larmor radius at rm, via the SAME 4-point
        ! Lagrange stencil the run-type uses to size its window.
        rhoL_rm = interp_rho_L(rm)
        if (.not. (rhoL_rm > 0.0_dp)) then
            print *, 'FAIL: rho_L(rm) not positive, = ', rhoL_rm
            error stop
        end if
        print *, 'rho_L(rm) [electron] = ', rhoL_rm

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
        print *, 'rho_L(rm) [electron]     = ', rhoL_rm

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
    !> With the r_g quadrature (n_rg = 192, converged in P3.2) and the window
    !> (dx_asis, dx_tr) FIXED, refine only the number of Fourier modes M (the
    !> basis has 2M+1 modes) and reconstruct delta-Phi on a FIXED diagnostic
    !> grid. The localized solution's Fourier coefficients decay for
    !> |k_m| >> 1/rho_L, so the Cauchy residual between successive M must fall
    !> and beat 1e-2 by M ~ 24-32. A failure to converge is a real finding about
    !> the basis truncation or the k-resolution, not a reason to loosen the
    !> tolerance (bump n_rg first: the exp(i(k_m - k_m')r) quadrature is exact
    !> for |m - m'| < n_rg by Nyquist, and 2M = 64 < 192, so n_rg = 192 resolves
    !> M up to 32).
    subroutine test_M_convergence()
        use config_m, only: profiles_in_memory, nml_config_path
        use species_m, only: set_profiles_from_arrays
        use kim_base_m, only: kim_t
        use kim_mod_m, only: from_kim_factory_get_kim
        use kim_resonances_m, only: r_res
        use rt_electrostatic_periodic_m, only: compute_periodic_delta_phi

        integer, parameter :: npts = 201
        integer, parameter :: n_rg = 192, ndiag = 41, nseq = 4
        integer, parameter :: M_seq(nseq) = [8, 16, 24, 32]

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
        print *, 'rho_L(rm) [electron]     = ', rhoL_rm

        ! FIXED quadrature and window (only M is refined here). n_rg = 192 is the
        ! P3.2-converged node count; it resolves the k-modes for M up to 32
        ! (2M = 64 < 192, Nyquist).
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
        print *, '--- M convergence (fixed n_rg, window, diagnostic grid) ---'
        print '(A)', '   M pair             res_L2            res_max'
        do k = 2, nseq
            print '(3X,I5," ->",I5,4X,ES16.8,2X,ES16.8)', &
                M_seq(k-1), M_seq(k), res_L2(k), res_max(k)
        end do
        print *, '-----------------------------------------------------------'

        ! Assert: the residual sequence DECREASES and the finest beats 1e-2.
        decreasing = (res_L2(4) < res_L2(3)) .and. (res_L2(3) < res_L2(2))
        if (.not. decreasing) then
            print *, 'FAIL: res_L2 sequence does not decrease monotonically:'
            print *, '   res_L2 =', res_L2(2), res_L2(3), res_L2(4)
            print *, '   -> real finding about basis truncation / k-resolution;'
            print *, '      investigate (bump n_rg first), do NOT loosen the'
            print *, '      threshold blindly.'
            error stop 'test_M_convergence: residual not decreasing'
        end if
        if (.not. (res_L2(4) < 1.0e-2_dp)) then
            print *, 'FAIL: finest res_L2(4) =', res_L2(4), ' not < 1.0e-2'
            print *, '   res_L2 =', res_L2(2), res_L2(3), res_L2(4)
            print *, '   -> Fourier basis has NOT converged; investigate (bump'
            print *, '      n_rg first), do NOT loosen the threshold blindly.'
            error stop 'test_M_convergence: basis not converged'
        end if

        print *, 'PASS: res_L2 decreasing and res_L2(finest) < 1.0e-2'
        print *, '=== test_M_convergence PASSED ==='
    end subroutine test_M_convergence

    !> 4-point Lagrange interpolation of the global electron Larmor radius
    !> plasma%spec(0)%rho_L at radius x, matching the run-type's interp_rho_L.
    real(dp) function interp_rho_L(x) result(rhoLx)
        use species_m, only: plasma
        real(dp), intent(in) :: x
        integer, parameter :: nlagr = 4, nder = 0
        real(dp) :: coef(0:nder, nlagr)
        integer :: gs, ir, ibeg, iend

        gs = plasma%grid_size
        call binsrc(plasma%r_grid, 1, gs, x, ir)
        ibeg = max(1, ir - nlagr/2)
        iend = ibeg + nlagr - 1
        if (iend > gs) then
            iend = gs
            ibeg = iend - nlagr + 1
        end if
        call plag_coeff(nlagr, nder, x, plasma%r_grid(ibeg:iend), coef)
        rhoLx = sum(coef(0, :) * plasma%spec(0)%rho_L(ibeg:iend))
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
