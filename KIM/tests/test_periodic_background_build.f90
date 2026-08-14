program test_periodic_background_build
    ! Validates build_periodic_plasma (Phase-2.3, THE CRUX): it must produce a
    ! fully-populated PERIODIC `plasma` on an equidistant window grid
    !   [rm - L/2, rm + L/2],  L = 2*(dx_asis + dx_tr),  n_rg boundary points,
    ! whose DERIVED quantities (rho_L, lambda_D, and the susceptibility
    ! functions I00..I21) are recomputed FROM the periodized primitives
    ! (n, Te, Ti, q, Er) and are therefore periodic.
    !
    ! Setup mirrors the in-memory-profile path of test_periodic_background.f90:
    !   write_test_namelist -> nml_config_path -> profiles_in_memory ->
    !   kim_init -> set_profiles_from_arrays -> electrostatic %init().
    ! The (m,n) = (-6, 2) mode makes q resonant at q = |m/n| = 3, so the
    ! q-profile is built to cross 3 inside the solver domain [r_min, r_plas].
    !
    ! The fixture contains D plus a Z=6 impurity with distinct density and
    ! temperature profiles. Assertions on the window:
    !   (1) rho_L, lambda_D finite and > 0 everywhere,
    !   (2) I00(j,0) finite everywhere (susceptibility recomputed),
    !   (3) as-is fidelity: inside the as-is core the BUILT derived quantities
    !       (rho_L, lambda_D) reproduce the TRUE global derived plasma there
    !       (the periodized primitives equal the true primitives in the core),
    !   (4) resonance: parallel wavenumber kp ~ 0 at the window point nearest rm
    !       (confirms q = 3 there and a correct kp recompute),
    !   (5) each ion profile survives independently and quasineutrality holds,
    !   (6) all finite (no NaN/Inf).

    use KIM_kinds_m, only: dp
    use periodic_background_m, only: build_periodic_plasma, &
                                     sample_periodic_species_profiles, &
                                     PERIODIC_BACKGROUND_NONPHYSICAL_PROFILE

    implicit none

    call test_build_periodic()

    print *, 'All tests PASSED'
    stop 0

contains

    subroutine test_build_periodic()
        use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
        use config_m, only: profiles_in_memory, nml_config_path
        use species_m, only: plasma, set_profiles_from_arrays
        use kim_base_m, only: kim_t
        use kim_mod_m, only: from_kim_factory_get_kim
        use kim_resonances_m, only: r_res
        use grid_m, only: rg_grid

        integer, parameter :: npts = 201

        real(dp) :: r_prof(npts), n_prof(npts), Te_prof(npts)
        real(dp) :: Ti_prof(npts), q_prof(npts), Er_prof(npts)
        class(kim_t), allocatable :: kim_instance

        real(dp) :: rm, dx_asis, dx_tr, rho_L_rm
        real(dp) :: r_asis, rho_L_true, lambda_D_true, rho_L_built, lambda_D_built
        real(dp) :: n_D_true, n_C_true, T_D_true, T_C_true
        real(dp) :: n_D_built, n_C_built, T_D_built, T_C_built
        real(dp) :: qn_residual, density_scale
        integer :: n_rg
        integer :: j, gs, j_rm, sp
        real(dp) :: kp_scale, tol_kp, tol_asis, rel
        logical :: any_bad

        call make_test_profiles(npts, r_prof, n_prof, Te_prof, Ti_prof, &
                                q_prof, Er_prof)

        call write_test_namelist('./KIM_config_periodic_build_test.nml')
        nml_config_path = './KIM_config_periodic_build_test.nml'

        profiles_in_memory = .true.
        call kim_init()
        call set_profiles_from_arrays(r_prof, n_prof, Te_prof, Ti_prof, &
                                      q_prof, Er_prof, npts)
        call set_distinct_two_ion_profiles()

        call from_kim_factory_get_kim('electrostatic', kim_instance)
        call kim_instance%init()

        ! Locate the resonance (q = |m/n| = 3) on the global plasma.
        call kim_prepare_resonances
        if (.not. (r_res > 0.0_dp)) then
            print *, 'FAIL: resonance not found, r_res = ', r_res
            error stop
        end if
        rm = r_res
        print *, 'resonant radius rm = ', rm

        call test_multispecies_cache_refresh(rm)

        ! Representative Larmor radius near rm (electrons) on the global plasma.
        rho_L_rm = rho_L_near(rm)
        print *, 'rho_L(rm) [electron] = ', rho_L_rm

        ! Window half-widths (task: 5 and 10 Larmor radii). Guard against a
        ! degenerate tiny window if rho_L is unexpectedly small.
        dx_asis =  5.0_dp * rho_L_rm
        dx_tr   = 10.0_dp * rho_L_rm
        if (dx_asis < 1.0e-3_dp) then
            dx_asis = 1.0_dp
            dx_tr   = 2.0_dp
            print *, 'NOTE: rho_L(rm) tiny; using fixed dx_asis/dx_tr [cm]'
        end if
        n_rg = 96

        ! Capture the TRUE (non-periodic) resonant-layer physics BEFORE the
        ! build. At this point the global `plasma` is the true derived plasma on
        ! the global rg_grid. Pick r_asis inside the as-is core (|r - rm| <=
        ! dx_asis), where the periodized primitives equal the true primitives so
        ! the derived quantities of the periodic build MUST reproduce the true
        ! global build. Interpolate the true rho_L / lambda_D at r_asis with the
        ! same 4-point Lagrange pattern KIM uses (binsrc + plag_coeff on
        ! plasma%r_grid). These are just scalars, safe across the mutation.
        r_asis = rm + 0.5_dp*dx_asis
        rho_L_true    = interp_lagrange4(plasma%r_grid, plasma%spec(0)%rho_L, &
                                         plasma%grid_size, r_asis)
        lambda_D_true = interp_lagrange4(plasma%r_grid, plasma%spec(0)%lambda_D, &
                                         plasma%grid_size, r_asis)
        n_D_true = interp_lagrange4(plasma%r_grid, plasma%spec(1)%n, &
                                    plasma%grid_size, r_asis)
        n_C_true = interp_lagrange4(plasma%r_grid, plasma%spec(2)%n, &
                                    plasma%grid_size, r_asis)
        T_D_true = interp_lagrange4(plasma%r_grid, plasma%spec(1)%T, &
                                    plasma%grid_size, r_asis)
        T_C_true = interp_lagrange4(plasma%r_grid, plasma%spec(2)%T, &
                                    plasma%grid_size, r_asis)

        call build_periodic_plasma(rm, dx_asis, dx_tr, n_rg)

        gs = plasma%grid_size
        if (gs /= n_rg) then
            print *, 'FAIL: grid_size /= n_rg after build; got ', gs
            error stop
        end if
        if (rg_grid%npts_b /= n_rg) then
            print *, 'FAIL: rg_grid%npts_b /= n_rg; got ', rg_grid%npts_b
            error stop
        end if

        ! (1) rho_L, lambda_D finite and strictly positive everywhere.
        any_bad = .false.
        do j = 1, gs
            if (.not. ieee_is_finite(plasma%spec(0)%rho_L(j)) .or. &
                .not. (plasma%spec(0)%rho_L(j) > 0.0_dp)) then
                print *, 'FAIL: rho_L not finite/positive at j=', j, &
                         ' val=', plasma%spec(0)%rho_L(j)
                any_bad = .true.
            end if
            if (.not. ieee_is_finite(plasma%spec(0)%lambda_D(j)) .or. &
                .not. (plasma%spec(0)%lambda_D(j) > 0.0_dp)) then
                print *, 'FAIL: lambda_D not finite/positive at j=', j, &
                         ' val=', plasma%spec(0)%lambda_D(j)
                any_bad = .true.
            end if
        end do
        if (any_bad) error stop

        ! (2) Susceptibility I00(:,0) recomputed and finite everywhere.
        do j = 1, gs
            if (.not. ieee_is_finite(real(plasma%spec(0)%I00(j, 0), dp)) .or. &
                .not. ieee_is_finite(aimag(plasma%spec(0)%I00(j, 0)))) then
                print *, 'FAIL: I00 non-finite at j=', j, ' val=', plasma%spec(0)%I00(j, 0)
                error stop
            end if
        end do
        print *, 'PASS: rho_L, lambda_D finite/positive; I00 finite everywhere'

        ! (3) As-is fidelity: the periodic build must reproduce the TRUE
        ! resonant-layer physics inside the as-is core. In |r - rm| <= dx_asis
        ! the periodized primitives equal the true primitives, so the DERIVED
        ! quantities recomputed by build_periodic_plasma must match the true
        ! global derived plasma there. We interpolate the BUILT rho_L / lambda_D
        ! at the SAME r_asis (now on the window grid) and compare to the true
        ! values captured before the build. The modest 1e-3 relative tolerance
        ! covers only the resolution difference between the global grid and the
        ! finer window grid -- it cannot pass trivially: a broken recompute
        ! (wrong B0/omega_c, wrong n/T redirection) would miss by orders of
        ! magnitude. This is the load-bearing correctness signal.
        rho_L_built    = interp_lagrange4(rg_grid%xb, plasma%spec(0)%rho_L, &
                                          gs, r_asis)
        lambda_D_built = interp_lagrange4(rg_grid%xb, plasma%spec(0)%lambda_D, &
                                          gs, r_asis)
        tol_asis = 1.0e-3_dp

        rel = abs(rho_L_built - rho_L_true) / max(abs(rho_L_true), tiny(1.0_dp))
        print *, 'as-is r_asis =', r_asis
        print *, '  rho_L  true =', rho_L_true, ' built =', rho_L_built, ' rel =', rel
        if (rel >= tol_asis) then
            print *, 'FAIL: built rho_L does not reproduce true value in as-is core'
            error stop
        end if

        rel = abs(lambda_D_built - lambda_D_true) / max(abs(lambda_D_true), tiny(1.0_dp))
        print *, '  lambda_D true =', lambda_D_true, ' built =', lambda_D_built, ' rel =', rel
        if (rel >= tol_asis) then
            print *, 'FAIL: built lambda_D does not reproduce true value in as-is core'
            error stop
        end if
        print *, 'PASS: periodic build reproduces true rho_L, lambda_D in as-is core'

        ! (5) Preserve every configured species independently. Species 1 is the
        ! dependent main ion, while the impurity density is preserved and the
        ! main-ion density enforces n_e = Z_D*n_D + Z_C*n_C. The fixture is
        ! exactly quasineutral before periodization, so all four profiles must
        ! reproduce the true as-is values to interpolation accuracy.
        n_D_built = interp_lagrange4(rg_grid%xb, plasma%spec(1)%n, gs, r_asis)
        n_C_built = interp_lagrange4(rg_grid%xb, plasma%spec(2)%n, gs, r_asis)
        T_D_built = interp_lagrange4(rg_grid%xb, plasma%spec(1)%T, gs, r_asis)
        T_C_built = interp_lagrange4(rg_grid%xb, plasma%spec(2)%T, gs, r_asis)

        call assert_relative_close('D density', n_D_built, n_D_true, tol_asis)
        call assert_relative_close('C density', n_C_built, n_C_true, tol_asis)
        call assert_relative_close('D temperature', T_D_built, T_D_true, tol_asis)
        call assert_relative_close('C temperature', T_C_built, T_C_true, tol_asis)

        density_scale = maxval(plasma%spec(0)%n)
        qn_residual = maxval(abs(plasma%spec(1)%n + &
                                 6.0_dp*plasma%spec(2)%n - plasma%spec(0)%n))
        if (qn_residual > 50.0_dp*epsilon(1.0_dp)*density_scale) then
            print *, 'FAIL: periodic multispecies background is not quasineutral'
            print *, '  max residual = ', qn_residual, ' density scale = ', density_scale
            error stop
        end if
        do j = 1, gs
            if (any([plasma%spec(0)%n(j), plasma%spec(1)%n(j), &
                     plasma%spec(2)%n(j)] <= 0.0_dp) .or. &
                any([plasma%spec(0)%T(j), plasma%spec(1)%T(j), &
                     plasma%spec(2)%T(j)] <= 0.0_dp)) then
                print *, 'FAIL: non-positive multispecies primitive at j=', j
                error stop
            end if
        end do
        print *, 'PASS: distinct two-ion profiles preserved and quasineutral'

        ! (4) Resonance: kp ~ 0 at the window point nearest rm (q = 3 there).
        j_rm = minloc(abs(rg_grid%xb - rm), dim=1)
        ! Scale for "small": kp away from resonance is O(n_mode/R0) ~ O(1e-2).
        kp_scale = maxval(abs(plasma%kp))
        tol_kp = 1.0e-3_dp * max(kp_scale, tiny(1.0_dp))
        print *, 'kp scale (max|kp|) = ', kp_scale
        print *, 'kp at j_rm=', j_rm, ' (r=', rg_grid%xb(j_rm), ') = ', plasma%kp(j_rm)
        if (abs(plasma%kp(j_rm)) >= tol_kp) then
            print *, 'FAIL: kp not ~0 at rm'
            print *, '  |kp(j_rm)| = ', abs(plasma%kp(j_rm)), ' tol = ', tol_kp
            error stop
        end if
        print *, 'PASS: kp ~ 0 at window point nearest rm (resonance q=3 confirmed)'

        ! (6) Global finiteness sweep of the load-bearing derived arrays.
        do j = 1, gs
            if (.not. ieee_is_finite(plasma%kp(j)) .or. &
                .not. ieee_is_finite(plasma%B0(j)) .or. &
                .not. ieee_is_finite(plasma%om_E(j)) .or. &
                .not. ieee_is_finite(plasma%spec(0)%vT(j)) .or. &
                .not. ieee_is_finite(plasma%spec(0)%nu(j))) then
                print *, 'FAIL: non-finite derived quantity at j=', j
                error stop
            end if
            if (.not. ieee_is_finite(plasma%spec(1)%vT(j)) .or. &
                .not. ieee_is_finite(plasma%spec(1)%nu(j)) .or. &
                .not. ieee_is_finite(plasma%spec(1)%rho_L(j)) .or. &
                .not. ieee_is_finite(plasma%spec(2)%vT(j)) .or. &
                .not. ieee_is_finite(plasma%spec(2)%nu(j)) .or. &
                .not. ieee_is_finite(plasma%spec(2)%rho_L(j))) then
                print *, 'FAIL: non-finite ion-derived quantity at j=', j
                error stop
            end if
            do sp = 0, 2
                if (.not. ieee_is_finite(plasma%spec(sp)%omega_c(j)) .or. &
                    .not. ieee_is_finite(plasma%spec(sp)%lambda_D(j)) .or. &
                    .not. ieee_is_finite(plasma%spec(sp)%A1(j)) .or. &
                    .not. ieee_is_finite(plasma%spec(sp)%A2(j)) .or. &
                    .not. ieee_is_finite(plasma%spec(sp)%x1(j)) .or. &
                    .not. ieee_is_finite(plasma%spec(sp)%x2(j, 0)) .or. &
                    .not. ieee_is_finite(real(plasma%spec(sp)%I00(j, 0), dp)) .or. &
                    .not. ieee_is_finite(aimag(plasma%spec(sp)%I00(j, 0)))) then
                    print *, 'FAIL: incomplete derived pipeline at j=', j, ' species=', sp
                    error stop
                end if
            end do
        end do
        print *, 'PASS: all derived quantities finite across the window'

        call test_multispecies_periodicity(rm, dx_asis, dx_tr)
        call test_negative_dependent_density_rejected(rm, dx_asis, dx_tr, n_rg)
    end subroutine test_build_periodic

    subroutine test_multispecies_periodicity(rm, dx_asis, dx_tr)
        real(dp), intent(in) :: rm, dx_asis, dx_tr
        real(dp) :: n_left(0:2), T_left(0:2), n_shift(0:2), T_shift(0:2)
        real(dp) :: q_left, Er_left, q_shift, Er_shift
        real(dp) :: x_values(5), x, period_length
        integer :: k, sp

        period_length = 2.0_dp*(dx_asis + dx_tr)
        x_values = [rm - dx_asis - dx_tr, rm - dx_asis, rm, &
                    rm + dx_asis, rm + dx_asis + dx_tr]
        do k = 1, size(x_values)
            x = x_values(k)
            call sample_periodic_species_profiles(rm, dx_asis, dx_tr, x, &
                                                  n_left, T_left, q_left, Er_left)
            call sample_periodic_species_profiles(rm, dx_asis, dx_tr, &
                x + period_length, n_shift, T_shift, q_shift, Er_shift)
            do sp = 0, 2
                call assert_relative_close('periodic species density', &
                                           n_shift(sp), n_left(sp), 1.0e-10_dp)
                call assert_relative_close('periodic species temperature', &
                                           T_shift(sp), T_left(sp), 1.0e-10_dp)
            end do
            call assert_relative_close('periodic q', q_shift, q_left, 1.0e-10_dp)
            call assert_relative_close('periodic Er', Er_shift, Er_left, 1.0e-10_dp)
        end do
        print *, 'PASS: all species primitives are periodic across the seam'
    end subroutine test_multispecies_periodicity

    subroutine test_negative_dependent_density_rejected(rm, dx_asis, dx_tr, n_rg)
        use species_m, only: plasma, profiles_generation
        use grid_m, only: rg_grid

        real(dp), intent(in) :: rm, dx_asis, dx_tr
        integer, intent(in) :: n_rg
        real(dp) :: old_grid_start, old_grid_end
        integer :: info, old_n_rg

        ! Z_C*n_C = 1.2*n_e leaves no positive density for the dependent
        ! deuterium species. The builder must report a profile error rather than
        ! advancing into logarithms, collision frequencies, or susceptibilities.
        old_n_rg = rg_grid%npts_b
        old_grid_start = rg_grid%xb(1)
        old_grid_end = rg_grid%xb(old_n_rg)
        plasma%spec(2)%n = 0.2_dp*plasma%spec(0)%n
        profiles_generation = profiles_generation + 1
        call build_periodic_plasma(rm, dx_asis, dx_tr, n_rg + 7, info)
        if (info /= PERIODIC_BACKGROUND_NONPHYSICAL_PROFILE) then
            print *, 'FAIL: negative dependent density status = ', info
            error stop
        end if
        if (rg_grid%npts_b /= old_n_rg .or. &
            abs(rg_grid%xb(1) - old_grid_start) > epsilon(1.0_dp) .or. &
            abs(rg_grid%xb(old_n_rg) - old_grid_end) > epsilon(1.0_dp)) then
            print *, 'FAIL: rejected profile mutated the global radial grid'
            error stop
        end if
        print *, 'PASS: negative dependent-ion density rejected as a profile error'
    end subroutine test_negative_dependent_density_rejected

    subroutine test_multispecies_cache_refresh(rm)
        ! Model two successive in-process profile cases. The generation bump is
        ! the same invalidation signal emitted by the profile adapter; case B
        ! keeps total ion pressure unchanged so the already-computed geometry
        ! remains a valid true background for the rest of this test.
        use species_m, only: plasma, profiles_generation

        real(dp), intent(in) :: rm
        real(dp) :: n_before(0:2), T_before(0:2), n_after(0:2), T_after(0:2)
        real(dp) :: q_sample, Er_sample, radius_fraction, impurity_fraction
        real(dp) :: ion_pressure, expected
        integer :: j

        call sample_periodic_species_profiles(rm, 1.0_dp, 2.0_dp, rm, &
                                              n_before, T_before, q_sample, Er_sample)

        do j = 1, plasma%grid_size
            radius_fraction = (plasma%r_grid(j) - plasma%r_grid(1)) / &
                              (plasma%r_grid(plasma%grid_size) - plasma%r_grid(1))
            ion_pressure = plasma%spec(1)%n(j)*plasma%spec(1)%T(j) + &
                           plasma%spec(2)%n(j)*plasma%spec(2)%T(j)
            impurity_fraction = 0.018_dp + 0.003_dp*radius_fraction
            plasma%spec(2)%n(j) = impurity_fraction*plasma%spec(0)%n(j)
            plasma%spec(1)%n(j) = plasma%spec(0)%n(j) - &
                                  6.0_dp*plasma%spec(2)%n(j)
            plasma%spec(2)%T(j) = 1050.0_dp + 100.0_dp*radius_fraction
            plasma%spec(1)%T(j) = (ion_pressure - &
                                    plasma%spec(2)%n(j)*plasma%spec(2)%T(j)) / &
                                   plasma%spec(1)%n(j)
        end do
        profiles_generation = profiles_generation + 1

        call sample_periodic_species_profiles(rm, 1.0_dp, 2.0_dp, rm, &
                                              n_after, T_after, q_sample, Er_sample)

        expected = interp_lagrange4(plasma%r_grid, plasma%spec(2)%n, &
                                    plasma%grid_size, rm)
        call assert_relative_close('refreshed C density', n_after(2), expected, 1.0e-10_dp)
        expected = interp_lagrange4(plasma%r_grid, plasma%spec(2)%T, &
                                    plasma%grid_size, rm)
        call assert_relative_close('refreshed C temperature', T_after(2), expected, 1.0e-10_dp)
        if (abs(n_after(2) - n_before(2)) <= 0.1_dp*abs(n_before(2)) .or. &
            abs(T_after(2) - T_before(2)) <= 0.1_dp*abs(T_before(2))) then
            print *, 'FAIL: multispecies cache retained the first profile case'
            error stop
        end if
        print *, 'PASS: multispecies cache refreshed for a second profile case'
    end subroutine test_multispecies_cache_refresh

    subroutine set_distinct_two_ion_profiles()
        ! Construct an exactly quasineutral D + C background. The carbon
        ! concentration and both ion temperatures have distinct radial shapes,
        ! so a single shared ion profile cannot satisfy this fixture.
        use species_m, only: plasma

        integer :: j
        real(dp) :: impurity_fraction, radius_fraction

        do j = 1, plasma%grid_size
            radius_fraction = (plasma%r_grid(j) - plasma%r_grid(1)) / &
                              (plasma%r_grid(plasma%grid_size) - plasma%r_grid(1))
            impurity_fraction = 0.010_dp + 0.005_dp*radius_fraction
            plasma%spec(2)%n(j) = impurity_fraction*plasma%spec(0)%n(j)
            plasma%spec(1)%n(j) = plasma%spec(0)%n(j) - &
                                  6.0_dp*plasma%spec(2)%n(j)
            plasma%spec(1)%T(j) = 700.0_dp + 80.0_dp*radius_fraction
            plasma%spec(2)%T(j) = 1400.0_dp - 120.0_dp*radius_fraction
        end do
    end subroutine set_distinct_two_ion_profiles

    subroutine assert_relative_close(label, actual, expected, tolerance)
        character(len=*), intent(in) :: label
        real(dp), intent(in) :: actual, expected, tolerance
        real(dp) :: relative_error

        relative_error = abs(actual - expected) / max(abs(expected), tiny(1.0_dp))
        print *, trim(label), ': true =', expected, ' built =', actual, &
                 ' rel =', relative_error
        if (relative_error >= tolerance) then
            print *, 'FAIL: ', trim(label), ' not preserved in as-is core'
            error stop
        end if
    end subroutine assert_relative_close

    real(dp) function rho_L_near(r) result(val)
        !> Return the electron Larmor radius of the global plasma at the grid
        !> point nearest radius r (the plasma is still on the global grid here).
        use species_m, only: plasma
        real(dp), intent(in) :: r
        integer :: jn
        jn = minloc(abs(plasma%r_grid - r), dim=1)
        val = plasma%spec(0)%rho_L(jn)
    end function rho_L_near

    real(dp) function interp_lagrange4(grid, vals, ngrid, x) result(fx)
        !> 4-point Lagrange interpolation of vals(:) on grid(:) at radius x.
        !> Mirrors the binsrc + plag_coeff stencil used throughout KIM
        !> (kim_aperfuns, interpolate_plasma_backs, calculate_equil).
        integer, intent(in) :: ngrid
        real(dp), intent(in) :: grid(ngrid), vals(ngrid), x
        integer, parameter :: nlagr = 4, nder = 0
        real(dp) :: coef(0:nder, nlagr)
        integer :: ir, ibeg, iend

        call binsrc(grid, 1, ngrid, x, ir)
        ibeg = max(1, ir - nlagr/2)
        iend = ibeg + nlagr - 1
        if (iend > ngrid) then
            iend = ngrid
            ibeg = iend - nlagr + 1
        end if
        call plag_coeff(nlagr, nder, x, grid(ibeg:iend), coef)
        fx = sum(coef(0, :) * vals(ibeg:iend))
    end function interp_lagrange4

    subroutine make_test_profiles(npts, r_prof, n_prof, Te_prof, Ti_prof, &
                                  q_prof, Er_prof)
        ! Analytic profiles on r in [2, 60] cm. q rises monotonically through 3
        ! inside the solver domain [r_min, r_plas] = [10, 50] so that the
        ! (m,n) = (-6, 2) mode (|m/n| = 3) has a resonance there.
        integer, intent(in) :: npts
        real(dp), intent(out) :: r_prof(npts), n_prof(npts), Te_prof(npts)
        real(dp), intent(out) :: Ti_prof(npts), q_prof(npts), Er_prof(npts)
        integer :: i

        do i = 1, npts
            r_prof(i) = 2.0_dp + 58.0_dp * real(i - 1, dp) / real(npts - 1, dp)
            n_prof(i) = 2.0e13_dp * (1.1_dp - r_prof(i) / 100.0_dp)
            Te_prof(i) = 1.0e3_dp * (1.2_dp - r_prof(i) / 100.0_dp)
            Ti_prof(i) = Te_prof(i)
            ! q: 1.5 at r=2 to 4.0 at r=60  ->  crosses 3 near r=36 cm.
            q_prof(i) = 1.5_dp + 2.5_dp * (r_prof(i) - 2.0_dp) / 58.0_dp
            Er_prof(i) = -0.5_dp
        end do
    end subroutine make_test_profiles

    subroutine write_test_namelist(path)
        ! Minimal electrostatic FokkerPlanck configuration on a coarse grid.
        ! m_mode = -6, n_mode = 2 makes q resonant at q = 3.
        character(len=*), intent(in) :: path
        integer :: iunit

        open(newunit=iunit, file=path, status='replace', action='write')
        write(iunit, '(A)') '&KIM_CONFIG'
        write(iunit, '(A)') ' number_of_ion_species = 2'
        write(iunit, '(A)') ' artificial_debye_case = 0'
        write(iunit, '(A)') " type_of_run = 'electrostatic'"
        write(iunit, '(A)') " collision_model = 'FokkerPlanck'"
        write(iunit, '(A)') ' read_species_from_namelist = .true.'
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
        write(iunit, '(A)') " output_path = './out_periodic_build_test/'"
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
        write(iunit, '(A)') '&KIM_species'
        write(iunit, '(A)') ' ai = 2, 12'
        write(iunit, '(A)') ' zi = 1, 6'
        write(iunit, '(A)') '/'
        close(iunit)
    end subroutine write_test_namelist

end program test_periodic_background_build
