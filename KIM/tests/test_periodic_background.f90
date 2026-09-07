program test_periodic_background
    ! Validates the KIM aperfuns callback and the periodic primitive sampler
    ! (periodic_background_m) against the vendored make_periodic utility.
    !
    ! Setup mirrors the in-memory-profile path of test_flr2_fourier_kernel.f90:
    !   write_test_namelist -> nml_config_path -> profiles_in_memory ->
    !   kim_init -> set_profiles_from_arrays -> electrostatic %init().
    ! After init the global `plasma` holds the TRUE primitive profiles
    ! (n_e, Te, Ti, q, Er) on plasma%r_grid, which kim_aperfuns interpolates.
    !
    ! Three properties are checked on a mid-domain period:
    !   (a) as-is fidelity: inside the as-is core the periodized value equals
    !       the direct interpolated (true) value to machine precision,
    !   (b) periodicity: sample(x) == sample(x + L),
    !   (c) finiteness + physicality: all components finite; n_e, Te, Ti > 0.

    use KIM_kinds_m, only: dp
    use periodic_background_m, only: kim_aperfuns, sample_periodic_primitives

    implicit none

    call test_periodic_sampling()
    call test_cache_invalidation()

    print *, 'All tests PASSED'
    stop 0

contains

    subroutine test_periodic_sampling()
        use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
        use config_m, only: profiles_in_memory, nml_config_path
        use species_m, only: plasma, set_profiles_from_arrays
        use kim_base_m, only: kim_t
        use kim_mod_m, only: from_kim_factory_get_kim

        integer, parameter :: npts = 101
        integer, parameter :: nfuns = 5

        real(dp) :: r_prof(npts), n_prof(npts), Te_prof(npts)
        real(dp) :: Ti_prof(npts), q_prof(npts), Er_prof(npts)
        class(kim_t), allocatable :: kim_instance

        real(dp) :: x_mid, dx_asis, dx_tr, period_len
        real(dp) :: x, x_arr(5)
        real(dp) :: funs_per(nfuns), funs_true(nfuns), funs_shift(nfuns)
        real(dp) :: r_min_grid, r_max_grid, tol
        integer :: i, k, gs

        call make_test_profiles(npts, r_prof, n_prof, Te_prof, Ti_prof, &
                                q_prof, Er_prof)

        call write_test_namelist('./KIM_config_periodic_bg_test.nml')
        nml_config_path = './KIM_config_periodic_bg_test.nml'

        profiles_in_memory = .true.
        call kim_init()
        call set_profiles_from_arrays(r_prof, n_prof, Te_prof, Ti_prof, &
                                      q_prof, Er_prof, npts)

        call from_kim_factory_get_kim('electrostatic', kim_instance)
        call kim_instance%init()

        gs = plasma%grid_size
        r_min_grid = plasma%r_grid(1)
        r_max_grid = plasma%r_grid(gs)

        ! Mid-domain period centre; small half-period so [x_mid +- L/2] stays
        ! well inside the profile domain. L = 2*(dx_asis + dx_tr) = 18 cm.
        x_mid = 0.5_dp*(r_min_grid + r_max_grid)
        dx_asis = 3.0_dp
        dx_tr = 6.0_dp
        period_len = 2.0_dp*(dx_asis + dx_tr)

        tol = 1.0e-10_dp

        ! (a) As-is fidelity: inside |x - x_mid| <= dx_asis the periodized value
        ! must equal the direct interpolated true value to machine precision.
        x_arr = [x_mid - dx_asis, x_mid - 0.5_dp*dx_asis, x_mid, &
                 x_mid + 0.5_dp*dx_asis, x_mid + dx_asis]
        do k = 1, size(x_arr)
            x = x_arr(k)
            call sample_periodic_primitives(x_mid, dx_asis, dx_tr, x, funs_per)
            call kim_aperfuns(nfuns, x, funs_true)
            do i = 1, nfuns
                if (abs(funs_per(i) - funs_true(i)) >= tol) then
                    print *, 'FAIL: as-is fidelity broken at x=', x, ' comp=', i
                    print *, '  periodized = ', funs_per(i)
                    print *, '  true       = ', funs_true(i)
                    print *, '  |diff|     = ', abs(funs_per(i) - funs_true(i))
                    error stop
                end if
            end do
            ! Physicality at as-is points: n_e, Te, Ti strictly positive.
            if (.not. (funs_per(1) > 0.0_dp)) then
                print *, 'FAIL: n_e not positive at x=', x, ' got ', funs_per(1)
                error stop
            end if
            if (.not. (funs_per(2) > 0.0_dp)) then
                print *, 'FAIL: Te not positive at x=', x, ' got ', funs_per(2)
                error stop
            end if
            if (.not. (funs_per(3) > 0.0_dp)) then
                print *, 'FAIL: Ti not positive at x=', x, ' got ', funs_per(3)
                error stop
            end if
        end do

        ! (b) Periodicity: sample(x) == sample(x + L) across the transition
        ! bands and the as-is core.
        x_arr = [x_mid - dx_asis - dx_tr, x_mid - dx_asis, x_mid, &
                 x_mid + dx_asis, x_mid + dx_asis + dx_tr]
        do k = 1, size(x_arr)
            x = x_arr(k)
            call sample_periodic_primitives(x_mid, dx_asis, dx_tr, x, funs_per)
            call sample_periodic_primitives(x_mid, dx_asis, dx_tr, &
                                            x + period_len, funs_shift)
            do i = 1, nfuns
                if (abs(funs_per(i) - funs_shift(i)) >= tol) then
                    print *, 'FAIL: periodicity broken at x=', x, ' comp=', i
                    print *, '  f(x)   = ', funs_per(i)
                    print *, '  f(x+L) = ', funs_shift(i)
                    print *, '  |diff| = ', abs(funs_per(i) - funs_shift(i))
                    error stop
                end if
            end do
        end do

        ! (c) Finiteness across the whole period (including transition bands).
        do k = 0, 20
            x = x_mid - 0.5_dp*period_len + period_len*real(k, dp)/20.0_dp
            call sample_periodic_primitives(x_mid, dx_asis, dx_tr, x, funs_per)
            do i = 1, nfuns
                if (.not. ieee_is_finite(funs_per(i))) then
                    print *, 'FAIL: non-finite component at x=', x, ' comp=', i
                    error stop
                end if
            end do
        end do

        print *, 'PASS: periodic sampling is as-is faithful, L-periodic, finite'
    end subroutine

    !> Prove the true-background cache AUTO-INVALIDATES when the profiles change.
    !> Two DIFFERENT cases in one process: init with profile set A (kim_aperfuns
    !> captures the cache at generation G), then swap in a DIFFERENT set B via
    !> set_profiles_from_arrays + re-init (which bumps
    !> species_m::profiles_generation to G+1 and rebuilds the true geometry),
    !> sample again, and assert the periodized as-is value now reflects B -- NOT
    !> the stale A snapshot. Without auto-invalidation kim_aperfuns would keep
    !> returning A's density. Mirrors the production QL-Balance / time-step path
    !> (set_profiles_from_arrays is always followed by an equilibrium re-init).
    subroutine test_cache_invalidation()
        use config_m, only: profiles_in_memory, nml_config_path
        use species_m, only: set_profiles_from_arrays, profiles_generation
        use kim_base_m, only: kim_t
        use kim_mod_m, only: from_kim_factory_get_kim

        integer, parameter :: npts = 101, nfuns = 5
        real(dp) :: r_prof(npts), n_prof(npts), Te_prof(npts)
        real(dp) :: Ti_prof(npts), q_prof(npts), Er_prof(npts)
        real(dp) :: x_core, funs_A(nfuns), funs_B(nfuns), n_expected_B
        class(kim_t), allocatable :: kim_instance
        integer :: gen_A, gen_B

        call write_test_namelist('./KIM_config_periodic_bg_test.nml')
        nml_config_path = './KIM_config_periodic_bg_test.nml'
        profiles_in_memory = .true.

        ! Case A: analytic profiles, init the equilibrium, sample n_e at a core
        ! radius. The direct kim_aperfuns callback captures the cache at gen A.
        call make_test_profiles(npts, r_prof, n_prof, Te_prof, Ti_prof, &
                                q_prof, Er_prof)
        call kim_init()
        call set_profiles_from_arrays(r_prof, n_prof, Te_prof, Ti_prof, &
                                      q_prof, Er_prof, npts)
        call from_kim_factory_get_kim('electrostatic', kim_instance)
        call kim_instance%init()
        gen_A = profiles_generation

        x_core = 30.0_dp
        call kim_aperfuns(nfuns, x_core, funs_A)

        ! Case B: DOUBLE the density everywhere (a clearly distinct case), swap it
        ! in and re-init the equilibrium -- this bumps profiles_generation.
        n_prof = 2.0_dp * n_prof
        call set_profiles_from_arrays(r_prof, n_prof, Te_prof, Ti_prof, &
                                      q_prof, Er_prof, npts)
        call kim_instance%init()
        gen_B = profiles_generation

        if (gen_B <= gen_A) then
            print *, 'FAIL: profiles_generation did not advance:', gen_A, gen_B
            error stop
        end if

        call kim_aperfuns(nfuns, x_core, funs_B)

        ! B's electron density at x_core must be ~2x A's (auto-invalidated cache).
        n_expected_B = 2.0_dp * funs_A(1)
        if (abs(funs_B(1) - n_expected_B) > 1.0e-6_dp * abs(n_expected_B)) then
            print *, 'FAIL: cache did NOT invalidate on profile change.'
            print *, '  n_e(A)          = ', funs_A(1)
            print *, '  n_e(B) returned = ', funs_B(1)
            print *, '  n_e(B) expected = ', n_expected_B, ' (2x A)'
            error stop
        end if

        print *, 'PASS: true-background cache auto-invalidates on profile change'
    end subroutine

    subroutine make_test_profiles(npts, r_prof, n_prof, Te_prof, Ti_prof, &
                                  q_prof, Er_prof)
        ! Analytic profiles on r in [2, 60] cm with q in [1.5, 2.5];
        ! solver domain [r_min, r_plas] = [10, 50].
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
        ! written in the exact group order kim_read_config reads them.
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
        write(iunit, '(A)') " output_path = './out_periodic_bg_test/'"
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
