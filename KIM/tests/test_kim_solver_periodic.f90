program test_kim_solver_periodic
    ! Integration test for the electrostatic_periodic run-type (Phase-2.6).
    !
    ! Drives the full lifecycle of the new run-type end-to-end:
    !   from_kim_factory_get_kim('electrostatic_periodic', kim)
    !   -> kim%init()   (global setup: grids + equilibrium + plasma quantities)
    !   -> kim%run()    (locate resonance rm, size the window from rho_L(rm),
    !                    build the periodic plasma, assemble the Fourier
    !                    matrices, solve, reconstruct dPhi on the window grid,
    !                    and pack it into fields_m::EBdat).
    !
    ! Setup mirrors test_periodic_assembly.f90 / test_periodic_solve.f90:
    ! in-memory analytic profiles, (m,n) = (-6, 2) so q is resonant at q = 3
    ! inside the solver domain, type_br_field = 12 (constant Br over the window).
    !
    ! This is the wiring/integration proof (localization/physics-vs-global is
    ! Phase-2.7). Assertions on the packed result in fields_m::EBdat:
    !   (Phi)    EBdat%Phi allocated, size == rg_grid%npts_b, all finite,
    !            and NOT all-zero (maxval(abs(EBdat%Phi)) > 0);
    !   (r_grid) EBdat%r_grid allocated, size == rg_grid%npts_b, spans the
    !            window [rm - L/2, rm + L/2].

    use KIM_kinds_m, only: dp

    implicit none

    call test_run_type_end_to_end()

    print *, 'All tests PASSED'
    stop 0

contains

    subroutine test_run_type_end_to_end()
        use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
        use config_m, only: profiles_in_memory, nml_config_path
        use species_m, only: plasma, set_profiles_from_arrays
        use kim_base_m, only: kim_t
        use kim_mod_m, only: from_kim_factory_get_kim
        use kim_resonances_m, only: r_res
        use grid_m, only: rg_grid
        use fields_m, only: EBdat

        integer, parameter :: npts = 201

        real(dp) :: r_prof(npts), n_prof(npts), Te_prof(npts)
        real(dp) :: Ti_prof(npts), q_prof(npts), Er_prof(npts)
        class(kim_t), allocatable :: kim_instance

        real(dp) :: rm, r_lo, r_hi, tol
        integer :: i, N

        call make_test_profiles(npts, r_prof, n_prof, Te_prof, Ti_prof, &
                                q_prof, Er_prof)

        call write_test_namelist('./KIM_config_periodic_run_test.nml')
        nml_config_path = './KIM_config_periodic_run_test.nml'

        profiles_in_memory = .true.
        call kim_init()
        call set_profiles_from_arrays(r_prof, n_prof, Te_prof, Ti_prof, &
                                      q_prof, Er_prof, npts)

        ! Full lifecycle through the new run-type.
        call from_kim_factory_get_kim('electrostatic_periodic', kim_instance)
        call kim_instance%init()
        call kim_instance%run()

        ! The resonance the run-type located (q = |m/n| = 3 on the global plasma).
        rm = r_res
        if (.not. (rm > 0.0_dp)) then
            print *, 'FAIL: resonance not found, r_res = ', rm
            error stop
        end if
        print *, 'resonant radius rm = ', rm

        N = rg_grid%npts_b

        ! (Phi) allocated, correct size, finite, non-zero.
        if (.not. allocated(EBdat%Phi)) then
            print *, 'FAIL: EBdat%Phi not allocated'
            error stop
        end if
        if (size(EBdat%Phi) /= N) then
            print *, 'FAIL: size(EBdat%Phi) /= rg_grid%npts_b; got ', &
                     size(EBdat%Phi), ' expected ', N
            error stop
        end if
        do i = 1, N
            if (.not. ieee_is_finite(real(EBdat%Phi(i), dp)) .or. &
                .not. ieee_is_finite(aimag(EBdat%Phi(i)))) then
                print *, 'FAIL: non-finite EBdat%Phi at', i
                error stop
            end if
        end do
        if (.not. (maxval(abs(EBdat%Phi)) > 0.0_dp)) then
            print *, 'FAIL: EBdat%Phi all zero (constant Br should drive response)'
            error stop
        end if
        print *, 'PASS: EBdat%Phi allocated, size ', N, &
                 ', finite & non-zero, max|Phi| =', maxval(abs(EBdat%Phi))

        ! (r_grid) allocated, correct size, spans the window [rm - L/2, rm + L/2].
        if (.not. allocated(EBdat%r_grid)) then
            print *, 'FAIL: EBdat%r_grid not allocated'
            error stop
        end if
        if (size(EBdat%r_grid) /= N) then
            print *, 'FAIL: size(EBdat%r_grid) /= rg_grid%npts_b; got ', &
                     size(EBdat%r_grid), ' expected ', N
            error stop
        end if

        ! The window is symmetric about rm; recover L from its span and check
        ! the endpoints match rm -+ L/2 to within a grid-spacing tolerance.
        r_lo = EBdat%r_grid(1)
        r_hi = EBdat%r_grid(N)
        tol = 1.0e-8_dp * (1.0_dp + abs(rm))
        if (abs(0.5_dp * (r_lo + r_hi) - rm) > (r_hi - r_lo) / real(N - 1, dp)) then
            print *, 'FAIL: window not centred on rm; midpoint =', &
                     0.5_dp * (r_lo + r_hi), ' rm =', rm
            error stop
        end if
        if (.not. (r_lo < rm .and. rm < r_hi)) then
            print *, 'FAIL: rm not inside window [', r_lo, ',', r_hi, ']'
            error stop
        end if
        print *, 'PASS: EBdat%r_grid spans window [', r_lo, ',', r_hi, &
                 '] about rm =', rm
    end subroutine test_run_type_end_to_end

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
        ! Minimal electrostatic-periodic FokkerPlanck configuration; m_mode = -6,
        ! n_mode = 2 makes q resonant at q = 3, type_br_field = 12 (constant Br).
        ! Deliberately OMITS the &KIM_PERIODIC group to prove the periodic_*
        ! config defaults apply when the optional namelist is absent.
        character(len=*), intent(in) :: path
        integer :: iunit

        open(newunit=iunit, file=path, status='replace', action='write')
        write(iunit, '(A)') '&KIM_CONFIG'
        write(iunit, '(A)') ' number_of_ion_species = 1'
        write(iunit, '(A)') ' artificial_debye_case = 0'
        write(iunit, '(A)') " type_of_run = 'electrostatic_periodic'"
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
        write(iunit, '(A)') " output_path = './out_periodic_run_test/'"
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

end program test_kim_solver_periodic
