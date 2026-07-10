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
    use flr2_fourier_kernel_m, only: hatG_rho_phi

    implicit none

    call test_populated_plasma_and_kernel_stub()

    print *, 'All tests PASSED'
    stop 0

contains

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
