program test_collisionless_ion_assembly

    use KIM_kinds_m, only: dp
    use config_m, only: profiles_in_memory, nml_config_path, ion_collision_model, &
        collisionless_kpar_epsilon, ion_fp_collision_scale, artificial_debye_case
    use species_m, only: plasma, scale_fp_collision_frequency, set_profiles_from_arrays
    use grid_m, only: xl_grid, kernel_taper_skip_threshold
    use kernel_m, only: kernel_spl_t, FP_fill_kernels, assemble_configured_fp_charge, &
        Krook_fill_kernel_phi_ions_collisionless, initialize_krook_mphi
    use integrands_gauss_m, only: integration_point_t
    use kim_base_m, only: kim_t
    use kim_mod_m, only: from_kim_factory_get_kim
    use, intrinsic :: ieee_arithmetic, only: ieee_is_finite

    implicit none

    integer, parameter :: npts = 101
    real(dp) :: r_prof(npts), n_prof(npts), Te_prof(npts)
    real(dp) :: Ti_prof(npts), q_prof(npts), Er_prof(npts)
    type(kernel_spl_t) :: Kphi, KB, Kjphi, KjB
    type(integration_point_t) :: krook_point
    class(kim_t), allocatable :: kim_instance
    complex(dp), allocatable :: ion_sum(:,:), electron_phi(:,:), electron_B(:,:)
    complex(dp), allocatable :: wide_band_ion(:,:), narrow_band_ion(:,:)
    real(dp) :: scale, offdiag_max
    integer :: i, j, wide_band_entries, narrow_band_entries

    krook_point%mphi = 17
    call initialize_krook_mphi(krook_point)
    if (krook_point%mphi /= 0) then
        print *, 'FAIL: collisionless Krook integration did not select mphi = 0'
        error stop
    end if

    if (scale_fp_collision_frequency(12.5_dp, 0, 0.25_dp) /= 12.5_dp) then
        print *, 'FAIL: FP ion collision scaling changed the electron frequency'
        error stop
    end if
    if (scale_fp_collision_frequency(12.5_dp, 1, 0.25_dp) /= 3.125_dp) then
        print *, 'FAIL: FP ion collision scaling was not applied to ions'
        error stop
    end if

    call make_test_profiles(r_prof, n_prof, Te_prof, Ti_prof, q_prof, Er_prof)
    call write_test_namelist('./KIM_config_collisionless_assembly.nml')
    nml_config_path = './KIM_config_collisionless_assembly.nml'
    profiles_in_memory = .true.
    call kim_init()
    if (trim(ion_collision_model) /= 'collisionless') then
        print *, 'FAIL: ion_collision_model was not read from KIM_CONFIG'
        error stop
    end if
    if (abs(collisionless_kpar_epsilon - 1.0e-5_dp) > 1.0e-15_dp) then
        print *, 'FAIL: collisionless_kpar_epsilon was not read from KIM_CONFIG'
        error stop
    end if
    if (abs(ion_fp_collision_scale - 0.25_dp) > 1.0e-15_dp) then
        print *, 'FAIL: ion_fp_collision_scale was not read from KIM_CONFIG'
        error stop
    end if
    call set_profiles_from_arrays(r_prof, n_prof, Te_prof, Ti_prof, q_prof, Er_prof, npts)
    call from_kim_factory_get_kim('electrostatic', kim_instance)
    call kim_instance%init()

    call Kphi%init_kernel(xl_grid%npts_b, xl_grid%npts_b)
    call KB%init_kernel(xl_grid%npts_b, xl_grid%npts_b)
    call Kjphi%init_kernel(xl_grid%npts_b, xl_grid%npts_b)
    call KjB%init_kernel(xl_grid%npts_b, xl_grid%npts_b)

    artificial_debye_case = 1
    call Krook_fill_kernel_phi_ions_collisionless(Kphi, KB)
    offdiag_max = 0.0_dp
    do i = 2, size(Kphi%Kllp_i, 1)
        offdiag_max = max(offdiag_max, abs(Kphi%Kllp_i(i, i-1, 1)))
    end do
    if (offdiag_max <= 0.0_dp) then
        print *, 'FAIL: collisionless F0 term omits adjacent overlapping basis functions'
        error stop
    end if
    artificial_debye_case = 0

    call FP_fill_kernels(Kphi, KB, Kjphi, KjB)
    electron_phi = Kphi%Kllp_e
    electron_B = KB%Kllp_e
    call assemble_configured_fp_charge(Kphi, KB)

    if (maxval(abs(Kphi%Kllp_e - electron_phi)) > 0.0_dp .or. &
        maxval(abs(KB%Kllp_e - electron_B)) > 0.0_dp) then
        print *, 'FAIL: hybrid fill modified FP electron matrices'
        error stop
    end if

    do i = 1, size(Kphi%Kllp, 1)
        do j = 1, size(Kphi%Kllp, 2)
            if (.not. ieee_is_finite(real(Kphi%Kllp(i,j), dp)) .or. &
                .not. ieee_is_finite(aimag(Kphi%Kllp(i,j))) .or. &
                .not. ieee_is_finite(real(KB%Kllp(i,j), dp)) .or. &
                .not. ieee_is_finite(aimag(KB%Kllp(i,j)))) then
                print *, 'FAIL: non-finite collisionless ion kernel at ', i, j
                error stop
            end if
        end do
    end do

    scale = max(1.0_dp, maxval(abs(Kphi%Kllp)))
    if (maxval(abs(Kphi%Kllp)) <= 0.0_dp) then
        print *, 'FAIL: collisionless ion rho-Phi matrix is identically zero'
        error stop
    end if
    if (maxval(abs(Kphi%Kllp - transpose(Kphi%Kllp))) > 1.0e-12_dp * scale) then
        print *, 'FAIL: collisionless ion rho-Phi matrix is not symmetric'
        error stop
    end if

    ion_sum = sum(Kphi%Kllp_i, dim=3)
    if (maxval(abs(Kphi%Kllp - electron_phi - ion_sum)) > 1.0e-12_dp * scale) then
        print *, 'FAIL: hybrid total differs from FP-electron plus collisionless-ion sum'
        error stop
    end if

    offdiag_max = 0.0_dp
    do i = 1, size(Kphi%Kllp, 1)
        do j = 1, size(Kphi%Kllp, 2)
            if (i /= j) offdiag_max = max(offdiag_max, abs(Kphi%Kllp(i,j)))
        end do
    end do
    if (offdiag_max <= 0.0_dp) then
        print *, 'FAIL: finite-ion-FLR collisionless matrix has no off-diagonal support'
        error stop
    end if

    if (minval(plasma%spec(0)%nu) <= 0.0_dp) then
        print *, 'FAIL: collisional electron frequency was not preserved'
        error stop
    end if

    wide_band_ion = Kphi%Kllp_i(:,:,1)
    scale = max(1.0_dp, maxval(abs(wide_band_ion)))
    wide_band_entries = count(abs(wide_band_ion) > 1.0e-12_dp * scale)

    kernel_taper_skip_threshold = exp(-1.0_dp)
    call Krook_fill_kernel_phi_ions_collisionless(Kphi, KB)
    narrow_band_ion = Kphi%Kllp_i(:,:,1)
    narrow_band_entries = count(abs(narrow_band_ion) > 1.0e-12_dp * scale)
    kernel_taper_skip_threshold = 1.0e-6_dp

    if (wide_band_entries <= narrow_band_entries) then
        print *, 'FAIL: collisionless ion kernel did not respond to the taper threshold'
        print *, '  wide-band resolved entries   = ', wide_band_entries
        print *, '  narrow-band resolved entries = ', narrow_band_entries
        error stop
    end if

    print *, 'PASS: ion-only collisionless assembly is finite, symmetric, and FLR-coupled'

contains

    subroutine make_test_profiles(r, n, Te, Ti, q, Er)
        real(dp), intent(out) :: r(npts), n(npts), Te(npts), Ti(npts), q(npts), Er(npts)
        integer :: k
        do k = 1, npts
            r(k) = 2.0_dp + 58.0_dp * real(k - 1, dp) / real(npts - 1, dp)
            n(k) = 2.0e13_dp * (1.1_dp - r(k) / 100.0_dp)
            Te(k) = 1.0e3_dp * (1.2_dp - r(k) / 100.0_dp)
            Ti(k) = Te(k)
            q(k) = 1.5_dp + (r(k) - 2.0_dp) / 58.0_dp
            Er(k) = -0.5_dp
        end do
    end subroutine make_test_profiles

    subroutine write_test_namelist(path)
        character(len=*), intent(in) :: path
        integer :: unit
        open(newunit=unit, file=path, status='replace', action='write')
        write(unit, '(A)') '&KIM_CONFIG'
        write(unit, '(A)') ' number_of_ion_species = 1'
        write(unit, '(A)') ' artificial_debye_case = 0'
        write(unit, '(A)') " type_of_run = 'electrostatic'"
        write(unit, '(A)') " collision_model = 'FokkerPlanck'"
        write(unit, '(A)') " ion_collision_model = 'collisionless'"
        write(unit, '(A)') ' collisionless_kpar_epsilon = 1.0e-5'
        write(unit, '(A)') ' ion_fp_collision_scale = 0.25'
        write(unit, '(A)') ' read_species_from_namelist = .false.'
        write(unit, '(A)') ' turn_off_ions = .false.'
        write(unit, '(A)') ' turn_off_electrons = .false.'
        write(unit, '(A)') " plasma_type = 'D'"
        write(unit, '(A)') ' rescale_density = .false.'
        write(unit, '(A)') ' number_density_rescale = 1.0'
        write(unit, '(A)') ' ion_flr_scale_factor = 1.0'
        write(unit, '(A)') '/'
        write(unit, '(A)') '&WKB_DISPERSION'
        write(unit, '(A)') '/'
        write(unit, '(A)') '&KIM_IO'
        write(unit, '(A)') " profile_location = './'"
        write(unit, '(A)') " output_path = './out_collisionless_assembly/'"
        write(unit, '(A)') ' hdf5_input = .false.'
        write(unit, '(A)') ' hdf5_output = .false.'
        write(unit, '(A)') ' log_level = 2'
        write(unit, '(A)') ' data_verbosity = 0'
        write(unit, '(A)') ' calculate_asymptotics = .false.'
        write(unit, '(A)') " h5_out_file = ''"
        write(unit, '(A)') ' write_diagnostics_dat = .false.'
        write(unit, '(A)') '/'
        write(unit, '(A)') '&KIM_SETUP'
        write(unit, '(A)') ' btor = -18000.0'
        write(unit, '(A)') ' R0 = 165.0'
        write(unit, '(A)') ' m_mode = -2'
        write(unit, '(A)') ' n_mode = 1'
        write(unit, '(A)') ' omega = 0.0'
        write(unit, '(A)') ' spline_base = 1'
        write(unit, '(A)') ' type_br_field = 12'
        write(unit, '(A)') ' collisions_off = .false.'
        write(unit, '(A)') ' set_profiles_constant = 0'
        write(unit, '(A)') ' bc_type = 3'
        write(unit, '(A)') ' mphi_max = 0'
        write(unit, '(A)') ' Br_boundary_re = 1.0'
        write(unit, '(A)') ' Br_boundary_im = 0.0'
        write(unit, '(A)') '/'
        write(unit, '(A)') '&KIM_GRID'
        write(unit, '(A)') " grid_spacing_rg = 'equidistant'"
        write(unit, '(A)') " grid_spacing_xl = 'equidistant'"
        write(unit, '(A)') ' l_space_dim = 32'
        write(unit, '(A)') ' rg_space_dim = 64'
        write(unit, '(A)') " theta_integration = 'GaussLegendre'"
        write(unit, '(A)') " theta_integration_method = ''"
        write(unit, '(A)') ' Larmor_skip_factor = 5'
        write(unit, '(A)') ' gauss_int_nodes_Ntheta = 7'
        write(unit, '(A)') ' gauss_int_nodes_Nx = 11'
        write(unit, '(A)') ' gauss_int_nodes_Nxp = 10'
        write(unit, '(A)') ' r_plas = 34.0'
        write(unit, '(A)') ' r_min = 28.0'
        write(unit, '(A)') ' width_res = 0.5'
        write(unit, '(A)') ' ampl_res = 15.0'
        write(unit, '(A)') ' hrmax_scaling = 1.0'
        write(unit, '(A)') ' rkf45_atol = 1.0e-9'
        write(unit, '(A)') ' rkf45_rtol = 1.0e-6'
        write(unit, '(A)') ' kernel_taper_skip_threshold = 1.0e-6'
        write(unit, '(A)') " quadpack_algorithm = 'QAG'"
        write(unit, '(A)') ' quadpack_key = 6'
        write(unit, '(A)') ' quadpack_limit = 500'
        write(unit, '(A)') ' quadpack_epsabs = 1.0e-10'
        write(unit, '(A)') ' quadpack_epsrel = 1.0e-10'
        write(unit, '(A)') ' quadpack_use_u_substitution = .true.'
        write(unit, '(A)') '/'
        write(unit, '(A)') '&KIM_PROFILES'
        write(unit, '(A)') '/'
        close(unit)
    end subroutine write_test_namelist

end program test_collisionless_ion_assembly
