program test_collisionless_ion_assembly

    use KIM_kinds_m, only: dp
    use config_m, only: profiles_in_memory, nml_config_path, ion_collision_model, &
        collisionless_kpar_epsilon, ion_fp_collision_scale, artificial_debye_case
    use species_m, only: plasma, scale_fp_collision_frequency, set_profiles_from_arrays
    use grid_m, only: xl_grid, kernel_taper_skip_threshold
    use kernel_m, only: kernel_spl_t, FP_fill_kernels, assemble_configured_fp_charge, &
        Krook_fill_kernel_phi_ions_collisionless, initialize_krook_mphi, &
        FP_calc_kernel_element_ions, Krook_calc_kernel_rho_term_by_term, &
        larmor_taper_band_limit
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
    complex(dp), allocatable :: weak_fp_ion_phi(:,:), weak_fp_ion_B(:,:)
    complex(dp) :: rho_phi_projection, rho_B_projection
    real(dp) :: scale, offdiag_max, rho_phi_error, rho_B_error
    real(dp) :: rho_phi_shape_error, rho_B_shape_error
    real(dp) :: requested_fp_scale, requested_kpar_epsilon
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
    requested_fp_scale = requested_ion_fp_scale()
    requested_kpar_epsilon = requested_collisionless_kpar_epsilon()
    call write_test_namelist(&
        './KIM_config_collisionless_assembly.nml', requested_fp_scale, &
        requested_kpar_epsilon)
    nml_config_path = './KIM_config_collisionless_assembly.nml'
    profiles_in_memory = .true.
    call kim_init()
    if (trim(ion_collision_model) /= 'collisionless') then
        print *, 'FAIL: ion_collision_model was not read from KIM_CONFIG'
        error stop
    end if
    if (abs(collisionless_kpar_epsilon - requested_kpar_epsilon) > 1.0e-15_dp) then
        print *, 'FAIL: collisionless_kpar_epsilon was not read from KIM_CONFIG'
        error stop
    end if
    if (abs(ion_fp_collision_scale - requested_fp_scale) > 1.0e-15_dp) then
        print *, 'FAIL: ion_fp_collision_scale was not read from KIM_CONFIG'
        error stop
    end if
    call set_profiles_from_arrays(r_prof, n_prof, Te_prof, Ti_prof, q_prof, Er_prof, npts)
    call from_kim_factory_get_kim('electrostatic', kim_instance)
    call kim_instance%init()
    print *, 'ion max |x1| = ', maxval(abs(plasma%spec(1)%x1_cc))
    print *, 'ion max |x2| = ', maxval(abs(plasma%spec(1)%x2_cc(:,0)))

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
    weak_fp_ion_phi = Kphi%Kllp_i(:,:,1)
    weak_fp_ion_B = KB%Kllp_i(:,:,1)
    call assemble_configured_fp_charge(Kphi, KB)

    rho_phi_error = relative_frobenius_error(&
        Kphi%Kllp_i(:,:,1), weak_fp_ion_phi)
    rho_B_error = relative_frobenius_error(&
        KB%Kllp_i(:,:,1), weak_fp_ion_B)
    rho_phi_projection = least_squares_projection(&
        Kphi%Kllp_i(:,:,1), weak_fp_ion_phi)
    rho_B_projection = least_squares_projection(&
        KB%Kllp_i(:,:,1), weak_fp_ion_B)
    rho_phi_shape_error = relative_frobenius_error(&
        Kphi%Kllp_i(:,:,1), rho_phi_projection * weak_fp_ion_phi)
    rho_B_shape_error = relative_frobenius_error(&
        KB%Kllp_i(:,:,1), rho_B_projection * weak_fp_ion_B)
    print *, 'weak-FP -> collisionless rho-Phi relative error = ', rho_phi_error
    print *, 'weak-FP -> collisionless rho-B relative error   = ', rho_B_error
    print *, 'rho-Phi collisionless/FP norm ratio             = ', &
        frobenius_norm(Kphi%Kllp_i(:,:,1)) / frobenius_norm(weak_fp_ion_phi)
    print *, 'rho-B collisionless/FP norm ratio               = ', &
        frobenius_norm(KB%Kllp_i(:,:,1)) / frobenius_norm(weak_fp_ion_B)
    print *, 'rho-Phi best-fit complex factor                 = ', rho_phi_projection
    print *, 'rho-B best-fit complex factor                   = ', rho_B_projection
    print *, 'rho-Phi residual after best-fit factor          = ', rho_phi_shape_error
    print *, 'rho-B residual after best-fit factor            = ', rho_B_shape_error
    call compare_rho_B_terms(weak_fp_ion_B, KB%Kllp_i(:,:,1))
    if (rho_phi_error > 0.01_dp .or. rho_B_error > 0.01_dp) then
        print *, 'FAIL: collisionless ion charge kernels do not converge to weak FP'
        error stop
    end if

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

    subroutine write_test_namelist(path, fp_scale, kpar_epsilon)
        character(len=*), intent(in) :: path
        real(dp), intent(in) :: fp_scale, kpar_epsilon
        integer :: unit
        open(newunit=unit, file=path, status='replace', action='write')
        write(unit, '(A)') '&KIM_CONFIG'
        write(unit, '(A)') ' number_of_ion_species = 1'
        write(unit, '(A)') ' artificial_debye_case = 0'
        write(unit, '(A)') " type_of_run = 'electrostatic'"
        write(unit, '(A)') " collision_model = 'FokkerPlanck'"
        write(unit, '(A)') " ion_collision_model = 'collisionless'"
        write(unit, '(A,ES15.8)') ' collisionless_kpar_epsilon = ', kpar_epsilon
        write(unit, '(A,ES15.8)') ' ion_fp_collision_scale = ', fp_scale
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

    real(dp) function requested_ion_fp_scale() result(scale)
        character(len=64) :: value
        integer :: status

        scale = 1.0e-5_dp
        call get_environment_variable(&
            'KIM_TEST_ION_FP_SCALE', value, status=status)
        if (status == 0 .and. len_trim(value) > 0) read(value, *) scale
    end function requested_ion_fp_scale

    real(dp) function requested_collisionless_kpar_epsilon() result(epsilon)
        character(len=64) :: value
        integer :: status

        epsilon = 1.0e-7_dp
        call get_environment_variable(&
            'KIM_TEST_KPAR_EPSILON', value, status=status)
        if (status == 0 .and. len_trim(value) > 0) read(value, *) epsilon
    end function requested_collisionless_kpar_epsilon

    real(dp) function relative_frobenius_error(actual, expected) result(error)
        complex(dp), intent(in) :: actual(:,:), expected(:,:)
        real(dp) :: expected_norm

        expected_norm = frobenius_norm(expected)
        error = frobenius_norm(actual - expected) &
            / max(expected_norm, tiny(1.0_dp))
    end function relative_frobenius_error

    real(dp) function frobenius_norm(matrix) result(norm)
        complex(dp), intent(in) :: matrix(:,:)

        norm = sqrt(sum(abs(matrix)**2))
    end function frobenius_norm

    complex(dp) function least_squares_projection(actual, expected) result(factor)
        complex(dp), intent(in) :: actual(:,:), expected(:,:)
        real(dp) :: expected_norm_squared

        expected_norm_squared = sum(abs(expected)**2)
        factor = sum(conjg(expected) * actual) &
            / max(expected_norm_squared, tiny(1.0_dp))
    end function least_squares_projection

    subroutine compare_rho_B_terms(fp_total, collisionless_total)
        use grid_m, only: gauss_int_nodes_Ntheta, gauss_int_nodes_Nx, &
            gauss_int_nodes_Nxp
        use integrals_gauss_m, only: gauss_config_t, init_gauss_int

        complex(dp), intent(in) :: fp_total(:,:), collisionless_total(:,:)
        type(gauss_config_t) :: gauss_conf
        complex(dp), allocatable :: fp_terms(:,:,:), collisionless_terms(:,:,:)
        complex(dp) :: fp_phi, fp_B, fp_j_phi, fp_j_B
        complex(dp) :: collisionless_phi, collisionless_B
        complex(dp) :: fp_element_terms(3), collisionless_element_terms(3)
        complex(dp) :: projection
        real(dp) :: dmax_global, error, shape_error
        integer :: l, lp, lp_lo, term

        allocate(fp_terms(size(fp_total,1), size(fp_total,2), 3))
        allocate(collisionless_terms(size(fp_total,1), size(fp_total,2), 3))
        fp_terms = (0.0_dp, 0.0_dp)
        collisionless_terms = (0.0_dp, 0.0_dp)

        gauss_conf%Nx = gauss_int_nodes_Nx
        gauss_conf%Nxp = gauss_int_nodes_Nxp
        gauss_conf%Ntheta = gauss_int_nodes_Ntheta
        call init_gauss_int(gauss_conf)
        dmax_global = larmor_taper_band_limit()

        do l = 1, size(fp_total, 1)
            lp_lo = l
            do
                if (lp_lo <= 1) exit
                if (abs(xl_grid%xb(lp_lo-1) - xl_grid%xb(l)) > dmax_global) exit
                lp_lo = lp_lo - 1
            end do
            lp_lo = min(lp_lo, l - 1)
            if (lp_lo < 1) lp_lo = 1

            do lp = lp_lo, l
                call FP_calc_kernel_element_ions(l, lp, fp_phi, fp_B, fp_j_phi, &
                    fp_j_B, gauss_conf, 1, k_rho_B_terms=fp_element_terms)
                call Krook_calc_kernel_rho_term_by_term(l, lp, collisionless_phi, &
                    collisionless_B, gauss_conf, species_first=1, species_last=1, &
                    collisionless=.true., k_rho_B_terms=collisionless_element_terms)
                fp_terms(l,lp,:) = fp_element_terms
                fp_terms(lp,l,:) = fp_element_terms
                collisionless_terms(l,lp,:) = collisionless_element_terms
                collisionless_terms(lp,l,:) = collisionless_element_terms
            end do
        end do

        print *, 'rho-B term reconstruction error (FP)          = ', &
            relative_frobenius_error(sum(fp_terms, dim=3), fp_total)
        print *, 'rho-B term reconstruction error (CL)          = ', &
            relative_frobenius_error(sum(collisionless_terms, dim=3), &
            collisionless_total)
        do term = 1, 3
            error = relative_frobenius_error(&
                collisionless_terms(:,:,term), fp_terms(:,:,term))
            projection = least_squares_projection(&
                collisionless_terms(:,:,term), fp_terms(:,:,term))
            shape_error = relative_frobenius_error(&
                collisionless_terms(:,:,term), projection * fp_terms(:,:,term))
            print *, 'rho-B F', term, 'G', term, ' relative error         = ', error
            print *, 'rho-B F', term, 'G', term, ' CL/FP norm ratio       = ', &
                frobenius_norm(collisionless_terms(:,:,term)) / &
                max(frobenius_norm(fp_terms(:,:,term)), tiny(1.0_dp))
            print *, 'rho-B F', term, 'G', term, ' FP fraction of total  = ', &
                frobenius_norm(fp_terms(:,:,term)) / &
                max(frobenius_norm(fp_total), tiny(1.0_dp))
            print *, 'rho-B F', term, 'G', term, ' CL fraction of total  = ', &
                frobenius_norm(collisionless_terms(:,:,term)) / &
                max(frobenius_norm(collisionless_total), tiny(1.0_dp))
            print *, 'rho-B F', term, 'G', term, ' best-fit factor        = ', &
                projection
            print *, 'rho-B F', term, 'G', term, ' shape residual         = ', &
                shape_error
        end do
    end subroutine compare_rho_B_terms

end program test_collisionless_ion_assembly
