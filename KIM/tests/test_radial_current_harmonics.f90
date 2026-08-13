program test_radial_current_harmonics
    use KIM_kinds_m, only: dp

    implicit none

    call run_harmonic_study()

    print *, 'Radial-current harmonic study PASSED'
    stop 0

contains

    subroutine run_harmonic_study()
        use config_m, only: nml_config_path, profiles_in_memory
        use grid_m, only: rg_grid
        use kim_base_m, only: kim_t
        use kim_mod_m, only: from_kim_factory_get_kim
        use kim_resonances_m, only: r_res
        use periodic_assembly_m, only: assemble_periodic_matrices
        use periodic_background_m, only: build_periodic_plasma
        use periodic_solve_m, only: reconstruct_jrad, solve_periodic
        use radial_current_fourier_kernel_m, only: hatG_jrad_bparallel, &
            hatG_jrad_br, hatG_jrad_phi
        use species_m, only: plasma, set_profiles_from_arrays

        integer, parameter :: max_cutoff = 6
        integer, parameter :: fourier_cutoff = 3
        integer, parameter :: n_profile = 201
        integer, parameter :: n_window = 64
        integer, parameter :: dim = 2 * fourier_cutoff + 1
        real(dp) :: r(n_profile), density(n_profile), te(n_profile), ti(n_profile)
        real(dp) :: q(n_profile), er(n_profile)
        class(kim_t), allocatable :: kim
        complex(dp), allocatable :: Kphi(:,:), KB(:,:), Kjphi(:,:), KjB(:,:)
        complex(dp), allocatable :: unused_jrphi(:,:), unused_jrB(:,:), phi_m(:)
        complex(dp) :: Kjrphi(dim, dim, 0:max_cutoff)
        complex(dp) :: KjrB(dim, dim, 0:max_cutoff)
        complex(dp) :: KjrBparallel(dim, dim, 0:max_cutoff)
        complex(dp) :: jrad(n_window, 0:max_cutoff)
        complex(dp) :: sample(3, 0:max_cutoff)
        real(dp) :: rm, rho_l_rm, dx_asis, dx_tr, period
        real(dp) :: full_norm(3), shell_rel(3), sample_rel(3)
        real(dp) :: jrad_full_norm, jrad_shell_rel, jrad_tail_rel
        integer :: cutoff, info, j, nearest

        call make_test_profiles(r, density, te, ti, q, er)
        call write_test_namelist('KIM_config_radial_harmonics_test.nml')
        nml_config_path = 'KIM_config_radial_harmonics_test.nml'
        profiles_in_memory = .true.
        call kim_init()
        call set_profiles_from_arrays(r, density, te, ti, q, er, n_profile)
        call from_kim_factory_get_kim('electrostatic', kim)
        call kim%init()

        call prepare_resonances
        if (.not. (r_res > 0.0_dp)) error stop 'radial harmonic study found no resonance'
        rm = r_res
        nearest = minloc(abs(plasma%r_grid - rm), dim=1)
        rho_l_rm = plasma%spec(1)%rho_L(nearest)
        dx_asis = 5.0_dp * rho_l_rm
        dx_tr = 10.0_dp * rho_l_rm
        period = 2.0_dp * (dx_asis + dx_tr)
        call build_periodic_plasma(rm, dx_asis, dx_tr, n_window)

        call assemble_periodic_matrices(plasma, period, fourier_cutoff, Kphi, KB, &
            Kjphi, KjB, unused_jrphi, unused_jrB)
        call solve_periodic(Kphi, KB, period, fourier_cutoff, (1.0_dp, 0.0_dp), &
            phi_m, info)
        if (info /= 0) error stop 'radial harmonic study periodic solve failed'

        j = rg_grid%npts_b / 2
        do cutoff = 0, max_cutoff
            call assemble_radial_matrices(period, fourier_cutoff, cutoff, &
                Kjrphi(:, :, cutoff), KjrB(:, :, cutoff), KjrBparallel(:, :, cutoff))
            jrad(:, cutoff) = reconstruct_jrad(Kjrphi(:, :, cutoff), &
                KjrB(:, :, cutoff), phi_m, (1.0_dp, 0.0_dp), period, &
                fourier_cutoff, rg_grid%xb)
            sample(1, cutoff) = hatG_jrad_phi(plasma, 2.0_dp / period, &
                -1.0_dp / period, j, cutoff)
            sample(2, cutoff) = hatG_jrad_br(plasma, 2.0_dp / period, &
                -1.0_dp / period, j, cutoff)
            sample(3, cutoff) = hatG_jrad_bparallel(plasma, 2.0_dp / period, &
                -1.0_dp / period, j, cutoff)
        end do

        call assert_finite_3d(Kjrphi, 'KjrPhi')
        call assert_finite_3d(KjrB, 'KjrBr')
        call assert_finite_3d(KjrBparallel, 'KjrBparallel')
        call assert_finite_2d(jrad, 'jrad')
        call assert_finite_2d(sample, 'fixed-sample kernels')

        full_norm = [frobenius_norm(Kjrphi(:, :, max_cutoff)), &
            frobenius_norm(KjrB(:, :, max_cutoff)), &
            frobenius_norm(KjrBparallel(:, :, max_cutoff))]
        if (any(full_norm <= tiny(1.0_dp))) then
            error stop 'a converged radial-current matrix is zero'
        end if
        jrad_full_norm = vector_norm(jrad(:, max_cutoff))
        if (jrad_full_norm <= tiny(1.0_dp)) error stop 'converged jrad is zero'

        print '(A)', 'RADIAL_HARMONIC_MATRIX_SHELLS'
        print '(A)', ' cutoff       KjrPhi         KjrBr  KjrBparallel'
        do cutoff = 1, max_cutoff
            shell_rel = [frobenius_norm(Kjrphi(:, :, cutoff) - Kjrphi(:, :, cutoff - 1)), &
                frobenius_norm(KjrB(:, :, cutoff) - KjrB(:, :, cutoff - 1)), &
                frobenius_norm(KjrBparallel(:, :, cutoff) &
                - KjrBparallel(:, :, cutoff - 1))] / full_norm
            print '(I7,3ES15.6)', cutoff, shell_rel
        end do

        print '(A)', 'RADIAL_HARMONIC_FIXED_SAMPLE_SHELLS'
        print '(A)', ' cutoff        jr-Phi          jr-Br   jr-Bparallel'
        do cutoff = 1, max_cutoff
            sample_rel = abs(sample(:, cutoff) - sample(:, cutoff - 1)) &
                / max(abs(sample(:, max_cutoff)), tiny(1.0_dp))
            print '(I7,3ES15.6)', cutoff, sample_rel
        end do

        print '(A)', 'RADIAL_HARMONIC_JRAD_CONVERGENCE'
        print '(A)', ' cutoff      shell/full        tail/full'
        do cutoff = 0, max_cutoff
            if (cutoff == 0) then
                jrad_shell_rel = vector_norm(jrad(:, 0)) / jrad_full_norm
            else
                jrad_shell_rel = vector_norm(jrad(:, cutoff) - jrad(:, cutoff - 1)) &
                    / jrad_full_norm
            end if
            jrad_tail_rel = vector_norm(jrad(:, max_cutoff) - jrad(:, cutoff)) &
                / jrad_full_norm
            print '(I7,2ES15.6)', cutoff, jrad_shell_rel, jrad_tail_rel
        end do
        print '(A,ES15.6)', 'RADIAL_HARMONIC_JRAD_FULL_NORM ', jrad_full_norm

        if (frobenius_norm(KjrBparallel(:, :, 1) - KjrBparallel(:, :, 0)) &
                <= tiny(1.0_dp)) then
            error stop 'ell=+/-1 double radial insertion unexpectedly vanished'
        end if
        if (vector_norm(jrad(:, 1) - jrad(:, 0)) <= tiny(1.0_dp)) then
            error stop 'ell=+/-1 did not change reconstructed jrad'
        end if
    end subroutine run_harmonic_study

    subroutine assemble_radial_matrices(period, M, cutoff, Kjrphi, KjrB, KjrBparallel)
        use constants_m, only: pi
        use grid_m, only: rg_grid
        use periodic_assembly_m, only: k_of_m
        use radial_current_fourier_kernel_m, only: hatG_jrad_bparallel, &
            hatG_jrad_br, hatG_jrad_phi
        use species_m, only: plasma

        real(dp), intent(in) :: period
        integer, intent(in) :: M, cutoff
        complex(dp), intent(out) :: Kjrphi(:, :), KjrB(:, :), KjrBparallel(:, :)
        complex(dp) :: acc_phi, acc_B, acc_Bparallel
        real(dp) :: kr_response, kr_source, weight
        integer :: im, imp, j, m_response, m_source

        weight = 2.0_dp * pi / real(rg_grid%npts_b, dp)
        do m_source = -M, M
            imp = m_source + M + 1
            kr_source = k_of_m(m_source, period)
            do m_response = -M, M
                im = m_response + M + 1
                kr_response = k_of_m(m_response, period)
                acc_phi = (0.0_dp, 0.0_dp)
                acc_B = (0.0_dp, 0.0_dp)
                acc_Bparallel = (0.0_dp, 0.0_dp)
                do j = 1, rg_grid%npts_b
                    acc_phi = acc_phi + hatG_jrad_phi(plasma, kr_response, kr_source, &
                        j, cutoff)
                    acc_B = acc_B + hatG_jrad_br(plasma, kr_response, kr_source, &
                        j, cutoff)
                    acc_Bparallel = acc_Bparallel + hatG_jrad_bparallel(plasma, &
                        kr_response, kr_source, j, cutoff)
                end do
                Kjrphi(im, imp) = weight * acc_phi
                KjrB(im, imp) = weight * acc_B
                KjrBparallel(im, imp) = weight * acc_Bparallel
            end do
        end do
    end subroutine assemble_radial_matrices

    real(dp) function frobenius_norm(matrix) result(norm)
        complex(dp), intent(in) :: matrix(:, :)
        norm = sqrt(sum(abs(matrix)**2))
    end function frobenius_norm

    real(dp) function vector_norm(vector) result(norm)
        complex(dp), intent(in) :: vector(:)
        norm = sqrt(sum(abs(vector)**2))
    end function vector_norm

    subroutine assert_finite_2d(values, label)
        use, intrinsic :: ieee_arithmetic, only: ieee_is_finite

        complex(dp), intent(in) :: values(:, :)
        character(*), intent(in) :: label

        if (.not. all(ieee_is_finite(real(values, dp))) .or. &
            .not. all(ieee_is_finite(aimag(values)))) error stop label//' is non-finite'
    end subroutine assert_finite_2d

    subroutine assert_finite_3d(values, label)
        use, intrinsic :: ieee_arithmetic, only: ieee_is_finite

        complex(dp), intent(in) :: values(:, :, :)
        character(*), intent(in) :: label

        if (.not. all(ieee_is_finite(real(values, dp))) .or. &
            .not. all(ieee_is_finite(aimag(values)))) error stop label//' is non-finite'
    end subroutine assert_finite_3d

    subroutine make_test_profiles(r, density, te, ti, q, er)
        real(dp), intent(out) :: r(:), density(:), te(:), ti(:), q(:), er(:)
        integer :: i
        real(dp) :: fraction

        do i = 1, size(r)
            fraction = real(i - 1, dp) / real(size(r) - 1, dp)
            r(i) = 2.0_dp + 58.0_dp * fraction
            density(i) = 2.0e13_dp * (1.1_dp - r(i) / 100.0_dp)
            te(i) = 1.0e3_dp * (1.2_dp - r(i) / 100.0_dp)
            ti(i) = te(i)
            q(i) = 1.5_dp + 2.5_dp * fraction
            er(i) = -0.5_dp
        end do
    end subroutine make_test_profiles

    subroutine write_test_namelist(path)
        character(*), intent(in) :: path
        integer :: unit

        open(newunit=unit, file=path, status='replace', action='write')
        write(unit, '(A)') '&KIM_CONFIG'
        write(unit, '(A)') ' number_of_ion_species = 1'
        write(unit, '(A)') ' artificial_debye_case = 0'
        write(unit, '(A)') " type_of_run = 'electrostatic_periodic'"
        write(unit, '(A)') " collision_model = 'FokkerPlanck'"
        write(unit, '(A)') " ion_collision_model = 'FokkerPlanck'"
        write(unit, '(A)') ' read_species_from_namelist = .false.'
        write(unit, '(A)') " plasma_type = 'D'"
        write(unit, '(A)') ' turn_off_electrons = .false.'
        write(unit, '(A)') ' turn_off_ions = .false.'
        write(unit, '(A)') ' rescale_density = .false.'
        write(unit, '(A)') ' number_density_rescale = 1.0'
        write(unit, '(A)') ' ion_flr_scale_factor = 1.0'
        write(unit, '(A)') '/'
        write(unit, '(A)') '&WKB_DISPERSION'
        write(unit, '(A)') '/'
        write(unit, '(A)') '&KIM_IO'
        write(unit, '(A)') " profile_location = './'"
        write(unit, '(A)') " output_path = './out_radial_harmonics/'"
        write(unit, '(A)') ' hdf5_input = .false.'
        write(unit, '(A)') ' hdf5_output = .false.'
        write(unit, '(A)') ' data_verbosity = 0'
        write(unit, '(A)') ' calculate_asymptotics = .false.'
        write(unit, '(A)') '/'
        write(unit, '(A)') '&KIM_SETUP'
        write(unit, '(A)') ' btor = -18000.0'
        write(unit, '(A)') ' R0 = 165.0'
        write(unit, '(A)') ' m_mode = -6'
        write(unit, '(A)') ' n_mode = 2'
        write(unit, '(A)') ' omega = 0.0'
        write(unit, '(A)') ' spline_base = 1'
        write(unit, '(A)') ' type_br_field = 12'
        write(unit, '(A)') ' collisions_off = .false.'
        write(unit, '(A)') ' set_profiles_constant = 0'
        write(unit, '(A)') ' bc_type = 3'
        write(unit, '(A)') ' mphi_max = 6'
        write(unit, '(A)') '/'
        write(unit, '(A)') '&KIM_GRID'
        write(unit, '(A)') " grid_spacing_rg = 'equidistant'"
        write(unit, '(A)') " grid_spacing_xl = 'equidistant'"
        write(unit, '(A)') ' l_space_dim = 32'
        write(unit, '(A)') ' rg_space_dim = 32'
        write(unit, '(A)') " theta_integration = 'GaussLegendre'"
        write(unit, '(A)') ' Larmor_skip_factor = 5'
        write(unit, '(A)') ' gauss_int_nodes_Ntheta = 17'
        write(unit, '(A)') ' gauss_int_nodes_Nx = 31'
        write(unit, '(A)') ' gauss_int_nodes_Nxp = 30'
        write(unit, '(A)') ' r_plas = 50.0'
        write(unit, '(A)') ' r_min = 10.0'
        write(unit, '(A)') ' width_res = 0.5'
        write(unit, '(A)') ' ampl_res = 15.0'
        write(unit, '(A)') ' hrmax_scaling = 1.0'
        write(unit, '(A)') ' kernel_taper_skip_threshold = 1.0e-6'
        write(unit, '(A)') '/'
        write(unit, '(A)') '&KIM_PROFILES'
        write(unit, '(A)') '/'
        write(unit, '(A)') '&KIM_PERIODIC'
        write(unit, '(A)') '/'
        close(unit)
    end subroutine write_test_namelist

end program test_radial_current_harmonics
