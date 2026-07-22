program test_periodic_assembly
    ! Validates assemble_periodic_matrices (Phase-2.4): assembly of the dense
    ! complex Fourier matrices K^{rhoPhi}_{m,m'} and K^{rhoB}_{m,m'} over one
    ! period, reusing the Phase-2.3 periodic plasma and the Phase-1 fused kernel
    ! (hatG_rho_phi / hatG_rho_B).
    !
    ! Discrete k-grid: k_m = 2*pi*m/L, m = -M..M (row/col index im = m + M + 1).
    ! Because the periodic background makes the integrand L-periodic and the
    ! Fourier phase exp(i(k_m - k_m') r_g) is exactly L-periodic
    ! ((k_m - k_m')*L = 2*pi*(m - m')), the periodic trapezoidal (equidistant)
    ! rule is spectrally accurate. The window grid rg_grid%xb is endpoint-
    ! exclusive over one period (h = L/N, xb(N) = rm + L/2 - h), so the N points
    ! are the N DISTINCT samples and all are summed with equal weight:
    !   K_{m,m'} = (2*pi/N) * sum_{j=1}^{N} hatG(plasma, k_m, k_m', j).
    !
    ! Setup mirrors test_periodic_background_build.f90 (in-memory profiles,
    ! (m,n) = (-6, 2) so q resonant at q = 3, build a periodic window at rm).
    !
    ! Assertions:
    !   (shape)          size(K,1) == size(K,2) == 2M+1 for both matrices;
    !   (characterization) for (m,m') = (0,0) and (1,-2), an INLINE brute-force
    !                    sum with the SAME formula reproduces the stored element
    !                    to < 1e-12*(1 + |elem|) -- pins k_m mapping, h factor,
    !                    index mapping, and seam handling;
    !   (diagonal nonzero) |Kphi(M+1, M+1)| > 0 (the m = m' = 0 element);
    !   (not Toeplitz)   elements with equal (m - m') but different m differ
    !                    (the FLR factors depend on both k_m and k_m');
    !   (finiteness)     all elements finite.

    use KIM_kinds_m, only: dp
    use periodic_assembly_m, only: assemble_periodic_matrices

    implicit none

    call test_assemble()

    print *, 'All tests PASSED'
    stop 0

contains

    subroutine test_assemble()
        use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
        use config_m, only: profiles_in_memory, nml_config_path
        use species_m, only: plasma, set_profiles_from_arrays
        use kim_base_m, only: kim_t
        use kim_mod_m, only: from_kim_factory_get_kim
        use kim_resonances_m, only: r_res
        use grid_m, only: rg_grid
        use periodic_background_m, only: build_periodic_plasma

        integer, parameter :: npts = 201
        integer, parameter :: M = 3        ! small: yields a 7x7 matrix

        real(dp) :: r_prof(npts), n_prof(npts), Te_prof(npts)
        real(dp) :: Ti_prof(npts), q_prof(npts), Er_prof(npts)
        class(kim_t), allocatable :: kim_instance

        complex(dp), allocatable :: Kphi(:,:), KB(:,:), Kjphi(:,:), KjB(:,:)
        real(dp) :: rm, dx_asis, dx_tr, rho_L_rm, L
        integer :: n_rg, N, dim, im, imp

        call make_test_profiles(npts, r_prof, n_prof, Te_prof, Ti_prof, &
                                q_prof, Er_prof)

        call write_test_namelist('./KIM_config_periodic_assembly_test.nml')
        nml_config_path = './KIM_config_periodic_assembly_test.nml'

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

        ! Representative ion Larmor radius near rm on the global plasma. The
        ! production kernel deliberately takes electrons in the zero-FLR limit,
        ! so ion FLR must set the test window and supply the non-Toeplitz signal.
        rho_L_rm = rho_L_near(rm)
        print *, 'rho_L(rm) [ion] = ', rho_L_rm

        ! Window half-widths: 5 and 10 ion Larmor radii (guard against a degenerate
        ! window if rho_L is unexpectedly small).
        dx_asis =  5.0_dp * rho_L_rm
        dx_tr   = 10.0_dp * rho_L_rm
        if (dx_asis < 1.0e-3_dp) then
            dx_asis = 1.0_dp
            dx_tr   = 2.0_dp
            print *, 'NOTE: rho_L(rm) tiny; using fixed dx_asis/dx_tr [cm]'
        end if
        n_rg = 96
        L = 2.0_dp * (dx_asis + dx_tr)

        call build_periodic_plasma(rm, dx_asis, dx_tr, n_rg)

        call test_configured_dispatch(plasma)

        N = rg_grid%npts_b
        if (N /= n_rg) then
            print *, 'FAIL: rg_grid%npts_b /= n_rg; got ', N
            error stop
        end if

        call assemble_periodic_matrices(plasma, L, M, Kphi, KB, Kjphi, KjB)

        ! (shape) all four matrices are (2M+1) x (2M+1).
        dim = 2 * M + 1
        if (size(Kphi, 1) /= dim .or. size(Kphi, 2) /= dim .or. &
            size(KB, 1)   /= dim .or. size(KB, 2)   /= dim .or. &
            size(Kjphi, 1) /= dim .or. size(Kjphi, 2) /= dim .or. &
            size(KjB, 1)  /= dim .or. size(KjB, 2)  /= dim) then
            print *, 'FAIL: matrix shape wrong'
            print *, '  Kphi:', size(Kphi, 1), size(Kphi, 2), &
                     '  KB:', size(KB, 1), size(KB, 2), &
                     '  Kjphi:', size(Kjphi, 1), size(Kjphi, 2), &
                     '  KjB:', size(KjB, 1), size(KjB, 2), ' expected ', dim
            error stop
        end if
        print *, 'PASS: shape (2M+1)x(2M+1) =', dim, 'x', dim

        ! (characterization) recompute two elements by an inline brute-force sum
        ! with the SAME quadrature formula and compare to the stored elements.
        call check_element(plasma, Kphi, KB, Kjphi, KjB, L, M, N,  0,  0)
        call check_element(plasma, Kphi, KB, Kjphi, KjB, L, M, N,  1, -2)
        print *, 'PASS: characterization -- brute-force sum reproduces elements'

        ! (diagonal nonzero) the m = m' = 0 element.
        im  = 0 + M + 1
        if (.not. (abs(Kphi(im, im)) > 0.0_dp)) then
            print *, 'FAIL: Kphi(0,0) element is zero'
            error stop
        end if
        print *, 'PASS: diagonal (m=m''=0) element nonzero, |Kphi| =', abs(Kphi(im, im))

        ! (not Toeplitz) elements sharing (m - m') but with different m must
        ! differ (the FLR factors depend on both k_m and k_m'). Try two
        ! independent super-diagonal pairs to guard against a coincidental
        ! equality; passing if EITHER pair differs above a tiny tolerance.
        if (.not. ( not_equal_pair(Kphi, M, -1, -2,  0, -1) .or. &
                    not_equal_pair(Kphi, M,  1,  0,  2,  1) )) then
            print *, 'FAIL: Kphi appears Toeplitz (constant along a diagonal)'
            error stop
        end if
        print *, 'PASS: Kphi is NOT Toeplitz (k_m/k_m'' dependence confirmed)'

        ! (finiteness) every element of all four matrices is finite.
        do imp = 1, dim
            do im = 1, dim
                if (.not. ieee_is_finite(real(Kphi(im, imp), dp)) .or. &
                    .not. ieee_is_finite(aimag(Kphi(im, imp))) .or. &
                    .not. ieee_is_finite(real(KB(im, imp), dp)) .or. &
                    .not. ieee_is_finite(aimag(KB(im, imp))) .or. &
                    .not. ieee_is_finite(real(Kjphi(im, imp), dp)) .or. &
                    .not. ieee_is_finite(aimag(Kjphi(im, imp))) .or. &
                    .not. ieee_is_finite(real(KjB(im, imp), dp)) .or. &
                    .not. ieee_is_finite(aimag(KjB(im, imp)))) then
                    print *, 'FAIL: non-finite matrix element at', im, imp
                    error stop
                end if
            end do
        end do
        print *, 'PASS: all matrix elements finite'
    end subroutine test_assemble

    subroutine test_configured_dispatch(plasma)
        use collisionless_fourier_kernel_m, only: collisionless_ion_cores, &
            configured_hatG_rho_phi, configured_hatG_rho_B, &
            configured_hatG_j_phi, configured_hatG_j_B
        use config_m, only: ion_collision_model, collisionless_kpar_epsilon, &
            turn_off_electrons, turn_off_ions
        use constants_m, only: com_unit, pi
        use flr2_fourier_kernel_m, only: hatG_rho_phi, hatG_rho_B, &
            hatG_j_phi, hatG_j_B, flr_arg_pair_sp, core_rho_phi_sp, &
            core_rho_B_sp, core_j_phi_sp, core_j_B_sp
        use grid_m, only: rg_grid
        use species_m, only: plasma_t

        type(plasma_t), intent(inout) :: plasma
        real(dp), parameter :: kr = 0.17_dp, krp = -0.09_dp
        integer :: j, sp
        real(dp) :: bplus, bcross, saved_nu
        complex(dp) :: phase, electron_phi, electron_B, electron_jphi, electron_jB
        complex(dp) :: ion_phi, ion_B, ion_jphi, ion_jB
        complex(dp) :: core_phi, core_B, core_jphi, core_jB
        complex(dp) :: actual_phi, actual_B, actual_jphi, actual_jB

        j = max(1, rg_grid%npts_b / 2)
        collisionless_kpar_epsilon = 1.0e-5_dp
        turn_off_electrons = .false.
        turn_off_ions = .false.

        ! The configured wrapper must be an exact delegation in the default FP
        ! mode so existing production results do not change.
        ion_collision_model = 'FokkerPlanck'
        call assert_dispatch_close('FP rho-Phi delegation', &
            configured_hatG_rho_phi(plasma, kr, krp, j), hatG_rho_phi(plasma, kr, krp, j))
        call assert_dispatch_close('FP rho-B delegation', &
            configured_hatG_rho_B(plasma, kr, krp, j), hatG_rho_B(plasma, kr, krp, j))
        call assert_dispatch_close('FP j-Phi delegation', &
            configured_hatG_j_phi(plasma, kr, krp, j), hatG_j_phi(plasma, kr, krp, j))
        call assert_dispatch_close('FP j-B delegation', &
            configured_hatG_j_B(plasma, kr, krp, j), hatG_j_B(plasma, kr, krp, j))

        phase = exp(-com_unit * (kr - krp) * rg_grid%xb(j))
        call flr_arg_pair_sp(plasma, 0, kr, krp, j, bplus, bcross)
        electron_phi = phase * core_rho_phi_sp(plasma, 0, bplus, bcross, j) / (8.0_dp * pi**2)
        electron_B = phase * core_rho_B_sp(plasma, 0, bplus, bcross, j) / (8.0_dp * pi**2)
        electron_jphi = phase * core_j_phi_sp(plasma, 0, bplus, bcross, j) / (8.0_dp * pi**2)
        electron_jB = phase * core_j_B_sp(plasma, 0, bplus, bcross, j) / (8.0_dp * pi**2)

        ion_phi = (0.0_dp, 0.0_dp)
        ion_B = (0.0_dp, 0.0_dp)
        ion_jphi = (0.0_dp, 0.0_dp)
        ion_jB = (0.0_dp, 0.0_dp)
        do sp = 1, plasma%n_species - 1
            call collisionless_ion_cores(plasma, sp, kr, krp, j, &
                collisionless_kpar_epsilon, core_phi, core_B, core_jphi, core_jB)
            ion_phi = ion_phi + phase * core_phi
            ion_B = ion_B + phase * core_B
            ion_jphi = ion_jphi + phase * core_jphi
            ion_jB = ion_jB + phase * core_jB
        end do

        ion_collision_model = 'collisionless'
        actual_phi = configured_hatG_rho_phi(plasma, kr, krp, j)
        actual_B = configured_hatG_rho_B(plasma, kr, krp, j)
        actual_jphi = configured_hatG_j_phi(plasma, kr, krp, j)
        actual_jB = configured_hatG_j_B(plasma, kr, krp, j)
        call assert_dispatch_close('hybrid rho-Phi composition', actual_phi, electron_phi + ion_phi)
        call assert_dispatch_close('hybrid rho-B composition', actual_B, electron_B + ion_B)
        call assert_dispatch_close('hybrid j-Phi composition', actual_jphi, electron_jphi + ion_jphi)
        call assert_dispatch_close('hybrid j-B composition', actual_jB, electron_jB + ion_jB)
        if (abs(ion_jphi) <= tiny(1.0_dp) .and. abs(ion_jB) <= tiny(1.0_dp)) then
            print *, 'FAIL: configured collisionless ion current is zero'
            error stop
        end if

        turn_off_electrons = .true.
        call assert_dispatch_close('electron-off rho-Phi', &
            configured_hatG_rho_phi(plasma, kr, krp, j), ion_phi)
        call assert_dispatch_close('electron-off j-Phi', &
            configured_hatG_j_phi(plasma, kr, krp, j), ion_jphi)
        saved_nu = plasma%spec(1)%nu(j)
        plasma%spec(1)%nu(j) = saved_nu * 17.0_dp
        call assert_dispatch_close('collisionless ion nu independence', &
            configured_hatG_rho_phi(plasma, kr, krp, j), ion_phi)
        plasma%spec(1)%nu(j) = saved_nu

        turn_off_electrons = .false.
        turn_off_ions = .true.
        call assert_dispatch_close('ion-off rho-Phi', &
            configured_hatG_rho_phi(plasma, kr, krp, j), electron_phi)
        call assert_dispatch_close('ion-off j-Phi', &
            configured_hatG_j_phi(plasma, kr, krp, j), electron_jphi)

        turn_off_ions = .false.
        ion_collision_model = 'FokkerPlanck'
        print *, 'PASS: configured FP/collisionless hybrid dispatch'
    end subroutine test_configured_dispatch

    subroutine assert_dispatch_close(label, actual, expected)
        character(*), intent(in) :: label
        complex(dp), intent(in) :: actual, expected
        real(dp) :: tolerance

        tolerance = 2.0e-12_dp * max(1.0_dp, abs(expected))
        if (abs(actual - expected) > tolerance) then
            print *, 'FAIL: ', label
            print *, '  actual   = ', actual
            print *, '  expected = ', expected
            print *, '  error    = ', abs(actual - expected), ' tolerance = ', tolerance
            error stop
        end if
    end subroutine assert_dispatch_close

    !> Recompute K^{rhoPhi}, K^{rhoB}, K^{jPhi} and K^{jB} at (m,m') by an INLINE
    !> brute-force periodic-trapezoidal sum using the SAME formula as the module,
    !> and assert every stored element matches to < 1e-12*(1 + |elem|).
    subroutine check_element(plasma, Kphi, KB, Kjphi, KjB, L, M, N, mm, mmp)
        use species_m, only: plasma_t
        use flr2_fourier_kernel_m, only: hatG_rho_phi, hatG_rho_B, hatG_j_phi, hatG_j_B
        use constants_m, only: pi

        type(plasma_t), intent(in) :: plasma
        complex(dp), intent(in) :: Kphi(:,:), KB(:,:), Kjphi(:,:), KjB(:,:)
        real(dp), intent(in) :: L
        integer, intent(in) :: M, N, mm, mmp

        complex(dp) :: acc_phi, acc_B, ref_phi, ref_B
        complex(dp) :: acc_jphi, acc_jB, ref_jphi, ref_jB
        real(dp) :: k_m, k_mp, tol
        integer :: im, imp, j

        k_m  = 2.0_dp * pi * real(mm,  dp) / L
        k_mp = 2.0_dp * pi * real(mmp, dp) / L

        ! SAME formula as assemble_periodic_matrices: endpoint-exclusive window,
        ! sum ALL N distinct samples with equal weight 2*pi/N.
        acc_phi  = (0.0_dp, 0.0_dp)
        acc_B    = (0.0_dp, 0.0_dp)
        acc_jphi = (0.0_dp, 0.0_dp)
        acc_jB   = (0.0_dp, 0.0_dp)
        do j = 1, N
            acc_phi  = acc_phi  + hatG_rho_phi(plasma, k_m, k_mp, j)
            acc_B    = acc_B    + hatG_rho_B(plasma, k_m, k_mp, j)
            acc_jphi = acc_jphi + hatG_j_phi(plasma, k_m, k_mp, j)
            acc_jB   = acc_jB   + hatG_j_B(plasma, k_m, k_mp, j)
        end do
        ref_phi  = (2.0_dp * pi / real(N, dp)) * acc_phi
        ref_B    = (2.0_dp * pi / real(N, dp)) * acc_B
        ref_jphi = (2.0_dp * pi / real(N, dp)) * acc_jphi
        ref_jB   = (2.0_dp * pi / real(N, dp)) * acc_jB

        im  = mm  + M + 1
        imp = mmp + M + 1

        tol = 1.0e-12_dp * (1.0_dp + abs(ref_phi))
        if (abs(Kphi(im, imp) - ref_phi) >= tol) then
            print *, 'FAIL: Kphi element mismatch at (m,mp)=', mm, mmp
            print *, '  stored =', Kphi(im, imp), ' ref =', ref_phi
            print *, '  |diff| =', abs(Kphi(im, imp) - ref_phi), ' tol =', tol
            error stop
        end if

        tol = 1.0e-12_dp * (1.0_dp + abs(ref_B))
        if (abs(KB(im, imp) - ref_B) >= tol) then
            print *, 'FAIL: KB element mismatch at (m,mp)=', mm, mmp
            print *, '  stored =', KB(im, imp), ' ref =', ref_B
            print *, '  |diff| =', abs(KB(im, imp) - ref_B), ' tol =', tol
            error stop
        end if

        tol = 1.0e-12_dp * (1.0_dp + abs(ref_jphi))
        if (abs(Kjphi(im, imp) - ref_jphi) >= tol) then
            print *, 'FAIL: Kjphi element mismatch at (m,mp)=', mm, mmp
            print *, '  stored =', Kjphi(im, imp), ' ref =', ref_jphi
            print *, '  |diff| =', abs(Kjphi(im, imp) - ref_jphi), ' tol =', tol
            error stop
        end if

        tol = 1.0e-12_dp * (1.0_dp + abs(ref_jB))
        if (abs(KjB(im, imp) - ref_jB) >= tol) then
            print *, 'FAIL: KjB element mismatch at (m,mp)=', mm, mmp
            print *, '  stored =', KjB(im, imp), ' ref =', ref_jB
            print *, '  |diff| =', abs(KjB(im, imp) - ref_jB), ' tol =', tol
            error stop
        end if
    end subroutine check_element

    !> True if the two elements at (m1,mp1) and (m2,mp2) differ by more than a
    !> tiny tolerance. Used to show K is not constant along a diagonal.
    logical function not_equal_pair(K, M, m1, mp1, m2, mp2) result(differ)
        complex(dp), intent(in) :: K(:,:)
        integer, intent(in) :: M, m1, mp1, m2, mp2
        real(dp) :: d
        d = abs(K(m1 + M + 1, mp1 + M + 1) - K(m2 + M + 1, mp2 + M + 1))
        differ = d > 1.0e3_dp * tiny(1.0_dp)
    end function not_equal_pair

    real(dp) function rho_L_near(r) result(val)
        !> Ion Larmor radius of the global plasma at the grid point nearest
        !> radius r (the plasma is still on the global grid here).
        use species_m, only: plasma
        real(dp), intent(in) :: r
        integer :: jn
        jn = minloc(abs(plasma%r_grid - r), dim=1)
        val = plasma%spec(1)%rho_L(jn)
    end function rho_L_near

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
        write(iunit, '(A)') " output_path = './out_periodic_assembly_test/'"
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

end program test_periodic_assembly
