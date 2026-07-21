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
    use flr2_fourier_kernel_m, only: hatG_rho_phi, hatG_rho_B
    use flr2_fourier_kernel_m, only: hatG_rho_phi_diag_sp, hatG_rho_B_diag_sp

    implicit none

    call test_populated_plasma_and_kernel_stub()
    call test_diagonal_matches_inline()
    call test_diagonal_j_matches_inline()
    call test_collapse_rho_phi()
    call test_phase_rho_phi()
    call test_collapse_rho_phi_ks2()
    call test_zero_wavenumber_limit()
    call test_collapse_rho_B()
    call test_collapse_rho_B_ks2()
    call test_zero_wavenumber_limit_rho_B()
    call test_scaled_bessel_large_arg()
    call test_large_bcross_kernel_finite()

    print *, 'All tests PASSED'
    stop 0

contains

    subroutine test_diagonal_matches_inline()
        ! Characterization test: the shared per-species diagonal integrand
        ! functions must reproduce, bit-for-bit within tolerance, the original
        ! inline expressions in calc_hatK_Phi_in_Fourier. The INLINE reference
        ! below is the permanent anchor (verbatim copy of the source-of-truth
        ! per-species rho-Phi Debye+FP term and signed rho-B FP term). It must
        ! not be refactored to call the functions it validates.
        use config_m, only: profiles_in_memory, nml_config_path
        use config_m, only: artificial_debye_case, turn_off_ions, turn_off_electrons
        use constants_m, only: com_unit, sol
        use fortnum_special, only: bessel_in
        use species_m, only: plasma, set_profiles_from_arrays, plasma_t
        use grid_m, only: rg_grid
        use kim_base_m, only: kim_t
        use kim_mod_m, only: from_kim_factory_get_kim

        integer, parameter :: npts = 101

        real(dp) :: r_prof(npts), n_prof(npts), Te_prof(npts)
        real(dp) :: Ti_prof(npts), q_prof(npts), Er_prof(npts)
        class(kim_t), allocatable :: kim_instance
        real(dp) :: kr_arr(3), kr, b, ks
        integer :: i, j, sp, jj_arr(3), k
        complex(dp) :: phi_inline, phi_func, B_inline, B_func
        real(dp) :: tol_phi, tol_B

        call make_test_profiles(npts, r_prof, n_prof, Te_prof, Ti_prof, &
                                q_prof, Er_prof)

        call write_test_namelist('./KIM_config_fourier_diag_test.nml')
        nml_config_path = './KIM_config_fourier_diag_test.nml'

        profiles_in_memory = .true.
        call kim_init()
        call set_profiles_from_arrays(r_prof, n_prof, Te_prof, Ti_prof, &
                                      q_prof, Er_prof, npts)

        call from_kim_factory_get_kim('electrostatic', kim_instance)
        call kim_instance%init()

        kr_arr = [0.1_dp, 1.0_dp, 5.0_dp]
        ! A few interior guiding-centre boundary-grid indices.
        jj_arr = [rg_grid%npts_b / 4, rg_grid%npts_b / 2, (3 * rg_grid%npts_b) / 4]

        do i = 1, size(kr_arr)
            kr = kr_arr(i)
            do k = 1, size(jj_arr)
                j = jj_arr(k)

                ! (i) INLINE reference: verbatim per-species expressions,
                ! summed over the non-turned-off species.
                phi_inline = (0.0_dp, 0.0_dp)
                B_inline = (0.0_dp, 0.0_dp)
                do sp = 0, plasma%n_species - 1
                    if (turn_off_ions .and. sp >= 1) cycle
                    if (turn_off_electrons .and. sp == 0) cycle

                    ks = plasma%ks(j)
                    b = (kr**2.0d0) * plasma%spec(sp)%rho_L(j)**2.0d0

                    if (artificial_debye_case <= 1) then
                        phi_inline = phi_inline - 1.0d0 / plasma%spec(sp)%lambda_D(j)**2.0d0
                    end if

                    if (artificial_debye_case == 0 .or. artificial_debye_case == 2) then
                        phi_inline = phi_inline + 1.0d0 / plasma%spec(sp)%lambda_D(j)**2.0d0 &
                            * com_unit * plasma%spec(sp)%vT(j)**2.0d0 * plasma%ks(j) &
                            / (plasma%spec(sp)%omega_c(j) * plasma%spec(sp)%nu(j)) * exp(-b) * &
                            (&
                                plasma%spec(sp)%I00(j, 0) * (&
                                    bessel_in(0, b) * (plasma%spec(sp)%A1(j) + plasma%spec(sp)%A2(j) * (1-b)) &
                                    + plasma%spec(sp)%A2(j) * b * bessel_in(-1, b) &
                                )&
                                + 0.5d0 * plasma%spec(sp)%I20(j, 0) * plasma%spec(sp)%A2(j) * bessel_in(0, b) &
                            )
                        B_inline = B_inline - 1.0d0 / plasma%spec(sp)%lambda_D(j)**2.0d0 * plasma%spec(sp)%vT(j)**3.0d0 &
                            / (plasma%spec(sp)%omega_c(j) * plasma%spec(sp)%nu(j) * sol) * exp(-b) * &
                            (&
                                plasma%spec(sp)%I01(j, 0) * (&
                                    bessel_in(0, b) * (plasma%spec(sp)%A1(j) + plasma%spec(sp)%A2(j) * (1-b)) &
                                    + plasma%spec(sp)%A2(j) * b * bessel_in(-1, b) &
                                )&
                                + 0.5d0 * plasma%spec(sp)%I21(j, 0) * plasma%spec(sp)%A2(j) * bessel_in(0, b) &
                            )
                    end if
                end do

                ! (ii) Shared functions: sum over the same non-turned-off species.
                phi_func = (0.0_dp, 0.0_dp)
                B_func = (0.0_dp, 0.0_dp)
                do sp = 0, plasma%n_species - 1
                    if (turn_off_ions .and. sp >= 1) cycle
                    if (turn_off_electrons .and. sp == 0) cycle
                    phi_func = phi_func + hatG_rho_phi_diag_sp(plasma, sp, kr, j)
                    B_func = B_func + hatG_rho_B_diag_sp(plasma, sp, kr, j)
                end do

                tol_phi = 1.0e-12_dp * (1.0_dp + abs(phi_inline))
                tol_B = 1.0e-12_dp * (1.0_dp + abs(B_inline))

                if (abs(phi_inline - phi_func) >= tol_phi) then
                    print *, 'FAIL: rho-Phi diagonal mismatch at kr=', kr, ' j=', j
                    print *, '  inline   = ', phi_inline
                    print *, '  fromfunc = ', phi_func
                    print *, '  |diff|   = ', abs(phi_inline - phi_func), ' tol = ', tol_phi
                    error stop
                end if
                if (abs(B_inline - B_func) >= tol_B) then
                    print *, 'FAIL: rho-B diagonal mismatch at kr=', kr, ' j=', j
                    print *, '  inline   = ', B_inline
                    print *, '  fromfunc = ', B_func
                    print *, '  |diff|   = ', abs(B_inline - B_func), ' tol = ', tol_B
                    error stop
                end if
            end do
        end do

        print *, 'PASS: shared diagonal integrand matches inline reference'
    end subroutine

    subroutine test_diagonal_j_matches_inline()
        ! Characterization test: the fused two-wavenumber j-kernels evaluated on
        ! the diagonal kr = krp must reproduce the diagonal j-kernel expressions
        ! of calc_hatK_Phi_in_Fourier (FLR2_asymptotics.f90:216-234), the
        ! source of truth for j-Phi (thesis 14.5) and j-B (thesis 14.6, reached
        ! via the identities I11 = I02 and I13 = I22).
        !
        ! The INLINE reference below is the permanent anchor (verbatim copy of
        ! the per-species source-of-truth expressions). It must not be
        ! refactored to call the functions it validates.
        !
        ! Note the deliberately absent /(4*pi): the rho kernels carry that factor
        ! only to cancel solve_periodic's Poisson 4*pi. Thesis (11.7) defines
        ! delta_j_par = K^{jPhi} delta_Phi + K^{jB} delta_Br with no 4*pi, so
        ! hatG_j_* omit it and the reference here is the raw per-species sum.
        use config_m, only: profiles_in_memory, nml_config_path
        use config_m, only: artificial_debye_case, turn_off_ions, turn_off_electrons
        use constants_m, only: com_unit, sol
        use fortnum_special, only: bessel_in
        use flr2_fourier_kernel_m, only: hatG_j_phi, hatG_j_B, kern_include_ks2
        use species_m, only: plasma, set_profiles_from_arrays
        use grid_m, only: rg_grid
        use kim_base_m, only: kim_t
        use kim_mod_m, only: from_kim_factory_get_kim

        integer, parameter :: npts = 101

        real(dp) :: r_prof(npts), n_prof(npts), Te_prof(npts)
        real(dp) :: Ti_prof(npts), q_prof(npts), Er_prof(npts)
        class(kim_t), allocatable :: kim_instance
        real(dp) :: kr_arr(3), kr, b, ks
        integer :: i, j, sp, jj_arr(3), k
        complex(dp) :: jphi_inline, jphi_func, jB_inline, jB_func
        real(dp) :: tol_jphi, tol_jB

        call make_test_profiles(npts, r_prof, n_prof, Te_prof, Ti_prof, &
                                q_prof, Er_prof)

        call write_test_namelist('./KIM_config_fourier_jdiag_test.nml')
        nml_config_path = './KIM_config_fourier_jdiag_test.nml'

        profiles_in_memory = .true.
        call kim_init()
        call set_profiles_from_arrays(r_prof, n_prof, Te_prof, Ti_prof, &
                                      q_prof, Er_prof, npts)

        call from_kim_factory_get_kim('electrostatic', kim_instance)
        call kim_instance%init()

        kern_include_ks2 = .false.

        kr_arr = [0.1_dp, 1.0_dp, 5.0_dp]
        jj_arr = [rg_grid%npts_b / 4, rg_grid%npts_b / 2, (3 * rg_grid%npts_b) / 4]

        do i = 1, size(kr_arr)
            kr = kr_arr(i)
            do k = 1, size(jj_arr)
                j = jj_arr(k)

                ! (i) INLINE reference: verbatim per-species expressions,
                ! summed over the non-turned-off species. NO /(4*pi).
                jphi_inline = (0.0_dp, 0.0_dp)
                jB_inline = (0.0_dp, 0.0_dp)
                do sp = 0, plasma%n_species - 1
                    if (turn_off_ions .and. sp >= 1) cycle
                    if (turn_off_electrons .and. sp == 0) cycle

                    ks = plasma%ks(j)
                    b = (kr**2.0d0) * plasma%spec(sp)%rho_L(j)**2.0d0

                    if (artificial_debye_case == 0 .or. artificial_debye_case == 2) then
                        jphi_inline = jphi_inline + 1.0d0 / plasma%spec(sp)%lambda_D(j)**2.0d0 &
                            * com_unit * plasma%spec(sp)%vT(j)**3.0d0 &
                            / (plasma%spec(sp)%omega_c(j) * plasma%spec(sp)%nu(j)) * ks * exp(-b) * &
                            (&
                                plasma%spec(sp)%I01(j, 0) * (&
                                    bessel_in(0, b) * (plasma%spec(sp)%A1(j) + plasma%spec(sp)%A2(j) * (1-b)) &
                                    + plasma%spec(sp)%A2(j) * b * bessel_in(-1, b) &
                                )&
                                + 0.5d0 * plasma%spec(sp)%I21(j, 0) * plasma%spec(sp)%A2(j) * bessel_in(0, b) &
                            )

                        jB_inline = jB_inline - 1.0d0 / plasma%spec(sp)%lambda_D(j)**2.0d0 &
                            * plasma%spec(sp)%vT(j)**4.0d0 &
                            / (plasma%spec(sp)%omega_c(j) * plasma%spec(sp)%nu(j) * sol) * exp(-b) * &
                            (&
                                plasma%spec(sp)%I11(j, 0) * (&
                                    bessel_in(0, b) * (plasma%spec(sp)%A1(j) + plasma%spec(sp)%A2(j) * (1-b)) &
                                    + plasma%spec(sp)%A2(j) * b * bessel_in(-1, b) &
                                )&
                                + 0.5d0 * plasma%spec(sp)%I13(j, 0) * plasma%spec(sp)%A2(j) * bessel_in(0, b) &
                            )
                    end if
                end do

                ! (ii) Fused two-wavenumber kernels on the diagonal: the phase
                ! exp(i*(kr-krp)*rg) is unity and (b_+, b_x) collapse to b.
                jphi_func = hatG_j_phi(plasma, kr, kr, j)
                jB_func = hatG_j_B(plasma, kr, kr, j)

                tol_jphi = 1.0e-12_dp * (1.0_dp + abs(jphi_inline))
                tol_jB = 1.0e-12_dp * (1.0_dp + abs(jB_inline))

                if (abs(jphi_inline - jphi_func) >= tol_jphi) then
                    print *, 'FAIL: j-Phi diagonal mismatch at kr=', kr, ' j=', j
                    print *, '  inline   = ', jphi_inline
                    print *, '  fromfunc = ', jphi_func
                    print *, '  |diff|   = ', abs(jphi_inline - jphi_func), ' tol = ', tol_jphi
                    error stop
                end if
                if (abs(jB_inline - jB_func) >= tol_jB) then
                    print *, 'FAIL: j-B diagonal mismatch at kr=', kr, ' j=', j
                    print *, '  inline   = ', jB_inline
                    print *, '  fromfunc = ', jB_func
                    print *, '  |diff|   = ', abs(jB_inline - jB_func), ' tol = ', tol_jB
                    error stop
                end if
            end do
        end do

        print *, 'PASS: fused j kernels match inline diagonal reference'
    end subroutine

    subroutine test_collapse_rho_phi()
        ! Collapse-by-construction: the fused two-wavenumber kernel evaluated on
        ! the diagonal k_r = k'_r must equal the per-species diagonal integrand
        ! summed over the non-turned-off species and divided by 4*pi (the
        ! source-of-truth normalization). Exact to machine precision because the
        ! diagonal and the fused kernel share the same per-species core.
        use config_m, only: profiles_in_memory, nml_config_path
        use config_m, only: turn_off_ions, turn_off_electrons
        use flr2_fourier_kernel_m, only: hatG_rho_phi, hatG_rho_phi_diag_sp, kern_include_ks2
        use constants_m, only: pi
        use species_m, only: plasma, set_profiles_from_arrays
        use grid_m, only: rg_grid
        use kim_base_m, only: kim_t
        use kim_mod_m, only: from_kim_factory_get_kim

        integer, parameter :: npts = 101

        real(dp) :: r_prof(npts), n_prof(npts), Te_prof(npts)
        real(dp) :: Ti_prof(npts), q_prof(npts), Er_prof(npts)
        class(kim_t), allocatable :: kim_instance
        real(dp) :: kr_arr(3), kr, tol
        integer :: i, j, sp, jj_arr(3), k
        complex(dp) :: dref, gfused

        call make_test_profiles(npts, r_prof, n_prof, Te_prof, Ti_prof, &
                                q_prof, Er_prof)

        call write_test_namelist('./KIM_config_fourier_collapse_test.nml')
        nml_config_path = './KIM_config_fourier_collapse_test.nml'

        profiles_in_memory = .true.
        call kim_init()
        call set_profiles_from_arrays(r_prof, n_prof, Te_prof, Ti_prof, &
                                      q_prof, Er_prof, npts)

        call from_kim_factory_get_kim('electrostatic', kim_instance)
        call kim_instance%init()

        kern_include_ks2 = .false.

        kr_arr = [0.1_dp, 1.0_dp, 5.0_dp]
        jj_arr = [rg_grid%npts_b / 4, rg_grid%npts_b / 2, (3 * rg_grid%npts_b) / 4]

        do i = 1, size(kr_arr)
            kr = kr_arr(i)
            do k = 1, size(jj_arr)
                j = jj_arr(k)

                dref = (0.0_dp, 0.0_dp)
                do sp = 0, plasma%n_species - 1
                    if (turn_off_ions .and. sp >= 1) cycle
                    if (turn_off_electrons .and. sp == 0) cycle
                    dref = dref + hatG_rho_phi_diag_sp(plasma, sp, kr, j)
                end do
                dref = dref / (4.0d0 * pi)

                gfused = hatG_rho_phi(plasma, kr, kr, j)

                tol = 1.0e-12_dp * (1.0_dp + abs(dref))
                if (abs(gfused - dref) >= tol) then
                    print *, 'FAIL: rho-Phi collapse mismatch at kr=', kr, ' j=', j
                    print *, '  diag/4pi = ', dref
                    print *, '  fused    = ', gfused
                    print *, '  |diff|   = ', abs(gfused - dref), ' tol = ', tol
                    error stop
                end if
            end do
        end do

        print *, 'PASS: fused rho-Phi kernel collapses to diagonal at kr=krp'
    end subroutine

    subroutine test_phase_rho_phi()
        ! Off-diagonal phase structure: for k_r /= k'_r the fused kernel carries
        ! the Fourier phase exp(i*(kr-krp)*rg). Because b_+ and b_x are symmetric
        ! under (kr,krp) -> (krp,kr), the per-species cores are equal, so the only
        ! difference between hatG(kr,krp) and hatG(krp,kr) is that phase factor.
        ! Stripping each of its own phase must leave identical (phase-free,
        ! still complex) core values. Tolerance is looser than the collapse
        ! test: this compares two independent complex-exp/Bessel evaluations,
        ! not a construction identity.
        use config_m, only: profiles_in_memory, nml_config_path
        use flr2_fourier_kernel_m, only: hatG_rho_phi, kern_include_ks2
        use constants_m, only: com_unit
        use species_m, only: plasma, set_profiles_from_arrays
        use grid_m, only: rg_grid
        use kim_base_m, only: kim_t
        use kim_mod_m, only: from_kim_factory_get_kim

        integer, parameter :: npts = 101

        real(dp) :: r_prof(npts), n_prof(npts), Te_prof(npts)
        real(dp) :: Ti_prof(npts), q_prof(npts), Er_prof(npts)
        class(kim_t), allocatable :: kim_instance
        real(dp) :: kr_arr(3), krp_arr(3), kr, krp, rg, tol
        integer :: i, j, jj_arr(3), k
        complex(dp) :: ci, g_fwd, g_rev, stripped_fwd, stripped_rev

        call make_test_profiles(npts, r_prof, n_prof, Te_prof, Ti_prof, &
                                q_prof, Er_prof)

        call write_test_namelist('./KIM_config_fourier_phase_test.nml')
        nml_config_path = './KIM_config_fourier_phase_test.nml'

        profiles_in_memory = .true.
        call kim_init()
        call set_profiles_from_arrays(r_prof, n_prof, Te_prof, Ti_prof, &
                                      q_prof, Er_prof, npts)

        call from_kim_factory_get_kim('electrostatic', kim_instance)
        call kim_instance%init()

        kern_include_ks2 = .false.
        ci = com_unit

        kr_arr = [0.1_dp, 1.0_dp, 2.0_dp]
        krp_arr = [1.0_dp, 5.0_dp, 0.5_dp]
        jj_arr = [rg_grid%npts_b / 4, rg_grid%npts_b / 2, (3 * rg_grid%npts_b) / 4]

        do i = 1, size(kr_arr)
            kr = kr_arr(i)
            krp = krp_arr(i)
            do k = 1, size(jj_arr)
                j = jj_arr(k)
                rg = rg_grid%xb(j)

                g_fwd = hatG_rho_phi(plasma, kr, krp, j)
                g_rev = hatG_rho_phi(plasma, krp, kr, j)

                stripped_fwd = g_fwd * exp(-ci * (kr - krp) * rg)
                stripped_rev = g_rev * exp(-ci * (krp - kr) * rg)

                tol = 1.0e-10_dp * (1.0_dp + abs(g_fwd))
                if (abs(stripped_fwd - stripped_rev) >= tol) then
                    print *, 'FAIL: rho-Phi phase mismatch at kr=', kr, ' krp=', krp, ' j=', j
                    print *, '  stripped_fwd = ', stripped_fwd
                    print *, '  stripped_rev = ', stripped_rev
                    print *, '  |diff|       = ', abs(stripped_fwd - stripped_rev), ' tol = ', tol
                    error stop
                end if
            end do
        end do

        print *, 'PASS: fused rho-Phi kernel carries the Fourier phase only'
    end subroutine

    subroutine test_collapse_rho_phi_ks2()
        ! Structural collapse for the k_s^2-inclusive path (kern_include_ks2 =
        ! .true.). This branch keeps k_s in the FLR arguments, so it does NOT
        ! collapse to the k_s-dropped diagonal source-of-truth. Instead, on the
        ! diagonal k_r = k'_r the two FLR arguments both reduce to the full
        ! perpendicular argument bf = (k_s^2 + k_r^2) * rho_L^2, so the fused
        ! kernel must equal the per-species core evaluated at (bf, bf), summed
        ! over the non-turned-off species and divided by 4*pi. This validates the
        ! sqrt/2 algebra of the .true. branch without needing a k_s = 0 setup.
        use config_m, only: profiles_in_memory, nml_config_path
        use config_m, only: turn_off_ions, turn_off_electrons
        use flr2_fourier_kernel_m, only: hatG_rho_phi, core_rho_phi_sp, kern_include_ks2
        use constants_m, only: pi
        use species_m, only: plasma, set_profiles_from_arrays
        use grid_m, only: rg_grid
        use kim_base_m, only: kim_t
        use kim_mod_m, only: from_kim_factory_get_kim

        integer, parameter :: npts = 101

        real(dp) :: r_prof(npts), n_prof(npts), Te_prof(npts)
        real(dp) :: Ti_prof(npts), q_prof(npts), Er_prof(npts)
        class(kim_t), allocatable :: kim_instance
        real(dp) :: kr_arr(3), kr, bf, tol
        integer :: i, j, sp, jj_arr(3), k
        complex(dp) :: ref, gfused

        call make_test_profiles(npts, r_prof, n_prof, Te_prof, Ti_prof, &
                                q_prof, Er_prof)

        call write_test_namelist('./KIM_config_fourier_collapse_ks2_test.nml')
        nml_config_path = './KIM_config_fourier_collapse_ks2_test.nml'

        profiles_in_memory = .true.
        call kim_init()
        call set_profiles_from_arrays(r_prof, n_prof, Te_prof, Ti_prof, &
                                      q_prof, Er_prof, npts)

        call from_kim_factory_get_kim('electrostatic', kim_instance)
        call kim_instance%init()

        kern_include_ks2 = .true.

        kr_arr = [0.1_dp, 1.0_dp, 5.0_dp]
        jj_arr = [rg_grid%npts_b / 4, rg_grid%npts_b / 2, (3 * rg_grid%npts_b) / 4]

        do i = 1, size(kr_arr)
            kr = kr_arr(i)
            do k = 1, size(jj_arr)
                j = jj_arr(k)

                ref = (0.0_dp, 0.0_dp)
                do sp = 0, plasma%n_species - 1
                    if (turn_off_ions .and. sp >= 1) cycle
                    if (turn_off_electrons .and. sp == 0) cycle
                    bf = (plasma%ks(j)**2 + kr**2) * plasma%spec(sp)%rho_L(j)**2
                    ref = ref + core_rho_phi_sp(plasma, sp, bf, bf, j)
                end do
                ref = ref / (4.0d0 * pi)

                gfused = hatG_rho_phi(plasma, kr, kr, j)

                tol = 1.0e-12_dp * (1.0_dp + abs(ref))
                if (abs(gfused - ref) >= tol) then
                    print *, 'FAIL: rho-Phi ks2 collapse mismatch at kr=', kr, ' j=', j
                    print *, '  core(bf,bf)/4pi = ', ref
                    print *, '  fused           = ', gfused
                    print *, '  |diff|          = ', abs(gfused - ref), ' tol = ', tol
                    error stop
                end if
            end do
        end do

        ! Defensive: restore the default so any later test sees the drop-k_s^2
        ! path in the shared module state.
        kern_include_ks2 = .false.

        print *, 'PASS: fused rho-Phi kernel (ks2 path) collapses to core(bf,bf) at kr=krp'
    end subroutine

    subroutine test_zero_wavenumber_limit()
        ! Analytic b -> 0 anchor: at kr = krp = 0 in the default drop-k_s^2 path,
        ! both FLR arguments vanish (b_+ = b_x = 0). The Bessel/exp assembly then
        ! reduces to a closed form INDEPENDENT of the Bessel routine, because
        ! exp(-b_+) = 1, (1 - b_+) = 1, I_0(0) = 1, I_{-1}(0) = 0, and the
        ! A2 * b_x * I_{-1}(b_x) term vanishes with b_x = 0. The Fourier phase
        ! exp(i*0*rg) = 1. So the per-species FP core collapses to
        !   1/lambda_D^2 * i * vT^2 * ks / (omega_c*nu) * ( I00*(A1+A2) + 0.5*I20*A2 )
        ! and the Debye term stays -1/lambda_D^2. We hard-code I_0(0)=1 and
        ! I_{-1}(0)=0 here rather than calling bessel_in, so this test validates
        ! the kernel's argument assembly against math, not against the Bessel
        ! routine it shares with the collapse/phase tests. Only valid for the
        ! .false. path (in the .true. path the arguments are rho_L^2*ks^2 /= 0).
        use config_m, only: profiles_in_memory, nml_config_path
        use config_m, only: artificial_debye_case, turn_off_ions, turn_off_electrons
        use flr2_fourier_kernel_m, only: hatG_rho_phi, kern_include_ks2
        use constants_m, only: com_unit, pi
        use species_m, only: plasma, set_profiles_from_arrays
        use grid_m, only: rg_grid
        use kim_base_m, only: kim_t
        use kim_mod_m, only: from_kim_factory_get_kim

        integer, parameter :: npts = 101

        real(dp) :: r_prof(npts), n_prof(npts), Te_prof(npts)
        real(dp) :: Ti_prof(npts), q_prof(npts), Er_prof(npts)
        class(kim_t), allocatable :: kim_instance
        real(dp) :: tol, ks
        integer :: j, sp, jj_arr(3), k
        complex(dp) :: closed, gfused

        call make_test_profiles(npts, r_prof, n_prof, Te_prof, Ti_prof, &
                                q_prof, Er_prof)

        call write_test_namelist('./KIM_config_fourier_zero_wavenumber_test.nml')
        nml_config_path = './KIM_config_fourier_zero_wavenumber_test.nml'

        profiles_in_memory = .true.
        call kim_init()
        call set_profiles_from_arrays(r_prof, n_prof, Te_prof, Ti_prof, &
                                      q_prof, Er_prof, npts)

        call from_kim_factory_get_kim('electrostatic', kim_instance)
        call kim_instance%init()

        ! Closed form derived above only holds for the drop-k_s^2 path; the
        ! toggle is deliberately fixed here (do NOT sweep it).
        kern_include_ks2 = .false.

        jj_arr = [rg_grid%npts_b / 4, rg_grid%npts_b / 2, (3 * rg_grid%npts_b) / 4]

        do k = 1, size(jj_arr)
            j = jj_arr(k)

            ! Per-species b -> 0 closed form, summed over the non-turned-off
            ! species, with I_0(0)=1 and I_{-1}(0)=0 folded in analytically.
            closed = (0.0_dp, 0.0_dp)
            do sp = 0, plasma%n_species - 1
                if (turn_off_ions .and. sp >= 1) cycle
                if (turn_off_electrons .and. sp == 0) cycle

                ks = plasma%ks(j)

                if (artificial_debye_case <= 1) then
                    closed = closed - 1.0d0 / plasma%spec(sp)%lambda_D(j)**2.0d0
                end if

                if (artificial_debye_case == 0 .or. artificial_debye_case == 2) then
                    closed = closed + 1.0d0 / plasma%spec(sp)%lambda_D(j)**2.0d0 &
                        * com_unit * plasma%spec(sp)%vT(j)**2.0d0 * ks &
                        / (plasma%spec(sp)%omega_c(j) * plasma%spec(sp)%nu(j)) * &
                        (&
                            plasma%spec(sp)%I00(j, 0) &
                                * (plasma%spec(sp)%A1(j) + plasma%spec(sp)%A2(j)) &
                            + 0.5d0 * plasma%spec(sp)%I20(j, 0) * plasma%spec(sp)%A2(j) &
                        )
                end if
            end do
            closed = closed / (4.0d0 * pi)

            gfused = hatG_rho_phi(plasma, 0.0_dp, 0.0_dp, j)

            tol = 1.0e-12_dp * (1.0_dp + abs(closed))
            if (abs(gfused - closed) >= tol) then
                print *, 'FAIL: rho-Phi b->0 limit mismatch at j=', j
                print *, '  closed = ', closed
                print *, '  fused  = ', gfused
                print *, '  |diff| = ', abs(gfused - closed), ' tol = ', tol
                error stop
            end if
        end do

        print *, 'PASS: fused rho-Phi kernel matches b->0 closed form'
    end subroutine

    subroutine test_collapse_rho_B()
        ! Collapse-by-construction for the rho-B kernel (drop-k_s^2 path): the
        ! fused two-wavenumber hatG_rho_B evaluated on the diagonal k_r = k'_r
        ! must equal the per-species signed diagonal integrand summed over the
        ! non-turned-off species and divided by 4*pi. Exact to machine precision
        ! because both share the same per-species core.
        use config_m, only: profiles_in_memory, nml_config_path
        use config_m, only: turn_off_ions, turn_off_electrons
        use flr2_fourier_kernel_m, only: hatG_rho_B, hatG_rho_B_diag_sp, kern_include_ks2
        use constants_m, only: pi
        use species_m, only: plasma, set_profiles_from_arrays
        use grid_m, only: rg_grid
        use kim_base_m, only: kim_t
        use kim_mod_m, only: from_kim_factory_get_kim

        integer, parameter :: npts = 101

        real(dp) :: r_prof(npts), n_prof(npts), Te_prof(npts)
        real(dp) :: Ti_prof(npts), q_prof(npts), Er_prof(npts)
        class(kim_t), allocatable :: kim_instance
        real(dp) :: kr_arr(3), kr, tol
        integer :: i, j, sp, jj_arr(3), k
        complex(dp) :: dref, gfused

        call make_test_profiles(npts, r_prof, n_prof, Te_prof, Ti_prof, &
                                q_prof, Er_prof)

        call write_test_namelist('./KIM_config_fourier_collapse_B_test.nml')
        nml_config_path = './KIM_config_fourier_collapse_B_test.nml'

        profiles_in_memory = .true.
        call kim_init()
        call set_profiles_from_arrays(r_prof, n_prof, Te_prof, Ti_prof, &
                                      q_prof, Er_prof, npts)

        call from_kim_factory_get_kim('electrostatic', kim_instance)
        call kim_instance%init()

        kern_include_ks2 = .false.

        kr_arr = [0.1_dp, 1.0_dp, 5.0_dp]
        jj_arr = [rg_grid%npts_b / 4, rg_grid%npts_b / 2, (3 * rg_grid%npts_b) / 4]

        do i = 1, size(kr_arr)
            kr = kr_arr(i)
            do k = 1, size(jj_arr)
                j = jj_arr(k)

                dref = (0.0_dp, 0.0_dp)
                do sp = 0, plasma%n_species - 1
                    if (turn_off_ions .and. sp >= 1) cycle
                    if (turn_off_electrons .and. sp == 0) cycle
                    dref = dref + hatG_rho_B_diag_sp(plasma, sp, kr, j)
                end do
                dref = dref / (4.0d0 * pi)

                gfused = hatG_rho_B(plasma, kr, kr, j)

                tol = 1.0e-12_dp * (1.0_dp + abs(dref))
                if (abs(gfused - dref) >= tol) then
                    print *, 'FAIL: rho-B collapse mismatch at kr=', kr, ' j=', j
                    print *, '  diag/4pi = ', dref
                    print *, '  fused    = ', gfused
                    print *, '  |diff|   = ', abs(gfused - dref), ' tol = ', tol
                    error stop
                end if
            end do
        end do

        print *, 'PASS: fused rho-B kernel collapses to diagonal at kr=krp'
    end subroutine

    subroutine test_collapse_rho_B_ks2()
        ! Structural collapse for the rho-B k_s^2-inclusive path
        ! (kern_include_ks2 = .true.). On the diagonal k_r = k'_r both FLR
        ! arguments reduce to the full perpendicular argument
        ! bf = (k_s^2 + k_r^2) * rho_L^2, so the fused kernel must equal the
        ! per-species core evaluated at (bf, bf), summed over the non-turned-off
        ! species and divided by 4*pi.
        use config_m, only: profiles_in_memory, nml_config_path
        use config_m, only: turn_off_ions, turn_off_electrons
        use flr2_fourier_kernel_m, only: hatG_rho_B, core_rho_B_sp, kern_include_ks2
        use constants_m, only: pi
        use species_m, only: plasma, set_profiles_from_arrays
        use grid_m, only: rg_grid
        use kim_base_m, only: kim_t
        use kim_mod_m, only: from_kim_factory_get_kim

        integer, parameter :: npts = 101

        real(dp) :: r_prof(npts), n_prof(npts), Te_prof(npts)
        real(dp) :: Ti_prof(npts), q_prof(npts), Er_prof(npts)
        class(kim_t), allocatable :: kim_instance
        real(dp) :: kr_arr(3), kr, bf, tol
        integer :: i, j, sp, jj_arr(3), k
        complex(dp) :: ref, gfused

        call make_test_profiles(npts, r_prof, n_prof, Te_prof, Ti_prof, &
                                q_prof, Er_prof)

        call write_test_namelist('./KIM_config_fourier_collapse_B_ks2_test.nml')
        nml_config_path = './KIM_config_fourier_collapse_B_ks2_test.nml'

        profiles_in_memory = .true.
        call kim_init()
        call set_profiles_from_arrays(r_prof, n_prof, Te_prof, Ti_prof, &
                                      q_prof, Er_prof, npts)

        call from_kim_factory_get_kim('electrostatic', kim_instance)
        call kim_instance%init()

        kern_include_ks2 = .true.

        kr_arr = [0.1_dp, 1.0_dp, 5.0_dp]
        jj_arr = [rg_grid%npts_b / 4, rg_grid%npts_b / 2, (3 * rg_grid%npts_b) / 4]

        do i = 1, size(kr_arr)
            kr = kr_arr(i)
            do k = 1, size(jj_arr)
                j = jj_arr(k)

                ref = (0.0_dp, 0.0_dp)
                do sp = 0, plasma%n_species - 1
                    if (turn_off_ions .and. sp >= 1) cycle
                    if (turn_off_electrons .and. sp == 0) cycle
                    bf = (plasma%ks(j)**2 + kr**2) * plasma%spec(sp)%rho_L(j)**2
                    ref = ref + core_rho_B_sp(plasma, sp, bf, bf, j)
                end do
                ref = ref / (4.0d0 * pi)

                gfused = hatG_rho_B(plasma, kr, kr, j)

                tol = 1.0e-12_dp * (1.0_dp + abs(ref))
                if (abs(gfused - ref) >= tol) then
                    print *, 'FAIL: rho-B ks2 collapse mismatch at kr=', kr, ' j=', j
                    print *, '  core(bf,bf)/4pi = ', ref
                    print *, '  fused           = ', gfused
                    print *, '  |diff|          = ', abs(gfused - ref), ' tol = ', tol
                    error stop
                end if
            end do
        end do

        ! Defensive: restore the default so any later test sees the drop-k_s^2
        ! path in the shared module state.
        kern_include_ks2 = .false.

        print *, 'PASS: fused rho-B kernel (ks2 path) collapses to core(bf,bf) at kr=krp'
    end subroutine

    subroutine test_zero_wavenumber_limit_rho_B()
        ! Analytic b -> 0 anchor for the rho-B kernel (drop-k_s^2 path): at
        ! kr = krp = 0 both FLR arguments vanish (b_+ = b_x = 0). With
        ! I_0(0) = 1, I_{-1}(0) = 0, exp(-b_+) = 1, (1 - b_+) = 1 and the
        ! A2 * b_x * I_{-1}(b_x) term vanishing, and the Fourier phase = 1, the
        ! per-species signed FP core collapses to
        !   -1/lambda_D^2 * vT^3/(omega_c*nu*sol) * ( I01*(A1+A2) + 0.5*I21*A2 ).
        ! We fold I_0(0)=1 and I_{-1}(0)=0 in analytically (no bessel_in call),
        ! so this validates the kernel's argument assembly against math, not the
        ! Bessel routine. rho-B has NO Debye term. Only valid for the .false. path.
        use config_m, only: profiles_in_memory, nml_config_path
        use config_m, only: artificial_debye_case, turn_off_ions, turn_off_electrons
        use flr2_fourier_kernel_m, only: hatG_rho_B, kern_include_ks2
        use constants_m, only: pi, sol
        use species_m, only: plasma, set_profiles_from_arrays
        use grid_m, only: rg_grid
        use kim_base_m, only: kim_t
        use kim_mod_m, only: from_kim_factory_get_kim

        integer, parameter :: npts = 101

        real(dp) :: r_prof(npts), n_prof(npts), Te_prof(npts)
        real(dp) :: Ti_prof(npts), q_prof(npts), Er_prof(npts)
        class(kim_t), allocatable :: kim_instance
        real(dp) :: tol
        integer :: j, sp, jj_arr(3), k
        complex(dp) :: closed, gfused

        call make_test_profiles(npts, r_prof, n_prof, Te_prof, Ti_prof, &
                                q_prof, Er_prof)

        call write_test_namelist('./KIM_config_fourier_zero_wavenumber_B_test.nml')
        nml_config_path = './KIM_config_fourier_zero_wavenumber_B_test.nml'

        profiles_in_memory = .true.
        call kim_init()
        call set_profiles_from_arrays(r_prof, n_prof, Te_prof, Ti_prof, &
                                      q_prof, Er_prof, npts)

        call from_kim_factory_get_kim('electrostatic', kim_instance)
        call kim_instance%init()

        ! Closed form derived above only holds for the drop-k_s^2 path.
        kern_include_ks2 = .false.

        jj_arr = [rg_grid%npts_b / 4, rg_grid%npts_b / 2, (3 * rg_grid%npts_b) / 4]

        do k = 1, size(jj_arr)
            j = jj_arr(k)

            ! Per-species b -> 0 closed form, summed over the non-turned-off
            ! species, with I_0(0)=1 and I_{-1}(0)=0 folded in analytically.
            ! rho-B carries no Debye term.
            closed = (0.0_dp, 0.0_dp)
            do sp = 0, plasma%n_species - 1
                if (turn_off_ions .and. sp >= 1) cycle
                if (turn_off_electrons .and. sp == 0) cycle

                if (artificial_debye_case == 0 .or. artificial_debye_case == 2) then
                    closed = closed - 1.0d0 / plasma%spec(sp)%lambda_D(j)**2.0d0 &
                        * plasma%spec(sp)%vT(j)**3.0d0 &
                        / (plasma%spec(sp)%omega_c(j) * plasma%spec(sp)%nu(j) * sol) * &
                        (&
                            plasma%spec(sp)%I01(j, 0) &
                                * (plasma%spec(sp)%A1(j) + plasma%spec(sp)%A2(j)) &
                            + 0.5d0 * plasma%spec(sp)%I21(j, 0) * plasma%spec(sp)%A2(j) &
                        )
                end if
            end do
            closed = closed / (4.0d0 * pi)

            gfused = hatG_rho_B(plasma, 0.0_dp, 0.0_dp, j)

            tol = 1.0e-12_dp * (1.0_dp + abs(closed))
            if (abs(gfused - closed) >= tol) then
                print *, 'FAIL: rho-B b->0 limit mismatch at j=', j
                print *, '  closed = ', closed
                print *, '  fused  = ', gfused
                print *, '  |diff| = ', abs(gfused - closed), ' tol = ', tol
                error stop
            end if
        end do

        print *, 'PASS: fused rho-B kernel matches b->0 closed form'
    end subroutine

    subroutine test_scaled_bessel_large_arg()
        ! Overflow-safe scaled Bessel products in the large-argument regime.
        ! At b_+ = b_x = 1000 the unscaled I_0(1000) overflows to +Inf (near
        ! the ~710 argument limit), so the naive exp(-b_+)*I_0(b_x) is Inf*0.
        ! The asymptotic branch must instead return finite values matching the
        ! analytic limits exp(-b)*I_0(b) -> 1/sqrt(2*pi*b) and I_{-1}/I_0 -> 1.
        use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
        use flr2_fourier_kernel_m, only: scaled_bessel_pair
        use constants_m, only: pi

        real(dp) :: b, sI0, sIm1, analytic_sI0, tol

        b = 1000.0_dp
        call scaled_bessel_pair(b, b, sI0, sIm1)

        if (.not. ieee_is_finite(sI0)) then
            print *, 'FAIL: scaled_bessel_pair sI0 not finite at b=1000: ', sI0
            error stop
        end if
        if (.not. ieee_is_finite(sIm1)) then
            print *, 'FAIL: scaled_bessel_pair sIm1 not finite at b=1000: ', sIm1
            error stop
        end if

        ! exp(-b)*I_0(b) -> 1/sqrt(2*pi*b) as b -> infinity.
        analytic_sI0 = 1.0_dp / sqrt(2.0_dp * pi * b)
        tol = 1.0e-3_dp * analytic_sI0
        if (abs(sI0 - analytic_sI0) >= tol) then
            print *, 'FAIL: scaled sI0 off analytic limit at b=1000'
            print *, '  sI0      = ', sI0
            print *, '  analytic = ', analytic_sI0
            print *, '  |diff|   = ', abs(sI0 - analytic_sI0), ' tol = ', tol
            error stop
        end if

        ! I_{-1}(b)/I_0(b) -> 1, so sIm1 is positive and within ~2% of sI0.
        if (.not. (sIm1 > 0.0_dp)) then
            print *, 'FAIL: scaled sIm1 must be positive at b=1000, got ', sIm1
            error stop
        end if
        if (abs(sIm1 - sI0) >= 0.02_dp * sI0) then
            print *, 'FAIL: scaled sIm1 not within 2% of sI0 at b=1000'
            print *, '  sI0    = ', sI0
            print *, '  sIm1   = ', sIm1
            print *, '  |diff| = ', abs(sIm1 - sI0), ' tol = ', 0.02_dp * sI0
            error stop
        end if

        print *, 'PASS: scaled_bessel_pair finite and matches large-arg asymptotics'
    end subroutine

    subroutine test_large_bcross_kernel_finite()
        ! End-to-end overflow guard: the per-species cores must stay finite when
        ! driven with a forced large FLR argument (b_x = 1000) that would
        ! overflow the naive exp(-b_+)*I_0(b_x) path. Uses the public
        ! core_*_sp entry points with an interior guiding-centre index and the
        ! default artificial_debye_case = 0 (FP term active).
        use config_m, only: profiles_in_memory, nml_config_path
        use config_m, only: artificial_debye_case
        use flr2_fourier_kernel_m, only: core_rho_phi_sp, core_rho_B_sp
        use species_m, only: plasma, set_profiles_from_arrays
        use grid_m, only: rg_grid
        use kim_base_m, only: kim_t
        use kim_mod_m, only: from_kim_factory_get_kim

        integer, parameter :: npts = 101

        real(dp) :: r_prof(npts), n_prof(npts), Te_prof(npts)
        real(dp) :: Ti_prof(npts), q_prof(npts), Er_prof(npts)
        class(kim_t), allocatable :: kim_instance
        real(dp) :: big
        integer :: j
        complex(dp) :: g_phi, g_B

        call make_test_profiles(npts, r_prof, n_prof, Te_prof, Ti_prof, &
                                q_prof, Er_prof)

        call write_test_namelist('./KIM_config_fourier_large_bcross_test.nml')
        nml_config_path = './KIM_config_fourier_large_bcross_test.nml'

        profiles_in_memory = .true.
        call kim_init()
        call set_profiles_from_arrays(r_prof, n_prof, Te_prof, Ti_prof, &
                                      q_prof, Er_prof, npts)

        call from_kim_factory_get_kim('electrostatic', kim_instance)
        call kim_instance%init()

        ! The namelist sets artificial_debye_case = 0, so the FP branch (which
        ! evaluates the Bessel products) is active.
        if (.not. (artificial_debye_case == 0 .or. artificial_debye_case == 2)) then
            print *, 'FAIL: FP branch inactive; artificial_debye_case = ', artificial_debye_case
            error stop
        end if

        j = rg_grid%npts_b / 2
        big = 1000.0_dp

        g_phi = core_rho_phi_sp(plasma, 0, big, big, j)
        g_B = core_rho_B_sp(plasma, 0, big, big, j)

        if (.not. is_finite_complex(g_phi)) then
            print *, 'FAIL: core_rho_phi_sp not finite at bcross=1000: ', g_phi
            error stop
        end if
        if (.not. is_finite_complex(g_B)) then
            print *, 'FAIL: core_rho_B_sp not finite at bcross=1000: ', g_B
            error stop
        end if

        print *, 'PASS: core_*_sp stay finite at large bcross (overflow guard)'
        print *, '  core_rho_phi_sp = ', g_phi
        print *, '  core_rho_B_sp   = ', g_B
    end subroutine

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
