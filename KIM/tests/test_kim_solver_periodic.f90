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
    call test_collisionless_ions_end_to_end()
    call test_multi_ion_order_independence()
    call test_global_approximation_enabled()
    call test_species_resolved_currents(.false.)
    call test_species_resolved_currents(.true.)

    print *, 'All tests PASSED'
    stop 0

contains

    subroutine test_species_resolved_currents(collisionless_ions)
        use config_m, only: profiles_in_memory, nml_config_path
        use fields_m, only: EBdat, EBdat_t
        use IO_collection_m, only: deinitialize_hdf5_output, h5id
        use KAMEL_hdf5_tools, only: h5_get, h5_obj_exists
        use kim_base_m, only: kim_t
        use kim_mod_m, only: from_kim_factory_get_kim
        use species_m, only: set_profiles_from_arrays

        logical, intent(in) :: collisionless_ions

        integer, parameter :: npts = 201
        integer, parameter :: ion_masses(2) = [2, 3]
        real(dp) :: r_prof(npts), n_prof(npts), Te_prof(npts)
        real(dp) :: Ti_prof(npts), q_prof(npts), Er_prof(npts)
        complex(dp), allocatable :: jpar(:), jpar_e(:), jpar_i(:)
        complex(dp), allocatable :: jpar_i1(:), jpar_i2(:)
        class(kim_t), allocatable :: kim_instance
        real(dp) :: scale
        integer :: N
        logical :: ex

        call make_test_profiles(npts, r_prof, n_prof, Te_prof, Ti_prof, &
                                q_prof, Er_prof)
        call write_test_namelist('./KIM_config_periodic_species_jpar_test.nml', &
            ion_masses=ion_masses, collisionless_ions=collisionless_ions, &
            hdf5_enabled=.true.)
        nml_config_path = './KIM_config_periodic_species_jpar_test.nml'

        profiles_in_memory = .true.
        call kim_init()
        call set_profiles_from_arrays(r_prof, n_prof, Te_prof, Ti_prof, &
                                      q_prof, Er_prof, npts)
        EBdat = EBdat_t()
        call from_kim_factory_get_kim('electrostatic_periodic', kim_instance)
        call kim_instance%init()
        call kim_instance%run()

        if (.not. allocated(EBdat%jpar_e)) then
            print *, 'FAIL: periodic EBdat%jpar_e not allocated'
            error stop
        end if
        if (.not. allocated(EBdat%jpar_i)) then
            print *, 'FAIL: periodic EBdat%jpar_i not allocated'
            error stop
        end if

        call h5_obj_exists(h5id, 'fields/jpar', ex)
        if (.not. ex) error stop 'periodic total jpar dataset missing'
        call h5_obj_exists(h5id, 'fields/jpar_e', ex)
        if (.not. ex) error stop 'periodic electron jpar dataset missing'
        call h5_obj_exists(h5id, 'fields/jpar_i', ex)
        if (.not. ex) error stop 'periodic ion-sum jpar dataset missing'
        call h5_obj_exists(h5id, 'fields/jpar_i1', ex)
        if (.not. ex) error stop 'periodic first-ion jpar dataset missing'
        call h5_obj_exists(h5id, 'fields/jpar_i2', ex)
        if (.not. ex) error stop 'periodic second-ion jpar dataset missing'

        N = size(EBdat%jpar)
        allocate(jpar(N), jpar_e(N), jpar_i(N), jpar_i1(N), jpar_i2(N))
        call h5_get(h5id, 'fields/jpar', jpar)
        call h5_get(h5id, 'fields/jpar_e', jpar_e)
        call h5_get(h5id, 'fields/jpar_i', jpar_i)
        call h5_get(h5id, 'fields/jpar_i1', jpar_i1)
        call h5_get(h5id, 'fields/jpar_i2', jpar_i2)

        scale = max(1.0_dp, maxval(abs(jpar)))
        if (maxval(abs(jpar_i - jpar_i1 - jpar_i2)) > 1.0e-12_dp * scale) then
            error stop 'periodic ion species currents do not sum to jpar_i'
        end if
        if (maxval(abs(jpar - jpar_e - jpar_i)) > 1.0e-12_dp * scale) then
            error stop 'periodic species currents do not sum to total jpar'
        end if
        if (maxval(abs(EBdat%jpar_e - jpar_e)) > 1.0e-12_dp * scale .or. &
                maxval(abs(EBdat%jpar_i - jpar_i)) > 1.0e-12_dp * scale) then
            error stop 'periodic EBdat species currents differ from HDF5 output'
        end if

        call deinitialize_hdf5_output()
        if (collisionless_ions) then
            print *, 'PASS: collisionless periodic species currents sum to total jpar'
        else
            print *, 'PASS: FokkerPlanck periodic species currents sum to total jpar'
        end if
    end subroutine test_species_resolved_currents

    subroutine test_collisionless_ions_end_to_end()
        use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
        use config_m, only: profiles_in_memory, nml_config_path, &
            ion_collision_model, collisionless_kpar_epsilon
        use fields_m, only: EBdat, EBdat_t
        use kim_base_m, only: kim_t
        use kim_mod_m, only: from_kim_factory_get_kim
        use species_m, only: set_profiles_from_arrays

        integer, parameter :: npts = 201
        real(dp) :: r_prof(npts), n_prof(npts), Te_prof(npts)
        real(dp) :: Ti_prof(npts), q_prof(npts), Er_prof(npts)
        complex(dp), allocatable :: phi_with_ions(:), jpar_with_ions(:)
        complex(dp), allocatable :: phi_electrons(:), jpar_electrons(:)
        class(kim_t), allocatable :: kim_instance
        real(dp) :: current_scale
        integer :: i

        call make_test_profiles(npts, r_prof, n_prof, Te_prof, Ti_prof, &
            q_prof, Er_prof)
        call write_test_namelist('./KIM_config_periodic_collisionless_test.nml', &
            collisionless_ions=.true.)
        nml_config_path = './KIM_config_periodic_collisionless_test.nml'
        profiles_in_memory = .true.
        call kim_init()
        call set_profiles_from_arrays(r_prof, n_prof, Te_prof, Ti_prof, &
            q_prof, Er_prof, npts)
        EBdat = EBdat_t()
        call from_kim_factory_get_kim('electrostatic_periodic', kim_instance)
        call kim_instance%init()
        call kim_instance%run()

        if (trim(ion_collision_model) /= 'collisionless' .or. &
                collisionless_kpar_epsilon <= 0.0_dp) then
            print *, 'FAIL: collisionless periodic configuration was not retained'
            error stop
        end if
        phi_with_ions = EBdat%Phi
        jpar_with_ions = EBdat%jpar
        do i = 1, size(phi_with_ions)
            if (.not. ieee_is_finite(real(phi_with_ions(i), dp)) .or. &
                .not. ieee_is_finite(aimag(phi_with_ions(i))) .or. &
                .not. ieee_is_finite(real(jpar_with_ions(i), dp)) .or. &
                .not. ieee_is_finite(aimag(jpar_with_ions(i)))) then
                print *, 'FAIL: non-finite collisionless periodic output at ', i
                error stop
            end if
        end do
        if (maxval(abs(phi_with_ions)) <= 0.0_dp .or. &
                maxval(abs(jpar_with_ions)) <= 0.0_dp) then
            print *, 'FAIL: collisionless periodic Phi or jpar is identically zero'
            error stop
        end if

        ! Repeat with ions disabled.  A changed current proves that the new ion
        ! current kernels reach reconstruction rather than only charge assembly.
        deallocate(kim_instance)
        call write_test_namelist('./KIM_config_periodic_collisionless_test.nml', &
            collisionless_ions=.true., ions_disabled=.true.)
        call kim_init()
        call set_profiles_from_arrays(r_prof, n_prof, Te_prof, Ti_prof, &
            q_prof, Er_prof, npts)
        EBdat = EBdat_t()
        call from_kim_factory_get_kim('electrostatic_periodic', kim_instance)
        call kim_instance%init()
        call kim_instance%run()
        phi_electrons = EBdat%Phi
        jpar_electrons = EBdat%jpar

        current_scale = max(1.0_dp, maxval(abs(jpar_with_ions)))
        if (maxval(abs(jpar_with_ions - jpar_electrons)) <= 1.0e-10_dp * current_scale) then
            print *, 'FAIL: disabling collisionless ions did not change jpar'
            error stop
        end if
        if (maxval(abs(phi_with_ions - phi_electrons)) <= &
                1.0e-10_dp * max(1.0_dp, maxval(abs(phi_with_ions)))) then
            print *, 'FAIL: disabling collisionless ions did not change Phi'
            error stop
        end if

        print *, 'PASS: collisionless periodic ions contribute finite Phi and jpar'
    end subroutine test_collisionless_ions_end_to_end

    subroutine test_run_type_end_to_end()
        use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
        use config_m, only: profiles_in_memory, nml_config_path, &
                            periodic_dr_asis_scale, periodic_dr_tr_scale, &
                            periodic_match_global_kernel_approximations
        use flr2_fourier_kernel_m, only: kern_include_ks2, kern_zero_flr_electrons
        use species_m, only: set_profiles_from_arrays
        use kim_base_m, only: kim_t
        use kim_mod_m, only: from_kim_factory_get_kim
        use kim_resonances_m, only: r_res
        use grid_m, only: rg_grid
        use fields_m, only: EBdat

        integer, parameter :: npts = 201

        real(dp) :: r_prof(npts), n_prof(npts), Te_prof(npts)
        real(dp) :: Ti_prof(npts), q_prof(npts), Er_prof(npts)
        class(kim_t), allocatable :: kim_instance

        real(dp) :: rm, r_lo, r_hi, tol, rho_i_rm, expected_r_lo
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

        ! Capture the reference scale from the TRUE global plasma before run()
        ! redirects it onto the periodic window. This is the same state used by
        ! run_electrostatic_periodic to derive and store the window parameters.
        call prepare_resonances
        rm = r_res
        if (.not. (rm > 0.0_dp)) then
            print *, 'FAIL: resonance not found, r_res = ', rm
            error stop
        end if
        rho_i_rm = interp_species_rho_L(1, rm)
        if (.not. (rho_i_rm > 0.0_dp)) then
            print *, 'FAIL: deuterium rho_L(rm) is not positive:', rho_i_rm
            error stop
        end if

        call kim_instance%run()
        print *, 'resonant radius rm = ', rm

        if (periodic_match_global_kernel_approximations) then
            print *, 'FAIL: global-kernel approximation must default to false'
            error stop
        end if
        if (.not. kern_include_ks2 .or. kern_zero_flr_electrons) then
            print *, 'FAIL: default periodic run did not use the full Fourier kernel'
            error stop
        end if
        print *, 'PASS: omitted periodic approximation flag selects full kernel'

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

        ! (jpar) allocated, correct size, all finite, and NOT all-zero: the same
        ! constant Br that drives Phi must drive a parallel current response via
        ! thesis (11.7), j_par = K^{jPhi} Phi + K^{jB} Br.
        if (.not. allocated(EBdat%jpar)) then
            print *, 'FAIL: EBdat%jpar not allocated'
            error stop
        end if
        if (size(EBdat%jpar) /= N) then
            print *, 'FAIL: size(EBdat%jpar) /= rg_grid%npts_b; got ', &
                     size(EBdat%jpar), ' expected ', N
            error stop
        end if
        do i = 1, N
            if (.not. ieee_is_finite(real(EBdat%jpar(i), dp)) .or. &
                .not. ieee_is_finite(aimag(EBdat%jpar(i)))) then
                print *, 'FAIL: non-finite EBdat%jpar at', i
                error stop
            end if
        end do
        if (.not. (maxval(abs(EBdat%jpar)) > 0.0_dp)) then
            print *, 'FAIL: EBdat%jpar all zero (constant Br should drive a current)'
            error stop
        end if
        print *, 'PASS: EBdat%jpar allocated, size ', N, &
                 ', finite & non-zero, max|jpar| =', maxval(abs(EBdat%jpar))

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
        expected_r_lo = rm - (periodic_dr_asis_scale + periodic_dr_tr_scale) * rho_i_rm
        if (abs(r_lo - expected_r_lo) > tol) then
            print *, 'FAIL: periodic window is not sized from deuterium rho_L'
            print *, '  actual lower edge =', r_lo
            print *, '  expected lower edge =', expected_r_lo
            print *, '  deuterium rho_L(rm) =', rho_i_rm
            error stop
        end if
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

        call assert_scale_metadata(rho_i_rm)
    end subroutine test_run_type_end_to_end

    subroutine assert_scale_metadata(expected_rho)
        use constants_m, only: pi
        use config_m, only: output_path, periodic_dr_asis_scale, &
                            periodic_dr_tr_scale, periodic_kmax_scale, periodic_n_rg

        real(dp), intent(in) :: expected_rho
        character(len=512) :: metadata_path, header
        integer :: io_unit, io_status, reference_species, charge, M, n_rg
        real(dp) :: mass, rho_ref, dx_asis, dx_tr, k_max, period, tol
        logical :: exists

        metadata_path = trim(output_path)//'setup/periodic_scale.dat'
        inquire(file=trim(metadata_path), exist=exists)
        if (.not. exists) then
            print *, 'FAIL: periodic scale metadata was not stored at ', trim(metadata_path)
            error stop
        end if

        open(newunit=io_unit, file=trim(metadata_path), status='old', action='read')
        read(io_unit, '(A)', iostat=io_status) header
        if (io_status /= 0 .or. index(header, 'reference_species') == 0) then
            print *, 'FAIL: periodic scale metadata header is missing or invalid'
            error stop
        end if
        read(io_unit, *, iostat=io_status) reference_species, charge, mass, rho_ref, &
            dx_asis, dx_tr, k_max, M, n_rg
        close(io_unit)
        if (io_status /= 0) then
            print *, 'FAIL: periodic scale metadata row could not be read'
            error stop
        end if

        tol = 1.0e-10_dp * max(1.0_dp, expected_rho)
        period = 2.0_dp * (periodic_dr_asis_scale + periodic_dr_tr_scale) * expected_rho
        if (reference_species /= 1 .or. charge /= 1 .or. .not. (mass > 0.0_dp) .or. &
            abs(rho_ref - expected_rho) > tol .or. &
            abs(dx_asis - periodic_dr_asis_scale * expected_rho) > tol .or. &
            abs(dx_tr - periodic_dr_tr_scale * expected_rho) > tol .or. &
            abs(k_max - periodic_kmax_scale / expected_rho) > tol .or. &
            M /= ceiling((periodic_kmax_scale / expected_rho) * period / (2.0_dp * pi)) .or. &
            n_rg /= periodic_n_rg) then
            print *, 'FAIL: stored periodic scale metadata does not match the run'
            print *, '  species/Z/mass =', reference_species, charge, mass
            print *, '  rho stored/expected/tol =', rho_ref, expected_rho, tol
            print *, '  dx_asis stored/expected =', dx_asis, &
                     periodic_dr_asis_scale * expected_rho
            print *, '  dx_tr stored/expected =', dx_tr, &
                     periodic_dr_tr_scale * expected_rho
            print *, '  k_max stored/expected =', k_max, &
                     periodic_kmax_scale / expected_rho
            print *, '  M/N_rg stored =', M, n_rg, ' expected N_rg =', periodic_n_rg
            error stop
        end if
        print *, 'PASS: periodic scale provenance stored in ', trim(metadata_path)
    end subroutine assert_scale_metadata

    subroutine test_multi_ion_order_independence()
        integer, parameter :: masses_forward(2) = [2, 12]
        integer, parameter :: masses_reverse(2) = [12, 2]
        real(dp) :: lower_forward, lower_reverse, expected_forward, expected_reverse
        real(dp) :: tol

        call run_multi_ion_case(masses_forward, lower_forward, expected_forward)
        call run_multi_ion_case(masses_reverse, lower_reverse, expected_reverse)

        tol = 1.0e-8_dp * (1.0_dp + abs(expected_forward))
        if (abs(lower_forward - expected_forward) > tol) then
            print *, 'FAIL: forward-order window does not cover the largest ion rho_L'
            print *, '  actual lower edge =', lower_forward
            print *, '  expected lower edge =', expected_forward
            error stop
        end if
        if (abs(lower_reverse - expected_reverse) > tol) then
            print *, 'FAIL: reverse-order window does not cover the largest ion rho_L'
            print *, '  actual lower edge =', lower_reverse
            print *, '  expected lower edge =', expected_reverse
            error stop
        end if
        if (abs(lower_forward - lower_reverse) > tol) then
            print *, 'FAIL: reference scale depends on ion storage order'
            print *, '  forward lower edge =', lower_forward
            print *, '  reverse lower edge =', lower_reverse
            error stop
        end if
        print *, 'PASS: largest active-ion rho_L is storage-order independent'
    end subroutine test_multi_ion_order_independence

    subroutine test_global_approximation_enabled()
        use config_m, only: profiles_in_memory, nml_config_path, &
                            periodic_match_global_kernel_approximations
        use flr2_fourier_kernel_m, only: kern_include_ks2, kern_zero_flr_electrons
        use species_m, only: set_profiles_from_arrays
        use kim_base_m, only: kim_t
        use kim_mod_m, only: from_kim_factory_get_kim

        integer, parameter :: npts = 201
        real(dp) :: r_prof(npts), n_prof(npts), Te_prof(npts)
        real(dp) :: Ti_prof(npts), q_prof(npts), Er_prof(npts)
        class(kim_t), allocatable :: kim_instance

        call make_test_profiles(npts, r_prof, n_prof, Te_prof, Ti_prof, &
                                q_prof, Er_prof)
        call write_test_namelist('./KIM_config_periodic_run_test.nml', &
                                 match_global_approximations=.true.)
        nml_config_path = './KIM_config_periodic_run_test.nml'

        profiles_in_memory = .true.
        call kim_init()
        call set_profiles_from_arrays(r_prof, n_prof, Te_prof, Ti_prof, &
                                      q_prof, Er_prof, npts)
        call from_kim_factory_get_kim('electrostatic_periodic', kim_instance)
        call kim_instance%init()
        call kim_instance%run()

        if (.not. periodic_match_global_kernel_approximations) then
            print *, 'FAIL: explicit global-kernel approximation flag was not read'
            error stop
        end if
        if (kern_include_ks2 .or. .not. kern_zero_flr_electrons) then
            print *, 'FAIL: enabled global approximation did not configure kernels'
            error stop
        end if
        print *, 'PASS: enabled global approximation matches global kernel assumptions'
    end subroutine test_global_approximation_enabled

    subroutine run_multi_ion_case(ion_masses, lower_edge, expected_lower_edge)
        use config_m, only: profiles_in_memory, nml_config_path, &
                            periodic_dr_asis_scale, periodic_dr_tr_scale
        use species_m, only: plasma, set_profiles_from_arrays
        use kim_base_m, only: kim_t
        use kim_mod_m, only: from_kim_factory_get_kim
        use kim_resonances_m, only: r_res
        use fields_m, only: EBdat, EBdat_t

        integer, intent(in) :: ion_masses(:)
        real(dp), intent(out) :: lower_edge, expected_lower_edge
        integer, parameter :: npts = 201
        real(dp) :: r_prof(npts), n_prof(npts), Te_prof(npts)
        real(dp) :: Ti_prof(npts), q_prof(npts), Er_prof(npts)
        real(dp) :: rho_max
        class(kim_t), allocatable :: kim_instance
        integer :: sp

        call make_test_profiles(npts, r_prof, n_prof, Te_prof, Ti_prof, &
                                q_prof, Er_prof)
        call write_test_namelist('./KIM_config_periodic_run_test.nml', ion_masses)
        nml_config_path = './KIM_config_periodic_run_test.nml'

        profiles_in_memory = .true.
        call kim_init()
        call set_profiles_from_arrays(r_prof, n_prof, Te_prof, Ti_prof, &
                                      q_prof, Er_prof, npts)
        EBdat = EBdat_t()
        call from_kim_factory_get_kim('electrostatic_periodic', kim_instance)
        call kim_instance%init()
        call kim_instance%run()

        rho_max = 0.0_dp
        do sp = 1, plasma%n_species - 1
            rho_max = max(rho_max, interp_species_rho_L(sp, r_res))
        end do
        lower_edge = EBdat%r_grid(1)
        expected_lower_edge = r_res &
            - (periodic_dr_asis_scale + periodic_dr_tr_scale) * rho_max
    end subroutine run_multi_ion_case

    real(dp) function interp_species_rho_L(sp, x) result(rhoLx)
        use species_m, only: plasma
        integer, intent(in) :: sp
        real(dp), intent(in) :: x
        integer, parameter :: nlagr = 4, nder = 0
        real(dp) :: coef(0:nder, nlagr)
        integer :: gs, ir, ibeg, iend

        gs = plasma%grid_size
        call binsrc(plasma%r_grid, 1, gs, x, ir)
        ibeg = max(1, ir - nlagr / 2)
        iend = ibeg + nlagr - 1
        if (iend > gs) then
            iend = gs
            ibeg = iend - nlagr + 1
        end if
        call plag_coeff(nlagr, nder, x, plasma%r_grid(ibeg:iend), coef)
        rhoLx = sum(coef(0, :) * plasma%spec(sp)%rho_L(ibeg:iend))
    end function interp_species_rho_L

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

    subroutine write_test_namelist(path, ion_masses, match_global_approximations, &
            collisionless_ions, ions_disabled, hdf5_enabled)
        ! Minimal electrostatic-periodic FokkerPlanck configuration; m_mode = -6,
        ! n_mode = 2 makes q resonant at q = 3, type_br_field = 12 (constant Br).
        ! Deliberately OMITS the &KIM_PERIODIC group to prove the periodic_*
        ! config defaults apply when the optional namelist is absent.
        character(len=*), intent(in) :: path
        integer, intent(in), optional :: ion_masses(:)
        logical, intent(in), optional :: match_global_approximations
        logical, intent(in), optional :: collisionless_ions, ions_disabled, hdf5_enabled
        integer :: iunit, i, nions
        logical :: custom_species, use_collisionless, disable_ions, write_hdf5

        custom_species = present(ion_masses)
        use_collisionless = .false.
        if (present(collisionless_ions)) use_collisionless = collisionless_ions
        disable_ions = .false.
        if (present(ions_disabled)) disable_ions = ions_disabled
        write_hdf5 = .false.
        if (present(hdf5_enabled)) write_hdf5 = hdf5_enabled
        nions = 1
        if (custom_species) nions = size(ion_masses)

        open(newunit=iunit, file=path, status='replace', action='write')
        write(iunit, '(A)') '&KIM_CONFIG'
        write(iunit, '(A,I0)') ' number_of_ion_species = ', nions
        write(iunit, '(A)') ' artificial_debye_case = 0'
        write(iunit, '(A)') " type_of_run = 'electrostatic_periodic'"
        write(iunit, '(A)') " collision_model = 'FokkerPlanck'"
        if (use_collisionless) then
            write(iunit, '(A)') " ion_collision_model = 'collisionless'"
            write(iunit, '(A)') ' collisionless_kpar_epsilon = 1.0e-5'
        else
            write(iunit, '(A)') " ion_collision_model = 'FokkerPlanck'"
        end if
        write(iunit, '(A,L1)') ' read_species_from_namelist = ', custom_species
        write(iunit, '(A,L1)') ' turn_off_ions = ', disable_ions
        write(iunit, '(A)') ' turn_off_electrons = .false.'
        write(iunit, '(A)') " plasma_type = 'D'"
        write(iunit, '(A)') ' rescale_density = .false.'
        write(iunit, '(A)') ' number_density_rescale = 1.0'
        write(iunit, '(A)') ' ion_flr_scale_factor = 1.0'
        write(iunit, '(A)') '/'
        if (custom_species) then
            write(iunit, '(A)') '&KIM_species'
            write(iunit, '(A)', advance='no') ' ai ='
            do i = 1, nions
                write(iunit, '(1X,I0)', advance='no') ion_masses(i)
            end do
            write(iunit, *)
            write(iunit, '(A)', advance='no') ' zi ='
            do i = 1, nions
                write(iunit, '(1X,I0)', advance='no') 1
            end do
            write(iunit, *)
            write(iunit, '(A)') '/'
        end if
        write(iunit, '(A)') '&WKB_DISPERSION'
        write(iunit, '(A)') '/'
        write(iunit, '(A)') '&KIM_IO'
        write(iunit, '(A)') " profile_location = './'"
        write(iunit, '(A)') " output_path = './out_periodic_run_test/'"
        write(iunit, '(A)') ' hdf5_input = .false.'
        write(iunit, '(A,L1)') ' hdf5_output = ', write_hdf5
        write(iunit, '(A)') ' log_level = 3'
        write(iunit, '(A)') ' data_verbosity = 0'
        write(iunit, '(A)') ' calculate_asymptotics = .false.'
        write(iunit, '(A)') " h5_out_file = 'KIM_species_jpar.h5'"
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
        if (present(match_global_approximations)) then
            write(iunit, '(A)') '&KIM_PERIODIC'
            write(iunit, '(A,L1)') &
                ' periodic_match_global_kernel_approximations = ', &
                match_global_approximations
            write(iunit, '(A)') '/'
        end if
        close(iunit)
    end subroutine write_test_namelist

end program test_kim_solver_periodic
