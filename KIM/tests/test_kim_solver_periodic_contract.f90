program test_kim_solver_periodic_contract
    ! End-to-end contract test for the periodic KIM response returned through
    ! kim_solver_t.  This is intentionally structural: it proves that the
    ! local periodic result contains the same physical field families needed
    ! by the QL-Balance adapter, without freezing a platform-dependent golden.
    use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
    use KIM_kinds_m, only: dp
    use kim_resonances_m, only: locate_periodic_resonance, &
        KIM_RESONANCE_OK, KIM_RESONANCE_NOT_FOUND, &
        KIM_RESONANCE_AMBIGUOUS
    use kim_solver_m, only: kim_solver_t, kim_results_t, kim_profiles_t, KIM_OK

    integer, parameter :: npts = 40
    integer, parameter :: m_mode = -6, n_mode = 2
    integer, parameter :: second_m_mode = -9, second_n_mode = 3
    complex(dp), parameter :: bparallel_only_drive = cmplx(0.25_dp, -0.5_dp, dp)
    complex(dp), parameter :: namelist_Br_drive = cmplx(1.0_dp, 0.0_dp, dp)
    complex(dp), parameter :: zero_drive = cmplx(0.0_dp, 0.0_dp, dp)
    real(dp), parameter :: forced_resonance = 37.0_dp
    real(dp), parameter :: profile_resonance = 36.8_dp
    type(kim_solver_t) :: kim
    type(kim_results_t) :: res
    type(kim_results_t) :: res_first
    type(kim_profiles_t) :: prof
    integer :: ierr, i
    real(dp) :: frac
    logical :: all_passed

    all_passed = .true.
    call test_signed_resonance_locator(all_passed)

    allocate(prof%r(npts), prof%n(npts), prof%Te(npts), prof%Ti(npts), &
             prof%q(npts), prof%Er(npts))
    do i = 1, npts
        prof%r(i) = 2.0_dp + 58.0_dp * real(i - 1, dp) / real(npts - 1, dp)
        frac = (prof%r(i) - 2.0_dp) / 58.0_dp
        prof%n(i) = 2.0e13_dp * (1.1_dp - prof%r(i) / 100.0_dp)
        prof%Te(i) = 1.0e3_dp * (1.2_dp - prof%r(i) / 100.0_dp)
        prof%Ti(i) = prof%Te(i)
        prof%q(i) = 1.5_dp + 2.5_dp * frac
        prof%Er(i) = -0.5_dp
    end do

    call kim%init('KIM_config_em_small.nml', run_type='electrostatic_periodic', &
                  profiles=prof, stat=ierr)
    call check('periodic init returns KIM_OK', ierr == KIM_OK, all_passed)
    if (ierr == KIM_OK) then
        call kim%solve(m=m_mode, n=n_mode, stat=ierr, Br_drive=zero_drive, &
                       Bparallel_drive=bparallel_only_drive, &
                       resonance_radius=forced_resonance)
        call check('first solve with scoped overrides returns KIM_OK', &
                   ierr == KIM_OK, all_passed)
    end if
    if (ierr == KIM_OK) then
        res = kim%results()
        call check('first solve honors prescribed resonance', &
                   res%r_resonance == forced_resonance, all_passed)
        call check_constant_field('first solve honors Br override', res%Br, &
                                  zero_drive, all_passed)
        call check_constant_field('first solve honors Bparallel override', &
                                  res%Bparallel, bparallel_only_drive, all_passed)
        call kim%set_profiles(prof, stat=ierr)
        call check('profile reset after first override returns KIM_OK', &
                   ierr == KIM_OK, all_passed)
    end if
    if (ierr == KIM_OK) then
        call kim%solve(m=m_mode, n=n_mode, stat=ierr)
        call check('periodic solve returns KIM_OK', ierr == KIM_OK, all_passed)
    end if

    if (ierr == KIM_OK) then
        res = kim%results()
        call check('local field grid populated', allocated(res%r_field), all_passed)
        if (allocated(res%r_field)) then
            call check('local field grid has useful size', size(res%r_field) > 1, all_passed)
        end if
        call check_complex_field('Phi', res%Phi, res%r_field, all_passed)
        call check_complex_field('Br', res%Br, res%r_field, all_passed)
        call check_complex_field('Bparallel', res%Bparallel, res%r_field, all_passed)
        call check_complex_field('Er', res%Er, res%r_field, all_passed)
        call check_complex_field('Etheta', res%Etheta, res%r_field, all_passed)
        call check_complex_field('Ez', res%Ez, res%r_field, all_passed)
        call check_complex_field('Es', res%Es, res%r_field, all_passed)
        call check_complex_field('Ep', res%Ep, res%r_field, all_passed)
        call check_complex_field('jpar', res%jpar, res%r_field, all_passed)
        call check_complex_field('jpar_e', res%jpar_e, res%r_field, all_passed)
        call check_complex_field('jpar_i', res%jpar_i, res%r_field, all_passed)
        call check('result mode matches request', res%m == m_mode .and. res%n == n_mode, all_passed)
        res_first = res
        call kim%set_profiles(prof, stat=ierr)
        call check('profile reset before second mode returns KIM_OK', &
                   ierr == KIM_OK, all_passed)
        call kim%solve(m=second_m_mode, n=second_n_mode, stat=ierr)
        call check('second periodic mode returns KIM_OK', ierr == KIM_OK, all_passed)
        call kim%set_profiles(prof, stat=ierr)
        call check('profile reset before repeated mode returns KIM_OK', &
                   ierr == KIM_OK, all_passed)
        call kim%solve(m=m_mode, n=n_mode, stat=ierr)
        call check('periodic mode can be repeated after another mode', ierr == KIM_OK, all_passed)
        res = kim%results()
        call check_finite_field('repeated periodic result remains finite', &
                                res%Phi, all_passed)
    end if

    if (ierr == KIM_OK) then
        call kim%solve(m=m_mode, n=n_mode, stat=ierr, Br_drive=zero_drive, &
                       Bparallel_drive=bparallel_only_drive, &
                       resonance_radius=forced_resonance)
        call check('Bparallel-only solve returns KIM_OK', ierr == KIM_OK, all_passed)
    end if
    if (ierr == KIM_OK) then
        res = kim%results()
        call check_constant_field('Bparallel-only Br is exactly zero', res%Br, &
                                  zero_drive, all_passed)
        call check_constant_field('prescribed Bparallel is returned exactly', &
                                  res%Bparallel, bparallel_only_drive, all_passed)
        call check('prescribed resonance is returned exactly', &
                   res%r_resonance == forced_resonance, all_passed)
        call check_nonzero_field('Bparallel-only potential response is nonzero', &
                                 res%Phi, all_passed)
        call check_nonzero_field('Bparallel-only current response is nonzero', &
                                 res%jpar, all_passed)

        ! Restore the independently specified profile before checking that the
        ! solve-scoped resonance override did not leak into the next solve.
        call kim%set_profiles(prof, stat=ierr)
        call check('profile reset before no-override solve returns KIM_OK', &
                   ierr == KIM_OK, all_passed)
        call kim%solve(m=m_mode, n=n_mode, stat=ierr)
        call check('no-override solve returns KIM_OK', ierr == KIM_OK, all_passed)
    end if
    if (ierr == KIM_OK) then
        res = kim%results()
        call check_constant_field('no-override Br restores the namelist drive', &
                                  res%Br, namelist_Br_drive, all_passed)
        call check_constant_field('no-override Bparallel restores exact zero', &
                                  res%Bparallel, zero_drive, all_passed)
        call check('no-override solve restores profile resonance selection', &
                   abs(res%r_resonance - profile_resonance) < 1.0e-12_dp, &
                   all_passed)
    end if

    call kim%finalize()
    call test_mode_order_isolation(prof, all_passed)
    call test_profile_update_grid_restore(prof, all_passed)
    call test_adaptive_grid_refresh(prof, all_passed)
    call test_file_profile_lifecycle(prof, all_passed)
    if (all_passed) then
        print *, 'All periodic KIM result-contract checks PASSED'
        stop 0
    else
        print *, 'Periodic KIM result-contract checks FAILED'
        stop 1
    end if

contains

    subroutine test_mode_order_isolation(profiles, passed)
        use setup_m, only: setup_m_mode => m_mode, setup_n_mode => n_mode

        type(kim_profiles_t), intent(in) :: profiles
        logical, intent(inout) :: passed
        integer, parameter :: other_m_mode = -7, other_n_mode = 2
        real(dp), parameter :: other_resonance = 48.4_dp
        type(kim_solver_t) :: ordered, fresh
        type(kim_results_t) :: ordered_first, ordered_second
        type(kim_results_t) :: fresh_first, fresh_second
        integer :: status

        setup_m_mode = m_mode
        setup_n_mode = n_mode
        call ordered%init('KIM_config_em_small.nml', &
            run_type='electrostatic_periodic', profiles=profiles, stat=status)
        call check('ordered-mode init returns KIM_OK', status == KIM_OK, passed)
        if (status == KIM_OK) then
            call ordered%solve(m_mode, n_mode, stat=status, Br_drive=zero_drive, &
                Bparallel_drive=bparallel_only_drive)
            call check('ordered first mode returns KIM_OK', status == KIM_OK, passed)
        end if
        if (status == KIM_OK) ordered_first = ordered%results()
        if (status == KIM_OK) then
            ! Deliberately do not call set_profiles here. The handle must restore
            ! the authoritative global profile before moving to a distant mode.
            call ordered%solve(other_m_mode, other_n_mode, stat=status, &
                Br_drive=zero_drive, Bparallel_drive=bparallel_only_drive, &
                resonance_radius=other_resonance)
            call check('ordered distant second mode returns KIM_OK', &
                status == KIM_OK, passed)
        end if
        if (status == KIM_OK) ordered_second = ordered%results()
        call ordered%finalize()

        setup_m_mode = m_mode
        setup_n_mode = n_mode
        call fresh%init('KIM_config_em_small.nml', &
            run_type='electrostatic_periodic', profiles=profiles, stat=status)
        call check('fresh first-mode init returns KIM_OK', status == KIM_OK, passed)
        if (status == KIM_OK) then
            call fresh%solve(m_mode, n_mode, stat=status, Br_drive=zero_drive, &
                Bparallel_drive=bparallel_only_drive)
            call check('fresh first mode returns KIM_OK', status == KIM_OK, passed)
        end if
        if (status == KIM_OK) then
            fresh_first = fresh%results()
            call check_matching_results('ordered first mode matches fresh solve', &
                ordered_first, fresh_first, passed)
        end if
        call fresh%finalize()

        setup_m_mode = other_m_mode
        setup_n_mode = other_n_mode
        call fresh%init('KIM_config_em_small.nml', &
            run_type='electrostatic_periodic', profiles=profiles, stat=status)
        call check('fresh second-mode init returns KIM_OK', status == KIM_OK, passed)
        if (status == KIM_OK) then
            call fresh%solve(other_m_mode, other_n_mode, stat=status, &
                Br_drive=zero_drive, Bparallel_drive=bparallel_only_drive, &
                resonance_radius=other_resonance)
            call check('fresh second mode returns KIM_OK', status == KIM_OK, passed)
        end if
        if (status == KIM_OK) then
            fresh_second = fresh%results()
            call check_matching_results('ordered second mode matches fresh solve', &
                ordered_second, fresh_second, passed)
        end if
        call fresh%finalize()

        setup_m_mode = m_mode
        setup_n_mode = n_mode
    end subroutine test_mode_order_isolation

    subroutine test_profile_update_grid_restore(profiles, passed)
        use setup_m, only: setup_m_mode => m_mode, setup_n_mode => n_mode

        type(kim_profiles_t), intent(in) :: profiles
        logical, intent(inout) :: passed
        type(kim_profiles_t) :: updated_profiles
        type(kim_solver_t) :: updated, fresh
        type(kim_results_t) :: updated_result, fresh_result
        integer :: status

        updated_profiles = profiles
        updated_profiles%n = 1.1_dp*profiles%n
        updated_profiles%Te = 0.9_dp*profiles%Te

        setup_m_mode = m_mode
        setup_n_mode = n_mode
        call updated%init('KIM_config_profile_update_large_rg.nml', &
            run_type='electrostatic_periodic', profiles=profiles, stat=status)
        call check('large-grid profile-update init returns KIM_OK', &
            status == KIM_OK, passed)
        if (status == KIM_OK) then
            call updated%solve(m_mode, n_mode, stat=status, Br_drive=zero_drive, &
                Bparallel_drive=bparallel_only_drive, &
                resonance_radius=profile_resonance)
            call check('large-grid pre-update solve returns KIM_OK', &
                status == KIM_OK, passed)
        end if
        if (status == KIM_OK) then
            call updated%set_profiles(updated_profiles, stat=status)
            call check('large-grid profile update returns KIM_OK', &
                status == KIM_OK, passed)
        end if
        if (status == KIM_OK) then
            call updated%solve(m_mode, n_mode, stat=status, Br_drive=zero_drive, &
                Bparallel_drive=bparallel_only_drive, &
                resonance_radius=profile_resonance)
            call check('large-grid post-update solve returns KIM_OK', &
                status == KIM_OK, passed)
        end if
        if (status == KIM_OK) updated_result = updated%results()
        call updated%finalize()

        setup_m_mode = m_mode
        setup_n_mode = n_mode
        call fresh%init('KIM_config_profile_update_large_rg.nml', &
            run_type='electrostatic_periodic', profiles=updated_profiles, &
            stat=status)
        call check('fresh large-grid updated init returns KIM_OK', &
            status == KIM_OK, passed)
        if (status == KIM_OK) then
            call fresh%solve(m_mode, n_mode, stat=status, Br_drive=zero_drive, &
                Bparallel_drive=bparallel_only_drive, &
                resonance_radius=profile_resonance)
            call check('fresh large-grid updated solve returns KIM_OK', &
                status == KIM_OK, passed)
        end if
        if (status == KIM_OK) then
            fresh_result = fresh%results()
            call check_matching_results('profile update matches fresh large-grid solve', &
                updated_result, fresh_result, passed)
        end if
        call fresh%finalize()
    end subroutine test_profile_update_grid_restore

    subroutine test_adaptive_grid_refresh(profiles, passed)
        type(kim_profiles_t), intent(in) :: profiles
        logical, intent(inout) :: passed
        integer, parameter :: other_m_mode = -7, other_n_mode = 2
        type(kim_profiles_t) :: shifted_profiles
        type(kim_solver_t) :: ordered, fresh
        type(kim_results_t) :: ordered_result, fresh_result
        real(dp) :: shifted_resonance
        integer :: locator_status, status

        call ordered%init('KIM_config_adaptive_mode_a.nml', &
            run_type='electrostatic_periodic', profiles=profiles, stat=status)
        call check('adaptive ordered init returns KIM_OK', &
            status == KIM_OK, passed)
        if (status == KIM_OK) then
            call ordered%solve(m_mode, n_mode, stat=status, Br_drive=zero_drive, &
                Bparallel_drive=bparallel_only_drive, &
                resonance_radius=profile_resonance)
            call check('adaptive ordered first mode returns KIM_OK', &
                status == KIM_OK, passed)
        end if
        if (status == KIM_OK) then
            call ordered%solve(other_m_mode, other_n_mode, stat=status, &
                Br_drive=zero_drive, Bparallel_drive=bparallel_only_drive)
            call check('adaptive ordered second mode returns KIM_OK', &
                status == KIM_OK, passed)
        end if
        if (status == KIM_OK) ordered_result = ordered%results()
        call ordered%finalize()

        call fresh%init('KIM_config_adaptive_mode_b.nml', &
            run_type='electrostatic_periodic', profiles=profiles, stat=status)
        call check('adaptive fresh second-mode init returns KIM_OK', &
            status == KIM_OK, passed)
        if (status == KIM_OK) then
            call fresh%solve(other_m_mode, other_n_mode, stat=status, &
                Br_drive=zero_drive, Bparallel_drive=bparallel_only_drive)
            call check('adaptive fresh second mode returns KIM_OK', &
                status == KIM_OK, passed)
        end if
        if (status == KIM_OK) then
            fresh_result = fresh%results()
            call check_matching_results(&
                'adaptive ordered second mode matches fresh mode-centered grid', &
                ordered_result, fresh_result, passed)
        end if
        call fresh%finalize()

        shifted_profiles = profiles
        shifted_profiles%q = profiles%q + 0.35_dp * &
            (profiles%r - profiles%r(1)) / &
            (profiles%r(size(profiles%r)) - profiles%r(1))
        call locate_periodic_resonance(shifted_profiles%r, shifted_profiles%q, &
            m_mode, n_mode, shifted_resonance, locator_status)
        call check('shifted adaptive profile has one signed resonance', &
            locator_status == KIM_RESONANCE_OK, passed)

        call ordered%init('KIM_config_adaptive_mode_a.nml', &
            run_type='electrostatic_periodic', profiles=profiles, stat=status)
        call check('adaptive profile-update init returns KIM_OK', &
            status == KIM_OK, passed)
        if (status == KIM_OK) then
            call ordered%solve(m_mode, n_mode, stat=status, Br_drive=zero_drive, &
                Bparallel_drive=bparallel_only_drive)
        end if
        if (status == KIM_OK) then
            call ordered%set_profiles(shifted_profiles, stat=status)
            call check('adaptive shifted profile update returns KIM_OK', &
                status == KIM_OK, passed)
        end if
        if (status == KIM_OK) then
            call ordered%solve(m_mode, n_mode, stat=status, Br_drive=zero_drive, &
                Bparallel_drive=bparallel_only_drive)
            call check('adaptive shifted-profile solve returns KIM_OK', &
                status == KIM_OK, passed)
        end if
        if (status == KIM_OK) ordered_result = ordered%results()
        call ordered%finalize()

        call fresh%init('KIM_config_adaptive_mode_a.nml', &
            run_type='electrostatic_periodic', profiles=shifted_profiles, &
            stat=status)
        call check('adaptive fresh shifted-profile init returns KIM_OK', &
            status == KIM_OK, passed)
        if (status == KIM_OK) then
            call fresh%solve(m_mode, n_mode, stat=status, Br_drive=zero_drive, &
                Bparallel_drive=bparallel_only_drive)
            call check('adaptive fresh shifted-profile solve returns KIM_OK', &
                status == KIM_OK, passed)
        end if
        if (status == KIM_OK) then
            fresh_result = fresh%results()
            call check_matching_results(&
                'adaptive shifted profile matches fresh resonance-centered grid', &
                ordered_result, fresh_result, passed)
        end if
        call fresh%finalize()
    end subroutine test_adaptive_grid_refresh

    subroutine test_file_profile_lifecycle(profiles, passed)
        use config_m, only: profiles_in_memory
        use kim_resonances_m, only: r_res
        use species_m, only: plasma

        type(kim_profiles_t), intent(in) :: profiles
        logical, intent(inout) :: passed
        type(kim_profiles_t) :: memory_profiles
        type(kim_solver_t) :: memory_handle, file_handle
        type(kim_results_t) :: file_result
        real(dp), parameter :: memory_resonance = 48.4_dp
        real(dp) :: located_resonance, radial_fraction
        integer :: locator_status, status

        memory_profiles = profiles
        radial_fraction = memory_resonance - memory_profiles%r(1)
        radial_fraction = radial_fraction / &
            (memory_profiles%r(size(memory_profiles%r)) - memory_profiles%r(1))
        memory_profiles%q = 1.0_dp + &
            (2.0_dp / radial_fraction) * &
            (memory_profiles%r - memory_profiles%r(1)) / &
            (memory_profiles%r(size(memory_profiles%r)) - memory_profiles%r(1))
        call locate_periodic_resonance(memory_profiles%r, memory_profiles%q, &
            m_mode, n_mode, located_resonance, locator_status)
        call check('in-memory adaptive fixture has its independent resonance', &
            locator_status == KIM_RESONANCE_OK .and. &
            abs(located_resonance - memory_resonance) < 1.0e-12_dp, passed)

        call memory_handle%init('KIM_config_adaptive_mode_a.nml', &
            run_type='electrostatic_periodic', profiles=memory_profiles, &
            stat=status)
        call check('memory-to-file lifecycle in-memory init returns KIM_OK', &
            status == KIM_OK, passed)
        call check('in-memory adaptive grid uses its profile resonance', &
            abs(r_res - memory_resonance) < 1.0e-12_dp, passed)
        call memory_handle%finalize()

        call file_handle%init('KIM_config_file_profiles_adaptive.nml', &
            run_type='electrostatic_periodic', stat=status)
        call check('file-backed init after in-memory handle returns KIM_OK', &
            status == KIM_OK, passed)
        call check('file-backed init selects the file profile source', &
            .not. profiles_in_memory, passed)
        if (status == KIM_OK) then
            call locate_periodic_resonance(plasma%r_grid, plasma%q, m_mode, &
                n_mode, located_resonance, locator_status)
            call check('file-backed init loads the independent q profile', &
                locator_status == KIM_RESONANCE_OK .and. &
                abs(located_resonance - profile_resonance) < 1.0e-4_dp, passed)
            call check('file-backed adaptive grid refreshes its cached resonance', &
                abs(r_res - profile_resonance) < 1.0e-4_dp, passed)
            call file_handle%solve(m_mode, n_mode, stat=status, &
                Br_drive=zero_drive, Bparallel_drive=zero_drive)
            call check('file-backed solve after in-memory handle returns KIM_OK', &
                status == KIM_OK, passed)
        end if
        if (status == KIM_OK) then
            file_result = file_handle%results()
            call check('file-backed solve uses its profile resonance', &
                abs(file_result%r_resonance - profile_resonance) < 1.0e-4_dp, &
                passed)
        end if
        call file_handle%finalize()

        call file_handle%init('KIM_config_file_profiles_adaptive.nml', &
            run_type='electrostatic_periodic', stat=status)
        call check('repeated file-backed init returns KIM_OK', &
            status == KIM_OK, passed)
        if (status == KIM_OK) then
            call locate_periodic_resonance(plasma%r_grid, plasma%q, m_mode, &
                n_mode, located_resonance, locator_status)
            call check('repeated file-backed init reloads the q profile', &
                locator_status == KIM_RESONANCE_OK .and. &
                abs(located_resonance - profile_resonance) < 1.0e-4_dp, passed)
            call check('repeated file-backed adaptive grid refreshes resonance', &
                abs(r_res - profile_resonance) < 1.0e-4_dp, passed)
        end if
        call file_handle%finalize()
    end subroutine test_file_profile_lifecycle

    subroutine check_matching_results(name, actual, expected, passed)
        character(*), intent(in) :: name
        type(kim_results_t), intent(in) :: actual, expected
        logical, intent(inout) :: passed
        logical :: ok
        real(dp) :: scale

        ok = actual%m == expected%m .and. actual%n == expected%n
        ok = ok .and. actual%r_resonance == expected%r_resonance
        ok = ok .and. allocated(actual%r_field) .and. allocated(expected%r_field)
        ok = ok .and. allocated(actual%kp) .and. allocated(expected%kp)
        ok = ok .and. allocated(actual%ks) .and. allocated(expected%ks)
        ok = ok .and. allocated(actual%Phi) .and. allocated(expected%Phi)
        ok = ok .and. allocated(actual%jpar) .and. allocated(expected%jpar)
        if (ok) then
            ok = size(actual%r_field) == size(expected%r_field) .and. &
                size(actual%kp) == size(expected%kp) .and. &
                size(actual%ks) == size(expected%ks) .and. &
                size(actual%Phi) == size(expected%Phi) .and. &
                size(actual%jpar) == size(expected%jpar)
        end if
        if (ok) then
            scale = max(1.0_dp, maxval(abs(expected%r_field)))
            ok = maxval(abs(actual%r_field - expected%r_field)) <= 1.0e-11_dp*scale
        end if
        if (ok) then
            scale = max(1.0_dp, maxval(abs(expected%kp)), maxval(abs(expected%ks)))
            ok = maxval(abs(actual%kp - expected%kp)) <= 1.0e-11_dp*scale .and. &
                maxval(abs(actual%ks - expected%ks)) <= 1.0e-11_dp*scale
        end if
        if (ok) then
            scale = max(1.0_dp, maxval(abs(expected%Phi)), maxval(abs(expected%jpar)))
            ok = maxval(abs(actual%Phi - expected%Phi)) <= 1.0e-10_dp*scale .and. &
                maxval(abs(actual%jpar - expected%jpar)) <= 1.0e-10_dp*scale
        end if
        call check(name, ok, passed)
    end subroutine check_matching_results

    subroutine test_signed_resonance_locator(passed)
        logical, intent(inout) :: passed
        real(dp), parameter :: radius(3) = [1.0_dp, 2.0_dp, 3.0_dp]
        real(dp), parameter :: q_increasing(3) = [2.0_dp, 3.0_dp, 4.0_dp]
        real(dp), parameter :: q_decreasing(3) = [4.0_dp, 3.0_dp, 2.0_dp]
        real(dp), parameter :: q_negative(3) = [-2.0_dp, -3.0_dp, -4.0_dp]
        real(dp), parameter :: q_multiple(3) = [2.0_dp, 4.0_dp, 2.0_dp]
        real(dp) :: resonance
        integer :: status

        call locate_periodic_resonance(radius, q_increasing, -6, 2, &
            resonance, status)
        call check('signed locator finds q=-m/n on increasing profile', &
                   status == KIM_RESONANCE_OK .and. resonance == 2.0_dp, passed)

        call locate_periodic_resonance(radius, q_decreasing, -6, 2, &
            resonance, status)
        call check('signed locator finds q=-m/n on decreasing profile', &
                   status == KIM_RESONANCE_OK .and. resonance == 2.0_dp, passed)

        call locate_periodic_resonance(radius, q_negative, 6, 2, &
            resonance, status)
        call check('signed locator retains a negative-q resonance', &
                   status == KIM_RESONANCE_OK .and. resonance == 2.0_dp, passed)

        call locate_periodic_resonance(radius, q_increasing, 6, 2, &
            resonance, status)
        call check('wrong-sign mode does not create an absolute-value resonance', &
                   status == KIM_RESONANCE_NOT_FOUND, passed)

        call locate_periodic_resonance(radius, q_multiple, -6, 2, &
            resonance, status)
        call check('multiple signed crossings are rejected as ambiguous', &
                   status == KIM_RESONANCE_AMBIGUOUS, passed)
    end subroutine test_signed_resonance_locator

    subroutine check(name, ok, passed)
        character(*), intent(in) :: name
        logical, intent(in) :: ok
        logical, intent(inout) :: passed
        if (ok) then
            print *, 'PASS: ', name
        else
            print *, 'FAIL: ', name
            passed = .false.
        end if
    end subroutine check

    subroutine check_complex_field(name, field, grid, passed)
        character(*), intent(in) :: name
        complex(dp), allocatable, intent(in) :: field(:)
        real(dp), allocatable, intent(in) :: grid(:)
        logical, intent(inout) :: passed
        logical :: ok
        ok = allocated(field) .and. allocated(grid)
        if (ok) ok = size(field) == size(grid)
        if (ok) ok = all(ieee_is_finite(real(field, dp))) .and. &
                     all(ieee_is_finite(aimag(field)))
        call check('periodic '//trim(name)//' populated, sized, finite', ok, passed)
    end subroutine check_complex_field

    subroutine check_finite_field(name, field, passed)
        character(*), intent(in) :: name
        complex(dp), allocatable, intent(in) :: field(:)
        logical, intent(inout) :: passed
        logical :: ok

        ok = allocated(field)
        if (ok) then
            ok = all(ieee_is_finite(real(field, dp))) .and. &
                 all(ieee_is_finite(aimag(field)))
        end if
        call check(name, ok, passed)
    end subroutine check_finite_field

    subroutine check_constant_field(name, field, expected, passed)
        character(*), intent(in) :: name
        complex(dp), allocatable, intent(in) :: field(:)
        complex(dp), intent(in) :: expected
        logical, intent(inout) :: passed
        logical :: ok

        ok = allocated(field)
        if (ok) ok = all(field == expected)
        call check(name, ok, passed)
    end subroutine check_constant_field

    subroutine check_nonzero_field(name, field, passed)
        character(*), intent(in) :: name
        complex(dp), allocatable, intent(in) :: field(:)
        logical, intent(inout) :: passed
        logical :: ok

        ok = allocated(field)
        if (ok) then
            ok = all(ieee_is_finite(real(field, dp))) .and. &
                 all(ieee_is_finite(aimag(field)))
        end if
        if (ok) ok = any(abs(field) > 0.0_dp)
        call check(name, ok, passed)
    end subroutine check_nonzero_field

end program test_kim_solver_periodic_contract
