program test_bparallel_linear_response
    use KIM_kinds_m, only: dp
    use constants_m, only: com_unit, pi, sol
    use config_m, only: artificial_debye_case, collision_model, &
        ion_collision_model, turn_off_electrons, turn_off_ions
    use species_m, only: plasma_t
    use flr2_fourier_kernel_m, only: core_rho_Bparallel_sp, &
        core_j_Bparallel_sp, kern_include_ks2
    use periodic_assembly_m, only: assemble_periodic_bparallel_matrices
    use grid_m, only: rg_grid
    use setup_m, only: collisions_off, mphi_max
    use rt_electrostatic_periodic_m, only: periodic_bparallel_support_status, &
        PERIODIC_BPAR_OK, PERIODIC_BPAR_UNSUPPORTED_COLLISION, &
        PERIODIC_BPAR_UNSUPPORTED_HARMONIC, PERIODIC_BPAR_UNSUPPORTED_APPROXIMATION
    use rt_electrostatic_periodic_m, only: PERIODIC_BPAR_UNSUPPORTED_DEBYE_MODEL

    implicit none

    integer, parameter :: ngrid = 3
    type(plasma_t) :: plasma
    character(len=32) :: test_mode

    call get_command_argument(1, test_mode)
    if (len_trim(test_mode) > 0) then
        call run_rejection_case(trim(test_mode))
        stop 0
    end if

    call populate_fixture(plasma)
    call test_direct_gyroaverage_oracle(plasma)
    call test_assembled_offdiagonal_oracle(plasma)
    call test_long_wavelength_sign_phase(plasma)
    call test_species_sign_and_moment_routes(plasma)
    call test_support_validation()

    print *, 'Bparallel linear-response tests PASSED'

contains

    subroutine populate_fixture(plasma)
        type(plasma_t), intent(out) :: plasma
        integer :: sp

        plasma%n_species = 2
        plasma%grid_size = ngrid
        allocate(plasma%ks(ngrid), plasma%spec(0:1))
        plasma%ks = [0.7_dp, 0.55_dp, 0.9_dp]

        do sp = 0, 1
            allocate(plasma%spec(sp)%lambda_D(ngrid), plasma%spec(sp)%vT(ngrid))
            allocate(plasma%spec(sp)%nu(ngrid), plasma%spec(sp)%omega_c(ngrid))
            allocate(plasma%spec(sp)%rho_L(ngrid), plasma%spec(sp)%A1(ngrid))
            allocate(plasma%spec(sp)%A2(ngrid))
            allocate(plasma%spec(sp)%I00(ngrid, 0:0), plasma%spec(sp)%I02(ngrid, 0:0))
            allocate(plasma%spec(sp)%I10(ngrid, 0:0), plasma%spec(sp)%I12(ngrid, 0:0))
            allocate(plasma%spec(sp)%I20(ngrid, 0:0), plasma%spec(sp)%I21(ngrid, 0:0))

            plasma%spec(sp)%lambda_D = 2.3_dp
            plasma%spec(sp)%vT = 3.1_dp
            plasma%spec(sp)%nu = 4.2_dp
            plasma%spec(sp)%rho_L = 0.37_dp
            plasma%spec(sp)%A1 = 0.12_dp
            plasma%spec(sp)%A2 = -0.07_dp
            plasma%spec(sp)%I00(:, 0) = cmplx(0.8_dp, -0.2_dp, dp)
            plasma%spec(sp)%I02(:, 0) = cmplx(-0.3_dp, 0.5_dp, dp)
            plasma%spec(sp)%I10(:, 0) = cmplx(0.4_dp, 0.1_dp, dp)
            plasma%spec(sp)%I12(:, 0) = cmplx(-0.2_dp, -0.6_dp, dp)
            ! Deliberately unrelated sentinels: a transposed source/output route
            ! would make the production value disagree with the direct oracle.
            plasma%spec(sp)%I20(:, 0) = cmplx(91.0_dp, -73.0_dp, dp)
            plasma%spec(sp)%I21(:, 0) = cmplx(-57.0_dp, 88.0_dp, dp)

            ! Break radial and species symmetries so phase, normalization, and
            ! observation/source mutations remain visible after assembly.
            plasma%spec(sp)%lambda_D(2:) = plasma%spec(sp)%lambda_D(2:) * &
                [1.1_dp, 0.9_dp]
            plasma%spec(sp)%vT(2:) = plasma%spec(sp)%vT(2:) * [0.85_dp, 1.2_dp]
            plasma%spec(sp)%nu(2:) = plasma%spec(sp)%nu(2:) * [1.25_dp, 0.8_dp]
            plasma%spec(sp)%rho_L(2:) = plasma%spec(sp)%rho_L(2:) * &
                [1.3_dp, 0.75_dp]
            plasma%spec(sp)%A1(2:) = plasma%spec(sp)%A1(2:) + &
                [0.04_dp, -0.03_dp]
            plasma%spec(sp)%A2(2:) = plasma%spec(sp)%A2(2:) + &
                [-0.02_dp, 0.05_dp]
            plasma%spec(sp)%I00(2:, 0) = plasma%spec(sp)%I00(2:, 0) * &
                [cmplx(0.8_dp, 0.2_dp, dp), cmplx(1.1_dp, -0.1_dp, dp)]
            plasma%spec(sp)%I02(2:, 0) = plasma%spec(sp)%I02(2:, 0) * &
                [cmplx(1.2_dp, -0.15_dp, dp), cmplx(0.7_dp, 0.25_dp, dp)]
            plasma%spec(sp)%I10(2:, 0) = plasma%spec(sp)%I10(2:, 0) * &
                [cmplx(0.9_dp, -0.3_dp, dp), cmplx(1.25_dp, 0.1_dp, dp)]
            plasma%spec(sp)%I12(2:, 0) = plasma%spec(sp)%I12(2:, 0) * &
                [cmplx(1.15_dp, 0.05_dp, dp), cmplx(0.65_dp, -0.2_dp, dp)]
        end do

        plasma%spec(0)%Zspec = -1
        plasma%spec(0)%omega_c = -6.0_dp
        plasma%spec(1)%Zspec = 1
        plasma%spec(1)%omega_c = 6.0_dp
    end subroutine populate_fixture

    subroutine run_rejection_case(mode)
        use collisionless_fourier_kernel_m, only: configured_hatG_Bparallel_all

        character(*), intent(in) :: mode
        type(plasma_t) :: empty_plasma
        complex(dp) :: rho_response, current_response

        collision_model = 'FokkerPlanck'
        ion_collision_model = 'FokkerPlanck'
        artificial_debye_case = 0
        collisions_off = .false.
        turn_off_electrons = .false.
        turn_off_ions = .false.
        empty_plasma%n_species = 0
        if (allocated(rg_grid%xb)) deallocate(rg_grid%xb)
        allocate(rg_grid%xb(1))
        rg_grid%npts_b = 1
        rg_grid%xb = 0.0_dp

        select case (mode)
        case ('harmonic')
            mphi_max = 1
            kern_include_ks2 = .true.
        case ('dropped-ks')
            mphi_max = 0
            kern_include_ks2 = .false.
        case default
            error stop 'unknown Bparallel rejection-test mode'
        end select

        ! Each CTest invocation is marked WILL_FAIL. Reaching STOP 0 proves
        ! that the public kernel seam failed to reject an unsupported model.
        call configured_hatG_Bparallel_all(empty_plasma, 0.2_dp, -0.4_dp, 1, &
            rho_response, current_response)
    end subroutine run_rejection_case

    subroutine test_direct_gyroaverage_oracle(plasma)
        type(plasma_t), intent(in) :: plasma
        real(dp), parameter :: kr_observation = -0.4_dp
        real(dp), parameter :: kr_source = 1.1_dp
        real(dp), parameter :: step = 2.0e-5_dp
        real(dp) :: ks, kp_observation
        complex(dp) :: direct_rho, direct_j, expected_rho, expected_j
        complex(dp) :: held_plus, held_minus, all_plus, all_minus
        complex(dp) :: held_derivative, all_derivative

        artificial_debye_case = 0
        kern_include_ks2 = .true.
        ks = plasma%ks(1)
        kp_observation = hypot(ks, kr_observation)

        held_plus = direct_gyro_bracket(ks + step, kr_source, kp_observation, &
            plasma%spec(0)%rho_L(1), plasma%spec(0)%A1(1), &
            plasma%spec(0)%A2(1), plasma%spec(0)%I00(1, 0), &
            plasma%spec(0)%I02(1, 0))
        held_minus = direct_gyro_bracket(ks - step, kr_source, kp_observation, &
            plasma%spec(0)%rho_L(1), plasma%spec(0)%A1(1), &
            plasma%spec(0)%A2(1), plasma%spec(0)%I00(1, 0), &
            plasma%spec(0)%I02(1, 0))
        held_derivative = (held_plus - held_minus) / (2.0_dp * step)
        expected_rho = -com_unit * plasma%spec(0)%vT(1)**2 * held_derivative / &
            (plasma%spec(0)%lambda_D(1)**2 * plasma%spec(0)%nu(1) * sol)

        held_plus = direct_gyro_bracket(ks + step, kr_source, kp_observation, &
            plasma%spec(0)%rho_L(1), plasma%spec(0)%A1(1), &
            plasma%spec(0)%A2(1), plasma%spec(0)%I10(1, 0), &
            plasma%spec(0)%I12(1, 0))
        held_minus = direct_gyro_bracket(ks - step, kr_source, kp_observation, &
            plasma%spec(0)%rho_L(1), plasma%spec(0)%A1(1), &
            plasma%spec(0)%A2(1), plasma%spec(0)%I10(1, 0), &
            plasma%spec(0)%I12(1, 0))
        held_derivative = (held_plus - held_minus) / (2.0_dp * step)
        expected_j = -com_unit * plasma%spec(0)%vT(1)**3 * held_derivative / &
            (plasma%spec(0)%lambda_D(1)**2 * plasma%spec(0)%nu(1) * sol)

        direct_rho = core_rho_Bparallel_sp(plasma, 0, kr_observation, kr_source, 1)
        direct_j = core_j_Bparallel_sp(plasma, 0, kr_observation, kr_source, 1)
        call require_relative_close(direct_rho, expected_rho, 2.0e-8_dp, &
            'rho-Bparallel agrees with direct gyrophase/finite-difference oracle')
        call require_relative_close(direct_j, expected_j, 2.0e-8_dp, &
            'j-Bparallel agrees with direct gyrophase/finite-difference oracle')

        ! Mutation oracle: changing the observation geometry together with the
        ! source is the forbidden full derivative and must be observably different.
        all_plus = direct_gyro_bracket(ks + step, kr_source, &
            hypot(ks + step, kr_observation), plasma%spec(0)%rho_L(1), &
            plasma%spec(0)%A1(1), plasma%spec(0)%A2(1), &
            plasma%spec(0)%I00(1, 0), plasma%spec(0)%I02(1, 0))
        all_minus = direct_gyro_bracket(ks - step, kr_source, &
            hypot(ks - step, kr_observation), plasma%spec(0)%rho_L(1), &
            plasma%spec(0)%A1(1), plasma%spec(0)%A2(1), &
            plasma%spec(0)%I00(1, 0), plasma%spec(0)%I02(1, 0))
        all_derivative = (all_plus - all_minus) / (2.0_dp * step)
        if (abs(all_derivative - &
                expected_rho * (com_unit * plasma%spec(0)%lambda_D(1)**2 * &
                plasma%spec(0)%nu(1) * sol / plasma%spec(0)%vT(1)**2)) &
                < 1.0e-5_dp) then
            error stop 'Bparallel source-only derivative mutation was not detected'
        end if
    end subroutine test_direct_gyroaverage_oracle

    subroutine test_assembled_offdiagonal_oracle(plasma)
        type(plasma_t), intent(in) :: plasma
        integer, parameter :: harmonic_cutoff = 1
        real(dp), parameter :: period = 5.3_dp
        real(dp), parameter :: step = 2.0e-5_dp
        complex(dp), allocatable :: assembled_rho(:,:), assembled_current(:,:)
        complex(dp) :: expected_rho, expected_current
        complex(dp) :: plus_value, minus_value, derivative, phase
        real(dp) :: k_observation, k_source, kp_observation, weight
        integer :: j, sp, observation_index, source_index

        collision_model = 'FokkerPlanck'
        ion_collision_model = 'FokkerPlanck'
        artificial_debye_case = 0
        collisions_off = .false.
        mphi_max = 0
        turn_off_electrons = .false.
        turn_off_ions = .false.
        kern_include_ks2 = .true.

        if (allocated(rg_grid%xb)) deallocate(rg_grid%xb)
        allocate(rg_grid%xb(ngrid))
        rg_grid%npts_b = ngrid
        rg_grid%xb = [-0.43_dp, 0.18_dp, 0.91_dp]

        call assemble_periodic_bparallel_matrices(plasma, period, &
            harmonic_cutoff, assembled_rho, assembled_current)

        ! Check row m=+1, source column m'=0. Unequal radial wavenumbers make
        ! observation/source swaps visible, and the nonsymmetric radial nodes
        ! make the sign of the Fourier phase visible.
        observation_index = harmonic_cutoff + 2
        source_index = harmonic_cutoff + 1
        k_observation = 2.0_dp * pi / period
        k_source = 0.0_dp
        weight = 2.0_dp * pi / real(ngrid, dp)
        expected_rho = (0.0_dp, 0.0_dp)
        expected_current = (0.0_dp, 0.0_dp)

        do j = 1, ngrid
            kp_observation = hypot(plasma%ks(j), k_observation)
            phase = exp(-com_unit * (k_observation - k_source) * rg_grid%xb(j))
            do sp = 0, plasma%n_species - 1
                plus_value = direct_gyro_bracket(plasma%ks(j) + step, &
                    k_source, kp_observation, plasma%spec(sp)%rho_L(j), &
                    plasma%spec(sp)%A1(j), plasma%spec(sp)%A2(j), &
                    plasma%spec(sp)%I00(j, 0), plasma%spec(sp)%I02(j, 0))
                minus_value = direct_gyro_bracket(plasma%ks(j) - step, &
                    k_source, kp_observation, plasma%spec(sp)%rho_L(j), &
                    plasma%spec(sp)%A1(j), plasma%spec(sp)%A2(j), &
                    plasma%spec(sp)%I00(j, 0), plasma%spec(sp)%I02(j, 0))
                derivative = (plus_value - minus_value) / (2.0_dp * step)
                expected_rho = expected_rho + phase * (-com_unit) * &
                    plasma%spec(sp)%vT(j)**2 * derivative / &
                    (plasma%spec(sp)%lambda_D(j)**2 * &
                    plasma%spec(sp)%nu(j) * sol * 8.0_dp * pi**2)

                plus_value = direct_gyro_bracket(plasma%ks(j) + step, &
                    k_source, kp_observation, plasma%spec(sp)%rho_L(j), &
                    plasma%spec(sp)%A1(j), plasma%spec(sp)%A2(j), &
                    plasma%spec(sp)%I10(j, 0), plasma%spec(sp)%I12(j, 0))
                minus_value = direct_gyro_bracket(plasma%ks(j) - step, &
                    k_source, kp_observation, plasma%spec(sp)%rho_L(j), &
                    plasma%spec(sp)%A1(j), plasma%spec(sp)%A2(j), &
                    plasma%spec(sp)%I10(j, 0), plasma%spec(sp)%I12(j, 0))
                derivative = (plus_value - minus_value) / (2.0_dp * step)
                expected_current = expected_current + phase * (-com_unit) * &
                    plasma%spec(sp)%vT(j)**3 * derivative / &
                    (plasma%spec(sp)%lambda_D(j)**2 * &
                    plasma%spec(sp)%nu(j) * sol * 8.0_dp * pi**2)
            end do
        end do
        expected_rho = weight * expected_rho
        expected_current = weight * expected_current

        call require_relative_close(assembled_rho(observation_index, source_index), &
            expected_rho, 3.0e-8_dp, &
            'assembled rho column has direct phase/normalization/source oracle')
        call require_relative_close(&
            assembled_current(observation_index, source_index), &
            expected_current, 3.0e-8_dp, &
            'assembled current column has direct phase/normalization/source oracle')
    end subroutine test_assembled_offdiagonal_oracle

    function direct_gyro_bracket(ks_source, kr_source, kp_observation, rho, &
            a1, a2, moment0, moment2) result(value)
        real(dp), intent(in) :: ks_source, kr_source, kp_observation
        real(dp), intent(in) :: rho, a1, a2
        complex(dp), intent(in) :: moment0, moment2
        complex(dp) :: value
        integer, parameter :: ntheta = 20000
        real(dp) :: bplus, bcross, kp_source, theta
        real(dp) :: average, cosine_average
        integer :: i

        kp_source = hypot(ks_source, kr_source)
        bplus = 0.5_dp * rho**2 * (kp_source**2 + kp_observation**2)
        bcross = rho**2 * kp_source * kp_observation
        average = 0.0_dp
        cosine_average = 0.0_dp
        do i = 0, ntheta - 1
            theta = 2.0_dp * acos(-1.0_dp) * real(i, dp) / real(ntheta, dp)
            average = average + exp(-bplus + bcross * cos(theta))
            cosine_average = cosine_average + &
                cos(theta) * exp(-bplus + bcross * cos(theta))
        end do
        average = average / real(ntheta, dp)
        cosine_average = cosine_average / real(ntheta, dp)

        value = ((a1 + a2 * (1.0_dp - bplus)) * average + &
            a2 * bcross * cosine_average) * moment0 + &
            0.5_dp * a2 * average * moment2
    end function direct_gyro_bracket

    subroutine test_long_wavelength_sign_phase(plasma)
        type(plasma_t), intent(inout) :: plasma
        real(dp), parameter :: rho_small = 1.0e-3_dp
        complex(dp), parameter :: drive = cmplx(0.3_dp, -0.4_dp, dp)
        real(dp) :: rho_saved, source_derivative_scale
        complex(dp) :: rho_prefactor, current_prefactor
        complex(dp) :: rho_coefficient, current_coefficient
        complex(dp) :: got_rho, got_current

        rho_saved = plasma%spec(0)%rho_L(1)
        plasma%spec(0)%rho_L(1) = rho_small
        source_derivative_scale = -rho_small**2 * plasma%ks(1)
        rho_prefactor = -com_unit * plasma%spec(0)%vT(1)**2 / &
            (plasma%spec(0)%lambda_D(1)**2 * plasma%spec(0)%nu(1) * sol)
        current_prefactor = -com_unit * plasma%spec(0)%vT(1)**3 / &
            (plasma%spec(0)%lambda_D(1)**2 * plasma%spec(0)%nu(1) * sol)
        rho_coefficient = (&
            (plasma%spec(0)%A1(1) + 2.0_dp * plasma%spec(0)%A2(1)) * &
                plasma%spec(0)%I00(1, 0) + &
            0.5_dp * plasma%spec(0)%A2(1) * plasma%spec(0)%I02(1, 0))
        current_coefficient = (&
            (plasma%spec(0)%A1(1) + 2.0_dp * plasma%spec(0)%A2(1)) * &
                plasma%spec(0)%I10(1, 0) + &
            0.5_dp * plasma%spec(0)%A2(1) * plasma%spec(0)%I12(1, 0))

        got_rho = drive * core_rho_Bparallel_sp(&
            plasma, 0, -0.4_dp, 1.1_dp, 1) / &
            (rho_prefactor * source_derivative_scale)
        got_current = drive * core_j_Bparallel_sp(&
            plasma, 0, -0.4_dp, 1.1_dp, 1) / &
            (current_prefactor * source_derivative_scale)
        call require_relative_close(got_rho, drive * rho_coefficient, &
            1.0e-5_dp, 'long-wave rho sign and complex drive phase')
        call require_relative_close(got_current, drive * current_coefficient, &
            1.0e-5_dp, 'long-wave current sign and complex drive phase')
        plasma%spec(0)%rho_L(1) = rho_saved
    end subroutine test_long_wavelength_sign_phase

    subroutine test_species_sign_and_moment_routes(plasma)
        type(plasma_t), intent(in) :: plasma
        complex(dp) :: rho_e, rho_i, j_e, j_i

        rho_e = core_rho_Bparallel_sp(plasma, 0, -0.4_dp, 1.1_dp, 1)
        rho_i = core_rho_Bparallel_sp(plasma, 1, -0.4_dp, 1.1_dp, 1)
        j_e = core_j_Bparallel_sp(plasma, 0, -0.4_dp, 1.1_dp, 1)
        j_i = core_j_Bparallel_sp(plasma, 1, -0.4_dp, 1.1_dp, 1)
        call require_close(rho_e, rho_i, 1.0e-14_dp, &
            'rho-Bparallel is even under signed cyclotron-frequency reversal')
        call require_close(j_e, j_i, 1.0e-14_dp, &
            'j-Bparallel is even under signed cyclotron-frequency reversal')
    end subroutine test_species_sign_and_moment_routes

    subroutine test_support_validation()
        complex(dp), parameter :: zero_drive = (0.0_dp, 0.0_dp)
        complex(dp), parameter :: active_drive = (0.3_dp, -0.2_dp)

        call require_status(periodic_bparallel_support_status(active_drive, &
            'FokkerPlanck', 'FokkerPlanck', .false., 0, 0, .false.), &
            PERIODIC_BPAR_OK, &
            'supported FP mphi=0 Bparallel drive')
        call require_status(periodic_bparallel_support_status(zero_drive, &
            'Krook', 'collisionless', .true., 1, 2, .true.), &
            PERIODIC_BPAR_OK, &
            'zero Bparallel preserves every pre-existing model')
        call require_status(periodic_bparallel_support_status(active_drive, &
            'FokkerPlanck', 'collisionless', .false., 0, 0, .false.), &
            PERIODIC_BPAR_UNSUPPORTED_COLLISION, &
            'collisionless nonzero Bparallel rejected')
        call require_status(periodic_bparallel_support_status(active_drive, &
            'Krook', 'FokkerPlanck', .false., 0, 0, .false.), &
            PERIODIC_BPAR_UNSUPPORTED_COLLISION, &
            'Krook nonzero Bparallel rejected')
        call require_status(periodic_bparallel_support_status(active_drive, &
            'FokkerPlanck', 'FokkerPlanck', .true., 0, 0, .false.), &
            PERIODIC_BPAR_UNSUPPORTED_COLLISION, &
            'collisions-off nonzero Bparallel rejected')
        call require_status(periodic_bparallel_support_status(active_drive, &
            'FokkerPlanck', 'FokkerPlanck', .false., 0, 1, .false.), &
            PERIODIC_BPAR_UNSUPPORTED_HARMONIC, &
            'nonzero cyclotron harmonic Bparallel rejected')
        call require_status(periodic_bparallel_support_status(active_drive, &
            'FokkerPlanck', 'FokkerPlanck', .false., 0, 0, .true.), &
            PERIODIC_BPAR_UNSUPPORTED_APPROXIMATION, &
            'ks-dropped Bparallel approximation rejected')
        call require_status(periodic_bparallel_support_status(active_drive, &
            'FokkerPlanck', 'FokkerPlanck', .false., 1, 0, .false.), &
            PERIODIC_BPAR_UNSUPPORTED_DEBYE_MODEL, &
            'Debye-only nonzero Bparallel rejected instead of zeroed')
    end subroutine test_support_validation

    subroutine require_close(got, want, tolerance, label)
        complex(dp), intent(in) :: got, want
        real(dp), intent(in) :: tolerance
        character(*), intent(in) :: label

        if (abs(got - want) > tolerance * max(1.0_dp, abs(want))) then
            print *, 'FAIL: ', label
            print *, '  got:  ', got
            print *, '  want: ', want
            error stop 'Bparallel oracle mismatch'
        end if
        print *, 'PASS: ', label
    end subroutine require_close

    subroutine require_relative_close(got, want, tolerance, label)
        complex(dp), intent(in) :: got, want
        real(dp), intent(in) :: tolerance
        character(*), intent(in) :: label

        if (abs(got - want) > tolerance * abs(want)) then
            print *, 'FAIL: ', label
            print *, '  got:  ', got
            print *, '  want: ', want
            error stop 'Bparallel long-wave oracle mismatch'
        end if
        print *, 'PASS: ', label
    end subroutine require_relative_close

    subroutine require_status(got, want, label)
        integer, intent(in) :: got, want
        character(*), intent(in) :: label

        if (got /= want) then
            print *, 'FAIL: ', label, ' got=', got, ' want=', want
            error stop 'Bparallel support status mismatch'
        end if
        print *, 'PASS: ', label
    end subroutine require_status

end program test_bparallel_linear_response
