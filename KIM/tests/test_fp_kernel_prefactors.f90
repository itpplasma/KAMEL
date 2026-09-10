program test_fp_kernel_prefactors
    use KIM_kinds_m, only: dp
    use species_m, only: plasma, species_t, init_electron_species
    use constants_m, only: com_unit, sol
    use FP_kernel_plasma_prefacs_m, only: FP_G1_rho_B, FP_G1_rho_phi, &
        FP_G3_j_phi, FP_kappa_j_phi

    implicit none

    type(species_t) :: spec
    complex(dp) :: got, want

    call init_electron_species(spec)
    call populate_species(spec)
    allocate(plasma%ks_cc(1))
    plasma%ks_cc = 0.7_dp

    got = FP_G3_j_phi(1, spec, 0)
    want = spec%I10_cc(1, 0) * spec%A2_cc(1) &
        * FP_kappa_j_phi(1, spec) * plasma%ks_cc(1)
    call require_close(got, want, 'FP G3 j_phi uses I10')

    got = FP_G1_rho_phi(1, spec, 0)
    want = (spec%I00_cc(1, 0) * (spec%A1_cc(1) + spec%A2_cc(1)) &
        + 0.5_dp * spec%A2_cc(1) * spec%I02_cc(1, 0)) &
        * com_unit * spec%vT_cc(1)**2 * plasma%ks_cc(1) &
        / (spec%lambda_D_cc(1)**2 * spec%omega_c_cc(1) * spec%nu_cc(1))
    call require_close(got, want, 'FP G1 rho_phi uses source-power I02')

    got = FP_G1_rho_B(1, spec, 0)
    want = (spec%I01_cc(1, 0) * (spec%A1_cc(1) + spec%A2_cc(1)) &
        + 0.5_dp * spec%A2_cc(1) * spec%I03_cc(1, 0)) &
        * (-spec%vT_cc(1)**3) &
        / (spec%lambda_D_cc(1)**2 * spec%omega_c_cc(1) &
        * spec%nu_cc(1) * sol)
    call require_close(got, want, 'FP G1 rho_B uses source-power I03')

    print *, 'FP kernel-prefactor tests PASSED'

contains

    subroutine populate_species(spec)
        type(species_t), intent(inout) :: spec

        allocate(spec%lambda_D_cc(1), spec%vT_cc(1), spec%nu_cc(1))
        allocate(spec%omega_c_cc(1), spec%A1_cc(1), spec%A2_cc(1))
        allocate(spec%I00_cc(1, 0:0), spec%I20_cc(1, 0:0))
        allocate(spec%I02_cc(1, 0:0), spec%I10_cc(1, 0:0))
        allocate(spec%I01_cc(1, 0:0), spec%I21_cc(1, 0:0))
        allocate(spec%I03_cc(1, 0:0))

        spec%lambda_D_cc = 2.0_dp
        spec%vT_cc = 3.0_dp
        spec%nu_cc = 5.0_dp
        spec%omega_c_cc = 7.0_dp
        spec%A1_cc = 0.11_dp
        spec%A2_cc = -0.04_dp
        spec%I00_cc(:, 0) = cmplx(0.2_dp, -0.7_dp, dp)
        spec%I20_cc(:, 0) = cmplx(1.4_dp, 0.9_dp, dp)
        spec%I02_cc(:, 0) = cmplx(-0.6_dp, 1.1_dp, dp)
        spec%I10_cc(:, 0) = cmplx(0.8_dp, -0.2_dp, dp)
        spec%I01_cc(:, 0) = cmplx(-0.3_dp, 0.6_dp, dp)
        spec%I21_cc(:, 0) = cmplx(1.7_dp, -1.3_dp, dp)
        spec%I03_cc(:, 0) = cmplx(-0.9_dp, 0.4_dp, dp)
    end subroutine populate_species

    subroutine require_close(got, want, label)
        complex(dp), intent(in) :: got, want
        character(*), intent(in) :: label
        real(dp) :: error

        error = abs(got - want) / max(abs(want), tiny(1.0_dp))
        if (error > 1.0e-12_dp) then
            print *, 'FAIL: ', label
            print *, '  got:  ', got
            print *, '  want: ', want
            error stop 'FP kernel-prefactor mismatch'
        end if
    end subroutine require_close

end program test_fp_kernel_prefactors
