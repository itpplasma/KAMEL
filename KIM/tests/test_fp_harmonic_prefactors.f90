program test_fp_harmonic_prefactors
    use KIM_kinds_m, only: dp
    use species_m, only: plasma, species_t, init_deuterium_species
    use FP_kernel_plasma_prefacs_m, only: FP_G1_j_phi, FP_G1_j_B, &
        FP_kappa_j_phi, FP_kappa_j_B

    implicit none

    type(species_t) :: spec
    complex(dp) :: got, want
    integer :: mphi

    call init_deuterium_species(spec)
    call populate_species(spec)
    allocate(plasma%ks_cc(1))
    plasma%ks_cc = 0.7_dp

    do mphi = -1, 1
        got = FP_G1_j_phi(1, spec, mphi)
        want = (spec%I10_cc(1, mphi) &
            * (spec%A1_cc(1) + spec%A2_cc(1) * (1.0_dp - mphi)) &
            + 0.5_dp * spec%A2_cc(1) * spec%I12_cc(1, mphi)) &
            * FP_kappa_j_phi(1, spec) * plasma%ks_cc(1)
        call require_close(got, want, 'FP G1 j_phi harmonic coefficient')

        got = FP_G1_j_B(1, spec, mphi)
        want = (spec%I11_cc(1, mphi) &
            * (spec%A1_cc(1) + spec%A2_cc(1) * (1.0_dp - mphi)) &
            + 0.5_dp * spec%A2_cc(1) * spec%I13_cc(1, mphi)) &
            * FP_kappa_j_B(1, spec)
        call require_close(got, want, 'FP G1 j_B harmonic coefficient')
    end do

    print *, 'FP harmonic-prefactor tests PASSED'

contains

    subroutine populate_species(spec)
        type(species_t), intent(inout) :: spec
        integer :: mphi

        allocate(spec%lambda_D_cc(1), spec%vT_cc(1), spec%nu_cc(1))
        allocate(spec%omega_c_cc(1), spec%A1_cc(1), spec%A2_cc(1))
        allocate(spec%I10_cc(1, -1:1), spec%I12_cc(1, -1:1))
        allocate(spec%I11_cc(1, -1:1), spec%I13_cc(1, -1:1))

        spec%lambda_D_cc = 2.0_dp
        spec%vT_cc = 3.0_dp
        spec%nu_cc = 5.0_dp
        spec%omega_c_cc = 7.0_dp
        spec%A1_cc = 0.12_dp
        spec%A2_cc = -0.04_dp

        do mphi = -1, 1
            spec%I10_cc(1, mphi) = cmplx(0.8_dp + 0.1_dp * mphi, &
                                         -0.2_dp + 0.05_dp * mphi, dp)
            spec%I12_cc(1, mphi) = cmplx(-0.3_dp + 0.07_dp * mphi, &
                                         0.6_dp - 0.02_dp * mphi, dp)
            spec%I11_cc(1, mphi) = cmplx(0.4_dp - 0.03_dp * mphi, &
                                         0.5_dp + 0.04_dp * mphi, dp)
            spec%I13_cc(1, mphi) = cmplx(-0.7_dp + 0.02_dp * mphi, &
                                         0.1_dp - 0.06_dp * mphi, dp)
        end do
    end subroutine populate_species

    subroutine require_close(got, want, label)
        complex(dp), intent(in) :: got, want
        character(*), intent(in) :: label
        real(dp) :: error

        error = abs(got - want) / max(abs(want), 1.0_dp)
        if (error > 1.0e-12_dp) then
            print *, 'FAIL: ', label
            print *, '  got:  ', got
            print *, '  want: ', want
            error stop 'FP harmonic-prefactor mismatch'
        end if
    end subroutine require_close

end program test_fp_harmonic_prefactors
