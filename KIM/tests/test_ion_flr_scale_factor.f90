program test_ion_flr_scale_factor
    use KIM_kinds_m, only: dp
    use species_m, only: plasma, species_t, init_electron_species, &
                         init_deuterium_species
    use config_m, only: ion_flr_scale_factor
    use FP_kernel_plasma_prefacs_m, only: FP_kappa_rho_B, FP_G1_rho_phi, &
                                          FP_G2_rho_phi, FP_G3_rho_phi, &
                                          FP_kappa_j_phi, FP_kappa_j_B
    use integrands_gauss_m, only: integration_point_t
    use kernel_m, only: set_fp_integration_radius

    implicit none

    real(dp), parameter :: scale = 4.0_dp
    type(species_t) :: electron, ion
    complex(dp) :: electron_base(6), electron_scaled(6)
    complex(dp) :: ion_base(6), ion_scaled(6)
    type(integration_point_t) :: point

    call init_electron_species(electron)
    call init_deuterium_species(ion)
    call populate_species(electron)
    call populate_species(ion)
    allocate(plasma%ks_cc(1))
    plasma%ks_cc = 0.7_dp

    ion_flr_scale_factor = 1.0_dp
    call evaluate_scaled_prefactors(electron, electron_base)
    call evaluate_scaled_prefactors(ion, ion_base)

    ion_flr_scale_factor = scale
    call evaluate_scaled_prefactors(electron, electron_scaled)
    call evaluate_scaled_prefactors(ion, ion_scaled)

    call require_close(electron_scaled, electron_base, &
                       'ion FLR scale must not change electron prefactors')
    call require_close(ion_scaled, scale * ion_base, &
                       'ion FLR scale must reach ion prefactors exactly once')
    call set_fp_integration_radius(point, ion, 1)
    call require_real_close(point%rhoT, ion%rho_L_cc(1), &
                            'global geometry must consume stored ion rho_L once')

    print *, 'Ion FLR scale-factor tests PASSED'

contains

    subroutine populate_species(spec)
        type(species_t), intent(inout) :: spec

        allocate(spec%lambda_D_cc(1), spec%vT_cc(1), spec%nu_cc(1))
        allocate(spec%omega_c_cc(1), spec%rho_L_cc(1))
        allocate(spec%A1_cc(1), spec%A2_cc(1))
        allocate(spec%I00_cc(1, 0:0), spec%I20_cc(1, 0:0))

        spec%lambda_D_cc = 2.0_dp
        spec%vT_cc = 3.0_dp
        spec%nu_cc = 5.0_dp
        spec%omega_c_cc = 7.0_dp
        spec%rho_L_cc = 0.37_dp
        spec%A1_cc = 0.11_dp
        spec%A2_cc = -0.04_dp
        spec%I00_cc(:, 0) = cmplx(0.8_dp, -0.2_dp, dp)
        spec%I20_cc(:, 0) = cmplx(-0.3_dp, 0.1_dp, dp)
    end subroutine populate_species

    subroutine evaluate_scaled_prefactors(spec, values)
        type(species_t), intent(in) :: spec
        complex(dp), intent(out) :: values(6)

        values(1) = FP_kappa_rho_B(1, spec)
        values(2) = FP_G1_rho_phi(1, spec, 0)
        values(3) = FP_G2_rho_phi(1, spec, 0)
        values(4) = FP_G3_rho_phi(1, spec, 0)
        values(5) = FP_kappa_j_phi(1, spec)
        values(6) = FP_kappa_j_B(1, spec)
    end subroutine evaluate_scaled_prefactors

    subroutine require_close(got, want, label)
        complex(dp), intent(in) :: got(:), want(:)
        character(*), intent(in) :: label
        real(dp) :: error

        error = maxval(abs(got - want) / max(abs(want), 1.0_dp))
        if (error > 1.0e-12_dp) then
            print *, 'FAIL: ', label
            print *, '  got:  ', got
            print *, '  want: ', want
            print *, '  max relative error: ', error
            error stop 'ion FLR scale-factor mismatch'
        end if
    end subroutine require_close

    subroutine require_real_close(got, want, label)
        real(dp), intent(in) :: got, want
        character(*), intent(in) :: label

        if (abs(got - want) > 1.0e-12_dp * max(abs(want), 1.0_dp)) then
            print *, 'FAIL: ', label
            print *, '  got:  ', got
            print *, '  want: ', want
            error stop 'ion FLR geometry mismatch'
        end if
    end subroutine require_real_close

end program test_ion_flr_scale_factor
