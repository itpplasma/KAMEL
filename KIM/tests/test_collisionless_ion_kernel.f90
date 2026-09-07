program test_collisionless_ion_kernel

    use KIM_kinds_m, only: dp
    use Krook_kernel_plasma_prefacs_m, only: Krook_collisionless_kpar, Krook_collisionless_kpar_magnitude, &
        Krook_collisionless_z0, Krook_z0_cc, &
        Krook_G1_rho_phi, Krook_G2_rho_phi, Krook_G3_rho_phi, &
        Krook_G1_rho_B, Krook_G2_rho_B, Krook_G3_rho_B, Krook_kappa_rho_B
    use species_m, only: plasma, species_t
    use setup_m, only: omega
    use, intrinsic :: ieee_arithmetic, only: ieee_is_finite

    implicit none

    complex(dp) :: actual, expected, explicit_c0, collisional
    complex(dp) :: pole_wide, pole_narrow, pole_cc, z0_cc, z0_finite, plasma_Z
    real(dp) :: magnitude_wide, magnitude_narrow, pole_magnitude
    type(species_t) :: spec

    pole_wide = Krook_collisionless_kpar(0.0_dp, 0.25_dp)
    expected = cmplx(0.0_dp, 0.25_dp, dp)
    if (abs(pole_wide - expected) > 1.0e-14_dp) then
        print *, 'FAIL: collisionless pole is not k_parallel + i epsilon'
        print *, '  actual   = ', pole_wide
        print *, '  expected = ', expected
        error stop
    end if

    pole_narrow = Krook_collisionless_kpar(0.0_dp, 0.125_dp)
    if (abs(1.0_dp / pole_narrow) <= abs(1.0_dp / pole_wide)) then
        print *, 'FAIL: causal pole does not grow as epsilon decreases'
        error stop
    end if

    pole_magnitude = Krook_collisionless_kpar_magnitude(3.0_dp, 4.0_dp)
    if (abs(pole_magnitude - 5.0_dp) > 1.0e-14_dp) then
        print *, 'FAIL: broadened collisionless |k_parallel| mismatch'
        error stop
    end if

    magnitude_wide = Krook_collisionless_kpar_magnitude(0.0_dp, 0.25_dp)
    magnitude_narrow = Krook_collisionless_kpar_magnitude(0.0_dp, 0.125_dp)
    if (1.0_dp / magnitude_narrow <= 1.0_dp / magnitude_wide) then
        print *, 'FAIL: broadened inverse magnitude does not recover the singular limit'
        error stop
    end if

    actual = Krook_collisionless_z0(2.0_dp, 0.5_dp, -3.0_dp, 4.0_dp, 0.25_dp)
    pole_magnitude = Krook_collisionless_kpar_magnitude(-3.0_dp, 0.25_dp)
    expected = cmplx(-1.5_dp / (pole_magnitude * sqrt(2.0_dp) * 4.0_dp), 0.0_dp, dp)

    if (abs(actual - expected) > 1.0e-14_dp) then
        print *, 'FAIL: collisionless z0 mismatch'
        print *, '  actual   = ', actual
        print *, '  expected = ', expected
        error stop
    end if
    if (aimag(actual) /= 0.0_dp) then
        print *, 'FAIL: magnitude-based collisionless z0 is not real: ', actual
        error stop
    end if

    actual = plasma_Z(Krook_collisionless_z0(-31334.6381_dp, 0.0_dp, &
        1.310539e-9_dp, 18950820.5_dp, 1.0e-5_dp))
    if (.not. ieee_is_finite(real(actual, dp)) .or. .not. ieee_is_finite(aimag(actual))) then
        print *, 'FAIL: collisionless z0 sends plasma Z into its overflowing lower half-plane'
        error stop
    end if

    allocate(plasma%om_E(2), plasma%kp(2), spec%vT(2), spec%z0(2))
    plasma%om_E = [2.0_dp, 4.0_dp]
    plasma%kp = [-3.0_dp, -5.0_dp]
    spec%vT = [4.0_dp, 8.0_dp]
    spec%z0 = [cmplx(1.0_dp, 2.0_dp, dp), cmplx(3.0_dp, 4.0_dp, dp)]
    omega = 0.5_dp

    actual = Krook_z0_cc(1, spec, collisionless=.false.)
    expected = cmplx(2.0_dp, 3.0_dp, dp)
    if (abs(actual - expected) > 1.0e-14_dp) then
        print *, 'FAIL: collisional cell-centred z0 did not preserve existing behavior'
        error stop
    end if

    pole_cc = Krook_collisionless_kpar(-4.0_dp, 0.25_dp)
    pole_magnitude = Krook_collisionless_kpar_magnitude(-4.0_dp, 0.25_dp)
    expected = Krook_collisionless_z0(3.0_dp, omega, -4.0_dp, 6.0_dp, 0.25_dp)
    actual = Krook_z0_cc(1, spec, collisionless=.true., epsilon=0.25_dp)
    if (abs(actual - expected) > 1.0e-14_dp) then
        print *, 'FAIL: explicit collisionless cell-centred z0 mismatch'
        print *, '  actual   = ', actual
        print *, '  expected = ', expected
        error stop
    end if
    if (any(spec%z0 /= [cmplx(1.0_dp, 2.0_dp, dp), cmplx(3.0_dp, 4.0_dp, dp)])) then
        print *, 'FAIL: collisionless z0 calculation mutated the species background'
        error stop
    end if

    allocate(plasma%ks(2), spec%A1(2), spec%A2(2), spec%rho_L(2), spec%lambda_D(2), spec%omega_c(2))
    plasma%ks = [0.2_dp, 0.3_dp]
    spec%A1 = [0.4_dp, 0.5_dp]
    spec%A2 = [-0.2_dp, -0.1_dp]
    spec%rho_L = [0.3_dp, 0.4_dp]
    spec%lambda_D = [0.7_dp, 0.9_dp]
    spec%omega_c = [2.0_dp, 2.4_dp]
    z0_cc = Krook_collisionless_z0(3.0_dp, omega, -4.0_dp, 6.0_dp, 0.25_dp)

    collisional = Krook_G1_rho_phi(1, spec)
    explicit_c0 = Krook_G1_rho_phi(1, spec, collisionless=.true., epsilon=0.25_dp)
    expected = 0.25_dp * 0.35_dp / (pole_magnitude * sqrt(2.0_dp)) * (&
        0.45_dp * plasma_Z(z0_cc) - 0.15_dp * plasma_Z(z0_cc) * (1.0_dp + z0_cc**2) &
        - 0.15_dp * z0_cc)
    call assert_prefactor('G1_rho_phi', collisional, explicit_c0, expected)

    collisional = Krook_G2_rho_phi(1, spec)
    explicit_c0 = Krook_G2_rho_phi(1, spec, collisionless=.true., epsilon=0.25_dp)
    expected = -0.25_dp * 0.35_dp * (-0.15_dp) * plasma_Z(z0_cc) / (pole_magnitude * sqrt(2.0_dp))
    call assert_prefactor('G2_rho_phi', collisional, explicit_c0, expected)

    collisional = Krook_G3_rho_phi(1, spec)
    explicit_c0 = Krook_G3_rho_phi(1, spec, collisionless=.true., epsilon=0.25_dp)
    expected = 0.25_dp * 0.35_dp * (-0.15_dp) * plasma_Z(z0_cc) / (&
        pole_magnitude * sqrt(2.0_dp))
    call assert_prefactor('G3_rho_phi', collisional, explicit_c0, expected)

    collisional = Krook_G1_rho_B(1, spec)
    z0_finite = Krook_z0_cc(1, spec, collisionless=.false.)
    expected = 0.45_dp * (z0_finite * plasma_Z(z0_finite) + 1.0_dp) - 0.15_dp * (&
        0.5_dp + (z0_finite * plasma_Z(z0_finite) + 1.0_dp) * (1.0_dp + z0_finite**2))
    call assert_close('finite Krook G1_rho_B', collisional, expected)

    explicit_c0 = Krook_G1_rho_B(1, spec, collisionless=.true., epsilon=0.25_dp)
    expected = 0.45_dp * (z0_cc * plasma_Z(z0_cc) + 1.0_dp) - 0.15_dp * (&
        0.5_dp + (z0_cc * plasma_Z(z0_cc) + 1.0_dp) * (1.0_dp + z0_cc**2))
    call assert_prefactor('G1_rho_B', collisional, explicit_c0, expected)

    collisional = Krook_G2_rho_B(1, spec)
    explicit_c0 = Krook_G2_rho_B(1, spec, collisionless=.true., epsilon=0.25_dp)
    expected = 0.15_dp * (z0_cc * plasma_Z(z0_cc) + 1.0_dp)
    call assert_prefactor('G2_rho_B', collisional, explicit_c0, expected)

    collisional = Krook_G3_rho_B(1, spec)
    explicit_c0 = Krook_G3_rho_B(1, spec, collisionless=.true., epsilon=0.25_dp)
    expected = -0.15_dp * (z0_cc * plasma_Z(z0_cc) + 1.0_dp)
    call assert_prefactor('G3_rho_B', collisional, explicit_c0, expected)

    collisional = Krook_kappa_rho_B(1, spec)
    explicit_c0 = Krook_kappa_rho_B(1, spec, collisionless=.true., epsilon=0.25_dp)
    expected = collisional * cmplx(-4.0_dp, 0.0_dp, dp) / pole_cc
    call assert_prefactor('kappa_rho_B', collisional, explicit_c0, expected)

    print *, 'PASS: collisionless magnitude terms and causal signed poles remain distinct'

contains

    subroutine assert_close(name, actual_value, reference_value)
        character(len=*), intent(in) :: name
        complex(dp), intent(in) :: actual_value, reference_value
        real(dp) :: tolerance

        tolerance = 1.0e-13_dp * (1.0_dp + abs(reference_value))
        if (abs(actual_value - reference_value) > tolerance) then
            print *, 'FAIL: prefactor mismatch for ', name
            print *, '  actual   = ', actual_value
            print *, '  expected = ', reference_value
            error stop
        end if
    end subroutine assert_close

    subroutine assert_prefactor(name, fp_value, c0_value, reference_value)
        character(len=*), intent(in) :: name
        complex(dp), intent(in) :: fp_value, c0_value, reference_value
        real(dp) :: tolerance

        tolerance = 1.0e-13_dp * (1.0_dp + abs(reference_value))
        if (abs(c0_value - reference_value) > tolerance) then
            print *, 'FAIL: explicit collisionless prefactor mismatch for ', name
            error stop
        end if
        if (abs(c0_value - fp_value) <= 1.0e-14_dp * (1.0_dp + abs(fp_value))) then
            print *, 'FAIL: collisionless selection had no effect for ', name
            error stop
        end if
    end subroutine assert_prefactor
end program test_collisionless_ion_kernel
