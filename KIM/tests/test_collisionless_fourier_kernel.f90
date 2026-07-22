program test_collisionless_fourier_kernel
    use KIM_kinds_m, only: dp
    use collisionless_fourier_kernel_m, only: collisionless_ion_cores
    use config_m, only: artificial_debye_case
    use constants_m, only: pi
    use flr2_fourier_kernel_m, only: set_global_kernel_approximations
    use setup_m, only: omega
    use species_m, only: plasma_t
    use, intrinsic :: ieee_arithmetic, only: ieee_is_finite

    implicit none

    type(plasma_t) :: test_plasma
    complex(dp) :: rho_phi, rho_B, j_phi, j_B
    complex(dp) :: rho_phi_ref, rho_B_ref, j_phi_ref, j_B_ref
    real(dp), parameter :: epsilon = 0.05_dp
    real(dp), parameter :: kr = 0.73_dp, krp = -0.41_dp

    call make_manufactured_plasma(test_plasma)
    call set_global_kernel_approximations(.false.)
    artificial_debye_case = 0
    omega = 0.17_dp

    call collisionless_ion_cores(test_plasma, 1, kr, krp, 1, epsilon, &
        rho_phi, rho_B, j_phi, j_B)

    rho_phi_ref = cmplx(-2.1836855175510750e-2_dp, -8.2515585946286270e-3_dp, dp)
    rho_B_ref = cmplx(-3.0215925066508310e-14_dp, 4.7259864203647080e-14_dp, dp)
    j_phi_ref = cmplx(3.3181741084489427e-3_dp, -2.0023934758877357e-2_dp, dp)
    j_B_ref = cmplx(3.0810039432056770e-14_dp, -4.8189101490788025e-14_dp, dp)

    call assert_close('rho-Phi Mathematica oracle', rho_phi, rho_phi_ref, 2.0e-13_dp)
    ! fortnum's integer-order Bessel approximation is accurate to about 1e-8
    ! for the I_{-1} value used by this manufactured case.
    call assert_close('rho-B Mathematica oracle', rho_B, rho_B_ref, 5.0e-8_dp)
    call assert_close('j-Phi Mathematica oracle', j_phi, j_phi_ref, 2.0e-13_dp)
    call assert_close('j-B Mathematica oracle', j_B, j_B_ref, 5.0e-8_dp)

    call test_magnetic_recurrence(test_plasma)
    call test_homogeneous_static_limit(test_plasma)
    call test_resonance_regularization(test_plasma)
    call test_collision_frequency_independence(test_plasma)

    print *, 'All collisionless Fourier-kernel tests PASSED'

contains

    subroutine make_manufactured_plasma(plasma)
        type(plasma_t), intent(out) :: plasma

        plasma%n_species = 2
        plasma%grid_size = 1
        allocate(plasma%spec(0:1))
        allocate(plasma%ks(1), plasma%kp(1), plasma%om_E(1))
        allocate(plasma%spec(1)%rho_L(1), plasma%spec(1)%vT(1))
        allocate(plasma%spec(1)%omega_c(1), plasma%spec(1)%lambda_D(1))
        allocate(plasma%spec(1)%A1(1), plasma%spec(1)%A2(1), plasma%spec(1)%nu(1))

        plasma%ks(1) = 0.37_dp
        plasma%kp(1) = -0.23_dp
        plasma%om_E(1) = 0.41_dp
        plasma%spec(1)%rho_L(1) = 0.62_dp
        plasma%spec(1)%vT(1) = 1.3_dp
        plasma%spec(1)%omega_c(1) = 2.1_dp
        plasma%spec(1)%lambda_D(1) = 0.8_dp
        plasma%spec(1)%A1(1) = 0.12_dp
        plasma%spec(1)%A2(1) = -0.07_dp
        plasma%spec(1)%nu(1) = 9.0e7_dp
    end subroutine make_manufactured_plasma

    subroutine test_magnetic_recurrence(plasma)
        use Krook_kernel_plasma_prefacs_m, only: Krook_collisionless_z0
        type(plasma_t), intent(in) :: plasma
        complex(dp) :: rp, rb, jp, jb, z0

        call collisionless_ion_cores(plasma, 1, kr, krp, 1, epsilon, rp, rb, jp, jb)
        z0 = Krook_collisionless_z0(plasma%om_E(1), omega, plasma%kp(1), &
            plasma%spec(1)%vT(1), epsilon)
        call assert_close('j-B/rho-B recurrence', jb, &
            sqrt(2.0_dp) * plasma%spec(1)%vT(1) * z0 * rb, 2.0e-13_dp)
    end subroutine test_magnetic_recurrence

    subroutine test_homogeneous_static_limit(plasma)
        type(plasma_t), intent(inout) :: plasma
        complex(dp) :: rp, rb, jp, jb, expected
        real(dp) :: saved_omega, saved_om_E, saved_ks, saved_rho_L
        real(dp) :: saved_A1, saved_A2

        saved_omega = omega
        saved_om_E = plasma%om_E(1)
        saved_ks = plasma%ks(1)
        saved_rho_L = plasma%spec(1)%rho_L(1)
        saved_A1 = plasma%spec(1)%A1(1)
        saved_A2 = plasma%spec(1)%A2(1)

        omega = 0.0_dp
        plasma%om_E(1) = 0.0_dp
        plasma%ks(1) = 0.0_dp
        plasma%spec(1)%rho_L(1) = 0.0_dp
        plasma%spec(1)%A1(1) = 0.0_dp
        plasma%spec(1)%A2(1) = 0.0_dp

        call collisionless_ion_cores(plasma, 1, kr, krp, 1, epsilon, rp, rb, jp, jb)
        expected = cmplx(-1.0_dp / (8.0_dp * pi**2 * plasma%spec(1)%lambda_D(1)**2), 0.0_dp, dp)
        call assert_close('homogeneous static Debye limit', rp, expected, 2.0e-13_dp)
        call assert_close('homogeneous rho-B vanishes', rb, cmplx(0.0_dp, 0.0_dp, dp), 1.0e-30_dp)
        call assert_close('homogeneous j-Phi vanishes', jp, cmplx(0.0_dp, 0.0_dp, dp), 1.0e-30_dp)
        call assert_close('homogeneous j-B vanishes', jb, cmplx(0.0_dp, 0.0_dp, dp), 1.0e-30_dp)

        omega = saved_omega
        plasma%om_E(1) = saved_om_E
        plasma%ks(1) = saved_ks
        plasma%spec(1)%rho_L(1) = saved_rho_L
        plasma%spec(1)%A1(1) = saved_A1
        plasma%spec(1)%A2(1) = saved_A2
    end subroutine test_homogeneous_static_limit

    subroutine test_resonance_regularization(plasma)
        type(plasma_t), intent(inout) :: plasma
        complex(dp) :: rp_wide, rb_wide, jp_wide, jb_wide
        complex(dp) :: rp_narrow, rb_narrow, jp_narrow, jb_narrow
        real(dp) :: saved_kp

        saved_kp = plasma%kp(1)
        plasma%kp(1) = 0.0_dp
        call collisionless_ion_cores(plasma, 1, kr, krp, 1, 0.10_dp, &
            rp_wide, rb_wide, jp_wide, jb_wide)
        call collisionless_ion_cores(plasma, 1, kr, krp, 1, 0.025_dp, &
            rp_narrow, rb_narrow, jp_narrow, jb_narrow)

        call assert_finite('rho-Phi at resonance', rp_narrow)
        call assert_finite('rho-B at resonance', rb_narrow)
        call assert_finite('j-Phi at resonance', jp_narrow)
        call assert_finite('j-B at resonance', jb_narrow)
        if (abs(rb_wide - rb_narrow) <= tiny(1.0_dp) .or. &
                abs(jp_wide - jp_narrow) <= tiny(1.0_dp)) then
            print *, 'FAIL: collisionless kernels are insensitive to epsilon at resonance'
            error stop
        end if
        plasma%kp(1) = saved_kp
    end subroutine test_resonance_regularization

    subroutine test_collision_frequency_independence(plasma)
        type(plasma_t), intent(inout) :: plasma
        complex(dp) :: rp_before, rb_before, jp_before, jb_before
        complex(dp) :: rp_after, rb_after, jp_after, jb_after

        call collisionless_ion_cores(plasma, 1, kr, krp, 1, epsilon, &
            rp_before, rb_before, jp_before, jb_before)
        plasma%spec(1)%nu(1) = 1.0e-99_dp
        call collisionless_ion_cores(plasma, 1, kr, krp, 1, epsilon, &
            rp_after, rb_after, jp_after, jb_after)
        call assert_close('rho-Phi collision-frequency independence', rp_after, rp_before, 0.0_dp)
        call assert_close('rho-B collision-frequency independence', rb_after, rb_before, 0.0_dp)
        call assert_close('j-Phi collision-frequency independence', jp_after, jp_before, 0.0_dp)
        call assert_close('j-B collision-frequency independence', jb_after, jb_before, 0.0_dp)
    end subroutine test_collision_frequency_independence

    subroutine assert_close(label, actual, expected, relative_tolerance)
        character(*), intent(in) :: label
        complex(dp), intent(in) :: actual, expected
        real(dp), intent(in) :: relative_tolerance
        real(dp) :: tolerance

        if (abs(expected) > 0.0_dp) then
            tolerance = relative_tolerance * abs(expected)
        else
            tolerance = relative_tolerance
        end if
        if (abs(actual - expected) > tolerance) then
            print *, 'FAIL: ', label
            print *, '  actual   = ', actual
            print *, '  expected = ', expected
            print *, '  error    = ', abs(actual - expected), ' tolerance = ', tolerance
            error stop
        end if
    end subroutine assert_close

    subroutine assert_finite(label, value)
        character(*), intent(in) :: label
        complex(dp), intent(in) :: value

        if (.not. ieee_is_finite(real(value, dp)) .or. .not. ieee_is_finite(aimag(value))) then
            print *, 'FAIL: ', label, ' is not finite: ', value
            error stop
        end if
    end subroutine assert_finite

end program test_collisionless_fourier_kernel
