module Krook_kernel_plasma_prefacs_m

    contains

    pure function Krook_collisionless_kpar(kpar, epsilon) result(val)

        use KIM_kinds_m, only: dp

        implicit none

        real(dp), intent(in) :: kpar, epsilon
        complex(dp) :: val

        val = cmplx(kpar, epsilon, dp)

    end function Krook_collisionless_kpar

    pure function Krook_collisionless_kpar_magnitude(kpar, epsilon) result(val)

        use KIM_kinds_m, only: dp

        implicit none

        real(dp), intent(in) :: kpar, epsilon
        real(dp) :: val

        ! Broaden analytical |k_parallel| terms without turning them into
        ! signed poles.  In particular, a complex z0 in the lower half-plane
        ! makes the analytic continuation of plasma_Z grow exponentially.
        val = sqrt(kpar**2 + epsilon**2)

    end function Krook_collisionless_kpar_magnitude

    pure function Krook_collisionless_z0(om_E, omega, kpar, vT, epsilon) result(val)

        use KIM_kinds_m, only: dp

        implicit none

        real(dp), intent(in) :: om_E, omega, kpar, vT
        real(dp), intent(in), optional :: epsilon
        complex(dp) :: val

        if (present(epsilon)) then
            val = cmplx(-(om_E - omega) / (&
                Krook_collisionless_kpar_magnitude(kpar, epsilon) * sqrt(2.0_dp) * vT), 0.0_dp, dp)
        else
            val = cmplx(-(om_E - omega) / (abs(kpar) * sqrt(2.0_dp) * vT), 0.0_dp, dp)
        end if

    end function Krook_collisionless_z0

    function Krook_z0_cc(j, spec, collisionless, epsilon) result(val)

        use KIM_kinds_m, only: dp
        use species_m, only: plasma, species_t
        use setup_m, only: omega

        implicit none

        integer, intent(in) :: j
        type(species_t), intent(in) :: spec
        logical, intent(in) :: collisionless
        real(dp), intent(in), optional :: epsilon
        complex(dp) :: val
        real(dp) :: kpar_cc, om_E_cc, vT_cc

        if (collisionless) then
            if (.not. present(epsilon)) then
                error stop 'Collisionless Krook z0 requires a causal-pole epsilon'
            end if
            kpar_cc = 0.5_dp * (plasma%kp(j) + plasma%kp(j+1))
            om_E_cc = 0.5_dp * (plasma%om_E(j) + plasma%om_E(j+1))
            vT_cc = 0.5_dp * (spec%vT(j) + spec%vT(j+1))
            if (vT_cc == 0.0_dp) then
                error stop 'Cannot form collisionless Krook z0 at vT = 0'
            end if
            val = Krook_collisionless_z0(om_E_cc, omega, kpar_cc, vT_cc, epsilon)
        else
            val = 0.5_dp * (spec%z0(j) + spec%z0(j+1))
        end if

    end function Krook_z0_cc

    function Krook_kappa_rho_phi(j, spec) result(val)

        use species_m, only: species_t
        use KIM_kinds_m, only: dp

        implicit none

        integer, intent(in) :: j
        type(species_t), intent(in) :: spec
        real(dp) :: val
        real(dp) :: lambda

        lambda = 0.5d0 * (spec%lambda_D(j) + spec%lambda_D(j+1))
        val = 1.0d0 / (lambda**2.0d0)  !/ sqrt(2.0d0)

    end function

    function Krook_kappa_rho_B(j, spec, collisionless, epsilon) result(val)

        use species_m, only: plasma, species_t
        use KIM_kinds_m, only: dp
        use constants_m, only: sol, com_unit

        implicit none

        integer, intent(in) :: j
        type(species_t), intent(in) :: spec
        logical, intent(in), optional :: collisionless
        real(dp), intent(in), optional :: epsilon
        complex(dp) :: val
        complex(dp) :: kpar_denom
        real(dp) :: kpar_cc
        logical :: use_collisionless

        kpar_cc = 0.5d0 * (plasma%kp(j) + plasma%kp(j+1))
        use_collisionless = .false.
        if (present(collisionless)) use_collisionless = collisionless
        if (use_collisionless) then
            if (.not. present(epsilon)) error stop 'Collisionless kappa_rho_B requires epsilon'
            ! kappa_rho_B contains a signed 1/k_parallel and therefore uses
            ! the causal pole rather than the broadened magnitude.
            kpar_denom = Krook_collisionless_kpar(kpar_cc, epsilon)
        else
            kpar_denom = cmplx(kpar_cc, 0.0_dp, dp)
        end if

        val = (0.5d0 * (spec%vT(j) + spec%vT(j+1)))**2.0d0 * com_unit &
            / (0.5d0 * (spec%lambda_D(j) + spec%lambda_D(j+1)))**2.0d0 &
            / (0.5d0 * (spec%omega_c(j) + spec%omega_c(j+1))) &
            / kpar_denom &
            / sol

    end function


    function Krook_G0_rho_phi(j, spec) result(val)

        use species_m, only: species_t
        use KIM_kinds_m, only: dp

        implicit none

        integer, intent(in) :: j
        type(species_t), intent(in) :: spec
        real(dp) :: val

        ! Unused arguments j and spec - required for interface consistency
        val = -1.0d0 !/ sqrt(2.0d0)

    end function

    function Krook_G1_rho_phi(j, spec, collisionless, epsilon) result(val)

        use KIM_kinds_m, only: dp
        use species_m, only: species_t, plasma

        implicit none

        integer, intent(in) :: j
        type(species_t), intent(in) :: spec
        logical, intent(in), optional :: collisionless
        real(dp), intent(in), optional :: epsilon
        complex(dp) :: val
        real(dp) :: ks_val, kpar, A1, A2, rhoT
        complex(dp) :: plasma_Z, z0, kpar_denom
        logical :: use_collisionless

        ks_val = 0.5d0 * (plasma%ks(j) + plasma%ks(j+1))
        kpar = 0.5d0 * (plasma%kp(j) + plasma%kp(j+1))
        A1 = 0.5d0 * (spec%A1(j) + spec%A1(j+1))
        A2 = 0.5d0 * (spec%A2(j) + spec%A2(j+1))
        use_collisionless = .false.
        if (present(collisionless)) use_collisionless = collisionless
        if (use_collisionless) then
            if (.not. present(epsilon)) error stop 'Collisionless G1_rho_phi requires epsilon'
            kpar_denom = cmplx(Krook_collisionless_kpar_magnitude(kpar, epsilon), 0.0_dp, dp)
            z0 = Krook_z0_cc(j, spec, .true., epsilon)
        else
            kpar_denom = cmplx(abs(kpar), 0.0_dp, dp)
            z0 = Krook_z0_cc(j, spec, .false.)
        end if
        rhoT = 0.5d0 * (spec%rho_L(j) + spec%rho_L(j+1))

        val = ks_val * rhoT /(kpar_denom * sqrt(2.0d0)) &
            * (&
                A1 * plasma_Z(z0) + A2 * plasma_Z(z0) * (1.0d0 + z0**2.0d0) + z0 * A2 &
            )

    end function

    function Krook_G2_rho_phi(j, spec, collisionless, epsilon) result(val)

        use KIM_kinds_m, only: dp
        use species_m, only: species_t, plasma

        implicit none

        integer, intent(in) :: j
        type(species_t), intent(in) :: spec
        logical, intent(in), optional :: collisionless
        real(dp), intent(in), optional :: epsilon
        complex(dp) :: val
        real(dp) :: ks_val, kpar, A2, rhoT
        complex(dp) :: plasma_Z, z0, kpar_denom
        logical :: use_collisionless

        ks_val = 0.5d0 * (plasma%ks(j) + plasma%ks(j+1))
        kpar = 0.5d0 * (plasma%kp(j) + plasma%kp(j+1))
        A2 = 0.5d0 * (spec%A2(j) + spec%A2(j+1))
        rhoT = 0.5d0 * (spec%rho_L(j) + spec%rho_L(j+1))
        use_collisionless = .false.
        if (present(collisionless)) use_collisionless = collisionless
        if (use_collisionless) then
            if (.not. present(epsilon)) error stop 'Collisionless G2_rho_phi requires epsilon'
            kpar_denom = cmplx(Krook_collisionless_kpar_magnitude(kpar, epsilon), 0.0_dp, dp)
            z0 = Krook_z0_cc(j, spec, .true., epsilon)
        else
            kpar_denom = cmplx(abs(kpar), 0.0_dp, dp)
            z0 = Krook_z0_cc(j, spec, .false.)
        end if

        val = - ks_val * rhoT * A2 * plasma_Z(z0) / (kpar_denom * sqrt(2.0d0))

    end function

    function Krook_G3_rho_phi(j, spec, collisionless, epsilon) result(val)

        use KIM_kinds_m, only: dp
        use species_m, only: species_t, plasma

        implicit none

        integer, intent(in) :: j
        type(species_t), intent(in) :: spec
        logical, intent(in), optional :: collisionless
        real(dp), intent(in), optional :: epsilon
        complex(dp) :: val
        real(dp) :: ks_val, kpar, A2, rhoT
        complex(dp) :: kpar_denom, plasma_Z, response, z0
        logical :: use_collisionless

        ks_val = 0.5d0 * (plasma%ks(j) + plasma%ks(j+1))
        kpar = 0.5d0 * (plasma%kp(j) + plasma%kp(j+1))
        A2 = 0.5d0 * (spec%A2(j) + spec%A2(j+1))
        rhoT = 0.5d0 * (spec%rho_L(j) + spec%rho_L(j+1))
        use_collisionless = .false.
        if (present(collisionless)) use_collisionless = collisionless
        response = cmplx(1.0_dp, 0.0_dp, dp)
        if (use_collisionless) then
            if (.not. present(epsilon)) error stop 'Collisionless G3_rho_phi requires epsilon'
            kpar_denom = cmplx(Krook_collisionless_kpar_magnitude(kpar, epsilon), 0.0_dp, dp)
            z0 = Krook_z0_cc(j, spec, .true., epsilon)
            response = plasma_Z(z0)
        else
            kpar_denom = cmplx(abs(kpar), 0.0_dp, dp)
        end if

        val = ks_val * rhoT /(kpar_denom * sqrt(2.0d0)) * A2 * response

    end function

    function Krook_G0_rho_B(j, spec) result(val)

        use species_m, only: species_t
        use KIM_kinds_m, only: dp

        implicit none

        integer, intent(in) :: j
        type(species_t), intent(in) :: spec
        real(dp) :: val

        ! Unused arguments j and spec - required for interface consistency
        val = 0.0d0

    end function

    function Krook_G1_rho_B(j, spec, collisionless, epsilon) result(val)

        use KIM_kinds_m, only: dp
        use species_m, only: species_t

        implicit none

        integer, intent(in) :: j
        type(species_t), intent(in) :: spec
        logical, intent(in), optional :: collisionless
        real(dp), intent(in), optional :: epsilon
        complex(dp) :: val
        real(dp) :: A1, A2
        complex(dp) :: plasma_Z, z0
        logical :: use_collisionless

        A1 = 0.5d0 * (spec%A1(j) + spec%A1(j+1))
        A2 = 0.5d0 * (spec%A2(j) + spec%A2(j+1))
        use_collisionless = .false.
        if (present(collisionless)) use_collisionless = collisionless
        if (use_collisionless) then
            if (.not. present(epsilon)) error stop 'Collisionless G1_rho_B requires epsilon'
            z0 = Krook_z0_cc(j, spec, .true., epsilon)
        else
            z0 = Krook_z0_cc(j, spec, .false.)
        end if

        val = A1 * (z0 * plasma_Z(z0) + 1.0d0) + A2 * &
            (&
                0.5d0 + (z0 * plasma_Z(z0) + 1.0d0) * (1.0d0 + z0**2.0d0) &
            )

    end function

    function Krook_G2_rho_B(j, spec, collisionless, epsilon) result(val)

        use KIM_kinds_m, only: dp
        use species_m, only: species_t

        implicit none

        integer, intent(in) :: j
        type(species_t),intent(in) :: spec
        logical, intent(in), optional :: collisionless
        real(dp), intent(in), optional :: epsilon
        complex(dp) :: val
        real(dp) :: A2
        complex(dp) :: plasma_Z, z0
        logical :: use_collisionless

        A2 = 0.5d0 * (spec%A2(j) + spec%A2(j+1))
        use_collisionless = .false.
        if (present(collisionless)) use_collisionless = collisionless
        if (use_collisionless) then
            if (.not. present(epsilon)) error stop 'Collisionless G2_rho_B requires epsilon'
            z0 = Krook_z0_cc(j, spec, .true., epsilon)
        else
            z0 = Krook_z0_cc(j, spec, .false.)
        end if

        val = - A2 * (z0 * plasma_Z(z0) + 1.0d0)

    end function

    function Krook_G3_rho_B(j, spec, collisionless, epsilon) result(val)

        use KIM_kinds_m, only: dp
        use species_m, only: species_t

        implicit none

        integer, intent(in) :: j
        type(species_t), intent(in) :: spec
        logical, intent(in), optional :: collisionless
        real(dp), intent(in), optional :: epsilon
        complex(dp) :: val
        real(dp) :: A2
        complex(dp) :: plasma_Z, z0
        logical :: use_collisionless

        A2 = 0.5d0 * (spec%A2(j) + spec%A2(j+1))
        use_collisionless = .false.
        if (present(collisionless)) use_collisionless = collisionless
        if (use_collisionless) then
            if (.not. present(epsilon)) error stop 'Collisionless G3_rho_B requires epsilon'
            z0 = Krook_z0_cc(j, spec, .true., epsilon)
            val = A2 * (z0 * plasma_Z(z0) + 1.0_dp)
        else
            val = cmplx(A2, 0.0_dp, dp)
        end if

    end function

end module
