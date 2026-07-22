module collisionless_fourier_kernel_m
    ! m_phi=0 collisionless ion response in the forced-periodicity Fourier
    ! basis.  The four cores implement thesis equations (13.53)--(13.59),
    ! (13.63), (13.67), (13.151), and (13.153), including the common
    ! continuous-Fourier normalization 1/(8*pi^2), but excluding the radial
    ! phase exp(-i*(kr-krp)*rg).  See docs/derivations for the reduction.
    use KIM_kinds_m, only: dp
    use species_m, only: plasma_t
    implicit none
    private

    public :: collisionless_ion_cores
    public :: configured_hatG_all
    public :: configured_hatG_rho_phi, configured_hatG_rho_B
    public :: configured_hatG_j_phi, configured_hatG_j_B

contains

    subroutine configured_hatG_all(plasma_in, kr, krp, j, rho_phi, rho_B, j_phi, j_B)
        ! Select the configured ion model without changing the electron model.
        ! In the default FP mode this delegates directly to the established
        ! functions.  In collisionless mode, electrons remain FP while every
        ! enabled ion species uses collisionless_ion_cores.
        use config_m, only: collisionless_kpar_epsilon, ion_collision_model, &
            turn_off_electrons, turn_off_ions
        use constants_m, only: com_unit, pi
        use flr2_fourier_kernel_m, only: hatG_rho_phi, hatG_rho_B, &
            hatG_j_phi, hatG_j_B, flr_arg_pair_sp, core_rho_phi_sp, &
            core_rho_B_sp, core_j_phi_sp, core_j_B_sp
        use grid_m, only: rg_grid

        type(plasma_t), intent(in) :: plasma_in
        real(dp), intent(in) :: kr, krp
        integer, intent(in) :: j
        complex(dp), intent(out) :: rho_phi, rho_B, j_phi, j_B

        integer :: sp
        real(dp) :: bplus, bcross
        complex(dp) :: phase, species_rho_phi, species_rho_B
        complex(dp) :: species_j_phi, species_j_B

        select case (trim(ion_collision_model))
        case ('FokkerPlanck')
            ! Preserve the original arithmetic path exactly by delegation.
            rho_phi = hatG_rho_phi(plasma_in, kr, krp, j)
            rho_B = hatG_rho_B(plasma_in, kr, krp, j)
            j_phi = hatG_j_phi(plasma_in, kr, krp, j)
            j_B = hatG_j_B(plasma_in, kr, krp, j)
            return
        case ('collisionless')
            continue
        case default
            error stop 'Periodic ion_collision_model must be FokkerPlanck or collisionless'
        end select

        if (collisionless_kpar_epsilon <= 0.0_dp) then
            error stop 'Collisionless periodic ions require collisionless_kpar_epsilon > 0'
        end if

        phase = exp(-com_unit * (kr - krp) * rg_grid%xb(j))
        rho_phi = (0.0_dp, 0.0_dp)
        rho_B = (0.0_dp, 0.0_dp)
        j_phi = (0.0_dp, 0.0_dp)
        j_B = (0.0_dp, 0.0_dp)

        if (.not. turn_off_electrons) then
            call flr_arg_pair_sp(plasma_in, 0, kr, krp, j, bplus, bcross)
            rho_phi = core_rho_phi_sp(plasma_in, 0, bplus, bcross, j) / (8.0_dp * pi**2)
            rho_B = core_rho_B_sp(plasma_in, 0, bplus, bcross, j) / (8.0_dp * pi**2)
            j_phi = core_j_phi_sp(plasma_in, 0, bplus, bcross, j) / (8.0_dp * pi**2)
            j_B = core_j_B_sp(plasma_in, 0, bplus, bcross, j) / (8.0_dp * pi**2)
        end if

        if (.not. turn_off_ions) then
            do sp = 1, plasma_in%n_species - 1
                call collisionless_ion_cores(plasma_in, sp, kr, krp, j, &
                    collisionless_kpar_epsilon, species_rho_phi, species_rho_B, &
                    species_j_phi, species_j_B)
                rho_phi = rho_phi + species_rho_phi
                rho_B = rho_B + species_rho_B
                j_phi = j_phi + species_j_phi
                j_B = j_B + species_j_B
            end do
        end if

        rho_phi = phase * rho_phi
        rho_B = phase * rho_B
        j_phi = phase * j_phi
        j_B = phase * j_B
    end subroutine configured_hatG_all

    complex(dp) function configured_hatG_rho_phi(plasma_in, kr, krp, j) result(value)
        type(plasma_t), intent(in) :: plasma_in
        real(dp), intent(in) :: kr, krp
        integer, intent(in) :: j
        complex(dp) :: unused_rho_B, unused_j_phi, unused_j_B

        call configured_hatG_all(plasma_in, kr, krp, j, value, unused_rho_B, &
            unused_j_phi, unused_j_B)
    end function configured_hatG_rho_phi

    complex(dp) function configured_hatG_rho_B(plasma_in, kr, krp, j) result(value)
        type(plasma_t), intent(in) :: plasma_in
        real(dp), intent(in) :: kr, krp
        integer, intent(in) :: j
        complex(dp) :: unused_rho_phi, unused_j_phi, unused_j_B

        call configured_hatG_all(plasma_in, kr, krp, j, unused_rho_phi, value, &
            unused_j_phi, unused_j_B)
    end function configured_hatG_rho_B

    complex(dp) function configured_hatG_j_phi(plasma_in, kr, krp, j) result(value)
        type(plasma_t), intent(in) :: plasma_in
        real(dp), intent(in) :: kr, krp
        integer, intent(in) :: j
        complex(dp) :: unused_rho_phi, unused_rho_B, unused_j_B

        call configured_hatG_all(plasma_in, kr, krp, j, unused_rho_phi, &
            unused_rho_B, value, unused_j_B)
    end function configured_hatG_j_phi

    complex(dp) function configured_hatG_j_B(plasma_in, kr, krp, j) result(value)
        type(plasma_t), intent(in) :: plasma_in
        real(dp), intent(in) :: kr, krp
        integer, intent(in) :: j
        complex(dp) :: unused_rho_phi, unused_rho_B, unused_j_phi

        call configured_hatG_all(plasma_in, kr, krp, j, unused_rho_phi, &
            unused_rho_B, unused_j_phi, value)
    end function configured_hatG_j_B

    subroutine collisionless_ion_cores(plasma_in, sp, kr, krp, j, epsilon, &
            rho_phi, rho_B, j_phi, j_B)
        use config_m, only: artificial_debye_case
        use constants_m, only: com_unit, pi, sol
        use flr2_fourier_kernel_m, only: flr_arg_pair_sp, scaled_bessel_pair
        use Krook_kernel_plasma_prefacs_m, only: Krook_collisionless_kpar, &
            Krook_collisionless_kpar_magnitude, Krook_collisionless_z0
        use setup_m, only: omega

        type(plasma_t), intent(in) :: plasma_in
        integer, intent(in) :: sp, j
        real(dp), intent(in) :: kr, krp, epsilon
        complex(dp), intent(out) :: rho_phi, rho_B, j_phi, j_B

        real(dp) :: grad_A1, grad_A2, bplus, bcross, ks, k_abs
        real(dp) :: lambda_D, omega_c, sI0, sIm1, vT
        real(dp) :: coeff0, coeff1, coeff2
        complex(dp) :: k_pole, zeta0, response_Z
        complex(dp) :: rho_phi_moment, j_phi_moment, magnetic_bracket
        complex(dp) :: plasma_Z

        call validate_inputs(plasma_in, sp, j, epsilon)
        if (artificial_debye_case /= 0) then
            error stop 'Collisionless periodic ions require artificial_debye_case=0'
        end if

        ks = plasma_in%ks(j)
        grad_A1 = plasma_in%spec(sp)%A1(j)
        grad_A2 = plasma_in%spec(sp)%A2(j)
        vT = plasma_in%spec(sp)%vT(j)
        omega_c = plasma_in%spec(sp)%omega_c(j)
        lambda_D = plasma_in%spec(sp)%lambda_D(j)

        call flr_arg_pair_sp(plasma_in, sp, kr, krp, j, bplus, bcross)
        call scaled_bessel_pair(bplus, bcross, sI0, sIm1)

        k_abs = Krook_collisionless_kpar_magnitude(plasma_in%kp(j), epsilon)
        k_pole = Krook_collisionless_kpar(plasma_in%kp(j), epsilon)
        zeta0 = Krook_collisionless_z0(plasma_in%om_E(j), omega, plasma_in%kp(j), vT, epsilon)
        response_Z = plasma_Z(zeta0)

        ! Thesis (13.53)--(13.55), m_phi=0.  Replacing k_parallel by k_abs
        ! in a1 is part of the even collisionless regularization: together
        ! with the rho-Phi prefactor it preserves the homogeneous static
        ! Debye limit independently of the sign of k_parallel.
        coeff0 = sI0 * (&
            -plasma_in%om_E(j) / omega_c &
            + ks * vT**2 / omega_c**2 * (grad_A1 + (1.0_dp - bplus) * grad_A2)) &
            + ks * vT**2 / omega_c**2 * grad_A2 * bcross * sIm1
        coeff1 = -k_abs / omega_c * sI0
        coeff2 = ks * grad_A2 / (2.0_dp * omega_c**2) * sI0

        ! Zeroth parallel-velocity moment, thesis (13.58), and the first
        ! moment obtained directly from the defining integral (13.56).  The
        ! printed expression (13.59) misses the common (1 + zeta0*Z) factor
        ! on the coeff1 and coeff2 terms and doubles the residual coeff2 term.
        rho_phi_moment = response_Z * (&
            coeff0 + sqrt(2.0_dp) * vT * zeta0 * coeff1 &
            + 2.0_dp * vT**2 * zeta0**2 * coeff2) &
            + sqrt(2.0_dp) * vT * (coeff1 + sqrt(2.0_dp) * vT * zeta0 * coeff2)
        j_phi_moment = (zeta0 * response_Z + 1.0_dp) * (&
            coeff0 + sqrt(2.0_dp) * vT * zeta0 * coeff1 &
            + 2.0_dp * vT**2 * zeta0**2 * coeff2) &
            + vT**2 * coeff2

        ! The original 2^(-7/2)/pi^2 and 1/(8*pi^2) prefactors are
        ! written with the same explicit 1/(8*pi^2) normalization.
        rho_phi = omega_c / (sqrt(2.0_dp) * lambda_D**2 * vT * k_abs) &
            * rho_phi_moment / (8.0_dp * pi**2)
        j_phi = omega_c / (lambda_D**2 * k_pole) &
            * j_phi_moment / (8.0_dp * pi**2)

        ! Mathematica-reduced common bracket in thesis (13.151),(13.153).
        magnetic_bracket = 0.5_dp * grad_A2 * sI0 &
            + (zeta0 * response_Z + 1.0_dp) * (&
                (grad_A1 + grad_A2 * (1.0_dp - bplus + zeta0**2)) * sI0 &
                + grad_A2 * bcross * sIm1)

        rho_B = com_unit * vT**2 / (sol * lambda_D**2 * omega_c * k_pole) &
            * magnetic_bracket / (8.0_dp * pi**2)
        j_B = com_unit * sqrt(2.0_dp) * vT**3 * zeta0 &
            / (sol * lambda_D**2 * omega_c * k_pole) &
            * magnetic_bracket / (8.0_dp * pi**2)
    end subroutine collisionless_ion_cores

    subroutine validate_inputs(plasma_in, sp, j, epsilon)
        type(plasma_t), intent(in) :: plasma_in
        integer, intent(in) :: sp, j
        real(dp), intent(in) :: epsilon

        if (sp < 1 .or. sp >= plasma_in%n_species) then
            error stop 'Collisionless periodic kernel requires an ion species index'
        end if
        if (j < 1 .or. j > plasma_in%grid_size) then
            error stop 'Collisionless periodic kernel radial index is out of range'
        end if
        if (epsilon <= 0.0_dp) then
            error stop 'Collisionless periodic kernel requires epsilon > 0'
        end if
        if (plasma_in%spec(sp)%vT(j) <= 0.0_dp) then
            error stop 'Collisionless periodic kernel requires vT > 0'
        end if
        if (plasma_in%spec(sp)%lambda_D(j) <= 0.0_dp) then
            error stop 'Collisionless periodic kernel requires lambda_D > 0'
        end if
        if (abs(plasma_in%spec(sp)%omega_c(j)) <= tiny(1.0_dp)) then
            error stop 'Collisionless periodic kernel requires nonzero omega_c'
        end if
    end subroutine validate_inputs

end module collisionless_fourier_kernel_m
