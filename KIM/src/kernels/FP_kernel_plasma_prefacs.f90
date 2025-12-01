! module for plasma factors in the Fokker-Planck case
module FP_kernel_plasma_prefacs_m

    implicit none

    contains

    function FP_kappa_rho_phi(j, spec) result(val)

        use species_m, only: species_t
        use KIM_kinds_m, only: dp

        implicit none

        integer, intent(in) :: j
        type(species_t), intent(in) :: spec
        real(dp) :: val

        val = 1.0d0 / (spec%lambda_D_cc(j)**2.0d0)

    end function FP_kappa_rho_phi

    function FP_kappa_rho_B(j, spec) result(val)

        use species_m, only: species_t
        use KIM_kinds_m, only: dp
        use constants_m, only: sol

        implicit none

        integer, intent(in) :: j
        type(species_t), intent(in) :: spec
        real(dp) :: val
        real(dp) :: lambda, vT, nu, omega_c

        lambda = spec%lambda_D_cc(j)
        vT = spec%vT_cc(j)
        nu = spec%nu_cc(j)
        omega_c = spec%omega_c_cc(j)

        val = - vT**3.0d0 / (lambda**2.0d0 * omega_c * nu * sol)

    end function FP_kappa_rho_B

    function FP_kappa_j_phi(j, spec) result(val)

        use species_m, only: species_t
        use KIM_kinds_m, only: dp
        use constants_m, only: com_unit

        implicit none

        integer, intent(in) :: j
        type(species_t), intent(in) :: spec
        complex(dp) :: val
        real(dp) :: lambda, vT, nu, omega_c

        lambda = spec%lambda_D_cc(j)
        vT = spec%vT_cc(j)
        nu = spec%nu_cc(j)
        omega_c = spec%omega_c_cc(j)

        val = com_unit * vT**3.0d0 / (lambda**2.0d0 * omega_c * nu)

    end function FP_kappa_j_phi

    function FP_kappa_j_B(j, spec) result(val)

        use species_m, only: species_t
        use KIM_kinds_m, only: dp
        use constants_m, only: sol

        implicit none

        integer, intent(in) :: j
        type(species_t), intent(in) :: spec
        real(dp) :: val
        real(dp) :: lambda, vT, nu, omega_c

        lambda = spec%lambda_D_cc(j)
        vT = spec%vT_cc(j)
        nu = spec%nu_cc(j)
        omega_c = spec%omega_c_cc(j)

        val = - vT**4.0d0 / (lambda**2.0d0 * omega_c * nu * sol)

    end function FP_kappa_j_B

    function FP_G0_rho_phi(j, spec) result(val)

        use species_m, only: species_t
        use KIM_kinds_m, only: dp

        implicit none

        integer, intent(in) :: j
        type(species_t), intent(in) :: spec
        real(dp) :: val

        val = - FP_kappa_rho_phi(j, spec)

    end function

    function FP_G1_rho_phi(j, spec, mphi) result(val)

        use species_m, only: species_t, plasma
        use KIM_kinds_m, only: dp
        use constants_m, only: com_unit

        implicit none

        integer, intent(in) :: j, mphi
        type(species_t), intent(in) :: spec
        complex(dp) :: val
        complex(dp) :: I00, I20
        real(dp) :: A1, A2
        complex(dp) :: prefactor

        A1 = spec%A1_cc(j)
        A2 = spec%A2_cc(j)
        I00 = spec%I00_cc(j, mphi)
        I20 = spec%I20_cc(j, mphi)

        prefactor = FP_kappa_rho_phi(j, spec) * com_unit * spec%vT_cc(j)**2.0d0 / &
            (spec%omega_c_cc(j) * spec%nu_cc(j)) * plasma%ks_cc(j)

        val = (I00 * (A1 + A2 * (1.0d0 - mphi)) + 0.5d0 * A2 * I20) * prefactor

    end function FP_G1_rho_phi


    function FP_G2_rho_phi(j, spec, mphi) result(val)

        use species_m, only: species_t, plasma
        use KIM_kinds_m, only: dp
        use constants_m, only: com_unit

        implicit none

        integer, intent(in) :: j, mphi
        type(species_t), intent(in) :: spec
        complex(dp) :: val
        complex(dp) :: I00
        real(dp) :: A2
        complex(dp) :: prefactor

        A2 = spec%A2_cc(j)
        I00 = spec%I00_cc(j, mphi)

        prefactor = com_unit * spec%vT_cc(j)**2.0d0 / &
            (spec%omega_c_cc(j) * spec%nu_cc(j)) * plasma%ks_cc(j)

        val = - I00 * A2 * prefactor * FP_kappa_rho_phi(j, spec)

    end function FP_G2_rho_phi

    function FP_G3_rho_phi(j, spec, mphi) result(val)

        use species_m, only: species_t, plasma
        use KIM_kinds_m, only: dp
        use constants_m, only: com_unit

        implicit none

        integer, intent(in) :: j, mphi
        type(species_t), intent(in) :: spec
        complex(dp) :: val
        complex(dp) :: I00
        real(dp) :: A2
        complex(dp) :: prefactor

        A2 = spec%A2_cc(j)
        I00 = spec%I00_cc(j, mphi)

        prefactor = com_unit * spec%vT_cc(j)**2.0d0 / &
            (spec%omega_c_cc(j) * spec%nu_cc(j)) * plasma%ks_cc(j)

        val = I00 * A2 * prefactor * FP_kappa_rho_phi(j, spec)

    end function FP_G3_rho_phi


    function FP_G1_rho_B(j, spec, mphi) result(val)

        use species_m, only: species_t
        use KIM_kinds_m, only: dp

        implicit none

        integer, intent(in) :: j, mphi
        type(species_t), intent(in) :: spec
        complex(dp) :: val
        complex(dp) :: I01, I21
        real(dp) :: A1, A2

        A1 = spec%A1_cc(j)
        A2 = spec%A2_cc(j)
        I01 = spec%I01_cc(j, mphi)
        I21 = spec%I21_cc(j, mphi)

        val = (I01 * (A1 + A2 * (1.0d0 - mphi)) + 0.5d0 * A2 * I21) * FP_kappa_rho_B(j, spec)

    end function


    function FP_G2_rho_B(j, spec, mphi) result(val)

        use species_m, only: species_t
        use KIM_kinds_m, only: dp

        implicit none

        integer, intent(in) :: j, mphi
        type(species_t), intent(in) :: spec
        complex(dp) :: val
        complex(dp) :: I01
        real(dp) :: A2

        A2 = spec%A2_cc(j)
        I01 = spec%I01_cc(j, mphi)

        val = - I01 * A2 * FP_kappa_rho_B(j, spec)

    end function


    function FP_G3_rho_B(j, spec, mphi) result(val)

        use species_m, only: species_t
        use KIM_kinds_m, only: dp

        implicit none

        integer, intent(in) :: j, mphi
        type(species_t), intent(in) :: spec
        complex(dp) :: val
        complex(dp) :: I01
        real(dp) :: A2

        A2 = spec%A2_cc(j)
        I01 = spec%I01_cc(j, mphi)

        val = I01 * A2 * FP_kappa_rho_B(j, spec)

    end function


    function FP_G1_j_phi(j, spec, mphi) result(val)

        use species_m, only: species_t, plasma
        use KIM_kinds_m, only: dp

        implicit none

        integer, intent(in) :: j, mphi
        type(species_t), intent(in) :: spec
        complex(dp) :: val
        complex(dp) :: I10, I12
        real(dp) :: A1, A2

        A1 = spec%A1_cc(j)
        A2 = spec%A2_cc(j)
        I10 = spec%I10_cc(j, mphi)
        I12 = spec%I12_cc(j, mphi)

        val = (I10 * (A1 + A2) + 0.5d0 * A2 * I12) * FP_kappa_j_phi(j, spec) * plasma%ks_cc(j)

    end function


    function FP_G2_j_phi(j, spec, mphi) result(val)

        use species_m, only: species_t, plasma
        use KIM_kinds_m, only: dp

        implicit none

        integer, intent(in) :: j, mphi
        type(species_t), intent(in) :: spec
        complex(dp) :: val
        complex(dp) :: I10
        real(dp) :: A2

        A2 = spec%A2_cc(j)
        I10 = spec%I10_cc(j, mphi)

        val = - I10 * A2 * FP_kappa_j_phi(j, spec) * plasma%ks_cc(j)

    end function


    function FP_G3_j_phi(j, spec, mphi) result(val)

        use species_m, only: species_t, plasma
        use KIM_kinds_m, only: dp

        implicit none

        integer, intent(in) :: j, mphi
        type(species_t), intent(in) :: spec
        complex(dp) :: val
        complex(dp) :: I10
        real(dp) :: A2

        A2 = spec%A2_cc(j)
        I10 = spec%I01_cc(j, mphi)

        val = I10 * A2 * FP_kappa_j_phi(j, spec) * plasma%ks_cc(j)

    end function


    function FP_G1_j_B(j, spec, mphi) result(val)

        use species_m, only: species_t
        use KIM_kinds_m, only: dp

        implicit none

        integer, intent(in) :: j, mphi
        type(species_t), intent(in) :: spec
        complex(dp) :: val
        complex(dp) :: I11, I13
        real(dp) :: A1, A2

        A1 = spec%A1_cc(j)
        A2 = spec%A2_cc(j)
        I11 = spec%I11_cc(j, mphi)
        I13 = spec%I13_cc(j, mphi)

        val = (I11 * (A1 + A2) + 0.5d0 * A2 * I13) * FP_kappa_j_B(j, spec)

    end function


    function FP_G2_j_B(j, spec, mphi) result(val)

        use species_m, only: species_t
        use KIM_kinds_m, only: dp

        implicit none

        integer, intent(in) :: j, mphi
        type(species_t), intent(in) :: spec
        complex(dp) :: val
        complex(dp) :: I11
        real(dp) :: A2

        A2 = spec%A2_cc(j)
        I11 = spec%I11_cc(j, mphi)

        val = - I11 * A2 * FP_kappa_j_B(j, spec)

    end function


    function FP_G3_j_B(j, spec, mphi) result(val)

        use species_m, only: species_t
        use KIM_kinds_m, only: dp

        implicit none

        integer, intent(in) :: j, mphi
        type(species_t), intent(in) :: spec
        complex(dp) :: val
        complex(dp) :: I11
        real(dp) :: A2

        A2 = spec%A2_cc(j)
        I11 = spec%I11_cc(j, mphi)

        val = I11 * A2 * FP_kappa_j_B(j, spec)

    end function

end module
