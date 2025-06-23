! module for plasma factors in the Fokker-Planck case
module FP_kernel_plasma_prefacs

    implicit none

    contains

    function FP_kappa_rho_phi(j, spec) result(val)

        use species, only: plasma, species_t
        use constants, only: pi
        use KIM_kinds, only: dp

        implicit none

        integer, intent(in) :: j
        type(species_t), intent(in) :: spec
        real(dp) :: val
        real(dp) :: lambda

        lambda = 0.5d0 * (spec%lambda_D(j) + spec%lambda_D(j+1))
        val = 1.0d0 / (lambda**2.0d0)  !/ sqrt(2.0d0)

    end function

    function FP_kappa_rho_B(j, spec) result(val)

        use species, only: plasma, species_t
        use KIM_kinds, only: dp
        use constants, only: com_unit, sol

        implicit none

        integer, intent(in) :: j
        type(species_t), intent(in) :: spec
        real(dp) :: val, lambda, vT, nu, omega_c

        lambda = 0.5d0 * (spec%lambda_D(j) + spec%lambda_D(j+1))
        vT = 0.5d0 * (spec%vT(j) + spec%vT(j+1))
        nu = 0.5d0 * (spec%nu(j) + spec%nu(j+1))
        omega_c = 0.5d0 * (spec%omega_c(j) + spec%omega_c(j+1))

        val = - vT**2.0d0 / (lambda**2.0d0 * omega_c * nu * sol)

    end function

    function FP_G0_rho_phi(j, spec) result(val)

        use species, only: plasma, species_t
        use constants, only: pi
        use KIM_kinds, only: dp

        implicit none

        integer, intent(in) :: j
        type(species_t), intent(in) :: spec
        real(dp) :: val

        val = -1.0d0 

    end function

    function FP_G1_rho_phi(j, spec) result(val)

        use species, only: plasma, species_t
        use KIM_kinds, only: dp
        use constants, only: com_unit

        implicit none

        integer, intent(in) :: j
        type(species_t), intent(in) :: spec
        complex(dp) :: val, I00, I20
        real(dp) :: A1, A2

        real(dp) :: lambda, vT, nu, omega_c, ks

        lambda = 0.5d0 * (spec%lambda_D(j) + spec%lambda_D(j+1))
        vT = 0.5d0 * (spec%vT(j) + spec%vT(j+1))
        nu = 0.5d0 * (spec%nu(j) + spec%nu(j+1))
        omega_c = 0.5d0 * (spec%omega_c(j) + spec%omega_c(j+1))
        ks = 0.5d0 * (plasma%ks(j) + plasma%ks(j+1))

        A1 = 0.5d0 * (spec%A1(j) + spec%A1(j+1))
        A2 = 0.5d0 * (spec%A2(j) + spec%A2(j+1))
        I00 = 0.5d0 * (spec%I00(j) + spec%I00(j+1))
        I20 = 0.5d0 * (spec%I20(j) + spec%I20(j+1))

        val = I00 * (A1 + A2) + 0.5d0 * A2 * I20
        val = val * com_unit * vT**2.0d0 * ks / (lambda**2.0d0 * omega_c * nu)

    end function


    function FP_G2_rho_phi(j, spec) result(val)

        use species, only: plasma, species_t
        use KIM_kinds, only: dp
        use constants, only: com_unit

        implicit none

        integer, intent(in) :: j
        type(species_t), intent(in) :: spec
        complex(dp) :: val, I00
        real(dp) :: A2
        real(dp) :: lambda, vT, nu, omega_c, ks

        lambda = 0.5d0 * (spec%lambda_D(j) + spec%lambda_D(j+1))
        vT = 0.5d0 * (spec%vT(j) + spec%vT(j+1))
        nu = 0.5d0 * (spec%nu(j) + spec%nu(j+1))
        omega_c = 0.5d0 * (spec%omega_c(j) + spec%omega_c(j+1))
        ks = 0.5d0 * (plasma%ks(j) + plasma%ks(j+1))

        A2 = 0.5d0 * (spec%A2(j) + spec%A2(j+1))
        I00 = 0.5d0 * (spec%I00(j) + spec%I00(j+1))

        val = I00 * A2 * com_unit * vT**2.0d0 * ks / (lambda**2.0d0 * omega_c * nu)

    end function


    function FP_G3_rho_phi(j, spec) result(val)

        use species, only: plasma, species_t
        use KIM_kinds, only: dp
        use constants, only: com_unit

        implicit none

        integer, intent(in) :: j
        type(species_t), intent(in) :: spec
        complex(dp) :: val, I00
        real(dp) :: A2

        real(dp) :: lambda, vT, nu, omega_c, ks

        lambda = 0.5d0 * (spec%lambda_D(j) + spec%lambda_D(j+1))
        vT = 0.5d0 * (spec%vT(j) + spec%vT(j+1))
        nu = 0.5d0 * (spec%nu(j) + spec%nu(j+1))
        omega_c = 0.5d0 * (spec%omega_c(j) + spec%omega_c(j+1))
        ks = 0.5d0 * (plasma%ks(j) + plasma%ks(j+1))

        A2 = 0.5d0 * (spec%A2(j) + spec%A2(j+1))
        I00 = 0.5d0 * (spec%I00(j) + spec%I00(j+1))

        val = I00 * A2 * com_unit * vT**2.0d0 * ks / (lambda**2.0d0 * omega_c * nu)

    end function


    function FP_G1_rho_B(j, spec) result(val)

        use species, only: plasma, species_t
        use KIM_kinds, only: dp

        implicit none

        integer, intent(in) :: j
        type(species_t), intent(in) :: spec
        complex(dp) :: val, I01, I21
        real(dp) :: A1, A2

        A1 = 0.5d0 * (spec%A1(j) + spec%A1(j+1))
        A2 = 0.5d0 * (spec%A2(j) + spec%A2(j+1))
        I01 = 0.5d0 * (spec%I01(j) + spec%I01(j+1))
        I21 = 0.5d0 * (spec%I21(j) + spec%I21(j+1))

        val = I01 * (A1 + A2) + 0.5d0 * A2 * I21

    end function


    function FP_G2_rho_B(j, spec) result(val)

        use species, only: plasma, species_t
        use KIM_kinds, only: dp

        implicit none

        integer, intent(in) :: j
        type(species_t), intent(in) :: spec
        complex(dp) :: val, I01
        real(dp) :: A2

        A2 = 0.5d0 * (spec%A2(j) + spec%A2(j+1))
        I01 = 0.5d0 * (spec%I01(j) + spec%I01(j+1))

        val = I01 * A2

    end function


    function FP_G3_rho_B(j, spec) result(val)

        use species, only: plasma, species_t
        use KIM_kinds, only: dp

        implicit none

        integer, intent(in) :: j
        type(species_t), intent(in) :: spec
        complex(dp) :: val
        real(dp) :: A2

        A2 = 0.5d0 * (spec%A2(j) + spec%A2(j+1))

        val = spec%symbI(0,1) * A2

    end function

end module