! module for plasma factors in the Fokker-Planck case
module FP_kernel_plasma_prefacs

    implicit none

    contains

    function FP_kappa_rho_phi(j, spec) result(val)

        use species, only: plasma, species_t
        use constants, only: pi, com_unit
        use KIM_kinds, only: dp

        implicit none

        integer, intent(in) :: j
        type(species_t), intent(in) :: spec
        complex(dp) :: val
        real(dp) :: lambda, vT, nu, omega_c, ks

        lambda = 0.5d0 * (spec%lambda_D(j) + spec%lambda_D(j+1))
        vT = 0.5d0 * (spec%vT(j) + spec%vT(j+1))
        nu = 0.5d0 * (spec%nu(j) + spec%nu(j+1))
        omega_c = 0.5d0 * (spec%omega_c(j) + spec%omega_c(j+1))
        ks = 0.5d0 * (plasma%ks(j) + plasma%ks(j+1))

        val = com_unit * vT**2.0d0 * ks / (lambda**2.0d0 * omega_c * nu)

    end function

    function FP_kappa_rho_B(j, spec) result(val)

        use species, only: plasma, species_t
        use KIM_kinds, only: dp

        implicit none

        integer, intent(in) :: j
        type(species_t), intent(in) :: spec
        real(dp) :: val, lambda, vT, nu, omega_c

        lambda = 0.5d0 * (spec%lambda_D(j) + spec%lambda_D(j+1))
        vT = 0.5d0 * (spec%vT(j) + spec%vT(j+1))
        nu = 0.5d0 * (spec%nu(j) + spec%nu(j+1))
        omega_c = 0.5d0 * (spec%omega_c(j) + spec%omega_c(j+1))

        val = - vT**2.0d0 / (lambda**2.0d0 * omega_c * nu)

    end function

    function FP_G1_rho_phi(j, spec) result(val)

        use species, only: plasma, species_t
        use KIM_kinds, only: dp

        implicit none

        integer, intent(in) :: j
        type(species_t), intent(in) :: spec
        complex(dp) :: val
        real(dp) :: A1, A2

        A1 = 0.5d0 * (spec%A1(j) + spec%A1(j+1))
        A2 = 0.5d0 * (spec%A2(j) + spec%A2(j+1))

        val = spec%symbI(0,0) * (A1 + A2) + 0.5d0 * A2 * spec%symbI(2,0)

    end function


    function FP_G2_rho_phi(j, spec) result(val)

        use species, only: plasma, species_t
        use KIM_kinds, only: dp

        implicit none

        integer, intent(in) :: j
        type(species_t), intent(in) :: spec
        complex(dp) :: val
        real(dp) :: A2

        A2 = 0.5d0 * (spec%A2(j) + spec%A2(j+1))

        val = spec%symbI(0,0) * A2

    end function


    function FP_G3_rho_phi(j, spec) result(val)

        use species, only: plasma, species_t
        use KIM_kinds, only: dp

        implicit none

        integer, intent(in) :: j
        type(species_t), intent(in) :: spec
        complex(dp) :: val
        real(dp) :: A2

        A2 = 0.5d0 * (spec%A2(j) + spec%A2(j+1))

        val = spec%symbI(0,0) * A2

    end function


    function FP_G1_rho_B(j, spec) result(val)

        use species, only: plasma, species_t
        use KIM_kinds, only: dp

        implicit none

        integer, intent(in) :: j
        type(species_t), intent(in) :: spec
        complex(dp) :: val
        real(dp) :: A1, A2

        A1 = 0.5d0 * (spec%A1(j) + spec%A1(j+1))
        A2 = 0.5d0 * (spec%A2(j) + spec%A2(j+1))

        val = spec%symbI(0,1) * (A1 + A2) + 0.5d0 * A2 * spec%symbI(2,1)

    end function


    function FP_G2_rho_B(j, spec) result(val)

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