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
        real(dp) :: val
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
        use constants, only: pi
        use KIM_kinds, only: dp

        implicit none

        integer, intent(in) :: j
        type(species_t), intent(in) :: spec
        real(dp) :: val
        !real(dp), dimension(dim), intent(in) :: vT, nu
        real(dp), dimension(:), allocatable :: x1, x2
        real(dp), dimension(:), allocatable :: epm2, brm2, epbr_re, epbr_im
        complex(dp), dimension(:, :, :), allocatable :: symbI

        val = -1.0d0 


    end function


    subroutine calc_susc_funcs(symbI, j, spec)

        use KIM_kinds, only: dp
        use grid_mod, only: rg_grid
        use resonances_mod, only: r_res, gg_width
        use species, only: plasma, species_t

        implicit none

        complex(dp), dimension(:, :, :), intent(inout) :: symbI
        type(species_t), intent(in) :: spec
        real(dp) :: x1, x2
        integer :: i

        symbI = 0.d0

        x1 = (plasma%kp(j+1) + plasma%kp(j)) * 0.5d0 &
            *(spec%vT(j+1) + spec%vT(j)) * 0.5d0 &
            /((spec%nu(j+1) + spec%nu(j)) * 0.5d0)
        x2 = -(plasma%om_E(j+1) + plasma%om_E(j)) * 0.5 &
            /((spec%nu(j+1) + spec%nu(j)) * 0.5d0)
    
        !if (rg_grid%xb(j) .lt. r_res - 2.d0*gg_width) cycle
        !if (rg_grid%xb(j) .gt. r_res + 2.d0*gg_width) cycle
        call getIfunc(x1, x2, symbI(:, :, i))

    end subroutine

end module