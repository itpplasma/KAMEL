module reduced_kernel

    use KIM_kinds, only: dp

    implicit none

    type :: kernel_spl_t
        integer :: npts_l, npts_lp
        complex(dp), allocatable :: Kllp(:,:)
        contains
            procedure :: init_kernel
    end type kernel_spl_t

    contains

    subroutine init_kernel(this, npts_l, npts_lp)

        implicit none

        class(kernel_spl_t), intent(inout) :: this
        integer, intent(in) :: npts_l, npts_lp

        this%npts_l = npts_l
        this%npts_lp = npts_lp
        allocate(this%Kllp(npts_l, npts_lp))
        this%Kllp = 0.0d0

    end subroutine init_kernel

    subroutine fill_kernel_phi(kernel_rho_phi_llp)

        use KIM_kinds, only: dp
        use gsl_mod, only: erf => gsl_sf_erf
        use gauss_quad, only: gauss_config_t, init_gauss_int

        implicit none

        type(kernel_spl_t), intent(inout) :: kernel_rho_phi_llp
        type(gauss_config_t) :: gauss_conf
        integer :: l, lp
        complex(dp) :: kernel_llp

        gauss_conf%n = 10
        call init_gauss_int(gauss_conf)

        !$omp parallel do collapse(2) private(l,lp, kernel_llp)
        do l = 2, kernel_rho_phi_llp%npts_l-1
            do lp = 2, kernel_rho_phi_llp%npts_lp-1
                if (abs(l - lp) > 1) cycle
                call calc_kernel_rho_phi(l, lp, kernel_llp, gauss_conf)
                kernel_rho_phi_llp%Kllp(l, lp) = kernel_llp
                if (isnan(real(kernel_llp))) then
                    print *, "kernel_llp is NaN for l = ", l, " lp = ", lp
                    print *, "kernel_llp = ", kernel_llp
                    stop
                end if
            end do
        end do
        !$omp end parallel do

    end subroutine fill_kernel_phi

    subroutine calc_kernel_rho_phi(l, lp, kernel_llp, gauss_conf)

        use KIM_kinds, only: dp
        use grid, only: xl_grid, rg_grid
        use gsl_mod, only: erf => gsl_sf_erf
        use gauss_quad, only: gauss_integrate_B0, gauss_integrate_B1, gauss_config_t
        use config, only: number_of_ion_species
        use species, only: plasma
        use constants, only: pi
        use reduced_integrands, only: int_B0_t, mathcal_A0, int_B1_t, mathcal_A1
        
        implicit none

        integer, intent(in) :: l, lp
        complex(dp) :: kernel_llp
        integer :: j, sigma
        type(gauss_config_t), intent(in) :: gauss_conf
        real(dp) :: integral_val
        type(int_B0_t) :: int_B0
        type(int_B1_t) :: int_B1
        
        kernel_llp = 0.0d0
        int_B0%l = l
        int_B0%lp = lp

        int_B1%l = l
        int_B1%lp = lp

        do sigma = 0, number_of_ion_species
            do j = 2, xl_grid%npts_b-1
                int_B0%j = j
                int_B0%rhoT = 0.5d0 * (plasma%spec(sigma)%rho_L(j) + plasma%spec(sigma)%rho_L(j+1))
                !rhoT = 0.5d0 * (plasma%spec(sigma)%rho_L(j) + plasma%spec(sigma)%rho_L(j+1))
                call gauss_integrate_B0(int_B0, rg_grid%xb(j-1), rg_grid%xb(j+1), integral_val, gauss_conf)
                kernel_llp = kernel_llp + integral_val * mathcal_A0(j, plasma%spec(sigma)) 

                int_B1%j = j
                int_B1%rhoT = 0.5d0 * (plasma%spec(sigma)%rho_L(j) + plasma%spec(sigma)%rho_L(j+1))
                call gauss_integrate_B1(int_B1, integral_val, gauss_conf)
                kernel_llp = kernel_llp + integral_val * mathcal_A1(j, plasma%spec(sigma))
            end do
        end do

        kernel_llp = kernel_llp / (8.0d0 * pi**3.0d0) 
            
    end subroutine

end module

