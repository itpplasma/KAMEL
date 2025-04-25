module electrostatic_kernel

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

    subroutine fill_kernel_phi_semi_analytic(kernel_rho_phi_llp, kernel_rho_B_llp)

        use KIM_kinds, only: dp
        use gsl_mod, only: erf => gsl_sf_erf
        use electrostatic_integrals, only: gauss_config_t, init_gauss_int

        implicit none

        type(kernel_spl_t), intent(inout) :: kernel_rho_phi_llp
        type(kernel_spl_t), intent(inout) :: kernel_rho_B_llp
        type(gauss_config_t) :: gauss_conf
        integer :: l, lp
        complex(dp) :: kernel_phi_llp, kernel_B_llp

        gauss_conf%n = 10
        call init_gauss_int(gauss_conf)

        !$omp parallel do collapse(2) private(l,lp, kernel_phi_llp, kernel_B_llp)
        do l = 1, kernel_rho_phi_llp%npts_l
            do lp = 1, kernel_rho_phi_llp%npts_lp
                if (abs(l - lp) > 1) cycle

                call calc_kernel_rho_semi_analytic(l, lp, kernel_phi_llp, kernel_B_llp, gauss_conf)
                kernel_rho_phi_llp%Kllp(l, lp) = kernel_phi_llp
                kernel_rho_B_llp%Kllp(l, lp) = kernel_B_llp

                if (isnan(real(kernel_phi_llp))) then
                    print *, "kernel_llp is NaN for l = ", l, " lp = ", lp
                    print *, "kernel_llp = ", kernel_phi_llp
                    stop
                end if
                if (isnan(real(kernel_B_llp))) then
                    print *, "kernel_B_llp is NaN for l = ", l, " lp = ", lp
                    print *, "kernel_B_llp = ", kernel_B_llp
                    stop
                end if
            end do
        end do
        !$omp end parallel do

    end subroutine

    subroutine calc_kernel_rho_semi_analytic(l, lp, kernel_phi_llp, kernel_B_llp, gauss_conf)

        use KIM_kinds, only: dp
        use gsl_mod, only: erf => gsl_sf_erf
        use electrostatic_integrals, only: gauss_integrate_B0, gauss_integrate_B1, gauss_config_t
        use species, only: plasma
        use constants, only: pi, sol, com_unit
        use electrostatic_integrands, only: int_B0_rho_phi_t, mathcal_A0_rho_phi, int_B1_rho_phi_t, mathcal_A1_rho_phi,&
            mathcal_A1_rho_B, mathcal_A2_rho_B
        use config, only: artificial_debye_case
        
        implicit none

        integer, intent(in) :: l, lp
        complex(dp) :: kernel_phi_llp, kernel_B_llp
        integer :: j, sigma
        type(gauss_config_t), intent(in) :: gauss_conf
        real(dp) :: integral_val
        type(int_B0_rho_phi_t) :: int_B0
        type(int_B1_rho_phi_t) :: int_B1
        
        kernel_phi_llp = 0.0d0
        kernel_B_llp = 0.0d0

        call set_xl_at_edge(l, lp, int_B0, int_B1)

        do sigma = 0, plasma%n_species - 1
            do j = 2, size(plasma%r_grid)-1
                int_B0%j = j
                int_B0%rhoT = 0.5d0 * (plasma%spec(sigma)%rho_L(j) + plasma%spec(sigma)%rho_L(j+1))
                call gauss_integrate_B0(int_B0, int_B0%xlm1, int_B0%xlp1, integral_val, gauss_conf)
                kernel_phi_llp = kernel_phi_llp + integral_val * mathcal_A0_rho_phi(j, plasma%spec(sigma)) 

                if (artificial_debye_case) cycle

                int_B1%j = j
                int_B1%rhoT = 0.5d0 * (plasma%spec(sigma)%rho_L(j) + plasma%spec(sigma)%rho_L(j+1))
                call gauss_integrate_B1(int_B1, integral_val, gauss_conf)
                kernel_phi_llp = kernel_phi_llp - integral_val * mathcal_A1_rho_phi(j, plasma%spec(sigma))
                kernel_B_llp = kernel_B_llp - integral_val * mathcal_A1_rho_B(j, plasma%spec(sigma)) &
                                * (0.5d0 * (plasma%spec(sigma)%vT(j) + plasma%spec(sigma)%vT(j+1)))**2.0d0 &
                                / (0.5d0 * (plasma%spec(sigma)%lambda_D(j) + plasma%spec(sigma)%lambda_D(j+1)))**2.0d0 &
                                / (0.5d0 * (plasma%spec(sigma)%omega_c(j) + plasma%spec(sigma)%omega_c(j+1))) &
                                / abs(0.5d0 * (plasma%kp(j) + plasma%kp(j+1)))

            end do
        end do

        kernel_phi_llp = kernel_phi_llp / (8.0d0 * pi**3.0d0) 
        kernel_B_llp = kernel_B_llp / (8.0d0 * pi**3.0d0 * sol) * com_unit
            
    end subroutine

    subroutine set_xl_at_edge(l, lp, intB0, intB1)
        
        use grid, only: xl_grid
        use KIM_kinds, only: dp
        use electrostatic_integrands, only: int_B0_rho_phi_t, int_B1_rho_phi_t

        implicit none

        integer, intent(in) :: l, lp
        type(int_B0_rho_phi_t), intent(inout) :: intB0
        type(int_B1_rho_phi_t), intent(inout) :: intB1   

        intB0%xl = xl_grid%xb(l)
        intB0%xlp = xl_grid%xb(lp)
        intB1%xl = xl_grid%xb(l)
        intB1%xlp = xl_grid%xb(lp)

        ! handle kernel edges
        if (l == 1) then
            intB0%xlm1 = xl_grid%xb(l)
            intB1%xlm1 = xl_grid%xb(l)
        else
            intB0%xlm1 = xl_grid%xb(l-1)
            intB1%xlm1 = xl_grid%xb(l-1)
        end if
        if (lp == 1) then
            intB0%xlpm1 = xl_grid%xb(lp)
            intB1%xlpm1 = xl_grid%xb(lp)
        else
            intB0%xlpm1 = xl_grid%xb(lp-1)
            intB1%xlpm1 = xl_grid%xb(lp-1)
        end if
        if (l == xl_grid%npts_b) then
            intB0%xlp1 = xl_grid%xb(l)
            intB1%xlp1 = xl_grid%xb(l)
        else
            intB0%xlp1 = xl_grid%xb(l+1)
            intB1%xlp1 = xl_grid%xb(l+1)
        end if
        if (lp == xl_grid%npts_b) then
            intB0%xlpp1 = xl_grid%xb(lp)
            intB1%xlpp1 = xl_grid%xb(lp)
        else
            intB0%xlpp1 = xl_grid%xb(lp+1)
            intB1%xlpp1 = xl_grid%xb(lp+1)
        end if

    end subroutine

end module

