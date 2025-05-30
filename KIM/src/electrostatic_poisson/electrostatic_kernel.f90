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

    subroutine fill_kernel_phi(kernel_rho_phi_llp, kernel_rho_B_llp)

        use KIM_kinds, only: dp
        use electrostatic_integrals, only: gauss_config_t, init_gauss_int
        use grid, only: delta_l_max

        implicit none

        type(kernel_spl_t), intent(inout) :: kernel_rho_phi_llp
        type(kernel_spl_t), intent(inout) :: kernel_rho_B_llp
        type(gauss_config_t) :: gauss_conf
        integer :: l, lp
        complex(dp) :: kernel_phi_llp, kernel_B_llp

        gauss_conf%n = 5
        call init_gauss_int(gauss_conf)

        !$omp parallel do collapse(2) private(l,lp, kernel_phi_llp, kernel_B_llp)
        do l = 1, kernel_rho_phi_llp%npts_l
            do lp = 1, kernel_rho_phi_llp%npts_lp
                if (abs(l - lp) > delta_l_max) cycle

                call calc_kernel_rho_all_terms_at_once(l, lp, kernel_phi_llp, kernel_B_llp, gauss_conf)
                !call calc_kernel_rho_term_by_term(l, lp, kernel_phi_llp, kernel_B_llp, gauss_conf)
                kernel_rho_phi_llp%Kllp(l, lp) = kernel_phi_llp
                kernel_rho_B_llp%Kllp(l, lp) = kernel_B_llp

                if (isnan(real(kernel_phi_llp))) then
                    print *, "semi analytical kernel_llp is NaN for l = ", l, " lp = ", lp
                    print *, "semi analytical kernel_llp = ", kernel_phi_llp
                    stop
                end if
                if (isnan(real(kernel_B_llp))) then
                    print *, "semi analytical kernel_B_llp is NaN for l = ", l, " lp = ", lp
                    print *, "semi analytical kernel_B_llp = ", kernel_B_llp
                    stop
                end if
            end do
        end do
        !$omp end parallel do

    end subroutine

    subroutine calc_kernel_rho_term_by_term(l, lp, kernel_phi_llp, kernel_B_llp, gauss_conf)

        use KIM_kinds, only: dp
        use electrostatic_integrals, only: gauss_integrate_F0, gauss_integrate_F1, gauss_integrate_F2, gauss_integrate_F3,&
            gauss_config_t
        use species, only: plasma
        use constants, only: pi, sol, com_unit
        use electrostatic_integrands, only: int_F0_rho_phi_t, int_F1_rho_phi_t, int_F2_rho_phi_t, int_F3_rho_phi_t, &
            integration_point_t
        use kernel_plasma_prefacs, only: G1_rho_phi, G1_rho_B, G2_rho_B, G3_rho_B, G0_rho_phi, G2_rho_phi, G3_rho_phi, &
            kappa_rho_phi, kappa_rho_B
        use config, only: artificial_debye_case
        
        implicit none

        integer, intent(in) :: l, lp
        complex(dp) :: kernel_phi_llp, kernel_B_llp
        integer :: j, sigma
        type(gauss_config_t), intent(in) :: gauss_conf
        real(dp) :: integral_val

        type(integration_point_t) :: int_point
        type(int_F0_rho_phi_t) :: int_F0
        type(int_F1_rho_phi_t) :: int_F1
        type(int_F2_rho_phi_t) :: int_F2
        type(int_F3_rho_phi_t) :: int_F3
        
        kernel_phi_llp = 0.0d0
        kernel_B_llp = 0.0d0

        call set_xl_at_edge(l, lp, int_point)

        do sigma = 0, plasma%n_species - 1
            do j = 2, size(plasma%r_grid)-1
                int_point%j = j
                int_point%rhoT = 0.5d0 * (plasma%spec(sigma)%rho_L(j) + plasma%spec(sigma)%rho_L(j+1))
                int_F0%int_point = int_point

                if (abs(l - lp) <= 1) then
                    call gauss_integrate_F0(int_F0, int_point%xlm1, int_point%xlp1, integral_val, gauss_conf)
                    kernel_phi_llp = kernel_phi_llp &
                                    + integral_val * G0_rho_phi(j, plasma%spec(sigma)) * kappa_rho_phi(j, plasma%spec(sigma))
                    integral_val = 0.0d0
                end if
                

                if (.not. artificial_debye_case) then
                    int_F1%int_point = int_point
                    int_F2%int_point = int_point
                    int_F3%int_point = int_point

                    call gauss_integrate_F1(int_F1, integral_val, gauss_conf)
                    kernel_phi_llp = kernel_phi_llp + integral_val * G1_rho_phi(j, plasma%spec(sigma)) * kappa_rho_phi(j, plasma%spec(sigma))
                    kernel_B_llp = kernel_B_llp + integral_val * G1_rho_B(j, plasma%spec(sigma)) * kappa_rho_B(j, plasma%spec(sigma))
                    integral_val = 0.0d0

                    call gauss_integrate_F2(int_F2, integral_val, gauss_conf)
                    kernel_phi_llp = kernel_phi_llp + integral_val * kappa_rho_phi(j, plasma%spec(sigma)) * G2_rho_phi(j, plasma%spec(sigma))
                    kernel_B_llp = kernel_B_llp + integral_val * G2_rho_B(j, plasma%spec(sigma)) * kappa_rho_B(j, plasma%spec(sigma))
                    integral_val = 0.0d0

                    call gauss_integrate_F3(int_F3, integral_val, gauss_conf)
                    kernel_phi_llp = kernel_phi_llp + integral_val * kappa_rho_phi(j, plasma%spec(sigma)) * G3_rho_phi(j, plasma%spec(sigma))
                    kernel_B_llp = kernel_B_llp + integral_val * G3_rho_B(j, plasma%spec(sigma)) * kappa_rho_B(j, plasma%spec(sigma))
                    integral_val = 0.0d0
                end if
                
            end do
        end do

        kernel_phi_llp = kernel_phi_llp / (8.0d0 * pi**3.0d0)  /2.0d0 ! factor 1/2 is somehow missing. Including this factor nicely reproduces the debye case.
        kernel_B_llp = kernel_B_llp / (8.0d0 * pi**3.0d0 * sol) * com_unit
            
    end subroutine

    subroutine calc_kernel_rho_all_terms_at_once(l, lp, kernel_phi_llp, kernel_B_llp, gauss_conf)

        use KIM_kinds, only: dp
        use electrostatic_integrals, only: gauss_integrate_rho_all_simultaneously, gauss_config_t
        use species, only: plasma
        use constants, only: pi, sol, com_unit
        use electrostatic_integrands, only: int_F0_rho_phi_t, int_F1_rho_phi_t, int_F2_rho_phi_t, int_F3_rho_phi_t, &
            integration_point_t
        use kernel_plasma_prefacs, only: G1_rho_phi, G1_rho_B, G2_rho_B, G3_rho_B, G0_rho_phi, G2_rho_phi, G3_rho_phi, &
            kappa_rho_phi, kappa_rho_B
        use config, only: artificial_debye_case
        
        implicit none

        integer, intent(in) :: l, lp
        complex(dp) :: kernel_phi_llp, kernel_B_llp
        integer :: j, sigma
        type(gauss_config_t), intent(in) :: gauss_conf
        complex(dp) :: integral_val_phi, integral_val_B

        type(integration_point_t) :: int_point
        
        kernel_phi_llp = 0.0d0
        kernel_B_llp = 0.0d0

        call set_xl_at_edge(l, lp, int_point)

        do sigma = 0, plasma%n_species - 1
            do j = 2, size(plasma%r_grid)-1
                int_point%j = j
                int_point%rhoT = 0.5d0 * (plasma%spec(sigma)%rho_L(j) + plasma%spec(sigma)%rho_L(j+1))

                call gauss_integrate_rho_all_simultaneously(int_point, plasma%spec(sigma), &
                    integral_val_phi, integral_val_B, gauss_conf)
                kernel_phi_llp = kernel_phi_llp + integral_val_phi * kappa_rho_phi(j, plasma%spec(sigma))
                kernel_B_llp = kernel_B_llp + integral_val_B * kappa_rho_B(j, plasma%spec(sigma))
                
            end do
        end do

        kernel_phi_llp = kernel_phi_llp / (8.0d0 * pi**3.0d0)  /2.0d0 ! factor 1/2 is somehow missing. Including this factor nicely reproduces the debye case.
        kernel_B_llp = kernel_B_llp / (8.0d0 * pi**3.0d0 * sol) * com_unit
            
    end subroutine


    subroutine set_xl_at_edge(l, lp, int_point)
        
        use grid, only: xl_grid
        use KIM_kinds, only: dp
        use electrostatic_integrands, only: integration_point_t

        implicit none

        integer, intent(in) :: l, lp
        type(integration_point_t), intent(inout) :: int_point

        int_point%xl = xl_grid%xb(l)
        int_point%xlp = xl_grid%xb(lp)

        ! handle kernel edges
        if (l == 1) then
            int_point%xlm1 = xl_grid%xb(l)
        else
            int_point%xlm1 = xl_grid%xb(l-1)
        end if
        if (lp == 1) then
            int_point%xlpm1 = xl_grid%xb(lp)
        else
            int_point%xlpm1 = xl_grid%xb(lp-1)
        end if
        if (l == xl_grid%npts_b) then
            int_point%xlp1 = xl_grid%xb(l)
        else
            int_point%xlp1 = xl_grid%xb(l+1)
        end if
        if (lp == xl_grid%npts_b) then
            int_point%xlpp1 = xl_grid%xb(lp)
        else
            int_point%xlpp1 = xl_grid%xb(lp+1)
        end if

    end subroutine

end module

