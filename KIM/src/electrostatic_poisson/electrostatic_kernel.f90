module electrostatic_kernel_m

    use KIM_kinds_m, only: dp

    implicit none
    
    ! Diagnostic variables to track maximum and minimum distances
    real(dp) :: max_distance_xl_xlp = 0.0d0
    real(dp) :: min_distance_xl_xlp = 0.0d0
    integer :: max_index_distance = 0
    integer :: min_index_distance = 0
    integer :: max_dist_l = 0, max_dist_lp = 0
    integer :: min_dist_l = 0, min_dist_lp = 0
    integer :: max_idx_l = 0, max_idx_lp = 0
    integer :: min_idx_l = 0, min_idx_lp = 0

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

    subroutine Krook_fill_kernel_phi(K_rho_phi_llp, K_rho_B_llp)

        use KIM_kinds_m, only: dp
        use electrostatic_integrals_gauss_mod, only: gauss_config_t, init_gauss_int
        use grid_m, only: gauss_int_nodes_Ntheta, gauss_int_nodes_Nx, gauss_int_nodes_Nxp

        implicit none

        type(kernel_spl_t), intent(inout) :: K_rho_phi_llp
        type(kernel_spl_t), intent(inout) :: K_rho_B_llp
        type(gauss_config_t) :: gauss_conf
        integer :: l, lp
        complex(dp) :: k_rho_phi, k_rho_B

        gauss_conf%Nx = gauss_int_nodes_Nx
        gauss_conf%Nxp = gauss_int_nodes_Nxp
        gauss_conf%Ntheta = gauss_int_nodes_Ntheta
        call init_gauss_int(gauss_conf)

        write(*,*) 'Filling Krook collision kernels...'

        !$omp parallel do collapse(1) private(l,lp, k_rho_phi, k_rho_B)
        do l = 1, K_rho_phi_llp%npts_l
            do lp = 1, l

                call Krook_calc_kernel_rho_term_by_term(l, lp, k_rho_phi, k_rho_B, gauss_conf)
                K_rho_phi_llp%Kllp(l, lp) = k_rho_phi
                K_rho_B_llp%Kllp(l, lp) = k_rho_B

                if (isnan(real(k_rho_phi))) then
                    print *, "semi analytical kernel_llp is NaN for l = ", l, " lp = ", lp
                    print *, "semi analytical kernel_llp = ", k_rho_phi
                    stop
                end if
                if (isnan(real(k_rho_B))) then
                    print *, "semi analytical k_rho_B is NaN for l = ", l, " lp = ", lp
                    print *, "semi analytical k_rho_B = ", k_rho_B
                    stop
                end if
            
                K_rho_phi_llp%Kllp(lp, l) = K_rho_phi_llp%Kllp(l, lp)
                K_rho_B_llp%Kllp(lp, l) = K_rho_B_llp%Kllp(l, lp)

            end do
        end do
        !$omp end parallel do
        
        write(*,*) '======== Kernel Distance Diagnostics (Krook) ========'
        write(*,'(A,F12.6)') ' Maximum |xl - xlp| distance: ', max_distance_xl_xlp
        write(*,'(A,I6,A,I6)') ' Occurred at l = ', max_dist_l, ', lp = ', max_dist_lp
        write(*,'(A,F12.6)') ' Minimum |xl - xlp| distance: ', min_distance_xl_xlp
        write(*,'(A,I6,A,I6)') ' Occurred at l = ', min_dist_l, ', lp = ', min_dist_lp
        write(*,'(A,I6)') ' Maximum index distance |l - lp|: ', max_index_distance
        write(*,'(A,I6,A,I6)') ' Occurred at l = ', max_idx_l, ', lp = ', max_idx_lp
        write(*,'(A,I6)') ' Minimum index distance |l - lp|: ', min_index_distance
        write(*,'(A,I6,A,I6)') ' Occurred at l = ', min_idx_l, ', lp = ', min_idx_lp
        write(*,*) '===================================================='

    end subroutine

    subroutine Krook_calc_kernel_rho_term_by_term(l, lp, k_rho_phi, k_rho_B, gauss_conf)

        use KIM_kinds_m, only: dp
        use electrostatic_integrals_gauss_mod, only: gauss_integrate_F0, gauss_integrate_F1, gauss_integrate_F2, gauss_integrate_F3,&
            gauss_config_t
        use species_m, only: plasma
        use constants_m, only: pi
        use electrostatic_integrands_gauss_mod, only: gauss_int_F0_rho_phi_t, gauss_int_F1_rho_phi_t, gauss_int_F2_rho_phi_t, gauss_int_F3_rho_phi_t, &
            integration_point_t
        use Krook_kernel_plasma_prefacs_m, only: Krook_G0_rho_phi, Krook_G1_rho_phi, Krook_G2_rho_phi, Krook_G3_rho_phi, &
            Krook_G1_rho_B, Krook_G2_rho_B, Krook_G3_rho_B, Krook_kappa_rho_phi, Krook_kappa_rho_B
        use config_m, only: artificial_debye_case, turn_off_ions
        use grid_m, only: Larmor_skip_factor
        
        implicit none

        integer, intent(in) :: l, lp
        complex(dp) :: k_rho_phi, k_rho_B
        integer :: j, sigma
        type(gauss_config_t), intent(in) :: gauss_conf
        real(dp) :: integral_val
        real(dp) :: current_distance
        integer :: current_idx_distance

        type(integration_point_t) :: int_point
        type(gauss_int_F0_rho_phi_t) :: int_F0
        type(gauss_int_F1_rho_phi_t) :: int_F1
        type(gauss_int_F2_rho_phi_t) :: int_F2
        type(gauss_int_F3_rho_phi_t) :: int_F3
        
        k_rho_phi = 0.0d0
        k_rho_B = 0.0d0

        call set_xl_at_edge(l, lp, int_point)

        do sigma = 0, plasma%n_species - 1
            if (turn_off_ions .and. sigma >= 1) cycle
            do j = 2, size(plasma%r_grid)-1
                int_point%j = j
                int_point%rhoT = 0.5d0 * (plasma%spec(sigma)%rho_L(j) + plasma%spec(sigma)%rho_L(j+1))
                int_F0%int_point = int_point

                if (l == lp) then
                    call gauss_integrate_F0(int_F0, int_point%xlm1, int_point%xlp1, integral_val, gauss_conf)
                    k_rho_phi = k_rho_phi &
                                    + integral_val * Krook_G0_rho_phi(j, plasma%spec(sigma)) * Krook_kappa_rho_phi(j, plasma%spec(sigma))
                end if

                ! Track maximum distances for diagnostics
                current_distance = abs(int_point%xl - int_point%xlp)

                if (current_distance > Larmor_skip_factor * int_point%rhoT) cycle

                current_idx_distance = abs(l - lp)
                
                !$omp critical
                if (current_distance > max_distance_xl_xlp) then
                    max_distance_xl_xlp = current_distance
                    max_dist_l = l
                    max_dist_lp = lp
                end if
                if (current_distance < min_distance_xl_xlp .and. current_distance > 0.0d0) then
                    min_distance_xl_xlp = current_distance
                    min_dist_l = l
                    min_dist_lp = lp
                end if
                if (current_idx_distance > max_index_distance) then
                    max_index_distance = current_idx_distance
                    max_idx_l = l
                    max_idx_lp = lp
                end if
                if (current_idx_distance < min_index_distance .and. current_idx_distance > 0) then
                    min_index_distance = current_idx_distance
                    min_idx_l = l
                    min_idx_lp = lp
                end if
                !$omp end critical
                
                if (.not. artificial_debye_case) then
                    int_F1%int_point = int_point
                    int_F2%int_point = int_point
                    int_F3%int_point = int_point

                    call gauss_integrate_F1(int_F1, integral_val, gauss_conf)
                    k_rho_phi = k_rho_phi + integral_val * Krook_G1_rho_phi(j, plasma%spec(sigma)) * Krook_kappa_rho_phi(j, plasma%spec(sigma))
                    k_rho_B = k_rho_B + integral_val * Krook_G1_rho_B(j, plasma%spec(sigma)) * Krook_kappa_rho_B(j, plasma%spec(sigma))

                    call gauss_integrate_F2(int_F2, integral_val, gauss_conf)
                    k_rho_phi = k_rho_phi + integral_val * Krook_kappa_rho_phi(j, plasma%spec(sigma)) * Krook_G2_rho_phi(j, plasma%spec(sigma))
                    k_rho_B = k_rho_B + integral_val * Krook_G2_rho_B(j, plasma%spec(sigma)) * Krook_kappa_rho_B(j, plasma%spec(sigma))

                    call gauss_integrate_F3(int_F3, integral_val, gauss_conf)
                    k_rho_phi = k_rho_phi + integral_val * Krook_kappa_rho_phi(j, plasma%spec(sigma)) * Krook_G3_rho_phi(j, plasma%spec(sigma))
                    k_rho_B = k_rho_B + integral_val * Krook_G3_rho_B(j, plasma%spec(sigma)) * Krook_kappa_rho_B(j, plasma%spec(sigma))
                end if
                
            end do
        end do

        k_rho_phi = k_rho_phi / (8.0d0 * pi**3.0d0) !/sqrt(2.0d0) ! factor sqrt(1/2) is somehow missing. Including this factor nicely reproduces the debye case.
        k_rho_B = k_rho_B / (8.0d0 * pi**3.0d0)
            
    end subroutine


    subroutine FP_fill_kernels(K_rho_phi_llp, K_rho_B_llp, K_j_phi_llp, K_j_B_llp)

        use KIM_kinds_m, only: dp
        use electrostatic_integrals_gauss_mod, only: gauss_config_t, init_gauss_int
        use grid_m, only: Larmor_skip_factor, gauss_int_nodes_Ntheta, gauss_int_nodes_Nx, gauss_int_nodes_Nxp

        implicit none

        type(kernel_spl_t), intent(inout) :: K_rho_phi_llp
        type(kernel_spl_t), intent(inout) :: K_rho_B_llp
        type(kernel_spl_t), intent(inout) :: K_j_phi_llp
        type(kernel_spl_t), intent(inout) :: K_j_B_llp
        type(gauss_config_t) :: gauss_conf
        integer :: l, lp

        gauss_conf%Nx = gauss_int_nodes_Nx
        gauss_conf%Nxp = gauss_int_nodes_Nxp
        gauss_conf%Ntheta = gauss_int_nodes_Ntheta

        call init_gauss_int(gauss_conf)

        write(*,*) 'Filling Fokker-Planck collision kernels...'

        !$omp parallel do collapse(1) private(l,lp)
        do l = 1, K_rho_phi_llp%npts_l
            do lp = 1, l

                call FP_calc_kernels(l, lp, K_rho_phi_llp%Kllp(l, lp),&
                                            K_rho_B_llp%Kllp(l, lp), &
                                            K_j_phi_llp%Kllp(l, lp), &
                                            K_j_B_llp%Kllp(l, lp), &
                                            gauss_conf)

                if (isnan(real(K_rho_phi_llp%Kllp(l,lp)))) then
                    print *, "K_rho_phi_llp is NaN for l = ", l, " lp = ", lp
                    stop
                end if
                if (isnan(real(K_rho_B_llp%Kllp(l,lp)))) then
                    print *, "K_rho_B_llp is NaN for l = ", l, " lp = ", lp
                    stop
                end if

                if (isnan(real(K_j_phi_llp%Kllp(l,lp)))) then
                    print *, "K_j_phi_llp is NaN for l = ", l, " lp = ", lp
                    stop
                end if
                if (isnan(real(K_j_B_llp%Kllp(l,lp)))) then
                    print *, "K_j_B_llp is NaN for l = ", l, " lp = ", lp
                    stop
                end if

                K_rho_phi_llp%Kllp(lp, l) = K_rho_phi_llp%Kllp(l, lp)
                K_rho_B_llp%Kllp(lp, l) = K_rho_B_llp%Kllp(l, lp)
                K_j_phi_llp%Kllp(lp, l) = K_j_phi_llp%Kllp(l, lp)
                K_j_B_llp%Kllp(lp, l) = K_j_B_llp%Kllp(l, lp)
            end do
        end do
        !$omp end parallel do
        
        ! Print diagnostic information
        write(*,*) '======== Kernel Distance Diagnostics (Fokker-Planck) ========'
        write(*,'(A,F12.6)') ' Maximum |xl - xlp| distance: ', max_distance_xl_xlp
        write(*,'(A,I6,A,I6)') ' Occurred at l = ', max_dist_l, ', lp = ', max_dist_lp
        write(*,'(A,F12.6)') ' Minimum |xl - xlp| distance: ', min_distance_xl_xlp
        write(*,'(A,I6,A,I6)') ' Occurred at l = ', min_dist_l, ', lp = ', min_dist_lp
        write(*,'(A,I6)') ' Maximum index distance |l - lp|: ', max_index_distance
        write(*,'(A,I6,A,I6)') ' Occurred at l = ', max_idx_l, ', lp = ', max_idx_lp
        write(*,'(A,I6)') ' Minimum index distance |l - lp|: ', min_index_distance
        write(*,'(A,I6,A,I6)') ' Occurred at l = ', min_idx_l, ', lp = ', min_idx_lp
        write(*,*) '============================================================='

    end subroutine

    
    subroutine FP_calc_kernels(l, lp, k_rho_phi, k_rho_B, k_j_phi, k_j_B, gauss_conf)

        use KIM_kinds_m, only: dp
        use electrostatic_integrals_gauss_mod, only: gauss_integrate_F0, gauss_integrate_F1, gauss_integrate_F2, gauss_integrate_F3,&
            gauss_config_t
        use species_m, only: plasma
        use constants_m, only: pi
        use electrostatic_integrands_gauss_mod, only: gauss_int_F0_rho_phi_t, gauss_int_F1_rho_phi_t, gauss_int_F2_rho_phi_t, gauss_int_F3_rho_phi_t, &
            integration_point_t
        use FP_kernel_plasma_prefacs_m, only: FP_G1_rho_phi, FP_G1_rho_B, FP_G2_rho_B, FP_G3_rho_B, &
            FP_G2_rho_phi, FP_G3_rho_phi, FP_kappa_rho_phi, FP_kappa_rho_B, FP_G0_rho_phi, &
            FP_kappa_j_phi, FP_kappa_j_B, FP_G1_j_phi, FP_G2_j_phi, FP_G3_j_phi, &
            FP_G1_j_B, FP_G2_j_B, FP_G3_j_B
        use grid_m, only: Larmor_skip_factor
        use config_m, only: turn_off_ions
        
        implicit none

        integer, intent(in) :: l, lp
        complex(dp) :: k_rho_phi, k_rho_B, k_j_phi, k_j_B
        integer :: j, sigma
        type(gauss_config_t), intent(in) :: gauss_conf
        real(dp) :: integral_val
        real(dp) :: current_distance
        integer :: current_idx_distance

        type(integration_point_t) :: int_point
        type(gauss_int_F0_rho_phi_t) :: int_F0
        type(gauss_int_F1_rho_phi_t) :: int_F1
        type(gauss_int_F2_rho_phi_t) :: int_F2
        type(gauss_int_F3_rho_phi_t) :: int_F3
        
        k_rho_phi = 0.0d0
        k_rho_B = 0.0d0
        k_j_phi = 0.0d0
        k_j_B = 0.0d0

        call set_xl_at_edge(l, lp, int_point)

        do sigma = 0, plasma%n_species - 1
            if (turn_off_ions .and. sigma >= 1) cycle
            do j = 2, size(plasma%r_grid)-1
                int_point%j = j
                int_point%rhoT = 0.5d0 * (plasma%spec(sigma)%rho_L(j) + plasma%spec(sigma)%rho_L(j+1))

                if (l == lp) then
                    int_F0%int_point = int_point
                    call gauss_integrate_F0(int_F0, int_point%xlm1, int_point%xlp1, integral_val, gauss_conf)
                    k_rho_phi = k_rho_phi &
                        + integral_val * FP_G0_rho_phi(j, plasma%spec(sigma)) * FP_kappa_rho_phi(j, plasma%spec(sigma))
                end if

                ! Track maximum distances for diagnostics
                current_distance = abs(int_point%xl - int_point%xlp)

                ! skip the kernel calculation for distances that are not connected
                if (current_distance > Larmor_skip_factor * int_point%rhoT) cycle

                current_idx_distance = abs(l - lp)
                
                !$omp critical
                if (current_distance > max_distance_xl_xlp) then
                    max_distance_xl_xlp = current_distance
                    max_dist_l = l
                    max_dist_lp = lp
                end if
                if (current_distance < min_distance_xl_xlp .and. current_distance > 0.0d0) then
                    min_distance_xl_xlp = current_distance
                    min_dist_l = l
                    min_dist_lp = lp
                end if
                if (current_idx_distance > max_index_distance) then
                    max_index_distance = current_idx_distance
                    max_idx_l = l
                    max_idx_lp = lp
                end if
                if (current_idx_distance < min_index_distance .and. current_idx_distance > 0) then
                    min_index_distance = current_idx_distance
                    min_idx_l = l
                    min_idx_lp = lp
                end if
                !$omp end critical

                int_F1%int_point = int_point
                int_F2%int_point = int_point
                int_F3%int_point = int_point

                call gauss_integrate_F1(int_F1, integral_val, gauss_conf)

                k_rho_phi = k_rho_phi + integral_val * FP_G1_rho_phi(j, plasma%spec(sigma)) * FP_kappa_rho_phi(j, plasma%spec(sigma))
                k_rho_B   = k_rho_B   + integral_val * FP_G1_rho_B(j, plasma%spec(sigma))   * FP_kappa_rho_B(j, plasma%spec(sigma))

                k_j_phi   = k_j_phi   + integral_val * FP_G1_j_phi(j, plasma%spec(sigma))   * FP_kappa_j_phi(j, plasma%spec(sigma))
                k_j_B     = k_j_B     + integral_val * FP_G1_j_B(j, plasma%spec(sigma))     * FP_kappa_j_B(j, plasma%spec(sigma))

                call gauss_integrate_F2(int_F2, integral_val, gauss_conf)

                k_rho_phi = k_rho_phi + integral_val * FP_G2_rho_phi(j, plasma%spec(sigma)) * FP_kappa_rho_phi(j, plasma%spec(sigma))
                k_rho_B   = k_rho_B   + integral_val * FP_G2_rho_B(j, plasma%spec(sigma))   * FP_kappa_rho_B(j, plasma%spec(sigma))

                k_j_phi   = k_j_phi   + integral_val * FP_G2_j_phi(j, plasma%spec(sigma))   * FP_kappa_j_phi(j, plasma%spec(sigma))
                k_j_B     = k_j_B     + integral_val * FP_G2_j_B(j, plasma%spec(sigma))     * FP_kappa_j_B(j, plasma%spec(sigma))

                call gauss_integrate_F3(int_F3, integral_val, gauss_conf)

                k_rho_phi = k_rho_phi + integral_val * FP_G3_rho_phi(j, plasma%spec(sigma)) * FP_kappa_rho_phi(j, plasma%spec(sigma))
                k_rho_B   = k_rho_B   + integral_val * FP_G3_rho_B(j, plasma%spec(sigma))   * FP_kappa_rho_B(j, plasma%spec(sigma))

                k_j_phi   = k_j_phi   + integral_val * FP_G3_j_phi(j, plasma%spec(sigma))   * FP_kappa_j_phi(j, plasma%spec(sigma))
                k_j_B     = k_j_B     + integral_val * FP_G3_j_B(j, plasma%spec(sigma))     * FP_kappa_j_B(j, plasma%spec(sigma))

            end do
        end do

        k_rho_phi = k_rho_phi / (8.0d0 * pi**3.0d0)
        k_rho_B = k_rho_B / (8.0d0 * pi**3.0d0)

        k_j_phi = k_j_phi / (8.0d0 * pi**3.0d0)
        k_j_B = k_j_B / (8.0d0 * pi**3.0d0)

            
    end subroutine
    
    subroutine set_xl_at_edge(l, lp, int_point)
        
        use grid_m, only: xl_grid
        use KIM_kinds_m, only: dp
        use electrostatic_integrands_gauss_mod, only: integration_point_t

        implicit none

        integer, intent(in) :: l, lp
        type(integration_point_t), intent(inout) :: int_point

        int_point%xl = xl_grid%xb(l)
        int_point%xlp = xl_grid%xb(lp)

        ! Handle lower boundary with symmetric extrapolation
        if (l == 1) then
            int_point%xlm1 = 2.0d0*xl_grid%xb(1) - xl_grid%xb(2)  ! Extrapolate
        else
            int_point%xlm1 = xl_grid%xb(l-1)
        end if
        
        if (lp == 1) then
            int_point%xlpm1 = 2.0d0*xl_grid%xb(1) - xl_grid%xb(2)  ! Fixed: symmetric extrapolation
        else
            int_point%xlpm1 = xl_grid%xb(lp-1)
        end if
        
        ! Handle upper boundary with symmetric extrapolation
        if (l == xl_grid%npts_b) then
            int_point%xlp1 = 2.0d0*xl_grid%xb(l) - xl_grid%xb(l-1)  ! Fixed: extrapolation
        else
            int_point%xlp1 = xl_grid%xb(l+1)
        end if
        
        if (lp == xl_grid%npts_b) then
            int_point%xlpp1 = 2.0d0*xl_grid%xb(lp) - xl_grid%xb(lp-1)  ! Fixed: extrapolation
        else
            int_point%xlpp1 = xl_grid%xb(lp+1)
        end if

    end subroutine

    subroutine fill_kernels_krook_fp(kernel_krook_rho_phi, kernel_krook_rho_B, &
                                      kernel_fp_rho_phi, kernel_fp_rho_B)
        !> Unified subroutine to fill both Krook and Fokker-Planck kernels
        !> Exploits shared Gaussian integration for efficiency
        
        use KIM_kinds_m, only: dp
        use electrostatic_integrals_gauss_mod, only: gauss_config_t, init_gauss_int, &
            gauss_integrate_F0, gauss_integrate_F1, gauss_integrate_F2, gauss_integrate_F3
        use grid_m, only: Larmor_skip_factor, gauss_int_nodes_Ntheta, gauss_int_nodes_Nx, gauss_int_nodes_Nxp
        use species_m, only: plasma
        use constants_m, only: pi
        use electrostatic_integrands_gauss_mod, only: gauss_int_F0_rho_phi_t, gauss_int_F1_rho_phi_t, &
            gauss_int_F2_rho_phi_t, gauss_int_F3_rho_phi_t, integration_point_t
        use Krook_kernel_plasma_prefacs_m, only: Krook_G0_rho_phi, Krook_G1_rho_phi, Krook_G2_rho_phi, Krook_G3_rho_phi, &
            Krook_G1_rho_B, Krook_G2_rho_B, Krook_G3_rho_B, Krook_kappa_rho_phi, Krook_kappa_rho_B
        use FP_kernel_plasma_prefacs_m, only: FP_G0_rho_phi, FP_G1_rho_phi, FP_G2_rho_phi, &
            FP_G3_rho_phi, FP_G1_rho_B, FP_G2_rho_B, FP_G3_rho_B, &
            FP_kappa_rho_phi, FP_kappa_rho_B
        use config_m, only: artificial_debye_case, turn_off_ions
        
        implicit none
        
        type(kernel_spl_t), intent(inout) :: kernel_krook_rho_phi, kernel_krook_rho_B
        type(kernel_spl_t), intent(inout) :: kernel_fp_rho_phi, kernel_fp_rho_B
        
        type(gauss_config_t) :: gauss_conf
        integer :: l, lp, j, sigma
        complex(dp) :: krook_phi_llp, krook_B_llp, fp_phi_llp, fp_B_llp
        real(dp) :: integral_F0, integral_F1, integral_F2, integral_F3
        real(dp) :: current_distance
        integer :: current_idx_distance
        
        type(integration_point_t) :: int_point
        type(gauss_int_F0_rho_phi_t) :: int_F0
        type(gauss_int_F1_rho_phi_t) :: int_F1
        type(gauss_int_F2_rho_phi_t) :: int_F2
        type(gauss_int_F3_rho_phi_t) :: int_F3
        
        ! Initialize Gaussian integration configuration
        gauss_conf%Nx = gauss_int_nodes_Nx
        gauss_conf%Nxp = gauss_int_nodes_Nxp
        gauss_conf%Ntheta = gauss_int_nodes_Ntheta
        call init_gauss_int(gauss_conf)

        write(*,*) 'Filling both Krook and Fokker-Planck kernels simultaneously...'
        
        !$omp parallel do collapse(1) private(l, lp, krook_phi_llp, krook_B_llp, &
        !$omp& fp_phi_llp, fp_B_llp, j, sigma, int_point, int_F0, int_F1, int_F2, int_F3, &
        !$omp& integral_F0, integral_F1, integral_F2, integral_F3, current_distance, current_idx_distance)
        do l = 1, kernel_krook_rho_phi%npts_l
            do lp = 1, l
                
                ! Initialize kernel values
                krook_phi_llp = 0.0d0
                krook_B_llp = 0.0d0
                fp_phi_llp = 0.0d0
                fp_B_llp = 0.0d0
                
                ! Set xl grid points at edges
                call set_xl_at_edge(l, lp, int_point)
                
                ! Loop over species and radial grid
                do sigma = 0, plasma%n_species - 1
                    if (turn_off_ions .and. sigma >= 1) cycle
                    do j = 2, size(plasma%r_grid)-1
                        int_point%j = j
                        int_point%rhoT = 0.5d0 * (plasma%spec(sigma)%rho_L(j) + plasma%spec(sigma)%rho_L(j+1))
                        
                        ! F0 integral (only for diagonal elements)
                        if (l == lp) then
                            int_F0%int_point = int_point
                            call gauss_integrate_F0(int_F0, int_point%xlm1, int_point%xlp1, integral_F0, gauss_conf)
                            
                            ! Krook contribution
                            krook_phi_llp = krook_phi_llp + integral_F0 * Krook_G0_rho_phi(j, plasma%spec(sigma)) * &
                                           Krook_kappa_rho_phi(j, plasma%spec(sigma))
                            
                            ! Fokker-Planck contribution
                            fp_phi_llp = fp_phi_llp + integral_F0 * FP_G0_rho_phi(j, plasma%spec(sigma)) * &
                                        FP_kappa_rho_phi(j, plasma%spec(sigma))
                        end if

                        ! Track maximum distances for diagnostics
                        current_distance = abs(int_point%xl - int_point%xlp)

                        if (current_distance > Larmor_skip_factor * int_point%rhoT) cycle

                        current_idx_distance = abs(l - lp)
                        
                        !$omp critical
                        if (current_distance > max_distance_xl_xlp) then
                            max_distance_xl_xlp = current_distance
                            max_dist_l = l
                            max_dist_lp = lp
                        end if
                        if (current_idx_distance > max_index_distance) then
                            max_index_distance = current_idx_distance
                            max_idx_l = l
                            max_idx_lp = lp
                        end if
                        !$omp end critical
                        
                        if (.not. artificial_debye_case) then
                            ! Set integration points for F1, F2, F3
                            int_F1%int_point = int_point
                            int_F2%int_point = int_point
                            int_F3%int_point = int_point
                            
                            ! Perform Gaussian integrations (shared between Krook and FP)
                            call gauss_integrate_F1(int_F1, integral_F1, gauss_conf)
                            call gauss_integrate_F2(int_F2, integral_F2, gauss_conf)
                            call gauss_integrate_F3(int_F3, integral_F3, gauss_conf)
                            
                            ! Krook contributions
                            krook_phi_llp = krook_phi_llp  &
                                + integral_F1 * Krook_G1_rho_phi(j, plasma%spec(sigma)) * Krook_kappa_rho_phi(j, plasma%spec(sigma))&
                                + integral_F2 * Krook_G2_rho_phi(j, plasma%spec(sigma)) * Krook_kappa_rho_phi(j, plasma%spec(sigma)) &
                                + integral_F3 * Krook_G3_rho_phi(j, plasma%spec(sigma)) * Krook_kappa_rho_phi(j, plasma%spec(sigma))
                            
                            krook_B_llp = krook_B_llp &
                                + integral_F1 * Krook_G1_rho_B(j, plasma%spec(sigma)) * Krook_kappa_rho_B(j, plasma%spec(sigma)) &
                                + integral_F2 * Krook_G2_rho_B(j, plasma%spec(sigma)) * Krook_kappa_rho_B(j, plasma%spec(sigma)) &
                                + integral_F3 * Krook_G3_rho_B(j, plasma%spec(sigma)) * Krook_kappa_rho_B(j, plasma%spec(sigma))
                            
                            ! Fokker-Planck contributions
                            fp_phi_llp = fp_phi_llp &
                                + integral_F1 * FP_G1_rho_phi(j, plasma%spec(sigma)) * FP_kappa_rho_phi(j, plasma%spec(sigma)) &
                                + integral_F2 * FP_G2_rho_phi(j, plasma%spec(sigma)) * FP_kappa_rho_phi(j, plasma%spec(sigma)) &
                                + integral_F3 * FP_G3_rho_phi(j, plasma%spec(sigma)) * FP_kappa_rho_phi(j, plasma%spec(sigma))
                            
                            fp_B_llp = fp_B_llp &
                                + integral_F1 * FP_G1_rho_B(j, plasma%spec(sigma)) * FP_kappa_rho_B(j, plasma%spec(sigma)) &
                                + integral_F2 * FP_G2_rho_B(j, plasma%spec(sigma)) * FP_kappa_rho_B(j, plasma%spec(sigma)) &
                                + integral_F3 * FP_G3_rho_B(j, plasma%spec(sigma)) * FP_kappa_rho_B(j, plasma%spec(sigma))
                        end if
                    end do
                end do
                
                ! Apply normalization factor
                krook_phi_llp = krook_phi_llp / (8.0d0 * pi**3.0d0)
                krook_B_llp = krook_B_llp / (8.0d0 * pi**3.0d0)
                fp_phi_llp = fp_phi_llp / (8.0d0 * pi**3.0d0)
                fp_B_llp = fp_B_llp / (8.0d0 * pi**3.0d0)
                
                ! Store results
                kernel_krook_rho_phi%Kllp(l, lp) = krook_phi_llp
                kernel_krook_rho_B%Kllp(l, lp) = krook_B_llp
                kernel_fp_rho_phi%Kllp(l, lp) = fp_phi_llp
                kernel_fp_rho_B%Kllp(l, lp) = fp_B_llp
                
                ! Check for NaN values
                if (isnan(real(krook_phi_llp)) .or. isnan(real(krook_B_llp)) .or. &
                    isnan(real(fp_phi_llp)) .or. isnan(real(fp_B_llp))) then
                    print *, "NaN detected in kernel calculation at l =", l, ", lp =", lp
                    print *, "Krook phi:", krook_phi_llp, "Krook B:", krook_B_llp
                    print *, "FP phi:", fp_phi_llp, "FP B:", fp_B_llp
                    stop
                end if

                ! exploit symmetry
                kernel_krook_rho_phi%Kllp(lp, l) = kernel_krook_rho_phi%Kllp(l, lp)
                kernel_krook_rho_B%Kllp(lp, l) = kernel_krook_rho_B%Kllp(l, lp)
                kernel_fp_rho_phi%Kllp(lp, l) = kernel_fp_rho_phi%Kllp(l, lp)
                kernel_fp_rho_B%Kllp(lp, l) = kernel_fp_rho_B%Kllp(l, lp)
            end do
        end do
        !$omp end parallel do

        write(*,*) ! New line after progress bar
        
        ! Print diagnostic information
        write(*,*) '======== Kernel Distance Diagnostics (Combined Krook+FP) ========'
        write(*,'(A,F12.6)') ' Maximum |xl - xlp| distance: ', max_distance_xl_xlp
        write(*,'(A,I6,A,I6)') ' Occurred at l = ', max_dist_l, ', lp = ', max_dist_lp
        write(*,'(A,F12.6)') ' Minimum |xl - xlp| distance: ', min_distance_xl_xlp
        write(*,'(A,I6,A,I6)') ' Occurred at l = ', min_dist_l, ', lp = ', min_dist_lp
        write(*,'(A,I6)') ' Maximum index distance |l - lp|: ', max_index_distance
        write(*,'(A,I6,A,I6)') ' Occurred at l = ', max_idx_l, ', lp = ', max_idx_lp
        write(*,'(A,I6)') ' Minimum index distance |l - lp|: ', min_index_distance
        write(*,'(A,I6,A,I6)') ' Occurred at l = ', min_idx_l, ', lp = ', min_idx_lp
        write(*,*) '=================================================================='
        
    end subroutine fill_kernels_krook_fp

end module

