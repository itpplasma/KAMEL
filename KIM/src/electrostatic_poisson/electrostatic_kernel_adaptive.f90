module electrostatic_kernel_adaptive_mod

    use KIM_kinds_m, only: dp
    use electrostatic_kernel_m, only: kernel_spl_t, max_distance_xl_xlp, min_distance_xl_xlp, &
        max_index_distance, min_index_distance, &
        max_dist_l, max_dist_lp, min_dist_l, min_dist_lp, &
        max_idx_l, max_idx_lp, min_idx_l, min_idx_lp

    implicit none

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

    subroutine FP_fill_kernels_adaptive(K_rho_phi_llp, K_rho_B_llp, K_j_phi_llp, K_j_B_llp)

        use KIM_kinds_m, only: dp
        use electrostatic_integrals_rkf45_mod, only: rkf45_config_t, init_rkf45_int
        use grid_m, only: Larmor_skip_factor, gauss_int_nodes_Ntheta, gauss_int_nodes_Nx, gauss_int_nodes_Nxp
        use loading_bar_m, only: updateLoadingBarWithETA

        implicit none

        type(kernel_spl_t), intent(inout) :: K_rho_phi_llp
        type(kernel_spl_t), intent(inout) :: K_rho_B_llp
        type(kernel_spl_t), intent(inout) :: K_j_phi_llp
        type(kernel_spl_t), intent(inout) :: K_j_B_llp
        type(rkf45_config_t) :: rkf45_conf
        integer :: l, lp
        integer :: total_iterations, current_iteration
        integer(kind=8) :: start_count, count_rate, count_max

        rkf45_conf%Nx = gauss_int_nodes_Nx
        rkf45_conf%Nxp = gauss_int_nodes_Nxp

        call init_rkf45_int(rkf45_conf)

        write(*,*) 'Filling Fokker-Planck collision kernels...'

        ! Calculate total number of iterations for the loading bar
        total_iterations = K_rho_phi_llp%npts_l * (K_rho_phi_llp%npts_l + 1) / 2
        current_iteration = 0

        ! Record start wall time for ETA calculation
        call system_clock(start_count, count_rate, count_max)

        do l = 1, K_rho_phi_llp%npts_l
            do lp = 1, l

                call FP_calc_kernels_adaptive(l, lp, K_rho_phi_llp%Kllp(l, lp),&
                                            K_rho_B_llp%Kllp(l, lp), &
                                            K_j_phi_llp%Kllp(l, lp), &
                                            K_j_B_llp%Kllp(l, lp), &
                                            rkf45_conf)

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

                current_iteration = current_iteration + 1
                call updateLoadingBarWithETA(current_iteration, total_iterations, start_count, count_rate)
            end do
        end do
        
        write(*,*)
        
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

    
    subroutine FP_calc_kernels_adaptive(l, lp, k_rho_phi, k_rho_B, k_j_phi, k_j_B, rkf45_conf)

        use KIM_kinds_m, only: dp
        use electrostatic_integrals_rkf45_mod, only: rkf45_integrate_F0, rkf45_integrate_F1, &
            rkf45_integrate_F2, rkf45_integrate_F3, rkf45_config_t
        use electrostatic_integrands_rkf45_mod, only: rkf45_integrand_context_t
        use species_m, only: plasma
        use constants_m, only: pi
        use FP_kernel_plasma_prefacs_m, only: FP_G1_rho_phi, FP_G1_rho_B, FP_G2_rho_B, FP_G3_rho_B, &
            FP_G2_rho_phi, FP_G3_rho_phi, FP_kappa_rho_phi, FP_kappa_rho_B, FP_G0_rho_phi, &
            FP_kappa_j_phi, FP_kappa_j_B, FP_G1_j_phi, FP_G2_j_phi, FP_G3_j_phi, &
            FP_G1_j_B, FP_G2_j_B, FP_G3_j_B
        use grid_m, only: Larmor_skip_factor, kernel_taper_skip_threshold
        use config_m, only: turn_off_ions
        
        implicit none

        integer, intent(in) :: l, lp
        complex(dp) :: k_rho_phi, k_rho_B, k_j_phi, k_j_B
        complex(dp) :: c_rho_phi, c_rho_B, c_j_phi, c_j_B  ! Kahan compensation terms
        integer :: j, sigma
        type(rkf45_config_t), intent(in) :: rkf45_conf
        real(dp) :: integral_val
        real(dp) :: current_distance
        integer :: current_idx_distance

        type(rkf45_integrand_context_t) :: context
        
        k_rho_phi = (0.0d0, 0.0d0)
        k_rho_B = (0.0d0, 0.0d0)
        k_j_phi = (0.0d0, 0.0d0)
        k_j_B = (0.0d0, 0.0d0)
        c_rho_phi = (0.0d0, 0.0d0)
        c_rho_B   = (0.0d0, 0.0d0)
        c_j_phi   = (0.0d0, 0.0d0)
        c_j_B     = (0.0d0, 0.0d0)

        call set_xl_at_edge(l, lp, context)

        do sigma = 0, plasma%n_species - 1
            if (turn_off_ions .and. sigma >= 1) cycle
            do j = 2, size(plasma%r_grid)-1
                context%j = j
                ! Use geometric mean to reduce bias when rho_L varies strongly across cell
                context%rhoT = sqrt(max(plasma%spec(sigma)%rho_L(j),0.0d0) * &
                                     max(plasma%spec(sigma)%rho_L(j+1),0.0d0))
                context%ks = 0.5d0 * (plasma%ks(context%j) + plasma%ks(context%j+1))

                if (l == lp) then
                    complex(dp) :: add
                    call rkf45_integrate_F0(integral_val, rkf45_conf, context)
                    add = integral_val * FP_G0_rho_phi(j, plasma%spec(sigma)) * FP_kappa_rho_phi(j, plasma%spec(sigma))
                    ! Kahan summation for k_rho_phi
                    block
                        complex(dp) :: y, t
                        y = add - c_rho_phi
                        t = k_rho_phi + y
                        c_rho_phi = (t - k_rho_phi) - y
                        k_rho_phi = t
                    end block
                end if

                !cycle

                ! Track maximum distances for diagnostics
                current_distance = abs(context%xl - context%xlp)

                ! Smoothly taper contributions for large separations to avoid discontinuities
                ! Weight per cell because rhoT varies with j
                ! w(d) = exp( - (d / (alpha * rhoT + eps))^p ) with p=2
                ! If the weight is below a small threshold, skip this (l,lp,j) contribution entirely.

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

                block
                    real(dp) :: eps_r, alpha, pexp, weight
                    eps_r = 1.0d-12
                    alpha = Larmor_skip_factor
                    pexp = 2.0d0
                    weight = exp( - ( current_distance / (alpha * max(context%rhoT, eps_r)) )**pexp )
                    if (weight < kernel_taper_skip_threshold) cycle

                    call rkf45_integrate_F1(integral_val, rkf45_conf, context)

                    block
                        complex(dp) :: add, y, t
                        add = weight * integral_val * FP_G1_rho_phi(j, plasma%spec(sigma)) * FP_kappa_rho_phi(j, plasma%spec(sigma))
                        y = add - c_rho_phi
                        t = k_rho_phi + y
                        c_rho_phi = (t - k_rho_phi) - y
                        k_rho_phi = t
                    end block

                    block
                        complex(dp) :: add, y, t
                        add = weight * integral_val * FP_G1_rho_B(j, plasma%spec(sigma)) * FP_kappa_rho_B(j, plasma%spec(sigma))
                        y = add - c_rho_B
                        t = k_rho_B + y
                        c_rho_B = (t - k_rho_B) - y
                        k_rho_B = t
                    end block

                    block
                        complex(dp) :: add, y, t
                        add = weight * integral_val * FP_G1_j_phi(j, plasma%spec(sigma)) * FP_kappa_j_phi(j, plasma%spec(sigma))
                        y = add - c_j_phi
                        t = k_j_phi + y
                        c_j_phi = (t - k_j_phi) - y
                        k_j_phi = t
                    end block

                    block
                        complex(dp) :: add, y, t
                        add = weight * integral_val * FP_G1_j_B(j, plasma%spec(sigma)) * FP_kappa_j_B(j, plasma%spec(sigma))
                        y = add - c_j_B
                        t = k_j_B + y
                        c_j_B = (t - k_j_B) - y
                        k_j_B = t
                    end block
                end block

                block
                    real(dp) :: eps_r, alpha, pexp, weight
                    eps_r = 1.0d-12
                    alpha = Larmor_skip_factor
                    pexp = 2.0d0
                    weight = exp( - ( current_distance / (alpha * max(context%rhoT, eps_r)) )**pexp )
                    if (weight < kernel_taper_skip_threshold) cycle

                    call rkf45_integrate_F2(integral_val, rkf45_conf, context)

                    block
                        complex(dp) :: add, y, t
                        add = weight * integral_val * FP_G2_rho_phi(j, plasma%spec(sigma)) * FP_kappa_rho_phi(j, plasma%spec(sigma))
                        y = add - c_rho_phi
                        t = k_rho_phi + y
                        c_rho_phi = (t - k_rho_phi) - y
                        k_rho_phi = t
                    end block

                    block
                        complex(dp) :: add, y, t
                        add = weight * integral_val * FP_G2_rho_B(j, plasma%spec(sigma)) * FP_kappa_rho_B(j, plasma%spec(sigma))
                        y = add - c_rho_B
                        t = k_rho_B + y
                        c_rho_B = (t - k_rho_B) - y
                        k_rho_B = t
                    end block

                    block
                        complex(dp) :: add, y, t
                        add = weight * integral_val * FP_G2_j_phi(j, plasma%spec(sigma)) * FP_kappa_j_phi(j, plasma%spec(sigma))
                        y = add - c_j_phi
                        t = k_j_phi + y
                        c_j_phi = (t - k_j_phi) - y
                        k_j_phi = t
                    end block

                    block
                        complex(dp) :: add, y, t
                        add = weight * integral_val * FP_G2_j_B(j, plasma%spec(sigma)) * FP_kappa_j_B(j, plasma%spec(sigma))
                        y = add - c_j_B
                        t = k_j_B + y
                        c_j_B = (t - k_j_B) - y
                        k_j_B = t
                    end block
                end block

                block
                    real(dp) :: eps_r, alpha, pexp, weight
                    eps_r = 1.0d-12
                    alpha = Larmor_skip_factor
                    pexp = 2.0d0
                    weight = exp( - ( current_distance / (alpha * max(context%rhoT, eps_r)) )**pexp )
                    if (weight < kernel_taper_skip_threshold) cycle

                    call rkf45_integrate_F3(integral_val, rkf45_conf, context)

                    block
                        complex(dp) :: add, y, t
                        add = weight * integral_val * FP_G3_rho_phi(j, plasma%spec(sigma)) * FP_kappa_rho_phi(j, plasma%spec(sigma))
                        y = add - c_rho_phi
                        t = k_rho_phi + y
                        c_rho_phi = (t - k_rho_phi) - y
                        k_rho_phi = t
                    end block

                    block
                        complex(dp) :: add, y, t
                        add = weight * integral_val * FP_G3_rho_B(j, plasma%spec(sigma)) * FP_kappa_rho_B(j, plasma%spec(sigma))
                        y = add - c_rho_B
                        t = k_rho_B + y
                        c_rho_B = (t - k_rho_B) - y
                        k_rho_B = t
                    end block

                    block
                        complex(dp) :: add, y, t
                        add = weight * integral_val * FP_G3_j_phi(j, plasma%spec(sigma)) * FP_kappa_j_phi(j, plasma%spec(sigma))
                        y = add - c_j_phi
                        t = k_j_phi + y
                        c_j_phi = (t - k_j_phi) - y
                        k_j_phi = t
                    end block

                    block
                        complex(dp) :: add, y, t
                        add = weight * integral_val * FP_G3_j_B(j, plasma%spec(sigma)) * FP_kappa_j_B(j, plasma%spec(sigma))
                        y = add - c_j_B
                        t = k_j_B + y
                        c_j_B = (t - k_j_B) - y
                        k_j_B = t
                    end block
                end block

            end do
        end do

        k_rho_phi = k_rho_phi / (8.0d0 * pi**3.0d0) !* 3.0d0 / 4.0d0 ! maybe pi instead of 3?
        k_rho_B = k_rho_B / (8.0d0 * pi**3.0d0)

        k_j_phi = k_j_phi / (8.0d0 * pi**3.0d0)
        k_j_B = k_j_B / (8.0d0 * pi**3.0d0)

            
    end subroutine
    
    subroutine set_xl_at_edge(l, lp, context)
        
        use grid_m, only: xl_grid
        use KIM_kinds_m, only: dp
        use electrostatic_integrands_rkf45_mod, only: rkf45_integrand_context_t

        implicit none

        integer, intent(in) :: l, lp
        type(rkf45_integrand_context_t), intent(inout) :: context

        context%xl = xl_grid%xb(l)
        context%xlp = xl_grid%xb(lp)

        ! Handle lower boundary with symmetric extrapolation
        if (l == 1) then
            context%xlm1 = 2.0d0*xl_grid%xb(1) - xl_grid%xb(2)
        else
            context%xlm1 = xl_grid%xb(l-1)
        end if
        
        if (lp == 1) then
            context%xlpm1 = 2.0d0*xl_grid%xb(1) - xl_grid%xb(2)  ! Fixed: symmetric extrapolation
        else
            context%xlpm1 = xl_grid%xb(lp-1)
        end if
        
        ! Handle upper boundary with symmetric extrapolation
        if (l == xl_grid%npts_b) then
            context%xlp1 = 2.0d0*xl_grid%xb(l) - xl_grid%xb(l-1)  ! Fixed: extrapolation
        else
            context%xlp1 = xl_grid%xb(l+1)
        end if
        
        if (lp == xl_grid%npts_b) then
            context%xlpp1 = 2.0d0*xl_grid%xb(lp) - xl_grid%xb(lp-1)  ! Fixed: extrapolation
        else
            context%xlpp1 = xl_grid%xb(lp+1)
        end if

    end subroutine

end module
