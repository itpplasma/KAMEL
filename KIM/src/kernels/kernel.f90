module kernel_m

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

    ! Precomputed per-(species, rg-cell) prefactors on cell centers
    complex(dp), allocatable, save :: pref_rho_phi_g1(:, :, :), &
                                    pref_rho_phi_g2(:, :, :), &
                                    pref_rho_phi_g3(:, :, :), &
                                    pref_rho_phi_g0(:, :)
    complex(dp), allocatable, save :: pref_rho_B_g1(:, :, :), &
                                    pref_rho_B_g2(:, :, :), &
                                    pref_rho_B_g3(:, :, :)
    complex(dp), allocatable, save :: pref_j_phi_g1(:,:, :), & 
                                    pref_j_phi_g2(:, :, :), &
                                    pref_j_phi_g3(:, :, :)
    complex(dp), allocatable, save :: pref_j_B_g1(:, :, :), &
                                    pref_j_B_g2(:, :, :), &
                                    pref_j_B_g3(:, :, :)
    logical, save :: pref_ready = .false.

    type :: kernel_spl_t
        integer :: npts_l, npts_lp
        complex(dp), allocatable :: Kllp(:,:)
        complex(dp), allocatable :: Kllp_e(:,:)
        complex(dp), allocatable :: Kllp_i(:,:,:)
        contains
            procedure :: init_kernel
    end type kernel_spl_t

    contains

    ! Thread-safe wrapper to update loading bar from within OpenMP
    subroutine update_bar(cur, tot, sc, cr)
        use loading_bar_m, only: updateLoadingBarWithETA
        implicit none
        integer, intent(in) :: cur, tot
        integer(kind=8), intent(in) :: sc, cr
        call updateLoadingBarWithETA(cur, tot, sc, cr)
    end subroutine update_bar

    subroutine init_kernel(this, npts_l, npts_lp)

        use species_m, only: plasma

        implicit none

        class(kernel_spl_t), intent(inout) :: this
        integer, intent(in) :: npts_l, npts_lp

        this%npts_l = npts_l
        this%npts_lp = npts_lp
        allocate(this%Kllp(npts_l, npts_lp))
        allocate(this%Kllp_e(npts_l, npts_lp))
        allocate(this%Kllp_i(npts_l, npts_lp, 1:plasma%n_species-1))
        this%Kllp = (0.0d0, 0.0d0)
        this%Kllp_e = (0.0d0, 0.0d0)
        this%Kllp_i = (0.0d0, 0.0d0)

    end subroutine init_kernel

    subroutine compute_cc_prefactors

        use species_m, only: plasma
        use grid_m, only: rg_grid
        use FP_kernel_plasma_prefacs_m, only: FP_G1_rho_phi, FP_G2_rho_phi, FP_G3_rho_phi, &
            FP_G1_rho_B, FP_G2_rho_B, FP_G3_rho_B, FP_G1_j_phi, FP_G2_j_phi, FP_G3_j_phi, &
            FP_G1_j_B, FP_G2_j_B, FP_G3_j_B, FP_G0_rho_phi
        use config_m, only: fdiagnostics
        use setup_m, only: mphi_max

        implicit none
        integer :: ns, j, sigma, mphi

        ns = plasma%n_species
        if (.not. allocated(pref_rho_phi_g1)) then
            allocate(pref_rho_phi_g0(ns, rg_grid%npts_c))
            allocate(pref_rho_phi_g1(ns, rg_grid%npts_c, -mphi_max:mphi_max))
            allocate(pref_rho_phi_g2(ns, rg_grid%npts_c, -mphi_max:mphi_max))
            allocate(pref_rho_phi_g3(ns, rg_grid%npts_c, -mphi_max:mphi_max))
            allocate(pref_rho_B_g1(ns,  rg_grid%npts_c, -mphi_max:mphi_max))
            allocate(pref_rho_B_g2(ns,  rg_grid%npts_c, -mphi_max:mphi_max))
            allocate(pref_rho_B_g3(ns,  rg_grid%npts_c, -mphi_max:mphi_max))
            allocate(pref_j_phi_g1(ns,  rg_grid%npts_c, -mphi_max:mphi_max))
            allocate(pref_j_phi_g2(ns,  rg_grid%npts_c, -mphi_max:mphi_max))
            allocate(pref_j_phi_g3(ns,  rg_grid%npts_c, -mphi_max:mphi_max))
            allocate(pref_j_B_g1(ns,    rg_grid%npts_c, -mphi_max:mphi_max))
            allocate(pref_j_B_g2(ns,    rg_grid%npts_c, -mphi_max:mphi_max))
            allocate(pref_j_B_g3(ns,    rg_grid%npts_c, -mphi_max:mphi_max))
        end if

        do sigma = 0, ns - 1
            do j = 1, rg_grid%npts_c
                pref_rho_phi_g0(sigma+1, j) = FP_G0_rho_phi(j, plasma%spec(sigma))

                do mphi = -mphi_max, mphi_max
                    pref_rho_phi_g1(sigma+1, j, mphi) = FP_G1_rho_phi(j, plasma%spec(sigma), mphi)
                    pref_rho_phi_g2(sigma+1, j, mphi) = FP_G2_rho_phi(j, plasma%spec(sigma), mphi)
                    pref_rho_phi_g3(sigma+1, j, mphi) = FP_G3_rho_phi(j, plasma%spec(sigma), mphi)

                    pref_rho_B_g1(sigma+1, j, mphi) = FP_G1_rho_B(j, plasma%spec(sigma), mphi)
                    pref_rho_B_g2(sigma+1, j, mphi) = FP_G2_rho_B(j, plasma%spec(sigma), mphi)
                    pref_rho_B_g3(sigma+1, j, mphi) = FP_G3_rho_B(j, plasma%spec(sigma), mphi)

                    pref_j_phi_g1(sigma+1, j, mphi) = FP_G1_j_phi(j, plasma%spec(sigma), mphi)
                    pref_j_phi_g2(sigma+1, j, mphi) = FP_G2_j_phi(j, plasma%spec(sigma), mphi)
                    pref_j_phi_g3(sigma+1, j, mphi) = FP_G3_j_phi(j, plasma%spec(sigma), mphi)

                    pref_j_B_g1(sigma+1,j, mphi) = FP_G1_j_B(j, plasma%spec(sigma), mphi)
                    pref_j_B_g2(sigma+1,j, mphi) = FP_G2_j_B(j, plasma%spec(sigma), mphi)
                    pref_j_B_g3(sigma+1,j, mphi) = FP_G3_j_B(j, plasma%spec(sigma), mphi)
                end do
            end do
        end do

        pref_ready = .true.

        if (fdiagnostics == 3) then
            call write_cc_prefactors
        end if

    end subroutine compute_cc_prefactors

    subroutine write_cc_prefactors

        use IO_collection_m, only: write_complex_profile, itoa
        use config_m, only: output_path
        use species_m, only: plasma
        use grid_m, only: rg_grid
        use KIM_kinds_m, only: dp
        use setup_m, only: mphi_max

        implicit none

        logical :: ex
        integer :: sp, mphi
        character(len=100) :: dirname
        complex(dp), allocatable :: temp_array(:)
        real(dp), allocatable :: r_temp(:)

        inquire(file=trim(output_path)//'diagnostics', exist=ex)
        if (.not. ex) then
            call system('mkdir -p '//trim(output_path)//'diagnostics')
        end if

        if (.not. allocated(temp_array)) then
            allocate(temp_array(rg_grid%npts_c))
            allocate(r_temp(rg_grid%npts_c))
        end if

        r_temp = rg_grid%xc

        do sp = 0, plasma%n_species-1
            write(dirname, '(A,A)') 'diagnostics/', plasma%spec(sp)%name
            inquire(file=trim(output_path)//trim(dirname), exist=ex)
            if (.not. ex) then
                call system('mkdir -p '//trim(output_path)//trim(dirname))
            end if

            temp_array = pref_rho_phi_g0(sp+1, :)
            call write_complex_profile(r_temp, temp_array, rg_grid%npts_c, trim(dirname)//'/pref_rho_phi_g0')

            do mphi = -mphi_max, mphi_max

                temp_array = pref_rho_phi_g1(sp+1, :, mphi)
                call write_complex_profile(r_temp, temp_array, rg_grid%npts_c, trim(dirname)//'/pref_rho_phi_g1_mphi_'//itoa(mphi))

                temp_array = pref_rho_phi_g2(sp+1, :, mphi)
                call write_complex_profile(r_temp, temp_array, rg_grid%npts_c, trim(dirname)//'/pref_rho_phi_g2_mphi_'//itoa(mphi))

                temp_array = pref_rho_phi_g3(sp+1, :, mphi)
                call write_complex_profile(r_temp, temp_array, rg_grid%npts_c, trim(dirname)//'/pref_rho_phi_g3_mphi_'//itoa(mphi))

                temp_array = pref_rho_B_g1(sp+1, :, mphi)
                call write_complex_profile(r_temp, temp_array, rg_grid%npts_c, trim(dirname)//'/pref_rho_B_g1_mphi_'//itoa(mphi))

                temp_array = pref_rho_B_g2(sp+1, :, mphi)
                call write_complex_profile(r_temp, temp_array, rg_grid%npts_c, trim(dirname)//'/pref_rho_B_g2_mphi_'//itoa(mphi))

                temp_array = pref_rho_B_g3(sp+1, :, mphi)
                call write_complex_profile(r_temp, temp_array, rg_grid%npts_c, trim(dirname)//'/pref_rho_B_g3_mphi_'//itoa(mphi))
            end do

        end do

    end subroutine

    ! TODO: Update Krook routines to the same procedure as FP routines
    subroutine Krook_fill_kernel_phi(K_rho_phi_llp, K_rho_B_llp)

        use KIM_kinds_m, only: dp
        use integrals_gauss_m, only: gauss_config_t, init_gauss_int
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
        use integrals_gauss_m, only: gauss_integrate_F0, gauss_integrate_F1, gauss_integrate_F2, gauss_integrate_F3,&
            gauss_config_t
        use species_m, only: plasma
        use constants_m, only: pi
        use integrands_gauss_m, only: gauss_int_F0_rho_phi_t, gauss_int_F1_rho_phi_t, gauss_int_F2_rho_phi_t, gauss_int_F3_rho_phi_t, &
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
                
                if (artificial_debye_case /=1) then
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

        k_rho_phi = k_rho_phi / (8.0d0 * pi**3.0d0) 
        k_rho_B = k_rho_B / (8.0d0 * pi**3.0d0)
            
    end subroutine


    subroutine FP_fill_kernels(K_rho_phi_llp, K_rho_B_llp, K_j_phi_llp, K_j_B_llp)

        use KIM_kinds_m, only: dp
        use integrals_gauss_m, only: gauss_config_t, init_gauss_int
        use grid_m, only: Larmor_skip_factor, gauss_int_nodes_Ntheta, gauss_int_nodes_Nx, gauss_int_nodes_Nxp, &
                        kernel_taper_skip_threshold, rg_grid, xl_grid
        use species_m, only: plasma
        use config_m, only: output_path, artificial_debye_case, fstatus, turn_off_ions, &
                            turn_off_electrons

        implicit none

        type(kernel_spl_t), intent(inout) :: K_rho_phi_llp
        type(kernel_spl_t), intent(inout) :: K_rho_B_llp
        type(kernel_spl_t), intent(inout) :: K_j_phi_llp
        type(kernel_spl_t), intent(inout) :: K_j_B_llp
        type(gauss_config_t) :: gauss_conf
        integer :: l, lp, sp
        real(dp) :: dmax_global, alpha, tau
        integer :: sigma
        integer :: total_iterations, current_iteration
        integer(kind=8) :: start_count, count_rate, count_max

        if (turn_off_electrons .and. turn_off_ions) then
            error stop 'Cannot turn off both electrons and ions!'
        end if

        gauss_conf%Nx = gauss_int_nodes_Nx
        gauss_conf%Nxp = gauss_int_nodes_Nxp
        gauss_conf%Ntheta = gauss_int_nodes_Ntheta

        call init_gauss_int(gauss_conf)

        if (.not. pref_ready) call compute_cc_prefactors

        ! Compute a global band-limit distance dmax using Larmor taper and skip threshold
        alpha = Larmor_skip_factor
        tau   = max(kernel_taper_skip_threshold, 1.0d-12)
        block
            real(dp) :: rhoT_max
            rhoT_max = 0.0d0
            do sigma = 0, plasma%n_species - 1
                if (allocated(plasma%spec(sigma)%rho_L_cc)) then
                    rhoT_max = max(rhoT_max, maxval(plasma%spec(sigma)%rho_L_cc))
                end if
            end do
            ! dmax_global = alpha * rhoT_max * sqrt(max(log(1.0d0/tau), 0.0d0))
            dmax_global = tau * rhoT_max 
        end block

        ! Calculate actual number of iterations accounting for band-limiting
        ! Also track maximum and minimum distances for diagnostics
        total_iterations = 0
        max_distance_xl_xlp = 0.0d0
        min_distance_xl_xlp = huge(1.0d0)
        max_index_distance = 0
        min_index_distance = huge(1)
        
        ! check how many lower triangle elements will be computed for the loading bar
        ! (upper triangle elements are given by symmetry)
        do l = 1, K_rho_phi_llp%npts_l
            block
                real(dp) :: xl_val, xlp_val, current_distance
                integer :: lp_lo, current_idx_distance
                xl_val = xl_grid%xb(l)
                lp_lo = l
                do
                    if (lp_lo <= 1) exit
                    if (abs(xl_grid%xb(lp_lo-1) - xl_val) > dmax_global) exit
                    lp_lo = lp_lo - 1
                end do

                ! ensure that at least one off diagonal element is computed for each l
                lp_lo = min(lp_lo, l - 1)
                if (lp_lo < 1) lp_lo = 1
                
                ! Track diagnostics for each (l,lp) pair that will be processed
                do lp = max(1,lp_lo), l
                    xlp_val = xl_grid%xb(lp)
                    current_distance = abs(xl_val - xlp_val)
                    current_idx_distance = abs(l - lp)
                    
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
                end do
                
                total_iterations = total_iterations + (l - max(1,lp_lo) + 1)
            end block
        end do
        current_iteration = 0
        if (fstatus >= 1) write(*,*) 'Total band-limited iterations: ', total_iterations
        if (fstatus >= 1) write(*,*) 'dmax_global: ', dmax_global, ' cm'
        if (fstatus >= 2 .and. artificial_debye_case /= 1) then
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
        end if


        write(*,*) 'Filling Fokker-Planck collision kernels (Gauss)...'
        call system_clock(start_count, count_rate, count_max)

        !$omp parallel do schedule(dynamic) default(shared) private(l,lp)
        do l = 1, K_rho_phi_llp%npts_l
            block
                real(dp) :: xl_val
                integer :: lp_lo
                xl_val = xl_grid%xb(l)
                lp_lo = l
                do
                    if (lp_lo <= 1) exit
                    if (abs(xl_grid%xb(lp_lo-1) - xl_val) > dmax_global) exit
                    lp_lo = lp_lo - 1
                end do

                ! ensure that at least one off diagonal element is computed for each l
                lp_lo = min(lp_lo, l - 1)
                if (lp_lo < 1) lp_lo = 1

                do lp = max(1,lp_lo), l

                    if (.not. turn_off_electrons) then
                        call FP_calc_kernel_element_electrons(l, lp, K_rho_phi_llp%Kllp_e(l, lp),&
                                                K_rho_B_llp%Kllp_e(l, lp), &
                                                K_j_phi_llp%Kllp_e(l, lp), &
                                                K_j_B_llp%Kllp_e(l, lp), &
                                                gauss_conf)

                        call check_is_nan(K_rho_phi_llp%Kllp_e(l,lp), 'K_rho_phi_llp electrons', l, lp)
                        call check_is_nan(K_rho_B_llp%Kllp_e(l,lp), 'K_rho_B_llp electrons', l, lp)
                        call check_is_nan(K_j_phi_llp%Kllp_e(l,lp), 'K_j_phi_llp electrons', l, lp)
                        call check_is_nan(K_j_B_llp%Kllp_e(l,lp), 'K_j_B_llp electrons', l, lp)

                        K_rho_phi_llp%Kllp_e(lp, l) = K_rho_phi_llp%Kllp_e(l, lp)
                        K_rho_B_llp%Kllp_e(lp, l) = K_rho_B_llp%Kllp_e(l, lp)
                        K_j_phi_llp%Kllp_e(lp, l) = K_j_phi_llp%Kllp_e(l, lp)
                        K_j_B_llp%Kllp_e(lp, l) = K_j_B_llp%Kllp_e(l, lp)
                    end if

                    if (.not. turn_off_ions) then
                        do sp = 1, plasma%n_species - 1
                            call FP_calc_kernel_element_ions(l, lp, K_rho_phi_llp%Kllp_i(l, lp, sp),&
                                                K_rho_B_llp%Kllp_i(l, lp, sp), &
                                                K_j_phi_llp%Kllp_i(l, lp, sp), &
                                                K_j_B_llp%Kllp_i(l, lp, sp), &
                                                gauss_conf, sp)

                            call check_is_nan(K_rho_phi_llp%Kllp_i(l,lp, sp), 'K_rho_phi_llp ions ', l, lp)
                            call check_is_nan(K_rho_B_llp%Kllp_i(l,lp, sp), 'K_rho_B_llp ions ', l, lp)
                            call check_is_nan(K_j_phi_llp%Kllp_i(l,lp, sp), 'K_j_phi_llp ions ', l, lp)
                            call check_is_nan(K_j_B_llp%Kllp_i(l,lp, sp), 'K_j_B_llp ions ', l, lp)

                            K_rho_phi_llp%Kllp_i(lp, l, sp) = K_rho_phi_llp%Kllp_i(l, lp, sp)
                            K_rho_B_llp%Kllp_i(lp, l, sp) = K_rho_B_llp%Kllp_i(l, lp, sp)
                            K_j_phi_llp%Kllp_i(lp, l, sp) = K_j_phi_llp%Kllp_i(l, lp, sp)
                            K_j_B_llp%Kllp_i(lp, l, sp) = K_j_B_llp%Kllp_i(l, lp, sp)

                        end do
                    end if

                    !$omp atomic
                    current_iteration = current_iteration + 1
                    !$omp critical(loading_bar)
                    if (mod(current_iteration, 32) == 0 .or. current_iteration == total_iterations) then
                        call update_bar(current_iteration, total_iterations, start_count, count_rate)
                    end if
                    !$omp end critical(loading_bar)
                end do
            end block
        end do
        !$omp end parallel do

        K_rho_phi_llp%Kllp = K_rho_phi_llp%Kllp_e
        K_rho_B_llp%Kllp = K_rho_B_llp%Kllp_e
        K_j_phi_llp%Kllp = K_j_phi_llp%Kllp_e
        K_j_B_llp%Kllp = K_j_B_llp%Kllp_e

        do sp = 1, plasma%n_species - 1
            K_rho_phi_llp%Kllp = K_rho_phi_llp%Kllp + K_rho_phi_llp%Kllp_i(:,:,sp)
            K_rho_B_llp%Kllp = K_rho_B_llp%Kllp + K_rho_B_llp%Kllp_i(:,:,sp)
            K_j_phi_llp%Kllp = K_j_phi_llp%Kllp + K_j_phi_llp%Kllp_i(:,:,sp)
            K_j_B_llp%Kllp = K_j_B_llp%Kllp + K_j_B_llp%Kllp_i(:,:,sp)
        end do

        write(*,*)
        write(*,*) 'Finished filling kernels.'

    end subroutine

    subroutine check_is_nan(value, name, l, lp)
        use KIM_kinds_m, only: dp
        implicit none
        complex(dp), intent(in) :: value
        character(len=*), intent(in) :: name
        integer, intent(in) :: l, lp

        if (isnan(real(value))) then
            print *, trim(name)//' is NaN for l = ', l, ' lp = ', lp
            stop
        end if
    end subroutine check_is_nan


    subroutine FP_calc_kernel_zero_FLR_limit_electrons(l, lp, k_rho_phi, k_rho_B, k_j_phi, k_j_B, gauss_conf)
        ! for benchmarking FLR2 and KIM against each other in the zero Larmor radius limit for electrons
        ! Debye term is omitted since quasineutrality is exploited in FLR2

        use KIM_kinds_m, only: dp
        use integrals_gauss_m, only: gauss_config_t
        use integrands_gauss_m, only: integration_point_t
        use species_m, only: plasma
        use constants_m, only: pi
        use grid_m, only: rg_grid
        use config_m, only: turn_off_ions, turn_off_electrons
        use functions_m, only: varphi_l

        implicit none

        integer, intent(in) :: l, lp
        complex(dp), intent(inout) :: k_rho_phi, k_rho_B, k_j_phi, k_j_B
        integer :: j
        real(dp) :: delta_rg_local

        type(gauss_config_t), intent(in) :: gauss_conf
        type(integration_point_t) :: int_point

        if (abs(l-lp)>2) then
            k_rho_phi = (0.0d0, 0.0d0)
            k_rho_B = (0.0d0, 0.0d0)
            k_j_phi = (0.0d0, 0.0d0)
            k_j_B = (0.0d0, 0.0d0)
            return
        end if

        k_rho_phi = (0.0d0, 0.0d0)
        k_rho_B = (0.0d0, 0.0d0)
        k_j_phi = (0.0d0, 0.0d0)
        k_j_B = (0.0d0, 0.0d0)

        call set_xl_at_edge(l, lp, int_point)

        ! Composite trapezoidal rule for non-equidistant grids
        ! All terms use the same grid structure: evaluate at cell centers but use boundary spacing
        do j = 1, rg_grid%npts_c-1

            int_point%j = j
            int_point%rhoT = plasma%spec(0)%rho_L_cc(j)

            ! Compute local spacing between adjacent boundary points
            ! This is the natural cell width for trapezoidal integration
            if (j == 1) then
                delta_rg_local = 0.5d0 * (rg_grid%xc(j+1) - rg_grid%xc(j))
            else if (j == rg_grid%npts_c) then
                delta_rg_local = 0.5d0 * (rg_grid%xc(j+1) - rg_grid%xc(j))
            else
                delta_rg_local = (rg_grid%xc(j+1) - rg_grid%xc(j))
            end if

            k_rho_phi = k_rho_phi + delta_rg_local &
                * ( & ! debye term plus first order
                    pref_rho_phi_g0(1, j) + pref_rho_phi_g1(1, j, 0) & ! zeroth mphi order is sufficient for electrons
                ) &
                * varphi_l(rg_grid%xc(j), int_point%xlm1, int_point%xl, int_point%xlp1) &
                * varphi_l(rg_grid%xc(j), int_point%xlpm1, int_point%xlp, int_point%xlpp1)

            k_rho_B = k_rho_B + delta_rg_local * pref_rho_B_g1(1, j, 0) &
                * varphi_l(rg_grid%xc(j), int_point%xlm1, int_point%xl, int_point%xlp1) &
                * varphi_l(rg_grid%xc(j), int_point%xlpm1, int_point%xlp, int_point%xlpp1)

            k_j_phi = k_j_phi + delta_rg_local * pref_j_phi_g1(1, j, 0) &
                * varphi_l(rg_grid%xc(j), int_point%xlm1, int_point%xl, int_point%xlp1) &
                * varphi_l(rg_grid%xc(j), int_point%xlpm1, int_point%xlp, int_point%xlpp1)

            k_j_B = k_j_B + delta_rg_local * pref_j_B_g1(1, j, 0) &
                * varphi_l(rg_grid%xc(j), int_point%xlm1, int_point%xl, int_point%xlp1) &
                * varphi_l(rg_grid%xc(j), int_point%xlpm1, int_point%xlp, int_point%xlpp1)

        end do

        ! Apply normalization factor
        k_rho_phi = k_rho_phi / (4.0d0 * pi)
        k_rho_B = k_rho_B / (4.0d0 * pi)
        k_j_phi = k_j_phi / (4.0d0 * pi)
        k_j_B = k_j_B / (4.0d0 * pi)

    end subroutine


    subroutine FP_calc_kernel_element_electrons(l, lp, k_rho_phi, k_rho_B, k_j_phi, k_j_B, gauss_conf)

        use KIM_kinds_m, only: dp
        use integrals_gauss_m, only: gauss_config_t

        implicit none

        integer, intent(in) :: l, lp
        complex(dp), intent(inout) :: k_rho_phi, k_rho_B, k_j_phi, k_j_B
        type(gauss_config_t), intent(in) :: gauss_conf

        k_rho_phi = (0.0d0, 0.0d0)
        k_rho_B = (0.0d0, 0.0d0)
        k_j_phi = (0.0d0, 0.0d0)
        k_j_B = (0.0d0, 0.0d0)

        ! use zero FLR limit for electrons
        call FP_calc_kernel_zero_FLR_limit_electrons(l, lp, k_rho_phi, k_rho_B, k_j_phi, k_j_B, gauss_conf)
            
    end subroutine

    subroutine FP_calc_kernel_element_ions(l, lp, k_rho_phi, k_rho_B, k_j_phi, k_j_B, gauss_conf, sigma)

        use KIM_kinds_m, only: dp
        use integrals_gauss_m, only: gauss_integrate_F0, gauss_integrate_F1, gauss_integrate_F2, gauss_integrate_F3,&
            gauss_config_t
        use species_m, only: plasma
        use constants_m, only: pi
        use integrands_gauss_m, only: gauss_int_F0_rho_phi_t, gauss_int_F1_rho_phi_t, gauss_int_F2_rho_phi_t, gauss_int_F3_rho_phi_t, &
            integration_point_t
        use FP_kernel_plasma_prefacs_m, only: FP_G0_rho_phi
        use grid_m, only: Larmor_skip_factor, rg_grid
        use constants_m, only: com_unit, sol
        use config_m, only: turn_off_ions, turn_off_electrons, artificial_debye_case
        use grid_m, only: xl_grid
        use resonances_mod, only: r_res
        use setup_m, only: mphi_max

        implicit none

        integer, intent(in) :: l, lp, sigma
        complex(dp) :: k_rho_phi, k_rho_B, k_j_phi, k_j_B
        integer :: j, mphi
        type(gauss_config_t), intent(in) :: gauss_conf
        real(dp) :: integral_val
        real(dp) :: current_distance

        type(integration_point_t) :: int_point
        type(gauss_int_F0_rho_phi_t) :: int_F0
        type(gauss_int_F1_rho_phi_t) :: int_F1
        type(gauss_int_F2_rho_phi_t) :: int_F2
        type(gauss_int_F3_rho_phi_t) :: int_F3

        k_rho_phi = (0.0d0, 0.0d0)
        k_rho_B = (0.0d0, 0.0d0)
        k_j_phi = (0.0d0, 0.0d0)
        k_j_B = (0.0d0, 0.0d0)

        call set_xl_at_edge(l, lp, int_point)

        do j = 1, rg_grid%npts_b-1

            int_point%j = j
            int_point%rhoT = max(plasma%spec(sigma)%rho_L_cc(j), 0.0d0)

            ! Debye term (F0) integration
            if (abs(l-lp)<=1 .and. artificial_debye_case /= 2 &
                .and. pref_rho_phi_G0(sigma+1, j) /= 0.0d0) then
                int_F0%int_point = int_point
                call gauss_integrate_F0(int_F0, int_point%xlm1, int_point%xlp1, integral_val, gauss_conf)
                k_rho_phi = k_rho_phi + integral_val * pref_rho_phi_g0(sigma+1, j)
            end if

            if (artificial_debye_case == 1) cycle

            ! skip term if species Larmor radius is too small to couple these grid points
            if (abs(l-lp) > 10 .and. abs(xl_grid%xb(l) - xl_grid%xb(lp))> Larmor_skip_factor * plasma%spec(sigma)%rho_L(j)) cycle
            if (abs(0.5d0 * (rg_grid%xb(j+1) + rg_grid%xb(j)) - 0.5d0 * (xl_grid%xb(l) + xl_grid%xb(lp))) &
                > Larmor_skip_factor * plasma%spec(sigma)%rho_L(j)) cycle

            int_F1%int_point = int_point
            int_F2%int_point = int_point
            int_F3%int_point = int_point

            do mphi = -mphi_max, mphi_max
                int_F1%int_point%mphi = mphi
                ! F1 integration
                call gauss_integrate_F1(int_F1, integral_val, gauss_conf)
                
                k_rho_phi = k_rho_phi + integral_val * pref_rho_phi_g1(sigma+1, j, mphi) * (-1.0d0)**mphi
                k_rho_B = k_rho_B + integral_val * pref_rho_B_g1(sigma+1, j, mphi) * (-1.0d0)**mphi
                k_j_phi = k_j_phi + integral_val * pref_j_phi_g1(sigma+1, j, mphi) * (-1.0d0)**mphi
                k_j_B = k_j_B + integral_val * pref_j_B_g1(sigma+1, j, mphi) * (-1.0d0)**mphi
            end do

            ! ignore FLR terms if resonance is too far from grid points
             if (abs(xl_grid%xb(l) - r_res) > 10.0d0 * int_point%rhoT .or. &
                abs(xl_grid%xb(lp) - r_res) > 10.0d0 * int_point%rhoT) then
                cycle
            end if

            do mphi = -mphi_max, mphi_max
                int_F2%int_point%mphi = mphi
                ! F2 integration
                call gauss_integrate_F2(int_F2, integral_val, gauss_conf)
                k_rho_phi = k_rho_phi + integral_val * pref_rho_phi_g2(sigma+1, j, mphi) * (-1.0d0)**mphi
                k_rho_B = k_rho_B + integral_val * pref_rho_B_g2(sigma+1, j, mphi) * (-1.0d0)**mphi
                k_j_phi = k_j_phi + integral_val * pref_j_phi_g2(sigma+1, j, mphi) * (-1.0d0)**mphi
                k_j_B = k_j_B + integral_val * pref_j_B_g2(sigma+1, j, mphi) * (-1.0d0)**mphi

                int_F3%int_point%mphi = mphi
                ! F3 integration
                call gauss_integrate_F3(int_F3, integral_val, gauss_conf)
                k_rho_phi = k_rho_phi + integral_val * pref_rho_phi_g3(sigma+1, j, mphi) * (-1.0d0)**mphi
                k_rho_B = k_rho_B + integral_val * pref_rho_B_g3(sigma+1, j, mphi) * (-1.0d0)**mphi
                k_j_phi = k_j_phi + integral_val * pref_j_phi_g3(sigma+1,j, mphi) * (-1.0d0)**mphi
                k_j_B = k_j_B + integral_val * pref_j_B_g3(sigma+1, j, mphi) * (-1.0d0)**mphi
            end do

        end do

        k_rho_phi = k_rho_phi / (8.0d0 * pi**3.0d0)
        k_rho_B = k_rho_B / (8.0d0 * pi**3.0d0)

        k_j_phi = k_j_phi / (8.0d0 * pi**3.0d0)
        k_j_B = k_j_B / (8.0d0 * pi**3.0d0)

    end subroutine


    subroutine FP_fill_kernels_flr2_benchmark(K_rho_phi_llp, K_rho_B_llp, K_j_phi_llp, K_j_B_llp)

        use KIM_kinds_m, only: dp
        use integrals_gauss_m, only: gauss_config_t, init_gauss_int
        use grid_m, only: Larmor_skip_factor, gauss_int_nodes_Ntheta, gauss_int_nodes_Nx, gauss_int_nodes_Nxp, &
                        kernel_taper_skip_threshold, rg_grid, xl_grid
        use species_m, only: plasma
        use config_m, only: output_path, artificial_debye_case, fstatus, turn_off_ions, &
                            turn_off_electrons
        use constants_m, only: pi

        implicit none

        type(kernel_spl_t), intent(inout) :: K_rho_phi_llp
        type(kernel_spl_t), intent(inout) :: K_rho_B_llp
        type(kernel_spl_t), intent(inout) :: K_j_phi_llp
        type(kernel_spl_t), intent(inout) :: K_j_B_llp
        type(gauss_config_t) :: gauss_conf
        integer :: l, lp, sp
        real(dp) :: dmax_global, alpha, tau
        integer :: sigma
        integer :: total_iterations, current_iteration
        integer(kind=8) :: start_count, count_rate, count_max

        if (turn_off_electrons .and. turn_off_ions) then
            error stop 'Cannot turn off both electrons and ions!'
        end if

        gauss_conf%Nx = gauss_int_nodes_Nx
        gauss_conf%Nxp = gauss_int_nodes_Nxp
        gauss_conf%Ntheta = gauss_int_nodes_Ntheta

        call init_gauss_int(gauss_conf)
        call rescale_susceptibility_functions()
        if (.not. pref_ready) call compute_cc_prefactors
        pref_rho_phi_G0 = 0.0d0 ! in FLR2 benchmark, 1/lambda_D^2 term is cancelled due to combination of thermodyanmic forces

        ! Compute a global band-limit distance dmax using Larmor taper and skip threshold
        alpha = Larmor_skip_factor
        tau   = max(kernel_taper_skip_threshold, 1.0d-12)
        block
            real(dp) :: rhoT_max
            rhoT_max = 0.0d0
            do sigma = 0, plasma%n_species - 1
                if (allocated(plasma%spec(sigma)%rho_L_cc)) then
                    rhoT_max = max(rhoT_max, maxval(plasma%spec(sigma)%rho_L_cc))
                end if
            end do
            ! dmax_global = alpha * rhoT_max * sqrt(max(log(1.0d0/tau), 0.0d0))
            dmax_global = tau * rhoT_max 
        end block

        ! Calculate actual number of iterations accounting for band-limiting
        ! Also track maximum and minimum distances for diagnostics
        total_iterations = 0
        max_distance_xl_xlp = 0.0d0
        min_distance_xl_xlp = huge(1.0d0)
        max_index_distance = 0
        min_index_distance = huge(1)
        
        do l = 1, K_rho_phi_llp%npts_l
            block
                real(dp) :: xl_val, xlp_val, current_distance
                integer :: lp_lo, current_idx_distance
                xl_val = xl_grid%xb(l)
                lp_lo = l
                do
                    if (lp_lo <= 1) exit
                    if (abs(xl_grid%xb(lp_lo-1) - xl_val) > dmax_global) exit
                    lp_lo = lp_lo - 1
                end do

                lp_lo = min(lp_lo, l - 1)
                if (lp_lo < 1) lp_lo = 1
                
                ! Track diagnostics for each (l,lp) pair that will be processed
                do lp = max(1,lp_lo), l
                    xlp_val = xl_grid%xb(lp)
                    current_distance = abs(xl_val - xlp_val)
                    current_idx_distance = abs(l - lp)
                    
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
                end do
                
                total_iterations = total_iterations + (l - max(1,lp_lo) + 1)
            end block
        end do
        current_iteration = 0
        if (fstatus >= 1) write(*,*) 'Total band-limited iterations: ', total_iterations
        if (fstatus >= 1) write(*,*) 'dmax_global: ', dmax_global, ' cm'
        if (fstatus >= 1 .and. artificial_debye_case /= 1) then
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
        end if


        write(*,*) 'Filling Fokker-Planck collision kernels (Gauss)...'
        call system_clock(start_count, count_rate, count_max)

        !$omp parallel do schedule(dynamic) default(shared) private(l,lp)
        do l = 1, K_rho_phi_llp%npts_l
            block
                real(dp) :: xl_val
                integer :: lp_lo
                xl_val = xl_grid%xb(l)
                lp_lo = l
                do
                    if (lp_lo <= 1) exit
                    if (abs(xl_grid%xb(lp_lo-1) - xl_val) > dmax_global) exit
                    lp_lo = lp_lo - 1
                end do

                ! Ensure at least one off-diagonal element (force lp_lo ≤ l-1)
                lp_lo = min(lp_lo, l - 1)
                if (lp_lo < 1) lp_lo = 1

                do lp = max(1,lp_lo), l

                    if (.not. turn_off_electrons) then
                        call FP_calc_kernel_electrons_FLR2_benchmark(l, lp, K_rho_phi_llp%Kllp_e(l, lp),&
                                                                    K_rho_B_llp%Kllp_e(l, lp), &
                                                                    K_j_phi_llp%Kllp_e(l, lp), &
                                                                    K_j_B_llp%Kllp_e(l, lp), &
                                                                    gauss_conf)

                        call check_is_nan(K_rho_phi_llp%Kllp_e(l,lp), 'K_rho_phi_llp electrons', l, lp)
                        call check_is_nan(K_rho_B_llp%Kllp_e(l,lp), 'K_rho_B_llp electrons', l, lp)
                        call check_is_nan(K_j_phi_llp%Kllp_e(l,lp), 'K_j_phi_llp electrons', l, lp)
                        call check_is_nan(K_j_B_llp%Kllp_e(l,lp), 'K_j_B_llp electrons', l, lp)

                        K_rho_phi_llp%Kllp_e(lp, l) = K_rho_phi_llp%Kllp_e(l, lp)
                        K_rho_B_llp%Kllp_e(lp, l) = K_rho_B_llp%Kllp_e(l, lp)
                        K_j_phi_llp%Kllp_e(lp, l) = K_j_phi_llp%Kllp_e(l, lp)
                        K_j_B_llp%Kllp_e(lp, l) = K_j_B_llp%Kllp_e(l, lp)
                    end if

                    if (.not. turn_off_ions) then
                        do sp = 1, plasma%n_species - 1
                            call FP_calc_kernel_ions_FLR2_benchmark(l, lp, K_rho_phi_llp%Kllp_i(l, lp, sp),&
                                                K_rho_B_llp%Kllp_i(l, lp, sp), &
                                                K_j_phi_llp%Kllp_i(l, lp, sp), &
                                                K_j_B_llp%Kllp_i(l, lp, sp), &
                                                gauss_conf, sp)

                            call check_is_nan(K_rho_phi_llp%Kllp_i(l,lp, sp), 'K_rho_phi_llp ions ', l, lp)
                            call check_is_nan(K_rho_B_llp%Kllp_i(l,lp, sp), 'K_rho_B_llp ions ', l, lp)
                            call check_is_nan(K_j_phi_llp%Kllp_i(l,lp, sp), 'K_j_phi_llp ions ', l, lp)
                            call check_is_nan(K_j_B_llp%Kllp_i(l,lp, sp), 'K_j_B_llp ions ', l, lp)

                            K_rho_phi_llp%Kllp_i(lp, l, sp) = K_rho_phi_llp%Kllp_i(l, lp, sp)
                            K_rho_B_llp%Kllp_i(lp, l, sp) = K_rho_B_llp%Kllp_i(l, lp, sp)
                            K_j_phi_llp%Kllp_i(lp, l, sp) = K_j_phi_llp%Kllp_i(l, lp, sp)
                            K_j_B_llp%Kllp_i(lp, l, sp) = K_j_B_llp%Kllp_i(l, lp, sp)

                        end do
                    end if

                    !$omp atomic
                    current_iteration = current_iteration + 1
                    !$omp critical(loading_bar)
                    if (mod(current_iteration, 32) == 0 .or. current_iteration == total_iterations) then
                        call update_bar(current_iteration, total_iterations, start_count, count_rate)
                    end if
                    !$omp end critical(loading_bar)
                end do
            end block
        end do
        !$omp end parallel do

        K_rho_phi_llp%Kllp = K_rho_phi_llp%Kllp_e
        K_rho_B_llp%Kllp = K_rho_B_llp%Kllp_e
        K_j_phi_llp%Kllp = K_j_phi_llp%Kllp_e
        K_j_B_llp%Kllp = K_j_B_llp%Kllp_e

        do sp = 1, plasma%n_species - 1
            K_rho_phi_llp%Kllp = K_rho_phi_llp%Kllp + K_rho_phi_llp%Kllp_i(:,:,sp)
            K_rho_B_llp%Kllp = K_rho_B_llp%Kllp + K_rho_B_llp%Kllp_i(:,:,sp)
            K_j_phi_llp%Kllp = K_j_phi_llp%Kllp + K_j_phi_llp%Kllp_i(:,:,sp)
            K_j_B_llp%Kllp = K_j_B_llp%Kllp + K_j_B_llp%Kllp_i(:,:,sp)
        end do

        write(*,*)
        write(*,*) 'Finished filling kernels.'

        contains
        
        subroutine rescale_susceptibility_functions()
            ! Rescale susceptibility functions for FLR2 benchmark
            ! has to be done in order to take quasineutrality into account as is done in FLR2
            
            use species_m, only: plasma
            use IO_collection_m, only: write_complex_profile, itoa
            use grid_m, only: rg_grid
            use config_m, only: output_path
            use setup_m, only: mphi_max

            implicit none

            integer :: sp
            integer :: j
            integer :: mphi

            do sp=0, plasma%n_species - 1
                ! Note: I00, I20, x1, x2 are now only available as _cc (cell center) versions
                ! They are computed after interpolation in calculate_thermodynamic_forces_and_susc
                do j = 1, rg_grid%npts_c
                    plasma%spec(sp)%I00_cc(j, :) = plasma%spec(sp)%I01_cc(j, :) * plasma%spec(sp)%x1_cc(j) / plasma%spec(sp)%x2_cc(j, :)
                    plasma%spec(sp)%I20_cc(j, :) = plasma%spec(sp)%I21_cc(j, :) * plasma%spec(sp)%x1_cc(j) / plasma%spec(sp)%x2_cc(j, :)
                end do

                do mphi = -mphi_max, mphi_max
                    call write_complex_profile(rg_grid%xc, plasma%spec(sp)%I01_cc(:, mphi), rg_grid%npts_c, &
                        'backs/'//trim(plasma%spec(sp)%name)//'/I01_cc_after_resc_mphi_0'//trim(adjustl(itoa(mphi))))
                    call write_complex_profile(rg_grid%xc, plasma%spec(sp)%I20_cc(:, mphi), rg_grid%npts_c, &
                        'backs/'//trim(plasma%spec(sp)%name)//'/I20_cc_after_resc_mphi_0'//trim(adjustl(itoa(mphi))))
                end do
            end do

            
        end subroutine

    end subroutine FP_fill_kernels_flr2_benchmark


    subroutine FP_calc_kernel_electrons_FLR2_benchmark(l, lp, k_rho_phi, k_rho_B, k_j_phi, k_j_B, gauss_conf)
        ! for benchmarking FLR2 and KIM against each other, FLR2 exploits quasineutrality which has to be taken into account for
        ! electrons and ions separately, most notably, by omitting the Debye shielding term and using different susceptibility functions

        use KIM_kinds_m, only: dp
        use integrals_gauss_m, only: gauss_config_t

        implicit none

        integer, intent(in) :: l, lp
        complex(dp), intent(inout) :: k_rho_phi, k_rho_B, k_j_phi, k_j_B
        type(gauss_config_t), intent(in) :: gauss_conf

        k_rho_phi = (0.0d0, 0.0d0)
        k_rho_B = (0.0d0, 0.0d0)
        k_j_phi = (0.0d0, 0.0d0)
        k_j_B = (0.0d0, 0.0d0)

        call FP_calc_kernel_zero_FLR_limit_electrons(l, lp, k_rho_phi, k_rho_B, k_j_phi, k_j_B, gauss_conf)

            
    end subroutine FP_calc_kernel_electrons_FLR2_benchmark

    subroutine FP_calc_kernel_ions_FLR2_benchmark(l, lp, k_rho_phi, k_rho_B, k_j_phi, k_j_B, gauss_conf, sigma)
        ! for benchmarking FLR2 and KIM against each other, FLR2 exploits quasineutrality which has to be taken into account for
        ! electrons and ions separately, most notably, by omitting the Debye shielding term and using different susceptibility functions

        use KIM_kinds_m, only: dp
        use integrals_gauss_m, only: gauss_config_t

        implicit none

        integer, intent(in) :: l, lp, sigma
        complex(dp), intent(inout) :: k_rho_phi, k_rho_B, k_j_phi, k_j_B
        type(gauss_config_t), intent(in) :: gauss_conf

        call FP_calc_kernel_element_ions(l, lp, k_rho_phi, k_rho_B, k_j_phi, k_j_B, gauss_conf, sigma)
        
    end subroutine FP_calc_kernel_ions_FLR2_benchmark

    
    subroutine set_xl_at_edge(l, lp, int_point)
        
        use grid_m, only: xl_grid
        use KIM_kinds_m, only: dp
        use integrands_gauss_m, only: integration_point_t

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

    
    subroutine write_kernels(kernel_rho_phi_llp, kernel_rho_B_llp, kernel_j_phi_llp, kernel_j_B_llp)

        use IO_collection_m, only: write_matrix
        use grid_m, only: xl_grid
        use species_m, only: plasma

        implicit none

        type(kernel_spl_t), intent(in) :: kernel_rho_phi_llp
        type(kernel_spl_t), intent(in) :: kernel_rho_B_llp
        type(kernel_spl_t), intent(in) :: kernel_j_phi_llp
        type(kernel_spl_t), intent(in) :: kernel_j_B_llp

        integer ::sp

        call write_matrix("kernel/K_rho_phi", real(kernel_rho_phi_llp%Kllp), xl_grid%npts_b, xl_grid%npts_b, &
            'Complex FLR2 benchmark kernel K_rho_phi', '1/cm^2')
        call write_matrix("kernel/K_rho_B", real(kernel_rho_B_llp%Kllp), xl_grid%npts_b, xl_grid%npts_b, &
            'Complex FLR2 benchmark kernel K_rho_B', '1/cm^2')
        call write_matrix("kernel/K_j_phi", real(kernel_j_phi_llp%Kllp), xl_grid%npts_b, xl_grid%npts_b, &
            'Complex FLR2 benchmark kernel K_j_phi', '1/cm^2')
        call write_matrix("kernel/K_j_B", real(kernel_j_B_llp%Kllp), xl_grid%npts_b, xl_grid%npts_b, &
            'Complex FLR2 benchmark kernel K_j_B', '1/cm^2')

        call write_matrix("kernel/K_rho_phi_e",   real(kernel_rho_phi_llp%Kllp_e), xl_grid%npts_b, xl_grid%npts_b, &
            'Complex FLR2 benchmark kernel K_rho_phi electrons', '1/cm^2')
        call write_matrix("kernel/K_rho_B_e",     real(kernel_rho_B_llp%Kllp_e), xl_grid%npts_b, xl_grid%npts_b, &
            'Complex FLR2 benchmark kernel K_rho_B electrons', '1/cm^2')
        call write_matrix("kernel/K_j_phi_e",     real(kernel_j_phi_llp%Kllp_e), xl_grid%npts_b, xl_grid%npts_b, &
            'Complex FLR2 benchmark kernel K_j_phi electrons', '1/cm^2')
        call write_matrix("kernel/K_j_B_e",       real(kernel_j_B_llp%Kllp_e), xl_grid%npts_b, xl_grid%npts_b, &
            'Complex FLR2 benchmark kernel K_j_B electrons', '1/cm^2')

        do sp = 1, plasma%n_species - 1
            call write_matrix("kernel/K_rho_phi_"//trim(plasma%spec(sp)%name),  real(kernel_rho_phi_llp%Kllp_i(:,:,sp)), xl_grid%npts_b, xl_grid%npts_b, &
                'Complex FLR2 benchmark kernel K_rho_phi electrons', '1/cm^2')
            call write_matrix("kernel/K_rho_B_"//trim(plasma%spec(sp)%name),     real(kernel_rho_B_llp%Kllp_i(:,:,sp)), xl_grid%npts_b, xl_grid%npts_b, &
                'Complex FLR2 benchmark kernel K_rho_B electrons', '1/cm^2')
            call write_matrix("kernel/K_j_phi_"//trim(plasma%spec(sp)%name),     real(kernel_j_phi_llp%Kllp_i(:,:,sp)), xl_grid%npts_b, xl_grid%npts_b, &
                'Complex FLR2 benchmark kernel K_j_phi electrons', '1/cm^2')
            call write_matrix("kernel/K_j_B_"//trim(plasma%spec(sp)%name),       real(kernel_j_B_llp%Kllp_i(:,:,sp)), xl_grid%npts_b, xl_grid%npts_b, &
                'Complex FLR2 benchmark kernel K_j_B electrons', '1/cm^2')
        end do


    end subroutine

    subroutine fill_kernels_krook_fp(kernel_krook_rho_phi, kernel_krook_rho_B, &
                                        kernel_fp_rho_phi, kernel_fp_rho_B)
        !> Unified subroutine to fill both Krook and Fokker-Planck kernels
        !> Exploits shared Gaussian integration for efficiency
        
        implicit none
        
        type(kernel_spl_t), intent(inout) :: kernel_krook_rho_phi, kernel_krook_rho_B
        type(kernel_spl_t), intent(inout) :: kernel_fp_rho_phi, kernel_fp_rho_B
        
        error stop 'fill_kernels_krook_fp currently not implemented'


    end subroutine fill_kernels_krook_fp
    

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! currently unused subroutines:

    ! subroutine FP_calc_kernel_element(l, lp, k_rho_phi, k_rho_B, k_j_phi, k_j_B, gauss_conf)

        ! use KIM_kinds_m, only: dp
        ! use integrals_gauss_m, only: gauss_integrate_F0, gauss_integrate_F1, gauss_integrate_F2, gauss_integrate_F3,&
            ! gauss_config_t, gauss_integrate_F1_electrons, gauss_integrate_F2_electrons
        ! use species_m, only: plasma
        ! use constants_m, only: pi
        ! use integrands_gauss_m, only: gauss_int_F0_rho_phi_t, gauss_int_F1_rho_phi_t, gauss_int_F2_rho_phi_t, gauss_int_F3_rho_phi_t, &
            ! integration_point_t, gauss_int_F1_rho_phi_electrons_t, gauss_int_F2_rho_phi_electrons_t
        ! use FP_kernel_plasma_prefacs_m, only: FP_G0_rho_phi
        ! use grid_m, only: Larmor_skip_factor, kernel_taper_skip_threshold, rg_grid, xl_grid
        ! use constants_m, only: com_unit, sol
        ! use config_m, only: turn_off_ions, turn_off_electrons, artificial_debye_case
        ! use functions_m, only: varphi_l

        ! implicit none

        ! integer, intent(in) :: l, lp
        ! complex(dp) :: k_rho_phi, k_rho_B, k_j_phi, k_j_B
        ! integer :: j, sigma
        ! type(gauss_config_t), intent(in) :: gauss_conf
        ! real(dp) :: integral_val
        ! real(dp) :: current_distance

        ! type(integration_point_t) :: int_point
        ! type(gauss_int_F0_rho_phi_t) :: int_F0
        ! type(gauss_int_F1_rho_phi_t) :: int_F1
        ! type(gauss_int_F2_rho_phi_t) :: int_F2
        ! type(gauss_int_F3_rho_phi_t) :: int_F3

        ! type(gauss_int_F1_rho_phi_electrons_t) :: int_F1_e
        ! type(gauss_int_F2_rho_phi_electrons_t) :: int_F2_e

        ! k_rho_phi = (0.0d0, 0.0d0)
        ! k_rho_B = (0.0d0, 0.0d0)
        ! k_j_phi = (0.0d0, 0.0d0)
        ! k_j_B = (0.0d0, 0.0d0)

        ! call set_xl_at_edge(l, lp, int_point)

        ! if (.not. turn_off_electrons) then
            ! ! zero FLR limit is sufficient for electrons
            ! call FP_calc_kernel_zero_FLR_limit_electrons(l, lp, k_rho_phi, k_rho_B, k_j_phi, k_j_B, gauss_conf)

            ! ! to counteract the division by (8 pi^3) later in the subroutine
            ! k_rho_phi = k_rho_phi * (8.0d0 * pi**3.0d0)
            ! k_rho_B = k_rho_B * (8.0d0 * pi**3.0d0)

            ! k_j_phi = k_j_phi * (8.0d0 * pi**3.0d0)
            ! k_j_B = k_j_B * (8.0d0 * pi**3.0d0)

        ! end if

        ! if (.not. turn_off_ions) then
            ! do sigma = 1, plasma%n_species - 1
                ! do j = 1, rg_grid%npts_b-1

                    ! int_point%j = j
                    ! int_point%rhoT = max(plasma%spec(sigma)%rho_L_cc(j), 0.0d0)

                    ! if (abs(l-lp)<=1 .and. artificial_debye_case /= 2) then
                        ! int_F0%int_point = int_point
                        ! call gauss_integrate_F0(int_F0, int_point%xlm1, int_point%xlp1, integral_val, gauss_conf)
                        ! k_rho_phi = k_rho_phi + integral_val * (-1.0d0) * (1.0d0 / (plasma%spec(sigma)%lambda_D_cc(j)**2.0d0))
                    ! end if

                    ! if (artificial_debye_case == 1) cycle

                    ! ! skip term if species Larmor radius is too small to couple these grid points
                    ! if (abs(l-lp) > 4 .and. abs(xl_grid%xb(l) - xl_grid%xb(lp))> 4.0d0 * plasma%spec(sigma)%rho_L(j)) cycle
                    ! if (abs(0.5d0 * (rg_grid%xb(j+1) + rg_grid%xb(j)) - 0.5d0 * (xl_grid%xb(l) + xl_grid%xb(lp))) > 4.0d0 * plasma%spec(sigma)%rho_L(j)) cycle

                    ! int_F1%int_point = int_point
                    ! int_F2%int_point = int_point
                    ! int_F3%int_point = int_point

                    ! ! F1 integration
                    ! call gauss_integrate_F1(int_F1, integral_val, gauss_conf)
                    ! k_rho_phi = k_rho_phi + integral_val * pref_rho_phi_g1(sigma+1,j)
                    ! k_rho_B = k_rho_B + integral_val * pref_rho_B_g1(sigma+1,j)
                    ! k_j_phi = k_j_phi + integral_val * pref_j_phi_g1(sigma+1,j)
                    ! k_j_B = k_j_B + integral_val * pref_j_B_g1(sigma+1,j)

                    ! ! cycle ! for testing (makes it cheaper)

                    ! ! F2 integration
                    ! call gauss_integrate_F2(int_F2, integral_val, gauss_conf)
                    ! k_rho_phi = k_rho_phi + integral_val * pref_rho_phi_g2(sigma+1,j)
                    ! k_rho_B = k_rho_B + integral_val * pref_rho_B_g2(sigma+1,j)
                    ! k_j_phi = k_j_phi + integral_val * pref_j_phi_g2(sigma+1,j)
                    ! k_j_B = k_j_B + integral_val * pref_j_B_g2(sigma+1,j)

                    ! ! F3 integration
                    ! call gauss_integrate_F3(int_F3, integral_val, gauss_conf)
                    ! k_rho_phi = k_rho_phi + integral_val * pref_rho_phi_g3(sigma+1,j)
                    ! k_rho_B = k_rho_B + integral_val * pref_rho_B_g3(sigma+1,j)
                    ! k_j_phi = k_j_phi + integral_val * pref_j_phi_g3(sigma+1,j)
                    ! k_j_B = k_j_B + integral_val * pref_j_B_g3(sigma+1,j)

                ! end do
            ! end do
        ! end if

        ! k_rho_phi = k_rho_phi / (8.0d0 * pi**3.0d0)
        ! k_rho_B = k_rho_B / (8.0d0 * pi**3.0d0)

        ! k_j_phi = k_j_phi / (8.0d0 * pi**3.0d0)
        ! k_j_B = k_j_B / (8.0d0 * pi**3.0d0)

            
    ! end subroutine


    ! ! TODO: Update combined Krook FP routine to use the same procedure as FP case
    ! subroutine fill_kernels_krook_fp(kernel_krook_rho_phi, kernel_krook_rho_B, &
                                      ! kernel_fp_rho_phi, kernel_fp_rho_B)
        ! !> Unified subroutine to fill both Krook and Fokker-Planck kernels
        ! !> Exploits shared Gaussian integration for efficiency
        
        ! use KIM_kinds_m, only: dp
        ! use integrals_gauss_m, only: gauss_config_t, init_gauss_int, &
            ! gauss_integrate_F0, gauss_integrate_F1, gauss_integrate_F2, gauss_integrate_F3
        ! use grid_m, only: Larmor_skip_factor, gauss_int_nodes_Ntheta, gauss_int_nodes_Nx, gauss_int_nodes_Nxp
        ! use species_m, only: plasma
        ! use constants_m, only: pi
        ! use integrands_gauss_m, only: gauss_int_F0_rho_phi_t, gauss_int_F1_rho_phi_t, &
            ! gauss_int_F2_rho_phi_t, gauss_int_F3_rho_phi_t, integration_point_t
        ! use Krook_kernel_plasma_prefacs_m, only: Krook_G0_rho_phi, Krook_G1_rho_phi, Krook_G2_rho_phi, Krook_G3_rho_phi, &
            ! Krook_G1_rho_B, Krook_G2_rho_B, Krook_G3_rho_B, Krook_kappa_rho_phi, Krook_kappa_rho_B
        ! use FP_kernel_plasma_prefacs_m, only: FP_G0_rho_phi, FP_G1_rho_phi, FP_G2_rho_phi, &
            ! FP_G3_rho_phi, FP_G1_rho_B, FP_G2_rho_B, FP_G3_rho_B, &
            ! FP_kappa_rho_phi, FP_kappa_rho_B
        ! use config_m, only: artificial_debye_case, turn_off_ions
        
        ! implicit none
        
        ! type(kernel_spl_t), intent(inout) :: kernel_krook_rho_phi, kernel_krook_rho_B
        ! type(kernel_spl_t), intent(inout) :: kernel_fp_rho_phi, kernel_fp_rho_B
        
        ! type(gauss_config_t) :: gauss_conf
        ! integer :: l, lp, j, sigma
        ! complex(dp) :: krook_phi_llp, krook_B_llp, fp_phi_llp, fp_B_llp
        ! real(dp) :: integral_F0, integral_F1, integral_F2, integral_F3
        ! real(dp) :: current_distance
        ! integer :: current_idx_distance
        
        ! type(integration_point_t) :: int_point
        ! type(gauss_int_F0_rho_phi_t) :: int_F0
        ! type(gauss_int_F1_rho_phi_t) :: int_F1
        ! type(gauss_int_F2_rho_phi_t) :: int_F2
        ! type(gauss_int_F3_rho_phi_t) :: int_F3
        
        ! ! Initialize Gaussian integration configuration
        ! gauss_conf%Nx = gauss_int_nodes_Nx
        ! gauss_conf%Nxp = gauss_int_nodes_Nxp
        ! gauss_conf%Ntheta = gauss_int_nodes_Ntheta
        ! call init_gauss_int(gauss_conf)

        ! write(*,*) 'Filling both Krook and Fokker-Planck kernels simultaneously...'
        
        ! !$omp parallel do collapse(1) private(l, lp, krook_phi_llp, krook_B_llp, &
        ! !$omp& fp_phi_llp, fp_B_llp, j, sigma, int_point, int_F0, int_F1, int_F2, int_F3, &
        ! !$omp& integral_F0, integral_F1, integral_F2, integral_F3, current_distance, current_idx_distance)
        ! do l = 1, kernel_krook_rho_phi%npts_l
            ! do lp = 1, l
                
                ! ! Initialize kernel values
                ! krook_phi_llp = 0.0d0
                ! krook_B_llp = 0.0d0
                ! fp_phi_llp = 0.0d0
                ! fp_B_llp = 0.0d0
                
                ! ! Set xl grid points at edges
                ! call set_xl_at_edge(l, lp, int_point)
                
                ! ! Loop over species and radial grid
                ! do sigma = 0, plasma%n_species - 1
                    ! if (turn_off_ions .and. sigma >= 1) cycle
                    ! do j = 2, size(plasma%r_grid)-1
                        ! int_point%j = j
                        ! int_point%rhoT = 0.5d0 * (plasma%spec(sigma)%rho_L(j) + plasma%spec(sigma)%rho_L(j+1))
                        
                        ! ! F0 integral (only for diagonal elements)
                        ! if (l == lp) then
                            ! int_F0%int_point = int_point
                            ! call gauss_integrate_F0(int_F0, int_point%xlm1, int_point%xlp1, integral_F0, gauss_conf)
                            
                            ! ! Krook contribution
                            ! krook_phi_llp = krook_phi_llp + integral_F0 * Krook_G0_rho_phi(j, plasma%spec(sigma)) * &
                                            ! Krook_kappa_rho_phi(j, plasma%spec(sigma))
                            
                            ! ! Fokker-Planck contribution
                            ! fp_phi_llp = fp_phi_llp + integral_F0 * FP_G0_rho_phi(j, plasma%spec(sigma)) * &
                                        ! FP_kappa_rho_phi(j, plasma%spec(sigma))
                        ! end if

                        ! ! Track maximum distances for diagnostics
                        ! current_distance = abs(int_point%xl - int_point%xlp)

                        ! if (current_distance > Larmor_skip_factor * int_point%rhoT) cycle

                        ! current_idx_distance = abs(l - lp)
                        
                        ! !$omp critical
                        ! if (current_distance > max_distance_xl_xlp) then
                            ! max_distance_xl_xlp = current_distance
                            ! max_dist_l = l
                            ! max_dist_lp = lp
                        ! end if
                        ! if (current_idx_distance > max_index_distance) then
                            ! max_index_distance = current_idx_distance
                            ! max_idx_l = l
                            ! max_idx_lp = lp
                        ! end if
                        ! !$omp end critical
                        
                        ! if (artificial_debye_case == 1) cycle

                        ! ! Set integration points for F1, F2, F3
                        ! int_F1%int_point = int_point
                        ! int_F2%int_point = int_point
                        ! int_F3%int_point = int_point
                        
                        ! ! Perform Gaussian integrations (shared between Krook and FP)
                        ! call gauss_integrate_F1(int_F1, integral_F1, gauss_conf)
                        ! call gauss_integrate_F2(int_F2, integral_F2, gauss_conf)
                        ! call gauss_integrate_F3(int_F3, integral_F3, gauss_conf)
                        
                        ! ! Krook contributions
                        ! krook_phi_llp = krook_phi_llp  &
                            ! + integral_F1 * Krook_G1_rho_phi(j, plasma%spec(sigma)) * Krook_kappa_rho_phi(j, plasma%spec(sigma))&
                            ! + integral_F2 * Krook_G2_rho_phi(j, plasma%spec(sigma)) * Krook_kappa_rho_phi(j, plasma%spec(sigma)) &
                            ! + integral_F3 * Krook_G3_rho_phi(j, plasma%spec(sigma)) * Krook_kappa_rho_phi(j, plasma%spec(sigma))
                        
                        ! krook_B_llp = krook_B_llp &
                            ! + integral_F1 * Krook_G1_rho_B(j, plasma%spec(sigma)) * Krook_kappa_rho_B(j, plasma%spec(sigma)) &
                            ! + integral_F2 * Krook_G2_rho_B(j, plasma%spec(sigma)) * Krook_kappa_rho_B(j, plasma%spec(sigma)) &
                            ! + integral_F3 * Krook_G3_rho_B(j, plasma%spec(sigma)) * Krook_kappa_rho_B(j, plasma%spec(sigma))
                        
                        ! ! Fokker-Planck contributions
                        ! fp_phi_llp = fp_phi_llp &
                            ! + integral_F1 * FP_G1_rho_phi(j, plasma%spec(sigma)) * FP_kappa_rho_phi(j, plasma%spec(sigma)) &
                            ! + integral_F2 * FP_G2_rho_phi(j, plasma%spec(sigma)) * FP_kappa_rho_phi(j, plasma%spec(sigma)) &
                            ! + integral_F3 * FP_G3_rho_phi(j, plasma%spec(sigma)) * FP_kappa_rho_phi(j, plasma%spec(sigma))
                        
                        ! fp_B_llp = fp_B_llp &
                            ! + integral_F1 * FP_G1_rho_B(j, plasma%spec(sigma)) * FP_kappa_rho_B(j, plasma%spec(sigma)) &
                            ! + integral_F2 * FP_G2_rho_B(j, plasma%spec(sigma)) * FP_kappa_rho_B(j, plasma%spec(sigma)) &
                            ! + integral_F3 * FP_G3_rho_B(j, plasma%spec(sigma)) * FP_kappa_rho_B(j, plasma%spec(sigma))

                    ! end do
                ! end do
                
                ! ! Apply normalization factor
                ! krook_phi_llp = krook_phi_llp / (8.0d0 * pi**3.0d0)
                ! krook_B_llp = krook_B_llp / (8.0d0 * pi**3.0d0)
                ! fp_phi_llp = fp_phi_llp / (8.0d0 * pi**3.0d0)
                ! fp_B_llp = fp_B_llp / (8.0d0 * pi**3.0d0)
                
                ! ! Store results
                ! kernel_krook_rho_phi%Kllp(l, lp) = krook_phi_llp
                ! kernel_krook_rho_B%Kllp(l, lp) = krook_B_llp
                ! kernel_fp_rho_phi%Kllp(l, lp) = fp_phi_llp
                ! kernel_fp_rho_B%Kllp(l, lp) = fp_B_llp
                
                ! ! Check for NaN values
                ! if (isnan(real(krook_phi_llp)) .or. isnan(real(krook_B_llp)) .or. &
                    ! isnan(real(fp_phi_llp)) .or. isnan(real(fp_B_llp))) then
                    ! print *, "NaN detected in kernel calculation at l =", l, ", lp =", lp
                    ! print *, "Krook phi:", krook_phi_llp, "Krook B:", krook_B_llp
                    ! print *, "FP phi:", fp_phi_llp, "FP B:", fp_B_llp
                    ! stop
                ! end if

                ! ! exploit symmetry
                ! kernel_krook_rho_phi%Kllp(lp, l) = kernel_krook_rho_phi%Kllp(l, lp)
                ! kernel_krook_rho_B%Kllp(lp, l) = kernel_krook_rho_B%Kllp(l, lp)
                ! kernel_fp_rho_phi%Kllp(lp, l) = kernel_fp_rho_phi%Kllp(l, lp)
                ! kernel_fp_rho_B%Kllp(lp, l) = kernel_fp_rho_B%Kllp(l, lp)
            ! end do
        ! end do
        ! !$omp end parallel do

        ! write(*,*) ! New line after progress bar
        
        ! ! Print diagnostic information
        ! write(*,*) '======== Kernel Distance Diagnostics (Combined Krook+FP) ========'
        ! write(*,'(A,F12.6)') ' Maximum |xl - xlp| distance: ', max_distance_xl_xlp
        ! write(*,'(A,I6,A,I6)') ' Occurred at l = ', max_dist_l, ', lp = ', max_dist_lp
        ! write(*,'(A,F12.6)') ' Minimum |xl - xlp| distance: ', min_distance_xl_xlp
        ! write(*,'(A,I6,A,I6)') ' Occurred at l = ', min_dist_l, ', lp = ', min_dist_lp
        ! write(*,'(A,I6)') ' Maximum index distance |l - lp|: ', max_index_distance
        ! write(*,'(A,I6,A,I6)') ' Occurred at l = ', max_idx_l, ', lp = ', max_idx_lp
        ! write(*,'(A,I6)') ' Minimum index distance |l - lp|: ', min_index_distance
        ! write(*,'(A,I6,A,I6)') ' Occurred at l = ', min_idx_l, ', lp = ', min_idx_lp
        ! write(*,*) '=================================================================='
        
    ! end subroutine fill_kernels_krook_fp


    ! subroutine FP_fill_kernels_flr2_benchmark_old(K_rho_phi_llp, K_rho_B_llp, K_j_phi_llp, K_j_B_llp, include_electrons, include_ions)

        ! use integrals_gauss_m, only: gauss_config_t, init_gauss_int
        ! use grid_m, only: Larmor_skip_factor, gauss_int_nodes_Ntheta, gauss_int_nodes_Nx, gauss_int_nodes_Nxp, &
                        ! kernel_taper_skip_threshold, xl_grid
        ! use species_m, only: plasma
        ! use config_m, only: fstatus

        ! implicit none

        ! type(kernel_spl_t), intent(inout) :: K_rho_phi_llp
        ! type(kernel_spl_t), intent(inout) :: K_rho_B_llp
        ! type(kernel_spl_t), intent(inout) :: K_j_phi_llp
        ! type(kernel_spl_t), intent(inout) :: K_j_B_llp
        ! logical, intent(in) :: include_electrons
        ! logical, intent(in) :: include_ions

        ! type(gauss_config_t) :: gauss_conf
        ! real(dp) :: dmax_global, alpha, tau
        ! integer :: total_iterations, current_iteration
        ! integer(kind=8) :: start_count, count_rate, count_max
        ! integer :: l, lp, lp_lo, lp_hi

        ! interface
            ! subroutine flr2_kernel_cb(l, lp, k_rho_phi, k_rho_B, k_j_phi, k_j_B, gauss_conf)
                ! use KIM_kinds_m, only: dp
                ! use integrals_gauss_m, only: gauss_config_t
                ! implicit none
                ! integer, intent(in) :: l, lp
                ! complex(dp), intent(inout) :: k_rho_phi, k_rho_B, k_j_phi, k_j_B
                ! type(gauss_config_t), intent(in) :: gauss_conf
            ! end subroutine flr2_kernel_cb
        ! end interface

        ! gauss_conf%Nx = gauss_int_nodes_Nx
        ! gauss_conf%Nxp = gauss_int_nodes_Nxp
        ! gauss_conf%Ntheta = gauss_int_nodes_Ntheta

        ! call init_gauss_int(gauss_conf)
        ! call rescale_susceptibility_functions()
        ! if (.not. pref_ready) call compute_cc_prefactors

        ! K_rho_phi_llp%Kllp = (0.0d0, 0.0d0)
        ! K_rho_B_llp%Kllp   = (0.0d0, 0.0d0)
        ! K_j_phi_llp%Kllp   = (0.0d0, 0.0d0)
        ! K_j_B_llp%Kllp     = (0.0d0, 0.0d0)

        ! alpha = Larmor_skip_factor
        ! tau   = max(kernel_taper_skip_threshold, 1.0d-12)
        ! block
            ! real(dp) :: rhoT_max
            ! integer :: sigma
            ! rhoT_max = 0.0d0
            ! do sigma = 0, plasma%n_species - 1
                ! if (allocated(plasma%spec(sigma)%rho_L_cc)) then
                    ! rhoT_max = max(rhoT_max, maxval(plasma%spec(sigma)%rho_L_cc))
                ! end if
            ! end do
            ! dmax_global = alpha * rhoT_max * sqrt(max(log(1.0d0 / tau), 0.0d0))
        ! end block

        ! total_iterations = 0
        ! block
            ! real(dp) :: xl_val
            ! integer :: lp_lo, lp_hi
            ! do l = 1, K_rho_phi_llp%npts_l

                ! xl_val = xl_grid%xb(l)

                ! lp_lo = l
                ! do
                    ! if (lp_lo <= 1) exit
                    ! if (abs(xl_grid%xb(lp_lo - 1) - xl_val) > dmax_global) exit
                    ! lp_lo = lp_lo - 1
                ! end do

                ! lp_hi = l
                ! do
                    ! if (lp_hi >= K_rho_phi_llp%npts_l) exit
                    ! if (abs(xl_grid%xb(lp_hi + 1) - xl_val) > dmax_global) exit
                    ! lp_hi = lp_hi + 1
                ! end do

                ! if (min(l, lp_hi) >= max(1, lp_lo)) then
                    ! total_iterations = total_iterations + (min(l, lp_hi) - max(1, lp_lo) + 1)
                ! end if
            ! end do
        ! end block

        ! if (fstatus == 1) write(*,*) 'FLR2 benchmark total band-limited iterations: ', total_iterations
        ! current_iteration = 0
        ! call system_clock(start_count, count_rate, count_max)

        ! if (include_electrons .and. .not. include_ions) then
            ! call fill_case(FP_calc_kernel_electrons_FLR2_benchmark)
        ! else if (include_ions .and. .not. include_electrons) then
            ! print *, 'FLR2 benchmark: ions only case is not implemented anymore. Need to adapt interface.'
            ! ! call fill_case(FP_calc_kernel_ions_FLR2_benchmark)
        ! else if (include_electrons .and. include_ions) then
            ! call fill_case(FP_calc_kernel_FLR2_benchmark)
        ! else
            ! write(*,*) 'FLR2 benchmark: both species disabled; kernels remain zero.'
        ! end if

        ! write(*,*)
        ! write(*,*) 'Finished filling kernels.'

    ! contains

        ! subroutine fill_case(calc_kernel)

            ! procedure(flr2_kernel_cb) :: calc_kernel
            ! integer :: l, lp, lp_lo, lp_hi
            ! real(dp) :: xl_val
            ! complex(dp) :: k_rho_phi, k_rho_B, k_j_phi, k_j_B

            ! !$omp parallel do default(shared) private(l, lp, lp_lo, lp_hi, xl_val, k_rho_phi, k_rho_B, k_j_phi, k_j_B) schedule(dynamic)
            ! do l = 1, K_rho_phi_llp%npts_l
                ! xl_val = xl_grid%xb(l)

                ! lp_lo = l
                ! do
                    ! if (lp_lo <= 2) exit
                    ! if (abs(xl_grid%xb(lp_lo - 1) - xl_val) > dmax_global) exit
                    ! lp_lo = lp_lo - 1
                ! end do

                ! lp_hi = l
                ! do
                    ! if (lp_hi >= K_rho_phi_llp%npts_l) exit
                    ! if (abs(xl_grid%xb(lp_hi + 1) - xl_val) > dmax_global) exit
                    ! lp_hi = lp_hi + 1
                ! end do

                ! do lp = max(1, lp_lo), min(l, lp_hi)
                    ! call calc_kernel(l, lp, k_rho_phi, k_rho_B, k_j_phi, k_j_B, gauss_conf)

                    ! K_rho_phi_llp%Kllp(l, lp) = k_rho_phi
                    ! K_rho_B_llp%Kllp(l, lp)   = k_rho_B
                    ! K_j_phi_llp%Kllp(l, lp)   = k_j_phi
                    ! K_j_B_llp%Kllp(l, lp)     = k_j_B

                    ! if (lp /= l) then
                        ! K_rho_phi_llp%Kllp(lp, l) = k_rho_phi
                        ! K_rho_B_llp%Kllp(lp, l)   = k_rho_B
                        ! K_j_phi_llp%Kllp(lp, l)   = k_j_phi
                        ! K_j_B_llp%Kllp(lp, l)     = k_j_B
                    ! end if

                    ! if (total_iterations > 0) then
                        ! !$omp atomic
                        ! current_iteration = current_iteration + 1
                        ! !$omp critical(loading_bar)
                        ! if (mod(current_iteration, 32) == 0 .or. current_iteration == total_iterations) then
                            ! call update_bar(current_iteration, total_iterations, start_count, count_rate)
                        ! end if
                        ! !$omp end critical(loading_bar)
                    ! end if
                ! end do
            ! end do
            ! !$omp end parallel do

        ! end subroutine fill_case

        ! subroutine rescale_susceptibility_functions()
            ! ! Rescale susceptibility functions for FLR2 benchmark
            ! ! has to be done in order to take quasineutrality into account as is done in FLR2
            
            ! use species_m, only: plasma
            ! use IO_collection_m, only: write_complex_profile
            ! use grid_m, only: rg_grid
            ! use config_m, only: output_path

            ! implicit none

            ! integer :: sp

            ! do sp=0, plasma%n_species - 1
                ! ! Note: I00, I20, x1, x2 are now only available as _cc (cell center) versions
                ! ! They are computed after interpolation in calculate_thermodynamic_forces_and_susc
                ! plasma%spec(sp)%I00_cc = plasma%spec(sp)%I01_cc * plasma%spec(sp)%x1_cc / plasma%spec(sp)%x2_cc
                ! plasma%spec(sp)%I20_cc = plasma%spec(sp)%I21_cc * plasma%spec(sp)%x1_cc / plasma%spec(sp)%x2_cc

                ! call write_complex_profile(rg_grid%xc, plasma%spec(sp)%I00_cc, rg_grid%npts_c, &
                    ! 'backs/'//trim(plasma%spec(sp)%name)//'/I00_cc_after_resc.dat')
                ! call write_complex_profile(rg_grid%xc, plasma%spec(sp)%I20_cc, rg_grid%npts_c, &
                    ! 'backs/'//trim(plasma%spec(sp)%name)//'/I20_cc_after_resc.dat')
            ! end do

            
        ! end subroutine

    ! end subroutine FP_fill_kernels_flr2_benchmark_old


    ! subroutine FP_calc_kernel_FLR2_benchmark(l, lp, k_rho_phi, k_rho_B, k_j_phi, k_j_B, gauss_conf)
        ! ! for benchmarking FLR2 and KIM against each other, FLR2 exploits quasineutrality which has to be taken into account for
        ! ! electrons and ions separately, most notably, by omitting the Debye shielding term and using different susceptibility functions

        ! use KIM_kinds_m, only: dp
        ! use integrals_gauss_m, only: gauss_integrate_F0, gauss_integrate_F1, gauss_integrate_F2, gauss_integrate_F3,&
            ! gauss_config_t, gauss_integrate_F1_electrons
        ! use species_m, only: plasma
        ! use constants_m, only: pi
        ! use integrands_gauss_m, only: gauss_int_F0_rho_phi_t, gauss_int_F1_rho_phi_t, gauss_int_F2_rho_phi_t, gauss_int_F3_rho_phi_t, &
            ! integration_point_t, gauss_int_F1_rho_phi_electrons_t
        ! use FP_kernel_plasma_prefacs_m, only: FP_G0_rho_phi
        ! use grid_m, only: Larmor_skip_factor, rg_grid
        ! use constants_m, only: com_unit, sol
        ! use config_m, only: turn_off_ions, turn_off_electrons, artificial_debye_case
        ! use grid_m, only: xl_grid

        ! implicit none

        ! integer, intent(in) :: l, lp
        ! complex(dp), intent(inout) :: k_rho_phi, k_rho_B, k_j_phi, k_j_B
        ! integer :: j, sigma
        ! type(gauss_config_t), intent(in) :: gauss_conf
        ! real(dp) :: integral_val
        ! real(dp) :: current_distance

        ! type(integration_point_t) :: int_point
        ! type(gauss_int_F0_rho_phi_t) :: int_F0
        ! type(gauss_int_F1_rho_phi_t) :: int_F1
        ! type(gauss_int_F2_rho_phi_t) :: int_F2
        ! type(gauss_int_F3_rho_phi_t) :: int_F3
        ! type(gauss_int_F1_rho_phi_electrons_t) :: int_F1_e

        ! k_rho_phi = (0.0d0, 0.0d0)
        ! k_rho_B = (0.0d0, 0.0d0)
        ! k_j_phi = (0.0d0, 0.0d0)
        ! k_j_B = (0.0d0, 0.0d0)

        ! call set_xl_at_edge(l, lp, int_point)

        ! do sigma = 0, plasma%n_species - 1
            ! do j = 1, rg_grid%npts_b-1

                ! int_point%j = j
                ! int_point%rhoT = max(plasma%spec(sigma)%rho_L_cc(j), 0.0d0)

                ! ! skip term if species Larmor radius is too small to couple these grid points
                ! if (abs(l-lp) > 5 .and. abs(xl_grid%xb(l) - xl_grid%xb(lp))> Larmor_skip_factor * plasma%spec(sigma)%rho_L(j)) cycle
                ! if (abs(0.5d0 * (rg_grid%xb(j+1) + rg_grid%xb(j)) - 0.5d0 * (xl_grid%xb(l) + xl_grid%xb(lp))) &
                    ! > Larmor_skip_factor * plasma%spec(sigma)%rho_L(j)) cycle

                ! if (sigma == 0) then
                    ! int_F1_e%int_point = int_point

                    ! ! F1 integration
                    ! integral_val = 0.0d0
                    ! call gauss_integrate_F1_electrons(int_F1_e, integral_val, gauss_conf)
                    ! integral_val = integral_val * pi ! pi because of missing Bessel function representation
                    ! k_rho_phi = k_rho_phi + integral_val * pref_rho_phi_g1(1,j)
                    ! k_rho_B = k_rho_B + integral_val * pref_rho_B_g1(1,j)
                    ! k_j_phi = k_j_phi + integral_val * pref_j_phi_g1(1,j)
                    ! k_j_B = k_j_B + integral_val * pref_j_B_g1(1,j)
                    ! cycle
                ! end if

                ! int_F1%int_point = int_point
                ! int_F2%int_point = int_point
                ! int_F3%int_point = int_point

                ! ! F1 integration
                ! call gauss_integrate_F1(int_F1, integral_val, gauss_conf)
                ! k_rho_phi = k_rho_phi + integral_val * pref_rho_phi_g1(sigma+1,j)
                ! k_rho_B = k_rho_B + integral_val * pref_rho_B_g1(sigma+1,j)
                ! k_j_phi = k_j_phi + integral_val * pref_j_phi_g1(sigma+1,j)
                ! k_j_B = k_j_B + integral_val * pref_j_B_g1(sigma+1,j)

                ! ! F2 integration
                ! call gauss_integrate_F2(int_F2, integral_val, gauss_conf)
                ! k_rho_phi = k_rho_phi + integral_val * pref_rho_phi_g2(sigma+1,j)
                ! k_rho_B = k_rho_B + integral_val * pref_rho_B_g2(sigma+1,j)
                ! k_j_phi = k_j_phi + integral_val * pref_j_phi_g2(sigma+1,j)
                ! k_j_B = k_j_B + integral_val * pref_j_B_g2(sigma+1,j)

                ! ! F3 integration
                ! call gauss_integrate_F3(int_F3, integral_val, gauss_conf)
                ! k_rho_phi = k_rho_phi + integral_val * pref_rho_phi_g3(sigma+1,j)
                ! k_rho_B = k_rho_B + integral_val * pref_rho_B_g3(sigma+1,j)
                ! k_j_phi = k_j_phi + integral_val * pref_j_phi_g3(sigma+1,j)
                ! k_j_B = k_j_B + integral_val * pref_j_B_g3(sigma+1,j)

            ! end do
        ! end do

        ! k_rho_phi = k_rho_phi / (8.0d0 * pi**3.0d0)
        ! k_rho_B = k_rho_B / (8.0d0 * pi**3.0d0)

        ! k_j_phi = k_j_phi / (8.0d0 * pi**3.0d0)
        ! k_j_B = k_j_B / (8.0d0 * pi**3.0d0)

    ! end subroutine FP_calc_kernel_FLR2_benchmark




    ! subroutine FP_calc_kernel_ions_FLR2_benchmark_old(l, lp, k_rho_phi, k_rho_B, k_j_phi, k_j_B, gauss_conf, sigma)
        ! ! for benchmarking FLR2 and KIM against each other, FLR2 exploits quasineutrality which has to be taken into account for
        ! ! electrons and ions separately, most notably, by omitting the Debye shielding term and using different susceptibility functions

        ! use KIM_kinds_m, only: dp
        ! use integrals_gauss_m, only: gauss_integrate_F0, gauss_integrate_F1, gauss_integrate_F2, gauss_integrate_F3,&
            ! gauss_config_t
        ! use species_m, only: plasma
        ! use constants_m, only: pi
        ! use integrands_gauss_m, only: gauss_int_F0_rho_phi_t, gauss_int_F1_rho_phi_t, gauss_int_F2_rho_phi_t, gauss_int_F3_rho_phi_t, &
            ! integration_point_t
        ! use FP_kernel_plasma_prefacs_m, only: FP_G0_rho_phi
        ! use grid_m, only: Larmor_skip_factor, kernel_taper_skip_threshold, rg_grid
        ! use constants_m, only: com_unit, sol
        ! use config_m, only: turn_off_ions, turn_off_electrons, artificial_debye_case
        ! use grid_m, only: xl_grid
        ! use resonances_mod, only: r_res

        ! implicit none

        ! integer, intent(in) :: l, lp, sigma
        ! complex(dp), intent(inout) :: k_rho_phi, k_rho_B, k_j_phi, k_j_B
        ! integer :: j
        ! type(gauss_config_t), intent(in) :: gauss_conf
        ! real(dp) :: integral_val
        ! real(dp) :: current_distance

        ! type(integration_point_t) :: int_point
        ! type(gauss_int_F0_rho_phi_t) :: int_F0
        ! type(gauss_int_F1_rho_phi_t) :: int_F1
        ! type(gauss_int_F2_rho_phi_t) :: int_F2
        ! type(gauss_int_F3_rho_phi_t) :: int_F3

        ! k_rho_phi = (0.0d0, 0.0d0)
        ! k_rho_B = (0.0d0, 0.0d0)
        ! k_j_phi = (0.0d0, 0.0d0)
        ! k_j_B = (0.0d0, 0.0d0)

        ! call set_xl_at_edge(l, lp, int_point)

        ! do j = 1, rg_grid%npts_b-1

            ! int_point%j = j
            ! int_point%rhoT = max(plasma%spec(sigma)%rho_L_cc(j), 0.0d0)

            ! ! skip term if species Larmor radius is too small to couple these grid points
            ! if (abs(l-lp) > 10 .and. abs(xl_grid%xb(l) - xl_grid%xb(lp))> Larmor_skip_factor * plasma%spec(sigma)%rho_L(j)) cycle
            ! if (abs(0.5d0 * (rg_grid%xb(j+1) + rg_grid%xb(j)) - 0.5d0 * (xl_grid%xb(l) + xl_grid%xb(lp))) &
                ! > Larmor_skip_factor * plasma%spec(sigma)%rho_L(j)) cycle
            
            ! int_F1%int_point = int_point
            ! int_F2%int_point = int_point
            ! int_F3%int_point = int_point

            ! ! F1 integration
            ! call gauss_integrate_F1(int_F1, integral_val, gauss_conf)
            ! k_rho_phi = k_rho_phi + integral_val * pref_rho_phi_g1(sigma+1,j)
            ! k_rho_B = k_rho_B + integral_val * pref_rho_B_g1(sigma+1,j)
            ! k_j_phi = k_j_phi + integral_val * pref_j_phi_g1(sigma+1,j)
            ! k_j_B = k_j_B + integral_val * pref_j_B_g1(sigma+1,j)

            ! ! ignore FLR terms if resonance is too far from grid points
            ! if (abs(xl_grid%xb(l) - r_res) > 10.0d0 * int_point%rhoT .or. &
                ! abs(xl_grid%xb(lp) - r_res) > 10.0d0 * int_point%rhoT) then
                ! cycle
            ! end if

            ! ! F2 integration
            ! call gauss_integrate_F2(int_F2, integral_val, gauss_conf)
            ! k_rho_phi = k_rho_phi + integral_val * pref_rho_phi_g2(sigma+1,j)
            ! k_rho_B = k_rho_B + integral_val * pref_rho_B_g2(sigma+1,j)
            ! k_j_phi = k_j_phi + integral_val * pref_j_phi_g2(sigma+1,j)
            ! k_j_B = k_j_B + integral_val * pref_j_B_g2(sigma+1,j)

            ! ! F3 integration
            ! call gauss_integrate_F3(int_F3, integral_val, gauss_conf)
            ! k_rho_phi = k_rho_phi + integral_val * pref_rho_phi_g3(sigma+1,j)
            ! k_rho_B = k_rho_B + integral_val * pref_rho_B_g3(sigma+1,j)
            ! k_j_phi = k_j_phi + integral_val * pref_j_phi_g3(sigma+1,j)
            ! k_j_B = k_j_B + integral_val * pref_j_B_g3(sigma+1,j)

        ! end do

        ! k_rho_phi = k_rho_phi / (8.0d0 * pi**3.0d0)
        ! k_rho_B = k_rho_B / (8.0d0 * pi**3.0d0)

        ! k_j_phi = k_j_phi / (8.0d0 * pi**3.0d0)
        ! k_j_B = k_j_B / (8.0d0 * pi**3.0d0)

    ! end subroutine FP_calc_kernel_ions_FLR2_benchmark_old

end module
