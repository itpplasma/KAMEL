module kernel

    implicit none

    !double complex, dimension(:,:), allocatable :: K_rho_B
    !double complex, dimension(:,:), allocatable :: K_rho_phi
    double complex, dimension(:,:,:), allocatable :: K_rho_phi_of_rg
    double complex, dimension(:,:,:), allocatable :: K_rho_B_of_rg
    !double complex, dimension(:,:,:), allocatable :: K_j_phi_of_rg
    !double complex, dimension(:,:,:), allocatable :: K_j_B_of_rg

    double complex, dimension(:,:), allocatable :: K_rho_phi_llp
    double complex, dimension(:,:), allocatable :: K_rho_B_llp

    double complex, dimension(:,:), allocatable :: K_j_phi_llp
    double complex, dimension(:,:), allocatable :: K_j_B_llp

    contains
        subroutine fill_all_kernels

            use kernel_functions, only: kernel_rho_phi_of_kr_krp_rg, kernel_rho_B_of_kr_krp_rg, kernel_j_B_of_kr_krp_rg, &
                kernel_j_phi_of_kr_krp_rg
            !use kernel_functions, only: K_rho_phi_of_rg, K_rho_B_of_rg, K_j_phi_of_rg, K_j_B_of_rg
            use config, only: fstatus
            use loading_bar
            use grid, only: l_space_dim, r_space_dim, varphi_lkr, npoib, rb
            use kr_grid, only: k_space_dim, kr, krp

            implicit none
            integer :: i_kr, i_krp, i_rg
            integer :: count_loading = 0
            integer :: total_count_loading = 0

            if (fstatus == 1) write(*,*) 'Status: Fill kernel rho phi'
            if (.not. allocated(K_rho_phi_of_rg)) allocate(K_rho_phi_of_rg(k_space_dim, k_space_dim, npoib))
            if (.not. allocated(K_rho_B_of_rg)) allocate(K_rho_B_of_rg(k_space_dim, k_space_dim, npoib))
            !if (.not. allocated(K_j_phi_of_rg)) allocate(K_j_phi_of_rg(k_space_dim, k_space_dim, npoib))
            !if (.not. allocated(K_j_B_of_rg)) allocate(K_j_B_of_rg(k_space_dim, k_space_dim, npoib))

            total_count_loading = k_space_dim * k_space_dim * npoib
        
            !$OMP PARALLEL DO collapse(3) default(none) schedule(guided) &
            !$OMP PRIVATE(i_krp, i_kr, i_rg, total_count_loading) &
            !$OMP SHARED(K_rho_phi_of_rg, kr, K_rho_B_of_rg, &
            !!$OMP K_j_phi_of_rg, K_j_B_of_rg, &
            !$OMP krp, rb, k_space_dim, npoib, count_loading)
            do i_krp = 1, k_space_dim
                do i_kr = 1, k_space_dim
                    do i_rg = 1, npoib
                        K_rho_phi_of_rg(i_krp, i_kr, i_rg) = kernel_rho_phi_of_kr_krp_rg(kr(i_kr), krp(i_krp), rb(i_rg)) 
                        K_rho_B_of_rg(i_krp, i_kr, i_rg) = kernel_rho_B_of_kr_krp_rg(kr(i_kr), krp(i_krp), rb(i_rg)) 
                        !K_j_phi_of_rg(i_krp, i_kr, i_rg) = kernel_j_phi_of_kr_krp_rg(kr(i_kr), krp(i_krp), rb(i_rg)) 
                        !K_j_B_of_rg(i_krp, i_kr, i_rg) = kernel_j_B_of_kr_krp_rg(kr(i_kr), krp(i_krp), rb(i_rg)) 
                    end do
                end do
            end do
            !$OMP END PARALLEL DO
            if (fstatus == 1) write(*,*) 'Status: Finished filling kernel rho phi'

        end subroutine fill_all_kernels


end module