module kernel

    implicit none

    double complex, dimension(:,:), allocatable :: K_rho_B
    double complex, dimension(:,:), allocatable :: K_rho_phi
    double complex, dimension(:,:,:), allocatable :: K_rho_phi_of_rg
    double complex, dimension(:,:), allocatable :: K_j_B
    double complex, dimension(:,:), allocatable :: K_j_phi

    double complex, dimension(:,:), allocatable :: K_rho_phi_llp
    double complex, dimension(:,:), allocatable :: K_rho_B_llp

end module