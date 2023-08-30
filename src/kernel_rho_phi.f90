subroutine kernel_rho_phi
! fill kernel correlating rho and phi

    use plas_parameter
    use constants
    use config
    use grid
    use back_quants

    implicit none

    double complex :: besselj ! complex bessel function from bessel.f90

    integer :: npoi_der, nder
    double precision, dimension(:,:), allocatable :: coef
    double precision, dimension(:), allocatable :: x

    double precision, dimension(:), allocatable :: kr ! radial wavenumber
    double precision, dimension(:), allocatable :: krp ! radial wavenumber prime

    integer :: sigma ! for loop over species

    if (fstatus == 1) write(*,*) 'Status: Generating kernel rho phi'

    npoi_der = 4 ! number of polynomials for derivative
    nder = 1 ! first derivative

    allocate(x(npoi_der), coef(0:nder, npoi_der))
    allocate(kr(k_space_dim), krp(k_space_dim))



end subroutine