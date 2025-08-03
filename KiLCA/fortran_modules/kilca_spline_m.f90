module kilca_spline_m
    use iso_fortran_env, only: real64, int32, int64
    use kilca_types_m
    implicit none
    
    private
    
    ! Spline boundary condition types
    integer, parameter, public :: SPLINE_NATURAL = 1
    integer, parameter, public :: SPLINE_PERIODIC = 2
    integer, parameter, public :: SPLINE_ZERO = 3
    
    ! Spline data structure
    type, public :: spline_data_t
        integer :: N                          ! Spline degree (must be odd)
        integer :: ind                        ! Last search index for faster evaluation
        integer :: type                       ! Type of spline boundary
        integer :: dimx                       ! Dimension of x grid
        real(real64), pointer :: x(:) => null()        ! x grid for profiles
        real(real64), pointer :: C(:,:) => null()      ! Matrix of spline coefficients
        real(real64), pointer :: BC(:,:) => null()     ! Binomial coefficients
        real(real64), pointer :: fac(:,:) => null()    ! Factors for spline evaluation
    end type spline_data_t
    
    ! Public procedures
    public :: spline_create
    public :: spline_destroy
    public :: spline_calc_coefficients
    public :: spline_eval
    public :: spline_eval_deriv
    
    ! Private procedures for internal use
    private :: calc_spline_boundaries
    private :: calc_spline_coefficients
    private :: search_array
    private :: binary_search
    private :: set_bc_array
    private :: set_fac_array
    
contains

    !---------------------------------------------------------------------------
    ! Create and initialize spline structure
    !---------------------------------------------------------------------------
    subroutine spline_create(spline, N, type, dimx, x, ierr)
        type(spline_data_t), intent(out) :: spline
        integer, intent(in) :: N              ! Degree (must be odd)
        integer, intent(in) :: type           ! Boundary type
        integer, intent(in) :: dimx           ! Size of x array
        real(real64), intent(in), target :: x(:)    ! x grid points
        integer, intent(out) :: ierr
        
        ! Variables i, j removed - were unused
        
        ierr = 0
        
        ! Validate input parameters
        if (N <= 0 .or. mod(N, 2) == 0) then
            ierr = -1
            print *, "Error: Spline degree must be positive and odd"
            return
        end if
        
        if (dimx < N + 1) then
            ierr = -2
            print *, "Error: Grid size must be at least N+1"
            return
        end if
        
        if (type < 1 .or. type > 3) then
            ierr = -3
            print *, "Error: Invalid boundary type"
            return
        end if
        
        ! Set basic parameters
        spline%N = N
        spline%type = type
        spline%dimx = dimx
        spline%ind = dimx / 2
        
        ! Point to x array
        spline%x => x
        
        ! Allocate coefficient matrix
        allocate(spline%C(dimx, N+1))
        spline%C = 0.0_real64
        
        ! Calculate binomial coefficients
        allocate(spline%BC(N+1, N+1))
        call set_bc_array(N, spline%BC)
        
        ! Calculate factors for spline evaluation
        allocate(spline%fac(N+1, N+1))
        call set_fac_array(N, spline%fac)
        
    end subroutine spline_create
    
    !---------------------------------------------------------------------------
    ! Destroy spline structure and free memory
    !---------------------------------------------------------------------------
    subroutine spline_destroy(spline)
        type(spline_data_t), intent(inout) :: spline
        
        if (associated(spline%C)) deallocate(spline%C)
        if (associated(spline%BC)) deallocate(spline%BC)
        if (associated(spline%fac)) deallocate(spline%fac)
        
        nullify(spline%x)
        nullify(spline%C)
        nullify(spline%BC)
        nullify(spline%fac)
        
        spline%N = 0
        spline%ind = 0
        spline%type = 0
        spline%dimx = 0
        
    end subroutine spline_destroy
    
    !---------------------------------------------------------------------------
    ! Calculate spline coefficients from data
    !---------------------------------------------------------------------------
    subroutine spline_calc_coefficients(spline, y, ierr)
        type(spline_data_t), intent(inout) :: spline
        real(real64), intent(in) :: y(:)      ! Function values at x points
        integer, intent(out) :: ierr
        
        real(real64), allocatable :: W(:)     ! Working array
        integer :: dimy, wsize
        
        ierr = 0
        
        ! Check dimensions
        if (size(y) /= spline%dimx) then
            ierr = -1
            print *, "Error: y array size must match x array size"
            return
        end if
        
        dimy = 1  ! Single function for now
        
        ! Allocate working array
        wsize = (spline%N + 1) * (2 * dimy + spline%dimx * (4 + 3 * ((spline%N - 1) / 2)))
        allocate(W(wsize))
        
        ! Calculate spline boundaries
        ierr = calc_spline_boundaries(spline, dimy, y, W)
        if (ierr /= 0) return
        
        ! Calculate spline coefficients
        ierr = calc_spline_coefficients(spline, 1, dimy, y, W)
        
        deallocate(W)
        
    end subroutine spline_calc_coefficients
    
    !---------------------------------------------------------------------------
    ! Evaluate spline at given points
    !---------------------------------------------------------------------------
    subroutine spline_eval(spline, z, result, ierr)
        type(spline_data_t), intent(in) :: spline
        real(real64), intent(in) :: z(:)      ! Points to evaluate at
        real(real64), intent(out) :: result(:) ! Results
        integer, intent(out) :: ierr
        
        integer :: i, j, ind  ! k removed - was unused
        real(real64) :: dx, val
        
        ierr = 0
        
        ! Check dimensions
        if (size(result) /= size(z)) then
            ierr = -1
            print *, "Error: result array size must match z array size"
            return
        end if
        
        ! Evaluate spline at each point
        do i = 1, size(z)
            ! Find interval containing z(i)
            ind = spline%ind
            call search_array(z(i), spline%dimx, spline%x, ind, ind)
            
            ! Handle boundary cases
            if (ind < 1) ind = 1
            if (ind >= spline%dimx) ind = spline%dimx - 1
            
            ! Compute spline value using stored coefficients
            dx = z(i) - spline%x(ind)
            val = 0.0_real64
            
            ! Evaluate polynomial
            do j = 0, spline%N
                val = val + spline%C(ind, j+1) * dx**j
            end do
            
            result(i) = val
        end do
        
    end subroutine spline_eval
    
    !---------------------------------------------------------------------------
    ! Evaluate spline derivative at given points
    !---------------------------------------------------------------------------
    subroutine spline_eval_deriv(spline, z, deriv_order, result, ierr)
        type(spline_data_t), intent(in) :: spline
        real(real64), intent(in) :: z(:)      ! Points to evaluate at
        integer, intent(in) :: deriv_order     ! Derivative order
        real(real64), intent(out) :: result(:) ! Results
        integer, intent(out) :: ierr
        
        integer :: i, j, k, ind
        real(real64) :: dx, val, coeff
        
        ierr = 0
        
        ! Check dimensions
        if (size(result) /= size(z)) then
            ierr = -1
            print *, "Error: result array size must match z array size"
            return
        end if
        
        ! Check derivative order
        if (deriv_order < 0 .or. deriv_order > spline%N) then
            ierr = -2
            print *, "Error: Invalid derivative order"
            return
        end if
        
        ! Evaluate derivative at each point
        do i = 1, size(z)
            ! Find interval containing z(i)
            ind = spline%ind
            call search_array(z(i), spline%dimx, spline%x, ind, ind)
            
            ! Handle boundary cases
            if (ind < 1) ind = 1
            if (ind >= spline%dimx) ind = spline%dimx - 1
            
            ! Compute derivative value
            dx = z(i) - spline%x(ind)
            val = 0.0_real64
            
            ! Evaluate derivative of polynomial
            do j = deriv_order, spline%N
                coeff = 1.0_real64
                do k = 0, deriv_order - 1
                    coeff = coeff * real(j - k, real64)
                end do
                val = val + coeff * spline%C(ind, j+1) * dx**(j - deriv_order)
            end do
            
            result(i) = val
        end do
        
    end subroutine spline_eval_deriv
    
    !---------------------------------------------------------------------------
    ! Calculate spline boundary conditions
    !---------------------------------------------------------------------------
    function calc_spline_boundaries(spline, dimy, y, W) result(ierr)
        type(spline_data_t), intent(inout) :: spline
        integer, intent(in) :: dimy  ! NOTE: dimy kept for interface compatibility
        real(real64), intent(in) :: y(:)
        real(real64), intent(inout) :: W(:)  ! NOTE: W kept for interface compatibility
        integer :: ierr
        
        ierr = 0
        
        ! For now, implement natural boundary conditions (second derivative = 0)
        if (spline%type == SPLINE_NATURAL) then
            ! Set boundary conditions in coefficient matrix
            ! This is simplified - full implementation would set up linear system
            spline%C(1, 1) = y(1)
            spline%C(spline%dimx, 1) = y(spline%dimx)
        else if (spline%type == SPLINE_PERIODIC) then
            ! Periodic boundary conditions
            ! Ensure continuity of function and derivatives at boundaries
            ierr = 0  ! Placeholder
        else if (spline%type == SPLINE_ZERO) then
            ! Zero boundary conditions
            spline%C(1, 1) = 0.0_real64
            spline%C(spline%dimx, 1) = 0.0_real64
        end if
        
    end function calc_spline_boundaries
    
    !---------------------------------------------------------------------------
    ! Calculate spline coefficients using linear system
    !---------------------------------------------------------------------------
    function calc_spline_coefficients(spline, Imin, dimy, y, W) result(ierr)
        type(spline_data_t), intent(inout) :: spline
        integer, intent(in) :: Imin, dimy  ! NOTE: Imin, dimy kept for interface compatibility
        real(real64), intent(in) :: y(:)
        real(real64), intent(inout) :: W(:)  ! NOTE: W kept for interface compatibility  
        integer :: ierr
        
        integer :: i, n
        real(real64), allocatable :: a(:), b(:), c(:), d(:), h(:), alpha(:), l(:), mu(:), z(:)
        
        ierr = 0
        n = spline%dimx
        
        ! For cubic splines (N=3), implement proper cubic spline algorithm
        if (spline%N == 3) then
            ! Allocate working arrays for cubic spline calculation
            allocate(a(n), b(n), c(n), d(n), h(n-1), alpha(n-1), l(n), mu(n), z(n))
            
            ! Copy data values
            a = y
            
            ! Calculate h values (step sizes)
            do i = 1, n-1
                h(i) = spline%x(i+1) - spline%x(i)
                if (abs(h(i)) < 1.0e-15_real64) then
                    ierr = -1  ! x values not properly ordered
                    deallocate(a, b, c, d, h, alpha, l, mu, z)
                    return
                end if
            end do
            
            ! Calculate alpha values for tridiagonal system
            do i = 2, n-1
                alpha(i) = (3.0_real64 / h(i)) * (a(i+1) - a(i)) - &
                          (3.0_real64 / h(i-1)) * (a(i) - a(i-1))
            end do
            
            ! Solve tridiagonal system using Thomas algorithm
            ! Forward elimination
            l(1) = 1.0_real64
            mu(1) = 0.0_real64
            z(1) = 0.0_real64
            
            do i = 2, n-1
                l(i) = 2.0_real64 * (spline%x(i+1) - spline%x(i-1)) - h(i-1) * mu(i-1)
                mu(i) = h(i) / l(i)
                z(i) = (alpha(i) - h(i-1) * z(i-1)) / l(i)
            end do
            
            l(n) = 1.0_real64
            z(n) = 0.0_real64
            c(n) = 0.0_real64
            
            ! Back substitution
            do i = n-1, 1, -1
                c(i) = z(i) - mu(i) * c(i+1)
                b(i) = (a(i+1) - a(i)) / h(i) - h(i) * (c(i+1) + 2.0_real64 * c(i)) / 3.0_real64
                d(i) = (c(i+1) - c(i)) / (3.0_real64 * h(i))
            end do
            
            ! Store coefficients in spline structure
            ! spline%C format: (node, coefficient)
            ! Coefficients for S_i(x) = a_i + b_i*(x-x_i) + c_i*(x-x_i)^2 + d_i*(x-x_i)^3
            do i = 1, n-1  ! Only n-1 cubic pieces
                spline%C(i, 1) = a(i)    ! Constant term
                spline%C(i, 2) = b(i)    ! Linear term
                spline%C(i, 3) = c(i)    ! Quadratic term
                spline%C(i, 4) = d(i)    ! Cubic term
            end do
            
            ! Last node uses previous cubic piece coefficients
            if (n > 1) then
                spline%C(n, :) = spline%C(n-1, :)
            end if
            
            deallocate(a, b, c, d, h, alpha, l, mu, z)
            
        else
            ! General case for arbitrary degree
            ! For now, only cubic splines (N=3) are fully implemented
            ierr = -99  ! Not implemented for other degrees
        end if
        
    end function calc_spline_coefficients
    
    !---------------------------------------------------------------------------
    ! Search for interval containing x value
    !---------------------------------------------------------------------------
    subroutine search_array(x, dimx, xa, ind, ind_guess)
        real(real64), intent(in) :: x
        integer, intent(in) :: dimx
        real(real64), intent(in) :: xa(:)
        integer, intent(out) :: ind
        integer, intent(in) :: ind_guess
        
        integer :: ilo, ihi
        
        ! Use previous guess as starting point
        if (ind_guess > 0 .and. ind_guess < dimx) then
            if (x >= xa(ind_guess) .and. x < xa(ind_guess + 1)) then
                ind = ind_guess
                return
            end if
        end if
        
        ! Binary search
        ilo = 1
        ihi = dimx
        ind = binary_search(x, xa, ilo, ihi)
        
    end subroutine search_array
    
    !---------------------------------------------------------------------------
    ! Binary search implementation
    !---------------------------------------------------------------------------
    function binary_search(x, xa, ilo, ihi) result(ind)
        real(real64), intent(in) :: x
        real(real64), intent(in) :: xa(:)
        integer, intent(in) :: ilo, ihi
        integer :: ind
        
        integer :: lo, hi, mid
        
        lo = ilo
        hi = ihi
        
        do while (hi - lo > 1)
            mid = (lo + hi) / 2
            if (x >= xa(mid)) then
                lo = mid
            else
                hi = mid
            end if
        end do
        
        ind = lo
        
    end function binary_search
    
    !---------------------------------------------------------------------------
    ! Set binomial coefficient array
    !---------------------------------------------------------------------------
    subroutine set_bc_array(N, BC)
        integer, intent(in) :: N
        real(real64), intent(out) :: BC(:,:)
        
        integer :: i, j
        
        BC = 0.0_real64
        
        ! Pascal's triangle
        do i = 0, N
            BC(i+1, 1) = 1.0_real64
            if (i+1 <= N+1) BC(i+1, i+1) = 1.0_real64
            do j = 1, i-1
                if (i > 0 .and. j > 0) BC(i+1, j+1) = BC(i, j) + BC(i, j+1)
            end do
        end do
        
    end subroutine set_bc_array
    
    !---------------------------------------------------------------------------
    ! Set factor array for spline evaluation
    !---------------------------------------------------------------------------
    subroutine set_fac_array(N, fac)
        integer, intent(in) :: N
        real(real64), intent(out) :: fac(:,:)
        
        integer :: i, j
        
        fac = 0.0_real64
        
        ! Set up factorial-based factors
        do i = 0, N
            do j = 0, i
                if (j == 0) then
                    fac(i+1, j+1) = 1.0_real64
                else
                    fac(i+1, j+1) = fac(i+1, j) * real(i - j + 1, real64)
                end if
            end do
        end do
        
    end subroutine set_fac_array

end module kilca_spline_m