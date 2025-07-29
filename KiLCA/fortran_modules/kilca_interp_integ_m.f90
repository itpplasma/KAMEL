module kilca_interp_integ_m
    use iso_fortran_env, only: real64, error_unit
    implicit none
    
    private
    
    ! Public procedures for interpolation
    public :: eval_neville_polynom
    public :: eval_lagrange_polynom
    public :: eval_complex_interp
    public :: search_array
    public :: binary_search
    public :: find_index_for_interp
    
    ! Public procedures for integration  
    public :: integrate_trapezoid
    public :: integrate_simpson
    public :: integrate_adaptive_quadrature
    public :: integrate_over_cylinder
    public :: integrate_gauss_legendre
    
    ! Public constants
    integer, parameter, public :: MAX_SUBDIVISIONS = 1000
    real(real64), parameter, public :: DEFAULT_TOLERANCE = 1.0e-8_real64
    
contains

    !---------------------------------------------------------------------------
    ! Neville polynomial interpolation with derivatives
    ! Translates eval_neville_polynom from interp.h
    !---------------------------------------------------------------------------
    subroutine eval_neville_polynom(dim, xg, yg, deg, x, Dmin, Dmax, ind, R, ierr)
        integer, intent(in) :: dim        ! Dimension of data arrays
        real(real64), intent(in) :: xg(dim)    ! x grid points
        real(real64), intent(in) :: yg(dim)    ! y values at grid points
        integer, intent(in) :: deg        ! Polynomial degree
        real(real64), intent(in) :: x     ! Evaluation point
        integer, intent(in) :: Dmin, Dmax ! Min and max derivative orders
        integer, intent(inout) :: ind     ! Starting index for search
        real(real64), intent(out) :: R(Dmax-Dmin+1)  ! Results (value + derivatives)
        integer, intent(out) :: ierr
        
        real(real64) :: p(0:Dmax+1, 0:deg, 0:deg)  ! Neville's table
        integer :: n, d, i, j
        
        ierr = 0
        
        ! Validate inputs
        if (deg < 0 .or. deg >= dim) then
            ierr = -1
            return
        end if
        
        if (Dmin < 0 .or. Dmax < Dmin) then
            ierr = -2
            return
        end if
        
        ! Check bounds and warn if outside
        if (x < xg(1) .or. x > xg(dim)) then
            write(error_unit, '(A,ES15.8,A,ES15.8,A,ES15.8)') &
                'Warning: eval_neville_polynom: x=', x, &
                ' is outside array bounds [', xg(1), ', ', xg(dim), ']'
        end if
        
        ! Find interpolation interval
        call find_index_for_interp(deg, x, dim, xg, ind)
        
        
        ! Initialize Neville's table
        p = 0.0_real64
        
        ! Set initial values (0th derivative) - following C++ code exactly
        do d = 0, deg
            p(1, d, d) = yg(ind + d)  ! C++ uses 0-based, Fortran uses ind+d for xa[d]
        end do
        
        ! Build the table for derivatives - following C++ algorithm
        do n = 0, Dmax
            do d = 1, deg
                do i = 0, deg - d
                    j = i + d
                    
                    ! Check for near-zero denominator (xa[i] - xa[j])
                    if (abs(xg(ind + i) - xg(ind + j)) < 1.0e-15_real64) then
                        ierr = -3  ! Duplicate x values
                        return
                    end if
                    
                    p(n+1, i, j) = real(n, real64) * (p(n, i, j-1) - p(n, i+1, j)) / &
                                   (xg(ind + i) - xg(ind + j)) + &
                                   ((x - xg(ind + j)) / (xg(ind + i) - xg(ind + j))) * &
                                   p(n+1, i, j-1) - &
                                   ((x - xg(ind + i)) / (xg(ind + i) - xg(ind + j))) * &
                                   p(n+1, i+1, j)
                end do
            end do
        end do
        
        ! Store results
        do n = Dmin, Dmax
            R(n - Dmin + 1) = p(n + 1, 0, deg)
        end do
        
    end subroutine eval_neville_polynom
    
    !---------------------------------------------------------------------------
    ! Lagrange polynomial interpolation
    ! Translates eval_lagrange_polynom from interp.h
    !---------------------------------------------------------------------------
    subroutine eval_lagrange_polynom(dim, xg, yg, deg, x, Dmin, Dmax, ind, result, ierr)
        integer, intent(in) :: dim, deg, Dmin, Dmax
        real(real64), intent(in) :: xg(dim), yg(dim), x
        integer, intent(inout) :: ind
        real(real64), intent(out) :: result
        integer, intent(out) :: ierr
        
        real(real64) :: fac
        integer :: n, i
        
        ierr = 0
        result = 0.0_real64
        
        ! Check for derivatives (not implemented yet)
        if (Dmax > 0) then
            write(error_unit, '(A)') &
                'Warning: eval_lagrange_polynom: derivatives not implemented yet!'
        end if
        
        ! Validate index
        if (ind < 0 .or. ind > dim - 1) then
            ind = dim / 2
        end if
        
        ! Find interpolation interval
        call find_index_for_interp(deg, x, dim, xg, ind)
        
        ! Lagrange interpolation formula
        do n = 0, deg
            fac = 1.0_real64
            do i = 0, deg
                if (i /= n) then
                    if (abs(xg(ind + n) - xg(ind + i)) < 1.0e-15_real64) then
                        ierr = -1  ! Duplicate x values
                        return
                    end if
                    fac = fac * (x - xg(ind + i)) / (xg(ind + n) - xg(ind + i))
                end if
            end do
            result = result + yg(ind + n) * fac
        end do
        
    end subroutine eval_lagrange_polynom
    
    !---------------------------------------------------------------------------
    ! Complex interpolation using Neville's method
    !---------------------------------------------------------------------------
    subroutine eval_complex_interp(dim, x, y, deg, xeval, ind, result, ierr)
        integer, intent(in) :: dim, deg
        real(real64), intent(in) :: x(dim), xeval
        complex(real64), intent(in) :: y(dim)
        integer, intent(inout) :: ind
        complex(real64), intent(out) :: result
        integer, intent(out) :: ierr
        
        real(real64) :: y_real(dim), y_imag(dim), result_real(1), result_imag(1)
        
        ierr = 0
        
        ! Split complex data into real and imaginary parts
        y_real = real(y, real64)
        y_imag = aimag(y)
        
        ! Interpolate real part
        call eval_neville_polynom(dim, x, y_real, deg, xeval, 0, 0, ind, result_real, ierr)
        if (ierr /= 0) return
        
        ! Interpolate imaginary part
        call eval_neville_polynom(dim, x, y_imag, deg, xeval, 0, 0, ind, result_imag, ierr)
        if (ierr /= 0) return
        
        ! Combine results
        result = cmplx(result_real(1), result_imag(1), real64)
        
    end subroutine eval_complex_interp
    
    !---------------------------------------------------------------------------
    ! Binary search in sorted array
    ! Translates binary_search from interp.h
    !---------------------------------------------------------------------------
    function binary_search(x, xa, ilo, ihi) result(index)
        real(real64), intent(in) :: x, xa(:)
        integer, intent(in) :: ilo, ihi
        integer :: index
        
        integer :: lo, hi, mid
        
        lo = ilo
        hi = ihi
        
        ! Binary search loop
        do while (hi > lo + 1)
            mid = (hi + lo) / 2
            if (xa(mid) > x) then
                hi = mid
            else
                lo = mid
            end if
        end do
        
        index = lo
        
    end function binary_search
    
    !---------------------------------------------------------------------------
    ! Search array for interpolation interval
    ! Translates search_array from interp.h
    !---------------------------------------------------------------------------
    subroutine search_array(x, dim, xa, ind)
        real(real64), intent(in) :: x
        integer, intent(in) :: dim
        real(real64), intent(in) :: xa(dim)
        integer, intent(inout) :: ind
        
        ! Bounds checking and correction
        if (ind < 1 .or. ind > dim - 1) then
            ind = dim / 2
        end if
        
        ! Search based on current position
        if (x < xa(ind)) then
            ind = binary_search(x, xa, 1, ind)
        else if (x >= xa(ind + 1)) then
            ind = binary_search(x, xa, ind, dim)
        end if
        
        ! Ensure valid range for interpolation
        ind = max(1, min(ind, dim - 1))
        
    end subroutine search_array
    
    !---------------------------------------------------------------------------
    ! Find index for interpolation
    ! Translates find_index_for_interp from interp.h
    !---------------------------------------------------------------------------
    subroutine find_index_for_interp(deg, x, dimx, xa, ind)
        integer, intent(in) :: deg, dimx
        real(real64), intent(in) :: x, xa(dimx)
        integer, intent(inout) :: ind
        
        ! Adjust index for polynomial degree
        ind = min(ind + (deg - 1) / 2, dimx - 2)
        
        ! Search for correct interval
        call search_array(x, dimx, xa, ind)
        
        ! Adjust for polynomial degree requirements
        ind = max(1, min(ind - (deg - 1) / 2, dimx - deg))
        
    end subroutine find_index_for_interp
    
    !---------------------------------------------------------------------------
    ! Trapezoidal rule integration
    !---------------------------------------------------------------------------
    subroutine integrate_trapezoid(n, x, y, result, ierr)
        integer, intent(in) :: n
        real(real64), intent(in) :: x(n), y(n)
        real(real64), intent(out) :: result
        integer, intent(out) :: ierr
        
        integer :: i
        real(real64) :: dx
        
        ierr = 0
        result = 0.0_real64
        
        if (n < 2) then
            ierr = -1
            return
        end if
        
        ! Trapezoidal rule: ∫f(x)dx ≈ Σ(h/2)[f(x_i) + f(x_i+1)]
        do i = 1, n - 1
            dx = x(i + 1) - x(i)
            if (dx <= 0.0_real64) then
                ierr = -2  ! Non-monotonic x array
                return
            end if
            result = result + 0.5_real64 * dx * (y(i) + y(i + 1))
        end do
        
    end subroutine integrate_trapezoid
    
    !---------------------------------------------------------------------------
    ! Simpson's rule integration
    !---------------------------------------------------------------------------
    subroutine integrate_simpson(n, x, y, result, ierr)
        integer, intent(in) :: n
        real(real64), intent(in) :: x(n), y(n)
        real(real64), intent(out) :: result
        integer, intent(out) :: ierr
        
        integer :: i
        real(real64) :: h, sum_odd, sum_even
        
        ierr = 0
        result = 0.0_real64
        
        if (n < 3) then
            ierr = -1
            return
        end if
        
        ! Check if we have evenly spaced points (simplified Simpson's rule)
        h = x(2) - x(1)
        do i = 2, n - 1
            if (abs((x(i + 1) - x(i)) - h) > 1.0e-12_real64) then
                ! Fall back to trapezoidal rule for non-uniform grid
                call integrate_trapezoid(n, x, y, result, ierr)
                return
            end if
        end do
        
        ! Simpson's 1/3 rule for uniform grid
        sum_odd = 0.0_real64
        sum_even = 0.0_real64
        
        do i = 2, n - 1, 2
            sum_even = sum_even + y(i)
        end do
        
        do i = 3, n - 2, 2
            sum_odd = sum_odd + y(i)
        end do
        
        result = (h / 3.0_real64) * (y(1) + 4.0_real64 * sum_even + 2.0_real64 * sum_odd + y(n))
        
    end subroutine integrate_simpson
    
    !---------------------------------------------------------------------------
    ! Adaptive quadrature integration (simplified GSL replacement)
    !---------------------------------------------------------------------------
    subroutine integrate_adaptive_quadrature(func, a, b, epsabs, epsrel, result, ierr)
        interface
            function func(x) result(y)
                import :: real64
                real(real64), intent(in) :: x
                real(real64) :: y
            end function func
        end interface
        real(real64), intent(in) :: a, b, epsabs, epsrel
        real(real64), intent(out) :: result
        integer, intent(out) :: ierr
        
        real(real64) :: error_estimate, tol
        integer :: max_levels
        
        ierr = 0
        max_levels = 20  ! Maximum recursion depth
        tol = max(epsabs, epsrel * abs(func(a) + func(b)) * abs(b - a) / 2.0_real64)
        
        call adaptive_simpson_recursive(func, a, b, tol, max_levels, result, error_estimate, ierr)
        
    end subroutine integrate_adaptive_quadrature
    
    !---------------------------------------------------------------------------
    ! Recursive adaptive Simpson integration
    !---------------------------------------------------------------------------
    recursive subroutine adaptive_simpson_recursive(func, a, b, tol, max_levels, result, error, ierr)
        interface
            function func(x) result(y)
                import :: real64
                real(real64), intent(in) :: x
                real(real64) :: y
            end function func
        end interface
        real(real64), intent(in) :: a, b, tol
        integer, intent(in) :: max_levels
        real(real64), intent(out) :: result, error
        integer, intent(out) :: ierr
        
        real(real64) :: c, h, fa, fb, fc, S, S1, S2, error1, error2
        
        ierr = 0
        
        ! Base case
        if (max_levels <= 0) then
            h = b - a
            result = h * (func(a) + 4.0_real64 * func((a + b) / 2.0_real64) + func(b)) / 6.0_real64
            error = abs(result) * 1.0e-6_real64  ! Rough estimate
            return
        end if
        
        c = (a + b) / 2.0_real64
        h = b - a
        
        fa = func(a)
        fb = func(b)
        fc = func(c)
        
        ! Simpson's rule for the whole interval
        S = h * (fa + 4.0_real64 * fc + fb) / 6.0_real64
        
        ! Simpson's rule for the two halves
        S1 = (h / 2.0_real64) * (fa + 4.0_real64 * func((a + c) / 2.0_real64) + fc) / 6.0_real64
        S2 = (h / 2.0_real64) * (fc + 4.0_real64 * func((c + b) / 2.0_real64) + fb) / 6.0_real64
        
        error = abs(S - (S1 + S2)) / 15.0_real64  ! Error estimate
        
        if (error < tol) then
            result = S1 + S2 + error  ! Richardson extrapolation
        else
            ! Recursive subdivision
            call adaptive_simpson_recursive(func, a, c, tol / 2.0_real64, max_levels - 1, S1, error1, ierr)
            if (ierr /= 0) return
            
            call adaptive_simpson_recursive(func, c, b, tol / 2.0_real64, max_levels - 1, S2, error2, ierr)
            if (ierr /= 0) return
            
            result = S1 + S2
            error = error1 + error2
        end if
        
    end subroutine adaptive_simpson_recursive
    
    !---------------------------------------------------------------------------
    ! Integration over cylindrical volume
    ! Translates integrate_over_cylinder from C++ code
    !---------------------------------------------------------------------------
    subroutine integrate_over_cylinder(n, r, f, vol_fac, result, ierr)
        integer, intent(in) :: n
        real(real64), intent(in) :: r(n), f(n), vol_fac
        real(real64), intent(out) :: result(n)
        integer, intent(out) :: ierr
        
        real(real64) :: integrand(n)
        integer :: i
        
        ierr = 0
        
        if (n < 2) then
            ierr = -1
            return
        end if
        
        ! Prepare integrand: f(r) * r for cylindrical coordinates
        integrand = f * r
        
        ! Integrate from r(1) to r(i) for each i
        result(1) = 0.0_real64
        
        do i = 2, n
            call integrate_trapezoid(i, r(1:i), integrand(1:i), result(i), ierr)
            if (ierr /= 0) return
            
            result(i) = result(i) * vol_fac
        end do
        
    end subroutine integrate_over_cylinder
    
    !---------------------------------------------------------------------------
    ! Gauss-Legendre quadrature (fixed 5-point rule)
    !---------------------------------------------------------------------------
    subroutine integrate_gauss_legendre(func, a, b, result, ierr)
        interface
            function func(x) result(y)
                import :: real64
                real(real64), intent(in) :: x
                real(real64) :: y
            end function func
        end interface
        real(real64), intent(in) :: a, b
        real(real64), intent(out) :: result
        integer, intent(out) :: ierr
        
        ! 5-point Gauss-Legendre quadrature nodes and weights
        real(real64), parameter :: nodes(5) = [ &
            -0.906179845938664_real64, &
            -0.538469310105683_real64, &
             0.000000000000000_real64, &
             0.538469310105683_real64, &
             0.906179845938664_real64 ]
        
        real(real64), parameter :: weights(5) = [ &
             0.236926885056189_real64, &
             0.478628670499366_real64, &
             0.568888888888889_real64, &
             0.478628670499366_real64, &
             0.236926885056189_real64 ]
        
        real(real64) :: x, scale
        integer :: i
        
        ierr = 0
        result = 0.0_real64
        scale = (b - a) / 2.0_real64
        
        do i = 1, 5
            x = (a + b) / 2.0_real64 + scale * nodes(i)
            result = result + weights(i) * func(x)
        end do
        
        result = result * scale
        
    end subroutine integrate_gauss_legendre

end module kilca_interp_integ_m