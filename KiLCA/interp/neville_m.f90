!> Neville polynomial interpolation with derivatives.
!>
!> Fortran port of the Fortran-facing wrappers eval_neville_polynom_ and
!> eval_neville_polynom_ready_ formerly defined in interp.cpp (the underlying
!> algorithm lived as inline functions in interp.h). The C symbols and their
!> by-reference int*/double* argument convention are preserved so the Fortran
!> callers (KiLCA/QL-Balance wave_code_data) link unchanged. The C++ inline
!> versions in interp.h are untouched and keep serving the C++ call sites.
module kilca_neville_m
    use, intrinsic :: iso_c_binding, only: c_int, c_double
    use, intrinsic :: iso_fortran_env, only: error_unit
    implicit none
    private

    public :: eval_neville_polynom, eval_neville_polynom_ready

contains

    subroutine eval_neville_polynom(dim, xg, yg, deg, x, Dmin, Dmax, ind, R) &
        bind(C, name="eval_neville_polynom_")
        integer(c_int), intent(in) :: dim, deg, Dmin, Dmax
        integer(c_int), intent(inout) :: ind
        real(c_double), intent(in) :: xg(0:*), yg(0:*), x
        real(c_double), intent(out) :: R(0:*)
        real(c_double), parameter :: tol = 1.0d-9
        integer(c_int) :: edeg

        if (x < xg(0) - tol .or. x > xg(dim - 1) + tol) &
            write (error_unit, '(a,es12.4)') &
            "warning: eval_neville_polynom: x is outside the array: x=", x

        ! Cap the moving-window degree to the points actually available. With
        ! enough points (dim > deg, the production case) edeg == deg, matching
        ! the C++ exactly; for a short grid the C++ formula dimx-deg-1 went
        ! negative and indexed out of bounds, so clamp to stay in range.
        edeg = min(deg, dim - 1)

        if (ind < 0 .or. ind > dim - 1) ind = dim/2
        call find_index_for_interp(edeg, x, dim, xg, ind)

        call neville_ready(xg(ind:ind + edeg), yg(ind:ind + edeg), edeg, x, Dmin, Dmax, R)
    end subroutine eval_neville_polynom

    subroutine eval_neville_polynom_ready(xa, ya, deg, x, Dmin, Dmax, R) &
        bind(C, name="eval_neville_polynom_ready_")
        integer(c_int), intent(in) :: deg, Dmin, Dmax
        real(c_double), intent(in) :: xa(0:*), ya(0:*), x
        real(c_double), intent(out) :: R(0:*)

        call neville_ready(xa(0:deg), ya(0:deg), deg, x, Dmin, Dmax, R)
    end subroutine eval_neville_polynom_ready

    subroutine neville_ready(xa, ya, deg, x, Dmin, Dmax, R)
        integer(c_int), intent(in) :: deg, Dmin, Dmax
        real(c_double), intent(in) :: xa(0:deg), ya(0:deg), x
        real(c_double), intent(out) :: R(0:*)
        real(c_double), parameter :: tol = 1.0d-9
        real(c_double), allocatable :: p(:, :, :)
        integer(c_int) :: d, i, j, n

        if (x < xa(0) - tol .or. x > xa(deg) + tol) &
            write (error_unit, '(a,3es24.16)') &
            "warning: eval_neville_polynom: x outside [a,b]: ", xa(0), xa(deg), x

        allocate (p(0:Dmax + 1, 0:deg, 0:deg))
        p = 0.0d0

        do d = 0, deg
            p(1, d, d) = ya(d)
        end do

        do n = 0, Dmax
            do d = 1, deg
                do i = 0, deg - d
                    j = i + d
                    p(n + 1, i, j) = n*(p(n, i, j - 1) - p(n, i + 1, j))/(xa(i) - xa(j)) + &
                                     ((x - xa(j))/(xa(i) - xa(j)))*p(n + 1, i, j - 1) - &
                                     ((x - xa(i))/(xa(i) - xa(j)))*p(n + 1, i + 1, j)
                end do
            end do
        end do

        do n = Dmin, Dmax
            R(n - Dmin) = p(n + 1, 0, deg)
        end do
    end subroutine neville_ready

    subroutine find_index_for_interp(deg, x, dimx, xa, ind)
        integer(c_int), intent(in) :: deg, dimx
        real(c_double), intent(in) :: x, xa(0:*)
        integer(c_int), intent(inout) :: ind

        ind = min(ind + (deg - 1)/2, dimx - 2)
        call search_array(x, dimx, xa, ind)
        ind = min(max(0, ind - (deg - 1)/2), dimx - deg - 1)
    end subroutine find_index_for_interp

    subroutine search_array(x, dimx, xa, ind)
        real(c_double), intent(in) :: x, xa(0:*)
        integer(c_int), intent(in) :: dimx
        integer(c_int), intent(inout) :: ind

        if (x < xa(ind)) then
            ind = binary_search(x, xa, 0, ind)
        else if (x >= xa(ind + 1)) then
            ind = binary_search(x, xa, ind, dimx - 1)
        end if
    end subroutine search_array

    integer(c_int) function binary_search(x, xa, ilo, ihi) result(res)
        real(c_double), intent(in) :: x, xa(0:*)
        integer(c_int), intent(in) :: ilo, ihi
        integer(c_int) :: lo, hi, i

        lo = ilo; hi = ihi
        do while (hi > lo + 1)
            i = (hi + lo)/2
            if (xa(i) > x) then
                hi = i
            else
                lo = i
            end if
        end do
        res = lo
    end function binary_search

end module kilca_neville_m
