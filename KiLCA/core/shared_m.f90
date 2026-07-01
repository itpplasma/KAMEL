!> The live subset of core/shared.{h,cpp}'s common helpers. Deliberately
!> partial: xmalloc/xrealloc/compare_doubles/Ckn/calc_interp_polynom/
!> eval_interp_polynom/localizator/localizator_4_derivs had zero live
!> callers anywhere in KiLCA or QL-Balance (verified by grep across all
!> C++, header and Fortran source files) and were dropped, matching this
!> port's established dead-code precedent.
!>
!> signum and sort_index_doubles keep C++ linkage (bind(C), no trailing
!> underscore) since their one remaining caller each is still C++
!> (progs/main_eig_param.cpp, math/adapt_grid/adaptive_grid.cpp - both S8
!> scope, not yet translated).
!>
!> binomial_coefficients is the live subset's odd one out: it is called
!> only by PRE-EXISTING legacy Fortran (flre/conductivity/kmatrices*.f90,
!> conduct_arrays*.f90) via implicit/external interface
!> (`call binomial_coefficients (%val(N), bico)`, no trailing underscore
!> in the source text) - gfortran mangles that call to the external
!> symbol "binomial_coefficients_". A free (non-module-contained)
!> subprogram below provides exactly that symbol, per this port's
!> established free-subprogram-vs-module-contained-procedure rule.
module kilca_shared_m
    use, intrinsic :: iso_c_binding, only: c_double, c_int, c_size_t
    use constants, only: dp
    use fortnum_multiroot, only: argsort
    implicit none
    private

    public :: signum, sort_index_doubles, strtol_int

contains

    !> Reproduces C's strtol(s, NULL, 10) exactly: skip leading whitespace,
    !> an optional sign, then consume a run of decimal digits, stopping at
    !> the first non-digit character (0 if no digits are found at all).
    !> This is NOT the same as a clean "parse the whole string as an
    !> integer, 0 on any failure" - the .in file readers throughout this
    !> port's oracle (read_line_2get_int_ in io/inout.cpp, and every
    !> zone-settings/eigmode-settings reader built on the same convention)
    !> use exactly this truncating strtol behavior, so genuinely malformed-
    !> looking input like "1.00e+05" for an integer field is not an error
    !> in the oracle - it silently becomes 1 (the leading digit run before
    !> the decimal point), and this port's Fortran translations must
    !> reproduce that silent truncation bit-for-bit rather than "properly"
    !> parsing the full numeric literal or defaulting to 0.
    function strtol_int(s) result(val)
        character(len=*), intent(in) :: s
        integer :: val
        integer :: i, n, sgn
        logical :: found_digit

        n = len(s)
        i = 1
        val = 0
        sgn = 1
        found_digit = .false.

        do while (i <= n)
            if (s(i:i) /= ' ' .and. s(i:i) /= achar(9)) exit
            i = i + 1
        end do

        if (i <= n) then
            if (s(i:i) == '-') then
                sgn = -1
                i = i + 1
            else if (s(i:i) == '+') then
                i = i + 1
            end if
        end if

        do while (i <= n)
            if (s(i:i) < '0' .or. s(i:i) > '9') exit
            val = val*10 + (iachar(s(i:i)) - iachar('0'))
            found_digit = .true.
            i = i + 1
        end do

        if (.not. found_digit) val = 0
        val = val*sgn
    end function strtol_int

    function signum(x) result(s) bind(C, name="signum")
        real(c_double), value :: x
        integer(c_int) :: s
        if (x < 0.0_c_double) then
            s = -1
        else if (x == 0.0_c_double) then
            s = 0
        else
            s = 1
        end if
    end function signum

    !> Ascending index sort: perm(k) receives the (0-based) index into x of
    !> the k-th smallest value. Uses fortnum_multiroot's native argsort
    !> (heapsort, 1-based) then converts to 0-based for the still-C++
    !> caller. argsort's heapsort is not stable; the conductivity grid
    !> shares zone-boundary nodes (duplicated x), so the sort must be
    !> deterministic for the spline to be well-defined - the oracle's own
    !> fixup (break exact-key ties by ascending original index) is
    !> reproduced verbatim below.
    subroutine sort_index_doubles(perm, xarr, n) bind(C, name="sort_index_doubles")
        integer(c_size_t), value :: n
        real(c_double), intent(in) :: xarr(0:n - 1)
        integer(c_size_t), intent(out) :: perm(0:n - 1)
        integer, allocatable :: fperm(:)
        integer :: nn, i
        integer(c_size_t) :: a, b

        nn = int(n)
        allocate (fperm(nn))
        call argsort(xarr, fperm)
        do i = 1, nn
            perm(i - 1) = int(fperm(i) - 1, c_size_t)
        end do
        deallocate (fperm)

        a = 0
        do while (a < n)
            b = a + 1
            do while (b < n)
                if (xarr(perm(b)) /= xarr(perm(a))) exit
                b = b + 1
            end do
            if (b - a > 1) call sort_ascending_size_t(perm(a:b - 1))
            a = b
        end do
    end subroutine sort_index_doubles

    !> Small-range insertion sort (std::sort's effect on a tie-block here
    !> is always small: at most a handful of coincident zone-boundary
    !> nodes), ascending by raw index value, matching std::sort(perm+a,
    !> perm+b) with the default operator< on size_t.
    subroutine sort_ascending_size_t(arr)
        integer(c_size_t), intent(inout) :: arr(:)
        integer :: i, j
        integer(c_size_t) :: key
        do i = 2, size(arr)
            key = arr(i)
            j = i - 1
            do while (j >= 1)
                if (arr(j) <= key) exit
                arr(j + 1) = arr(j)
                j = j - 1
            end do
            arr(j + 1) = key
        end do
    end subroutine sort_ascending_size_t

end module kilca_shared_m

!> Free (non-module-contained) subprogram: see the module banner above for
!> why this cannot live inside kilca_shared_m's `contains` block.
subroutine binomial_coefficients(N, BC)
    use, intrinsic :: iso_c_binding, only: c_int, c_double
    implicit none
    integer(c_int), value :: N
    real(c_double), intent(out) :: BC(0:N, 0:N)
    integer :: k, n_idx
    real(c_double) :: tmp

    do n_idx = 0, N
        tmp = 1.0_c_double
        BC(n_idx, 0) = tmp
        do k = 1, n_idx
            tmp = tmp*real(n_idx - k + 1, c_double)/real(k, c_double)
            BC(n_idx, k) = tmp
        end do
    end do
end subroutine binomial_coefficients
