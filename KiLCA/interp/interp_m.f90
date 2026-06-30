!> Adaptive grid thinning by polynomial-interpolation error, formerly
!> KiLCA/math/adapt_grid/adaptive_grid_pol.cpp's `sparse_grid_polynom` -- the
!> 12-argument, multi-component overload declared in adaptive_grid_pol.h and
!> actually reachable from flre_zone.cpp (which includes adaptive_grid_pol.h,
!> not interp.h). There is an unrelated, scalar 11-argument function with the
!> same name in KiLCA/interp/interp.cpp; flre_zone.cpp never includes
!> interp.h, so that one is not its caller's target and is left untranslated.
!>
!> find_index_for_interp/search_array/binary_search/eval_interp_polynom/
!> func_interp below are private helpers ported from the same
!> adaptive_grid_pol.{h,cpp} file. They are intentionally NOT kilca_neville_m
!> (a different module, ported earlier for a different call site): this
!> algorithm evaluates a plain Lagrange interpolation value only (no
!> derivatives, no out-of-range tolerance warning), with an explicit `dimy`
!> stride for vector-valued y data.
!>
!> The C++ pointer-swapping between (xold,yold,iold) and (xnew,ynew,inew)
!> each iteration (and the `if (x1 != xold)` skip-copy optimization at the
!> end) is a memory-reuse detail with no numerical effect; here the swap is
!> done with `move_alloc` and the final result is always copied into the
!> caller's x1/y1/ind1, unconditionally.
module kilca_interp_m
    use constants, only: dp
    implicit none
    private

    public :: sparse_grid_polynom

contains

    !> x(0:dimx_in-1)/y(dimy,0:dimx_in-1): the fine input grid/data.
    !> eps_a/eps_r: in/out -- read as the requested tolerances, overwritten
    !> on return with the achieved maximum error (matching the oracle, which
    !> stores merr into both regardless of their distinct input meanings).
    !> x1/y1/ind1(0:dimx_in-1): caller-allocated output buffers (thinned
    !> grid, thinned data, source-index map), valid for 0:dimout-1 on return.
    subroutine sparse_grid_polynom(x, dimx_in, y, dimy, deg, eps_a, eps_r, step, &
                                    dimout, x1, y1, ind1)
        integer, intent(in) :: dimx_in
        integer, intent(in) :: dimy
        real(dp), intent(in) :: x(0:dimx_in - 1)
        real(dp), intent(in) :: y(dimy, 0:dimx_in - 1)
        integer, intent(in) :: deg
        real(dp), intent(inout) :: eps_a, eps_r
        real(dp), intent(in) :: step
        integer, intent(out) :: dimout
        real(dp), intent(out) :: x1(0:dimx_in - 1)
        real(dp), intent(out) :: y1(dimy, 0:dimx_in - 1)
        integer, intent(out) :: ind1(0:dimx_in - 1)

        real(dp), allocatable :: xold(:), xnew(:), xtmp(:)
        real(dp), allocatable :: yold(:, :), ynew(:, :), ytmp(:, :)
        integer, allocatable :: iold(:), inew(:), itmp(:)
        real(dp), allocatable :: ym(:)
        integer :: xdim, dimx, k, ind, ic, node, j
        real(dp) :: xc, yc, errr, erra, errt, merr, xmerr, ymerr
        real(dp), parameter :: pi_local = 3.141592653589793238_dp
        logical :: flag

        xdim = dimx_in
        dimx = deg + 1

        allocate (xold(0:xdim - 1), yold(dimy, 0:xdim - 1), iold(0:xdim - 1))
        allocate (xnew(0:xdim - 1), ynew(dimy, 0:xdim - 1), inew(0:xdim - 1))
        allocate (ym(dimy))

        iold(0) = 0
        xold(0) = x(iold(0))

        iold(dimx - 1) = xdim - 1
        xold(dimx - 1) = x(iold(dimx - 1))

        do k = 1, deg - 1
            iold(k) = int(0.5_dp*xdim* &
                          (1.0_dp - cos(real(2*k - 1, dp)*pi_local/real(2*(deg - 1), dp))))
            xold(k) = x(iold(k))
        end do

        do k = 0, dimx - 1
            call func_interp(iold(k), yold(:, k), y, dimy)
        end do

        merr = 0.0_dp

        do
            ind = 0
            node = 0
            merr = 0.0_dp
            xmerr = 0.0_dp
            ymerr = 0.0_dp

            do k = 0, dimx - 2
                xnew(node) = xold(k)
                inew(node) = iold(k)
                ynew(:, node) = yold(:, k)
                node = node + 1

                ic = iold(k) + (iold(k + 1) - iold(k))/2

                if (ic == iold(k)) cycle

                xc = x(ic)

                call find_index_for_interp(deg, xc, dimx, xold, ind)

                call eval_interp_polynom(deg, xold(ind:), yold(:, ind:), dimy, xc, ym)

                flag = .false.
                do j = 1, dimy
                    yc = y(j, ic)

                    erra = abs(yc - ym(j))
                    errr = eps_r*abs(yc)

                    errt = erra
                    if (abs(yc) > 1.0e-16_dp) errt = min(errt, errt/abs(yc))

                    if (errt > merr) then
                        merr = errt
                        xmerr = xc
                        ymerr = yc
                    end if

                    if (.not. (erra < eps_a .or. erra < errr) .or. xc - xold(k) > step) then
                        flag = .true.
                    end if
                end do

                if (.not. flag) cycle

                xnew(node) = xc
                inew(node) = ic
                call func_interp(ic, ynew(:, node), y, dimy)
                node = node + 1
            end do

            xnew(node) = xold(dimx - 1)
            inew(node) = iold(dimx - 1)
            ynew(:, node) = yold(:, dimx - 1)
            node = node + 1

            call move_alloc(xold, xtmp)
            call move_alloc(xnew, xold)
            call move_alloc(xtmp, xnew)

            call move_alloc(yold, ytmp)
            call move_alloc(ynew, yold)
            call move_alloc(ytmp, ynew)

            call move_alloc(iold, itmp)
            call move_alloc(inew, iold)
            call move_alloc(itmp, inew)

            if (node == dimx) exit

            dimx = node
        end do

        x1(0:dimx - 1) = xold(0:dimx - 1)
        ind1(0:dimx - 1) = iold(0:dimx - 1)
        y1(:, 0:dimx - 1) = yold(:, 0:dimx - 1)

        dimout = dimx
        eps_a = merr
        eps_r = merr

        if (dimx == xdim) then
            write (*, '(a,i0,a)') &
                'sparse_grid_polynom: warning: no points removed from grid: dim = ', dimx, '.'
        end if
    end subroutine sparse_grid_polynom

    subroutine find_index_for_interp(deg, xc, dimx, xa, ind)
        integer, intent(in) :: deg, dimx
        real(dp), intent(in) :: xc, xa(0:*)
        integer, intent(inout) :: ind
        ind = min(ind + (deg - 1)/2, dimx - 2)
        call search_array(xc, dimx, xa, ind)
        ind = min(max(0, ind - (deg - 1)/2), dimx - deg - 1)
    end subroutine find_index_for_interp

    subroutine search_array(xc, dimx, xa, ind)
        real(dp), intent(in) :: xc, xa(0:*)
        integer, intent(in) :: dimx
        integer, intent(inout) :: ind
        if (xc < xa(ind)) then
            ind = binary_search(xc, xa, 0, ind)
        else if (xc >= xa(ind + 1)) then
            ind = binary_search(xc, xa, ind, dimx - 1)
        end if
    end subroutine search_array

    integer function binary_search(xc, xa, ilo, ihi) result(res)
        real(dp), intent(in) :: xc, xa(0:*)
        integer, intent(in) :: ilo, ihi
        integer :: lo, hi, i
        lo = ilo
        hi = ihi
        do while (hi > lo + 1)
            i = (hi + lo)/2
            if (xa(i) > xc) then
                hi = i
            else
                lo = i
            end if
        end do
        res = lo
    end function binary_search

    !> Plain (non-Neville) Lagrange evaluation at xc of the degree-`deg`
    !> polynomial through (xg(0:deg), yg(:,0:deg)), vectorized over dimy
    !> independent components.
    subroutine eval_interp_polynom(deg, xg, yg, dimy, xc, yout)
        integer, intent(in) :: deg, dimy
        real(dp), intent(in) :: xg(0:*), yg(dimy, 0:*), xc
        real(dp), intent(out) :: yout(dimy)
        integer :: j, n, i
        real(dp) :: fac
        do j = 1, dimy
            yout(j) = 0.0_dp
            do n = 0, deg
                fac = 1.0_dp
                do i = 0, deg
                    if (i /= n) fac = fac*(xc - xg(i))/(xg(n) - xg(i))
                end do
                yout(j) = yout(j) + yg(j, n)*fac
            end do
        end do
    end subroutine eval_interp_polynom

    !> Oracle signature also takes (x, dimx); dropped here since the body
    !> only ever reads y(:,ind) (the x/dimx parameters are unused in the
    !> oracle's own definition too).
    subroutine func_interp(ind, yout, y, dimy)
        integer, intent(in) :: ind, dimy
        real(dp), intent(in) :: y(dimy, 0:*)
        real(dp), intent(out) :: yout(dimy)
        yout = y(:, ind)
    end subroutine func_interp

end module kilca_interp_m
