!> Vector-valued adaptive 1D grid refinement by moving polynomial
!> interpolation, formerly KiLCA/math/adapt_grid/adaptive_grid_pol.cpp. The two
!> live entry points are adaptive_grid_polynom_res and adaptive_grid_polynom_err,
!> both called by the Fortran conductivity-grid generator (kilca_cond_profiles_m)
!> through bind(C) interface blocks, so the C symbol names are preserved.
!>
!> Dropped as dead (verified by grep across all C++, header and Fortran
!> sources):
!>   - the two-overload adaptive_grid_polynom: the dimy/no-ind_err overload is
!>     reached only from the orphaned demo main_agp.cpp (in no CMake target);
!>     the dim_err/ind_err overload had exactly one caller, the thin extern "C"
!>     adaptive_grid_polynom_err forwarder, whose algorithm body is therefore
!>     folded directly into adaptive_grid_polynom_err here.
!>   - sparse_grid_polynom and its helpers eval_interp_polynom/func_interp were
!>     already ported (kilca_interp_m) and have no caller of this file's copies.
!>   - sign was never called.
!>   - the DEBUG_FLAG (compile-time 0) and local debug (always 0) warning prints
!>     never executed, so their dead branches and the xmerr/ymerr/jmerr trace
!>     variables that only fed them are omitted.
!>
!> find_index_for_interp/search_array/binary_search and this file's own
!> yshift-strided eval_neville_polynom (distinct from kilca_neville_m's
!> Dmin/Dmax variant) are kept as private helpers. The C++ pointer-swap between
!> (xold,yold) and (xnew,ynew) plus the "if (x1 != xold)" skip-copy is a
!> memory-reuse detail with no numerical effect; here the swap uses move_alloc
!> and the result is always copied into the caller's x1/y1.
module adaptive_grid_pol_m
    use, intrinsic :: iso_c_binding, only: c_double, c_int, c_ptr, c_funptr, &
                                           c_f_procpointer
    implicit none
    private

    public :: adaptive_grid_polynom_res, adaptive_grid_polynom_err
    public :: find_index_for_interp, eval_neville_polynom

    abstract interface
        subroutine sample_cb(r, fval, ctx) bind(C)
            import :: c_double, c_ptr
            real(c_double), intent(in) :: r
            real(c_double), intent(out) :: fval(*)
            type(c_ptr), value :: ctx
        end subroutine sample_cb
    end interface

contains

    !> Refines a Chebyshev-seeded grid until, for each component listed in
    !> ind_err, the error of the central moving polynomial against its left and
    !> right neighbours falls below eps. Returns 0 on convergence, 1 if the
    !> buffer dimension xdim was reached first (then the previous grid and its
    !> error are restored). Body folded from the C++ dim_err overload of
    !> adaptive_grid_polynom that the extern "C" forwarder called.
    function adaptive_grid_polynom_err(f, p, a, b, dimy, deg, xdim, eps, &
                                       dim_err, ind_err, x1, y1) result(stat) &
        bind(C, name="adaptive_grid_polynom_err")
        type(c_funptr), value :: f
        type(c_ptr), value :: p
        real(c_double), value :: a, b
        integer(c_int), value :: dimy, deg, dim_err
        integer(c_int), intent(inout) :: xdim
        real(c_double), intent(inout) :: eps
        integer(c_int), intent(in) :: ind_err(0:*)
        real(c_double), intent(inout) :: x1(0:*), y1(0:*)
        integer(c_int) :: stat

        procedure(sample_cb), pointer :: cb
        real(c_double), allocatable :: xold(:), xnew(:), yold(:), ynew(:)
        real(c_double), allocatable :: yl(:), ym(:), yr(:)
        integer(c_int) :: maxdim, dimx, k, j, l, ind, indl, indr, node, m, jj, flag
        real(c_double) :: pi, xc, errl, errr, errt, merr, peps, err

        call c_f_procpointer(f, cb)

        maxdim = xdim

        allocate (xold(0:maxdim - 1), xnew(0:maxdim - 1))
        allocate (yold(0:maxdim*dimy - 1), ynew(0:maxdim*dimy - 1))
        allocate (yl(0:dim_err - 1), ym(0:dim_err - 1), yr(0:dim_err - 1))

        dimx = deg + 3

        xold(0) = a
        xold(dimx - 1) = b

        pi = 3.141592653589793238_c_double

        do k = 1, deg + 1
            xold(k) = 0.5_c_double*(a + b) - 0.5_c_double*(b - a) &
                      *cos(real(2*k - 1, c_double)*pi/real(2*(deg + 1), c_double))
        end do

        do k = 0, dimx - 1
            call cb(xold(k), yold(k*dimy), p)
        end do

        merr = 0.0_c_double

        do
            ind = 0
            node = 0

            peps = merr
            merr = 0.0_c_double

            do k = 0, dimx - 2
                xnew(node) = xold(k)
                do j = 0, dimy - 1
                    ynew(node*dimy + j) = yold(k*dimy + j)
                end do
                node = node + 1

                xc = xold(k) + 0.5_c_double*(xold(k + 1) - xold(k))

                call find_index_for_interp(deg, xc, dimx, xold, ind)

                indl = max(0, ind - 1)
                indr = min(ind + 1, dimx - deg - 1)

                flag = 0

                do l = 0, dim_err - 1
                    j = ind_err(l)

                    call eval_neville_polynom(xold(ind), yold(ind*dimy + j), &
                                              dimy, deg, xc, ym(l))
                    call eval_neville_polynom(xold(indl), yold(indl*dimy + j), &
                                              dimy, deg, xc, yl(l))
                    call eval_neville_polynom(xold(indr), yold(indr*dimy + j), &
                                              dimy, deg, xc, yr(l))

                    errl = abs(yl(l) - ym(l))
                    errr = abs(yr(l) - ym(l))

                    errt = max(errl, errr)

                    if (abs(ym(l)) > 1.0_c_double) errt = errt/abs(ym(l))

                    if (errt > merr) merr = errt

                    err = eps

                    if (errt > err) flag = 1
                end do

                if (node >= maxdim - 2) then
                    do m = 0, dimx - 1
                        x1(m) = xold(m)
                        do jj = 0, dimy - 1
                            y1(m*dimy + jj) = yold(m*dimy + jj)
                        end do
                    end do
                    xdim = dimx
                    eps = peps
                    stat = 1
                    return
                end if

                if (flag == 0) cycle

                xnew(node) = xc
                call cb(xnew(node), ynew(node*dimy), p)
                node = node + 1
            end do

            xnew(node) = xold(dimx - 1)
            do j = 0, dimy - 1
                ynew(node*dimy + j) = yold((dimx - 1)*dimy + j)
            end do
            node = node + 1

            call swap_real(xold, xnew)
            call swap_real(yold, ynew)

            if (node == dimx) exit

            dimx = node
        end do

        do m = 0, dimx - 1
            x1(m) = xold(m)
            do jj = 0, dimy - 1
                y1(m*dimy + jj) = yold(m*dimy + jj)
            end do
        end do

        xdim = dimx
        eps = merr
        stat = 0
    end function adaptive_grid_polynom_err

    !> Same refinement as adaptive_grid_polynom_err, but the per-point error
    !> threshold is relaxed away from the resonance radius r_res by a Gaussian
    !> of width D, dropping to eps_res at r_res and rising to eps far from it.
    function adaptive_grid_polynom_res(f, p, a, b, dimy, deg, xdim, eps, &
                                       r_res, D, eps_res, dim_err, ind_err, x1, y1) &
        result(stat) bind(C, name="adaptive_grid_polynom_res")
        type(c_funptr), value :: f
        type(c_ptr), value :: p
        real(c_double), value :: a, b, r_res, D, eps_res
        real(c_double), intent(inout) :: eps
        integer(c_int), value :: dimy, deg, dim_err
        integer(c_int), intent(inout) :: xdim
        integer(c_int), intent(in) :: ind_err(0:*)
        real(c_double), intent(inout) :: x1(0:*), y1(0:*)
        integer(c_int) :: stat

        procedure(sample_cb), pointer :: cb
        real(c_double), allocatable :: xold(:), xnew(:), yold(:), ynew(:)
        real(c_double), allocatable :: yl(:), ym(:), yr(:)
        integer(c_int) :: maxdim, dimx, k, j, l, ind, indl, indr, node, m, jj, flag
        real(c_double) :: pi, xc, errl, errr, errt, merr, peps, err

        call c_f_procpointer(f, cb)

        maxdim = xdim

        allocate (xold(0:maxdim - 1), xnew(0:maxdim - 1))
        allocate (yold(0:maxdim*dimy - 1), ynew(0:maxdim*dimy - 1))
        allocate (yl(0:dim_err - 1), ym(0:dim_err - 1), yr(0:dim_err - 1))

        dimx = deg + 3

        xold(0) = a
        xold(dimx - 1) = b

        pi = 3.141592653589793238_c_double

        do k = 1, deg + 1
            xold(k) = 0.5_c_double*(a + b) - 0.5_c_double*(b - a) &
                      *cos(real(2*k - 1, c_double)*pi/real(2*(deg + 1), c_double))
        end do

        do k = 0, dimx - 1
            call cb(xold(k), yold(k*dimy), p)
        end do

        merr = 0.0_c_double

        do
            ind = 0
            node = 0

            peps = merr
            merr = 0.0_c_double

            do k = 0, dimx - 2
                xnew(node) = xold(k)
                do j = 0, dimy - 1
                    ynew(node*dimy + j) = yold(k*dimy + j)
                end do
                node = node + 1

                xc = xold(k) + 0.5_c_double*(xold(k + 1) - xold(k))

                call find_index_for_interp(deg, xc, dimx, xold, ind)

                indl = max(0, ind - 1)
                indr = min(ind + 1, dimx - deg - 1)

                flag = 0

                do l = 0, dim_err - 1
                    j = ind_err(l)

                    call eval_neville_polynom(xold(ind), yold(ind*dimy + j), &
                                              dimy, deg, xc, ym(l))
                    call eval_neville_polynom(xold(indl), yold(indl*dimy + j), &
                                              dimy, deg, xc, yl(l))
                    call eval_neville_polynom(xold(indr), yold(indr*dimy + j), &
                                              dimy, deg, xc, yr(l))

                    errl = abs(yl(l) - ym(l))
                    errr = abs(yr(l) - ym(l))

                    errt = max(errl, errr)

                    if (abs(ym(l)) > 1.0_c_double) errt = errt/abs(ym(l))

                    if (errt > merr) merr = errt

                    err = (eps - eps_res)*(1.0_c_double &
                          - exp(-(xc - r_res)*(xc - r_res)/D/D)) + eps_res

                    if (errt > err) flag = 1
                end do

                if (node >= maxdim - 2) then
                    do m = 0, dimx - 1
                        x1(m) = xold(m)
                        do jj = 0, dimy - 1
                            y1(m*dimy + jj) = yold(m*dimy + jj)
                        end do
                    end do
                    xdim = dimx
                    eps = peps
                    stat = 1
                    return
                end if

                if (flag == 0) cycle

                xnew(node) = xc
                call cb(xnew(node), ynew(node*dimy), p)
                node = node + 1
            end do

            xnew(node) = xold(dimx - 1)
            do j = 0, dimy - 1
                ynew(node*dimy + j) = yold((dimx - 1)*dimy + j)
            end do
            node = node + 1

            call swap_real(xold, xnew)
            call swap_real(yold, ynew)

            if (node == dimx) exit

            dimx = node
        end do

        do m = 0, dimx - 1
            x1(m) = xold(m)
            do jj = 0, dimy - 1
                y1(m*dimy + jj) = yold(m*dimy + jj)
            end do
        end do

        xdim = dimx
        eps = merr
        stat = 0
    end function adaptive_grid_polynom_res

    subroutine swap_real(u, v)
        real(c_double), allocatable, intent(inout) :: u(:), v(:)
        real(c_double), allocatable :: t(:)
        call move_alloc(u, t)
        call move_alloc(v, u)
        call move_alloc(t, v)
    end subroutine swap_real

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

    function binary_search(x, xa, ilo, ihi) result(res)
        real(c_double), intent(in) :: x, xa(0:*)
        integer(c_int), intent(in) :: ilo, ihi
        integer(c_int) :: res, lo, hi, i
        lo = ilo
        hi = ihi
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

    !> Neville interpolation of one component sampled with stride yshift:
    !> R = value at x of the degree-deg polynomial through xa(0:deg) and
    !> ya(0), ya(yshift), ..., ya(deg*yshift).
    subroutine eval_neville_polynom(xa, ya, yshift, deg, x, R)
        real(c_double), intent(in) :: xa(0:*), ya(0:*)
        integer(c_int), intent(in) :: yshift, deg
        real(c_double), intent(in) :: x
        real(c_double), intent(out) :: R
        real(c_double) :: pp(0:deg, 0:deg)
        integer(c_int) :: d, i, j
        do d = 0, deg
            pp(d, d) = ya(d*yshift)
        end do
        do d = 1, deg
            do i = 0, deg - d
                j = i + d
                pp(i, j) = ((x - xa(j))/(xa(i) - xa(j)))*pp(i, j - 1) &
                           - ((x - xa(i))/(xa(i) - xa(j)))*pp(i + 1, j)
            end do
        end do
        R = pp(0, deg)
    end subroutine eval_neville_polynom

end module adaptive_grid_pol_m
