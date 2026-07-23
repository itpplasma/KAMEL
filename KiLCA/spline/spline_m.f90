!> Splines of arbitrary odd degree for arrays of data.
!>
!> Fortran port of the former spline.cpp. The C ABI (spline_alloc_,
!> spline_calc_, spline_eval_, spline_eval_d_, spline_free_) is preserved so the
!> KiLCA C++ callers link unchanged. The opaque handle sid is the C address of
!> the spline_data_t allocated here; x and the coefficient array C stay owned by
!> the caller and are referenced through stored C pointers.
module kilca_spline_m
    use, intrinsic :: iso_c_binding, only: c_int, c_double, c_intptr_t, c_ptr, &
        c_loc, c_f_pointer, c_associated, c_null_ptr
    implicit none
    private

    public :: spline_alloc, spline_calc, spline_eval, spline_eval_d, spline_free

    type :: spline_data_t
        integer(c_int) :: N
        integer(c_int) :: ind
        integer(c_int) :: stype
        integer(c_int) :: dimx
        type(c_ptr) :: x_ptr = c_null_ptr
        type(c_ptr) :: C_ptr = c_null_ptr
        real(c_double), allocatable :: BC(:)
        real(c_double), allocatable :: fac(:)
    end type spline_data_t

    interface
        subroutine dgesv(n, nrhs, a, lda, ipiv, b, ldb, info) bind(C, name="dgesv_")
            import :: c_int, c_double
            integer(c_int) :: n, nrhs, lda, ldb, info
            integer(c_int) :: ipiv(*)
            real(c_double) :: a(*), b(*)
        end subroutine dgesv
        subroutine dgbsv(n, kl, ku, nrhs, ab, ldab, ipiv, b, ldb, info) &
            bind(C, name="dgbsv_")
            import :: c_int, c_double
            integer(c_int) :: n, kl, ku, nrhs, ldab, ldb, info
            integer(c_int) :: ipiv(*)
            real(c_double) :: ab(*), b(*)
        end subroutine dgbsv
    end interface

contains

    subroutine spline_alloc(N, stype, dimx, x, Carr, sid) bind(C, name="spline_alloc_")
        integer(c_int), value :: N, stype, dimx
        type(c_ptr), value :: x, Carr
        integer(c_intptr_t), intent(out) :: sid
        type(spline_data_t), pointer :: sd

        if (.not. (N > 0 .and. mod(N, 2) == 1)) then
            sid = 0_c_intptr_t
            call alloc_fail()
            return
        end if
        if (.not. (dimx >= N .and. stype == 1)) then
            sid = 0_c_intptr_t
            call alloc_fail()
            return
        end if

        allocate (sd)
        sd%N = N
        sd%stype = stype
        sd%dimx = dimx
        sd%x_ptr = x
        sd%C_ptr = Carr
        sd%ind = int(0.5d0*dimx)

        allocate (sd%BC(0:(N + 1)*(N + 1) - 1))
        call set_bc_array(N, sd%BC)

        allocate (sd%fac(0:(N + 1)*(N + 1) - 1))
        call set_fac_array(N, sd%fac)

        sid = transfer(c_loc(sd), sid)
    end subroutine spline_alloc

    subroutine alloc_fail()
        write (0, '(a)') &
            "error: check of input parameters failed in spline_alloc_."
    end subroutine alloc_fail

    subroutine spline_calc(sid, y, Imin, Imax, W, ierr) bind(C, name="spline_calc_")
        integer(c_intptr_t), value :: sid
        type(c_ptr), value :: y, W
        integer(c_int), value :: Imin, Imax
        integer(c_int), intent(out) :: ierr
        type(spline_data_t), pointer :: sd
        real(c_double), pointer :: yp(:), Wp(:)
        real(c_double), allocatable, target :: Wloc(:)
        integer(c_int) :: dimy, Wn

        call handle_to_sd(sid, sd)
        dimy = Imax - Imin + 1
        Wn = (sd%N + 1)*(2*dimy + sd%dimx*(4 + 3*((sd%N - 1)/2)))

        call cptr0(y, yp, sd%dimx*dimy)
        if (c_associated(W)) then
            call cptr0(W, Wp, Wn)
        else
            allocate (Wloc(0:Wn - 1))
            Wp(0:Wn - 1) => Wloc
        end if

        ierr = 0
        ierr = calc_spline_boundaries(sd, dimy, yp, Wp)
        if (ierr /= 0) write (0, '(a)') &
            "error: spline_calc: calc_spline_boundaries failed!"

        ierr = calc_spline_coefficients(sd, Imin, dimy, yp, Wp)
        if (ierr /= 0) write (0, '(a)') &
            "error: spline_calc: calc_spline_coefficients failed!"
    end subroutine spline_calc

    integer(c_int) function calc_spline_boundaries(sd, dimy, y, W) result(info)
        type(spline_data_t), intent(in) :: sd
        integer(c_int), intent(in) :: dimy
        real(c_double), intent(in) :: y(0:)
        real(c_double), intent(inout) :: W(0:)
        real(c_double), pointer :: x(:)
        integer(c_int) :: N, dimx, D, NRHS, LDA, LDB
        integer(c_int) :: k, p, i, j, idx, boff, am_off
        integer(c_int), allocatable :: ipiv(:)

        N = sd%N
        dimx = sd%dimx
        call cptr0(sd%x_ptr, x, dimx)

        D = N + 1; NRHS = dimy; LDA = D; LDB = D; info = 0
        allocate (ipiv(D))
        am_off = 2*(N + 1)*dimy

        do k = 0, dimx - 1, dimx - 1
            boff = (k/(dimx - 1))*(N + 1)*dimy
            do i = 0, N
                idx = k + i*sign_int(1 - k)
                W(am_off + i) = 1.0d0
                W(am_off + i + N + 1) = x(idx) - x(k)
                do p = 2, N
                    W(am_off + i + p*(N + 1)) = &
                        W(am_off + i + (p - 1)*(N + 1))*W(am_off + i + N + 1)
                end do
                do j = 0, dimy - 1
                    W(boff + i + j*(N + 1)) = y(idx + j*dimx)
                end do
            end do
            call dgesv(D, NRHS, W(am_off:), LDA, ipiv, W(boff:), LDB, info)
            if (info /= 0) write (0, '(a,i0)') &
                "error: calc_spline_boundaries: INFO=", info
        end do
    end function calc_spline_boundaries

    integer(c_int) function calc_spline_coefficients(sd, Imin, dimy, y, W) result(info)
        type(spline_data_t), intent(in) :: sd
        integer(c_int), intent(in) :: Imin, dimy
        real(c_double), intent(in) :: y(0:)
        real(c_double), intent(inout) :: W(0:)
        real(c_double), pointer :: x(:), Cf(:)
        integer(c_int) :: N, dimx, s, KL, KU, LDAB, len, ieqn, iunk, n_, KLKU
        integer(c_int) :: j, k, p, idx, sc_off, moff, D, NRHS, LDB
        integer(c_int), allocatable :: ipiv(:)
        real(c_double) :: dx

        N = sd%N
        dimx = sd%dimx
        call cptr0(sd%x_ptr, x, dimx)
        call cptr0(sd%C_ptr, Cf, (N + 1)*dimx*(Imin + dimy))

        s = (N - 1)/2
        KL = s + 1; KU = s + 1; LDAB = 2*KL + KU + 1
        len = (N + 1)*dimx
        moff = 2*(N + 1)*dimy

        do j = 0, LDAB*len - 1
            W(moff + j) = 0.0d0
        end do

        ieqn = 0
        sc_off = (N + 1)*dimx*Imin
        KLKU = KL + KU

        do n_ = 0, s
            iunk = n_
            W(moff + KLKU + ieqn + iunk*(LDAB - 1)) = 1.0d0
            do j = 0, dimy - 1
                Cf(sc_off + ieqn + j*len) = W(n_ + j*(N + 1))
            end do
            ieqn = ieqn + 1
        end do

        do k = 1, dimx - 2
            do n_ = 0, N - 1
                p = n_
                iunk = (k - 1)*(N + 1) + p
                dx = 1.0d0
                idx = n_*(N + 1)
                W(moff + KLKU + ieqn + iunk*(LDAB - 1)) = sd%BC(p + idx)*dx
                do p = n_ + 1, N
                    iunk = (k - 1)*(N + 1) + p
                    dx = dx*(x(k) - x(k - 1))
                    W(moff + KLKU + ieqn + iunk*(LDAB - 1)) = sd%BC(p + idx)*dx
                end do
                iunk = k*(N + 1) + n_
                W(moff + KLKU + ieqn + iunk*(LDAB - 1)) = -1.0d0
                do j = 0, dimy - 1
                    Cf(sc_off + ieqn + j*len) = 0.0d0
                end do
                ieqn = ieqn + 1
            end do

            iunk = k*(N + 1)
            W(moff + KLKU + ieqn + iunk*(LDAB - 1)) = 1.0d0
            do j = 0, dimy - 1
                Cf(sc_off + ieqn + j*len) = y(k + j*dimx)
            end do
            ieqn = ieqn + 1
        end do

        k = dimx - 1
        do n_ = 0, s
            p = n_
            iunk = (k - 1)*(N + 1) + p
            dx = 1.0d0
            idx = n_*(N + 1)
            W(moff + KLKU + ieqn + iunk*(LDAB - 1)) = sd%BC(p + idx)*dx
            do p = n_ + 1, N
                iunk = (k - 1)*(N + 1) + p
                dx = dx*(x(k) - x(k - 1))
                W(moff + KLKU + ieqn + iunk*(LDAB - 1)) = sd%BC(p + idx)*dx
            end do
            iunk = k*(N + 1) + n_
            W(moff + KLKU + ieqn + iunk*(LDAB - 1)) = -1.0d0
            do j = 0, dimy - 1
                Cf(sc_off + ieqn + j*len) = 0.0d0
            end do
            ieqn = ieqn + 1
        end do

        do n_ = 0, N
            iunk = k*(N + 1) + n_
            W(moff + KLKU + ieqn + iunk*(LDAB - 1)) = 1.0d0
            do j = 0, dimy - 1
                Cf(sc_off + ieqn + j*len) = W((N + 1)*dimy + n_ + j*(N + 1))
            end do
            ieqn = ieqn + 1
        end do

        D = len; NRHS = dimy; LDB = D; info = 0
        allocate (ipiv(D))
        call dgbsv(D, KL, KU, NRHS, W(moff:), LDAB, ipiv, Cf(sc_off:), LDB, info)
        if (info /= 0) write (0, '(a,i0)') &
            "error: calc_spline_coefficients: INFO=", info
    end function calc_spline_coefficients

    subroutine spline_eval(sid, dimz, z, Dmin, Dmax, Imin, Imax, R) &
        bind(C, name="spline_eval_")
        integer(c_intptr_t), value :: sid
        integer(c_int), value :: dimz, Dmin, Dmax, Imin, Imax
        type(c_ptr), value :: z, R
        type(spline_data_t), pointer :: sd
        real(c_double), pointer :: zp(:), Rp(:), x(:), Cf(:), fac(:)
        integer(c_int) :: N, p, n_, j, k, ic0, ic1, ir0, ir1, idx, len, D1, D2
        real(c_double) :: tmp

        call handle_to_sd(sid, sd)
        N = sd%N
        call cptr0(z, zp, dimz)
        call cptr0(R, Rp, dimz*(Dmax - Dmin + 1)*(Imax - Imin + 1))
        call cptr0(sd%x_ptr, x, sd%dimx)
        call cptr0(sd%C_ptr, Cf, (N + 1)*sd%dimx*(Imax + 1))
        fac => sd%fac

        len = (N + 1)*sd%dimx
        D1 = Dmax - Dmin + 1; D2 = D1*(Imax - Imin + 1)

        do k = 0, dimz - 1
            call search_array(zp(k), sd%dimx, x, sd%ind)
            ic0 = sd%ind*(N + 1)
            ir0 = k*D2 - Dmin
            do j = Imin, Imax
                ic1 = ic0 + j*len
                ir1 = ir0 + (j - Imin)*D1
                do n_ = Dmin, Dmax
                    idx = n_*(N + 1)
                    tmp = fac(N + idx)*Cf(N + ic1)
                    do p = N - n_, 1, -1
                        tmp = fac(p + n_ - 1 + idx)*Cf(p + n_ - 1 + ic1) + &
                              (zp(k) - x(sd%ind))*tmp
                    end do
                    Rp(n_ + ir1) = tmp
                end do
            end do
        end do
    end subroutine spline_eval

    subroutine spline_eval_d(sid, dimz, z, Dmin, Dmax, Imin, Imax, R) &
        bind(C, name="spline_eval_d_")
        integer(c_intptr_t), value :: sid
        integer(c_int), value :: dimz, Dmin, Dmax, Imin, Imax
        type(c_ptr), value :: z, R
        type(spline_data_t), pointer :: sd
        real(c_double), pointer :: zp(:), Rp(:), x(:), Cf(:), fac(:)
        integer(c_int) :: N, p, n_, j, k, ic0, ic1, ir0, ir1, idx, len, D1, D2
        real(c_double) :: tmp

        call handle_to_sd(sid, sd)
        N = sd%N
        call cptr0(z, zp, dimz)
        call cptr0(R, Rp, dimz*(Dmax - Dmin + 1)*(Imax - Imin + 1))
        call cptr0(sd%x_ptr, x, sd%dimx)
        call cptr0(sd%C_ptr, Cf, (N + 1)*sd%dimx*(Imax + 1))
        fac => sd%fac

        len = (N + 1)*sd%dimx
        D1 = Imax - Imin + 1; D2 = D1*(Dmax - Dmin + 1)

        do k = 0, dimz - 1
            call search_array(zp(k), sd%dimx, x, sd%ind)
            ic0 = sd%ind*(N + 1)
            ir0 = k*D2 - (Imin + D1*Dmin)
            do j = Imin, Imax
                ic1 = ic0 + j*len
                ir1 = ir0 + j
                do n_ = Dmin, Dmax
                    idx = n_*(N + 1)
                    tmp = fac(N + idx)*Cf(N + ic1)
                    do p = N - n_, 1, -1
                        tmp = fac(p + n_ - 1 + idx)*Cf(p + n_ - 1 + ic1) + &
                              (zp(k) - x(sd%ind))*tmp
                    end do
                    Rp(ir1 + D1*n_) = tmp
                end do
            end do
        end do
    end subroutine spline_eval_d

    subroutine spline_free(sid) bind(C, name="spline_free_")
        integer(c_intptr_t), value :: sid
        type(spline_data_t), pointer :: sd

        if (sid == 0_c_intptr_t) return
        call handle_to_sd(sid, sd)
        if (associated(sd)) deallocate (sd)
    end subroutine spline_free

    subroutine handle_to_sd(sid, sd)
        integer(c_intptr_t), value :: sid
        type(spline_data_t), pointer, intent(out) :: sd
        type(c_ptr) :: cp
        cp = transfer(sid, cp)
        call c_f_pointer(cp, sd)
    end subroutine handle_to_sd

    subroutine cptr0(cp, p, n)
        type(c_ptr), value :: cp
        real(c_double), pointer, intent(out) :: p(:)
        integer(c_int), value :: n
        real(c_double), pointer :: tmp(:)
        call c_f_pointer(cp, tmp, [n])
        p(0:n - 1) => tmp
    end subroutine cptr0

    subroutine search_array(xv, dimx, xa, ind)
        real(c_double), intent(in) :: xv
        integer(c_int), intent(in) :: dimx
        real(c_double), intent(in) :: xa(0:)
        integer(c_int), intent(inout) :: ind
        if (xv < xa(ind)) then
            ind = binary_search(xv, xa, 0, ind)
        else if (xv >= xa(ind + 1)) then
            ind = binary_search(xv, xa, ind, dimx - 1)
        end if
    end subroutine search_array

    integer(c_int) function binary_search(xv, xa, ilo, ihi) result(res)
        real(c_double), intent(in) :: xv
        real(c_double), intent(in) :: xa(0:)
        integer(c_int), intent(in) :: ilo, ihi
        integer(c_int) :: lo, hi, i
        lo = ilo; hi = ihi
        do while (hi > lo + 1)
            i = (hi + lo)/2
            if (xa(i) > xv) then
                hi = i
            else
                lo = i
            end if
        end do
        res = lo
    end function binary_search

    integer(c_int) function sign_int(x) result(res)
        integer(c_int), intent(in) :: x
        if (x < 0) then
            res = -1
        else if (x == 0) then
            res = 0
        else
            res = 1
        end if
    end function sign_int

    subroutine set_bc_array(N, BC)
        integer(c_int), intent(in) :: N
        real(c_double), intent(out) :: BC(0:)
        integer(c_int) :: k, n_
        real(c_double) :: tmp
        do n_ = 0, N
            tmp = 1.0d0
            BC(n_) = tmp
            do k = 1, n_
                tmp = tmp*dble(n_ - k + 1)/dble(k)
                BC(n_ + k*(N + 1)) = tmp
            end do
        end do
    end subroutine set_bc_array

    subroutine set_fac_array(N, fac)
        integer(c_int), intent(in) :: N
        real(c_double), intent(out) :: fac(0:)
        integer(c_int) :: p, n_
        real(c_double) :: tmp
        do p = 0, N
            tmp = 1.0d0
            fac(p) = tmp
            do n_ = 1, p
                tmp = tmp*dble(p - n_ + 1)
                fac(p + n_*(N + 1)) = tmp
            end do
        end do
    end subroutine set_fac_array

end module kilca_spline_m
