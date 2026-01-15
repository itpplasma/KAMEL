module muller_root_finding
    ! Finds multiple complex roots using Müller's method with deflation
    use KIM_kinds_m, only: dp

    implicit none

    private

    public :: roots

    contains

    subroutine roots(calcF, known, nmore, km, realrt, ep1, ep2, &
                    guess, maxits, rts, nfound, fnv, ifail)
        ! Main driver: finds nmore roots starting from known roots using deflation
        interface
        subroutine calcF(x, f)
            import dp
            complex(dp), intent(in)  :: x
            complex(dp), intent(out) :: f
        end subroutine calcF
        end interface

        integer, intent(inout) :: known, nmore, km, guess, maxits
        logical, intent(in)    :: realrt
        real(dp), intent(inout) :: ep1, ep2
        complex(dp), intent(inout) :: rts(km), fnv(km)
        integer, intent(out) :: nfound, ifail

        ! Three points for quadratic interpolation
        complex(dp) :: x1, x2, x3, f1, f2, f3
        complex(dp) :: x21, x32, f21, f32, f31, f321
        complex(dp) :: b, xnew, fnew, denom, fsave, radicl, rt
        complex(dp) :: dif, cstep, fslast

        real(dp), parameter :: ep1def = 0.5d-6, ep2def = 1.0d-7
        integer, parameter :: itsdef = 100 ! minimum number of max iterations
        real(dp), parameter :: zero = 0.0d0, half = 0.5d0, one = 1.0d0, &
                                two = 2.0d0, four = 4.0d0, hundrd = 100.0d0

        integer :: i, ii, loop1, loop2, ilo, new, newm1, kount

        ! Initialization
        ! ep1 = max(ep1, ep1def)
        ! ep2 = max(ep2, ep2def)
        maxits = max(maxits, itsdef)
        ifail = 0
        nfound = 0

        ! Verify known roots are actual roots
        if (known >= 1) then
            do i=1, known
                ii = i
                call calcF(rts(i), fnv(i))
                if (abs(fnv(ii)) >= ep2) then
                ifail = 1
                guess = guess + known - ii + 1
                nmore = nmore + known - ii + 1
                known = ii - 1
                exit
                end if
            end do
        end if

        loop1 = known + 1
        loop2 = known + nmore
        if (loop1 > loop2 .or. loop1 < 1) then
        ifail = 2
        return
        end if

        if (guess < nmore) then
        ilo = guess + 1
        do i = ilo, loop2
            rts(i) = (0.0_dp, 0.0_dp)
        end do
        else
        guess = nmore
        end if

        ! Find each requested root
        do new = loop1, loop2
        kount = 3
        newm1 = new - 1
        rt = rts(new)
        x1 = rt - half
        x2 = rt + half
        x3 = rt

        call calcF(x1, f1)
        call calcF(x2, f2)
        call calcF(x3, f3)
        fslast = f3

        if (new > 1) call test(x1, f1, fsave, rts, newm1, ep2, kount, calcF)
        if (new > 1) call test(x2, f2, fsave, rts, newm1, ep2, kount, calcF)
        if (new > 1) call test(x3, f3, fslast, rts, newm1, ep2, kount, calcF)

        f21 = (f2-f1)/(x2-x1)

        ! Müller iteration
        do
            ! Quadratic interpolation coefficients
            x32 = x3 - x2
            f32 = (f3-f2)/x32
            f321 = (f32-f21)/(x3-x1)
            b = f32 + x32*f321
            radicl = b*b - four*f321*f3
            if (realrt .and. real(radicl) < zero) radicl = (0.0_dp, 0.0_dp)
            radicl = sqrt(radicl)
            ! Choose sign for maximum denominator
            if (real(b)*real(radicl) + aimag(b)*aimag(radicl) < zero) radicl = -radicl
            denom = b + radicl

            if (abs(denom) /= zero) then
            cstep = two*f3/denom
            if (realrt .and. abs(f3) /= zero .and. abs(f32) /= zero) then
                ! Keep step real for real root mode
                cstep = f32/abs(f32) * f3/abs(f3) * abs(cstep)
            end if
            xnew = x3 - cstep
            else
            if (abs(f3) < ep2) then
                xnew = x3
            else
                xnew = x3 + x32
            end if
            end if

            call calcF(xnew, fnew)
            fsave = fnew
            if (new > 1) call test(xnew, fnew, fsave, rts, newm1, ep2, kount, calcF)

            kount = kount + 1

            ! Convergence check
            if (kount > maxits) then
            rts(new) = xnew
            fnv(new) = fsave
            ifail = 3
            return
            end if

            dif = xnew - x3
            if (abs(dif) < ep1*abs(xnew)) exit
            if (abs(fsave) < ep2) exit
            if (.not. realrt .and. abs(fsave) >= hundrd*abs(fslast)) then
            ! Diverging - reduce step
            cstep = cstep*half
            xnew = xnew + cstep
            call calcF(xnew, fnew)
            fsave = fnew
            else
            x1 = x2; x2 = x3; x3 = xnew
            f1 = f2; f2 = f3; f3 = fnew
            fslast = fsave
            f21 = f32
            end if
        end do

        ! Store converged root
        rts(new) = xnew
        fnv(new) = fsave
        nfound = new

        end do

    end subroutine roots

    subroutine test(x, f, fsave, rts, k, eps, kount, calcF)
        ! Deflation: divide f by (x-r) for each known root r
        interface
        subroutine calcF(x, f)
            import dp
            complex(dp), intent(in)  :: x
            complex(dp), intent(out) :: f
        end subroutine calcF
        end interface

        complex(dp), intent(inout) :: x, f, fsave
        complex(dp), intent(in)    :: rts(:)
        integer, intent(in)        :: k
        real(dp), intent(in)       :: eps
        integer, intent(inout)     :: kount

        complex(dp) :: d
        real(dp), parameter :: pertb = 0.01d0
        integer :: i

        !do
        do i=1,k

            d = x - rts(i)

            if (abs(d) <= eps) then
            ! Too close to known root - perturb
                x = x + pertb
                call calcF(x, f)
                fsave = f
                kount = kount + 1
                cycle
            end if

            f = f/d

        end do
        !exit
        !end do

    end subroutine test

end module muller_root_finding

