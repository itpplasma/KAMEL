!<W plasma function: for huge z asymptotic expansion is used.

subroutine wfunc (kp, z, W)

implicit none;

real(8), intent(in) :: kp;
complex(8), intent(in) :: z;
complex(8), intent(out) :: W;

real(8), parameter :: pi = 3.141592653589793238462643383279502884197;
real(8), parameter :: srpi = sqrt(pi);
complex(8), parameter :: I = cmplx(0.0d0, 1.0d0, 8), U = cmplx(1.0d0, 0.0d0, 8);

real(8), parameter :: zabs = 1.0d5 !1.0d1; !switch to asymptotics

real(8) :: zr, zi, Wr, Wi, sgn;
logical :: flag;

complex(8) :: x, t;

sgn = sign(1.0d0, kp);

if (abs(z) < zabs) then !expansion works fine but less accurate
    zr = sgn * real (z);
    zi = sgn * aimag(z);

    call wofz (zr, zi, Wr, Wi, flag);
    if (flag .eqv. .true.) then
        if (-real(z**2) > 708.503061461606d0) then
            print *, 'warning: wfunc: exp overflow!'
        end if
        print *, 'wfunc: failed to compute W function for z =', z
        print *, 'kp =', kp, 'W = (', sgn*Wr, ',', sgn*Wi, ')'
    end if

    W = sgn*(Wr + I*Wi);

else !expansion works fine but less accurate (it is not clear how accuracy depends on z)
    x = U/z;
    t = x**2;

    W = x*(4096.0d0+t*(2048.0d0+t*(3072.0d0+t*(7680.0d0+t*(26880.0d0+ &
            t*(120960.0d0+t*(665280.0d0+t*(4324320.0d0+t*(32432400.0d0+ &
            t*(275675400.0d0+t*(2618916300.0d0+ &
            t*(27498621150.0d0+316234143225.0d0*t))))))))))))/4096.0d0;

    W = W*I/srpi;

    if (sgn*aimag(z) < 0.0d0) then

        W = W + 2.0d0*sgn*exp(-z**2);

    else if (aimag(z) == 0.d0) then

        W = W + 1.0d0*sgn*exp(-z**2);

    end if

    if (-real(z**2) > 708.503061461606d0) then
        print *, 'warning: wfunc: exp overflow:'
        print *, 'z=', z, 'W=', W
    end if

end if

end subroutine
