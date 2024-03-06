program test_wfunc

implicit none;

external :: wfunc

real(8), parameter :: pi = 3.141592653589793238462643383279502884197
complex(8), parameter :: im = (0.0d0, 1.0d0)

complex(8) :: z, W, x, t1, Wexp, WN

real(8) :: kp = 1.0d0

integer :: k;

do k=1,100

print *, 'input z='
read *, z

!print *, 'input kp='
!read *, kp

kp = 1.0;

call wfunc (kp, z, W)

!expansion:
    x = 1.0d0/z;
    t1 = x**2;

    Wexp = x*(4096.0d0+t1*(2048.0d0+t1*(3072.0d0+t1*(7680.0d0+t1*(26880.0d0+ &
            t1*(120960.0d0+t1*(665280.0d0+t1*(4324320.0d0+t1*(32432400.0d0+ &
            t1*(275675400.0d0+t1*(2618916300.0d0+ &
            t1*(27498621150.0d0+316234143225.0d0*t1))))))))))))/4096.0d0;

    Wexp = Wexp*im/sqrt(pi);

    if (sign(1.0d0, kp)*aimag(z) < 0.0d0) then
        Wexp = Wexp + 2.0d0*exp(-z**2)*sign(1.0d0, kp);
    else if (aimag(z) == 0.d0) then
        Wexp = Wexp + exp(-z**2)*sign(1.0d0, kp);
    end if


!Nikolays:
call wofz2(z, WN);

print *, 'z =', z
print *, 'kp =', kp
print *, 'WE =',Wexp
print *, 'WI =',W
print *, 'WN =',WN
print *, 'err_rel =', abs(WN-W), 'err_abs =', abs(WN-W)/abs(W);

print *

end do

end program
