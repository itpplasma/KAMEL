program test_hyper

implicit none;

INTERFACE
    COMPLEX*16 FUNCTION CONHYP (A,B,Z,LNCHF,IP)
        INTEGER LNCHF,IP
        COMPLEX*16 A,B,Z,CONHYP
    END FUNCTION
END INTERFACE

real(8), parameter :: pi = 3.141592653589793238462643383279502884197
complex(8), parameter :: im = cmplx(0.0d0, 1.0d0, 8), one = cmplx(1.0d0, 0.0d0, 8)

complex(8) :: z, a, H
complex(8) :: t1, t2;
real(8) :: kp, Vt, nu, omega, omega0, F_re, F_im, dom;

integer :: k, LNCHF = 0, IP = 0;

kp = 1.0d-4;
Vt = 1.0d8;
! omega = 1.0d4;
! omega0 = 1.0d5;

do k=1,100

print *, 'input nu='
read *, nu
print *, 'input dom='
read *, dom


t1 = kp**2 * Vt**2 / nu**2;

!t2 = -im*(omega - omega0)/nu + t1;
t2 = -im*(dom)/nu + t1;

print *, 't1=', t1, 't2=', t2;

!
! print *, 'input z='
! read *, z

H = CONHYP (one, one + t2, t1, LNCHF, IP);

print *
print *, 'CONHYP=', H;

call Hypergeometric1F1_kummer (real(one+t2), aimag(one+t2), real(t1), aimag(t1), F_re, F_im);

print *
print *, 'kummer: Fre=',F_re,'Fim=',F_im;

call Hypergeometric1F1_quad (real(one+t2), aimag(one+t2), real(t1), aimag(t1), F_re, F_im);

print *
print *, 'quad: Fre=',F_re,'Fim=',F_im;

!print *, sqrt(2*pi) * Vt * H /nu /t2;

end do

end program
