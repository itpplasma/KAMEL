!Here is an automatically generated fortran subroutine:

subroutine wm_asympt (zz, res)

implicit none;

complex(8), intent(in) :: zz
complex(8), intent(out) :: res
complex(8) :: x, t1

real(8), parameter :: pi = 3.141592653589793238462643383279502884197;

x = cmplx(1.0d0, 0.0d0, 8)/zz;
t1 = x**2;
res = 3.0d0/4.0d0+3.0d0/4096.0d0*t1*(2560.0d0+t1*(8960.0d0+t1* &
    (40320.0d0+t1*(221760.0d0+ t1*(1441440.0d0+t1*(10810800.0d0+t1* &
    (91891800.0d0+t1*(872972100.0d0+(9166207050.0d0+ &
    105411381075.0d0*t1)*t1))))))));

res = res*cmplx(0.0d0, 1.0d0, 8)/sqrt(pi);

end subroutine
