!<The modified W plasma function: fake singularity at kp=0 is removed. Definition: Wm = - z2*(I*srpi*W*z2*z + z2 + 0.5d0).

subroutine wmfunc (kp, z, Wm)

implicit none;

real(8), intent(in) :: kp;
complex(8), intent(in) :: z;
complex(8), intent(out) :: Wm;

real(8), parameter :: pi = 3.141592653589793238462643383279502884197;
real(8), parameter :: srpi = sqrt(pi);
complex(8), parameter :: I = cmplx(0.0d0, 1.0d0, 8), U = cmplx(1.0d0, 0.0d0, 8);

real(8), parameter :: zabs = 1.0d1; !switch to asymptotics

complex(8) :: W, z2, t2, pole, row;
real(8) :: sgn, sigma;

z2 = z*z;

if (abs(z) < zabs) then !use normal definition of Wm

    call wfunc (kp, z, W);

    Wm = - z2*(I*srpi*W*z2*z + z2 + 0.5d0);

else !use asymptotic expansion

    sgn = aimag(kp*z);

    if (sgn > 0.0d0) then

        sigma = 0.0d0;

    else if (sgn == 0.0d0) then

        sigma = sign(1.0d0, kp);

    else

        sigma = 2.0d0*sign(1.0d0, kp);

    end if

    if (sigma /= 0.0d0) then
        pole = - sigma*I*srpi*z*z2*z2*exp(-z2); !pole contribution
    else
        pole = cmplx(0.0d0, 0.0d0, 8);
    end if

    t2 = U/z2;
    row = 3.0d0/4.0d0+3.0d0/4096.0d0*t2*(2560.0d0+t2*(8960.0d0+t2* &
          (40320.0d0+t2*(221760.0d0+ t2*(1441440.0d0+t2*(10810800.0d0+t2* &
          (91891800.0d0+t2*(872972100.0d0+(9166207050.0d0+ &
           105411381075.0d0*t2)*t2))))))));

    Wm = pole + row;
end if

end subroutine
