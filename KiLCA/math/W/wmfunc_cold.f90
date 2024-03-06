subroutine wmfunc (kp, z, Wm)

implicit none;

real(8), intent(in) :: kp;
complex(8), intent(in) :: z;
complex(8), intent(out) :: Wm;

Wm = cmplx(0.0d0, 0.0d0, 8);

end subroutine
