!<Calculates array of W2 special functions in collisionless case.

!------------------------------------------------------------------------------

subroutine calc_W2_array (spec)

use constants, only: dpc, dp, sr2, srpi;
use conduct_arrays, only: vE, W2;
use conduct_parameters, only: omega_, omc_, ks_, kp_, vT_, Vp_;
use flre_sett, only: Nmax;

implicit none;

integer, intent(in) :: spec;

real(dp), dimension(1:5) :: Vt;
real(dp) :: omega0;

complex(8), parameter :: I = cmplx(0.0d0, 1.0d0, 8), U = cmplx(1.0d0, 0.0d0, 8);
complex(8), parameter :: isrpi = I*srpi;

complex(dpc) :: dom, zv, Wm;
complex(dpc), dimension(-5:-1) :: z;

integer :: l;

Vt(1) = vT_;
Vt(2) = Vt(1) * Vt(1);
Vt(3) = Vt(2) * Vt(1);
Vt(4) = Vt(3) * Vt(1);
Vt(5) = Vt(4) * Vt(1);

do l = -Nmax,Nmax

    omega0 = dble(l) * omc_ + ks_ * vE + kp_ * Vp_;

    dom = omega_ - omega0;

    zv = dom / (sr2 * kp_* vT_);

    z(-1) = U / zv;
    z(-2) = z(-1) / zv;
    z(-3) = z(-2) / zv;
    z(-4) = z(-3) / zv;
    z(-5) = z(-4) / zv;

    !W -> I/Sqrt[Pi](1/z + 1/2/z^3 + Wm/z^5)
    !Wm = - (I * Sqrt[Pi] * Wv + z(-1) + 0.5d0 * z(-3)) / z(-5);
    call wmfunc (kp_, zv, Wm);

!complex(dpc), dimension(0:1, 0:3, -Nmax:Nmax) :: W2

!begin collisionless model:
W2(0,0,l) = (isrpi/sr2)*Vt(1)*(2 + 2*Wm*z(-4) + z(-2))/dom

W2(0,1,l) = isrpi*Vt(2)*(2*Wm*z(-3) + z(-1))/dom

W2(0,2,l) = isrpi*sr2*Vt(3)*(1 + 2*Wm*z(-2))/dom

W2(0,3,l) = 4.0d0*isrpi*Wm*Vt(4)*z(-1)/dom

W2(1,0,l) = W2(0,1,l)

W2(1,1,l) = W2(0,2,l)

W2(1,2,l) = W2(0,3,l)

W2(1,3,l) = 4.0d0*isrpi*sr2*Wm*Vt(5)/dom
!end collisionless model:

end do

end subroutine

!------------------------------------------------------------------------------
