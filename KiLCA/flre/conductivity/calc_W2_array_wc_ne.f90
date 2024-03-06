!<Calculates array of W2 special functions.

!------------------------------------------------------------------------------

subroutine calc_W2_array (spec)

use constants, only: dpc, dp, sqrt2p;
use conduct_arrays, only: vE, W2;
use conduct_parameters, only: omega_, omc_, kp_, ks_, nu_, vT_, Vp_, r_, q_;
use flre_sett, only: Nmax;

implicit none;

integer, intent(in) :: spec;

real(dp), dimension(1:7) :: Vt;
real(dp) :: omega0;

integer :: l, m, n;

complex(dpc) :: x1, x2;

complex(dpc), dimension(0:3,0:3) :: Imn;

complex(dpc) :: den, Ip;

Vt(1) = vT_;

do l = 2,7
    Vt(l) = Vt(l-1) * Vt(1);
end do

do l = -Nmax,Nmax

    omega0 = dble(l) * omc_ + ks_ * vE + kp_ * Vp_;

    x1 = cmplx(kp_ * vT_ / nu_, 0.0d0, 8);
    x2 = (omega_ - omega0) / nu_;

    call calc_Imn_array(x1, x2, Imn);

    !FP collisions: particles conservation
    !do m = 0,1
    !    do n = 0,3
    !        W2(m,n,l) = sqrt2p / nu_ * Vt(m+n+1) * Imn(m,n);
    !    end do
    !end do

    !FP collisions: particles and energy conservation
    den = cmplx(1.0d0, 0.0d0, 8) - Imn(0,0) + 2.0d0 * Imn(2,0) - Imn(2,2);

    do m = 0,1
        do n = 0,3

            Ip = Imn(m,n) + ((Imn(m,0) - Imn(m,2)) * (Imn(n,0) - Imn(n,2))) / den;

            W2(m,n,l) = sqrt2p / nu_ * Vt(m+n+1) * Ip;

        end do
    end do

end do

! Values for (0,2), (0,3), (1,2), (1,3) are somewhat different for l ~= 0. Probbaly due to roundoff errors
! l = 0;
!
! write(1000,*) r_, real(W2(0,0,l)), aimag(W2(0,0,l))
! write(1000,*) r_, real(W2(0,1,l)), aimag(W2(0,1,l))
! write(1000,*) r_, real(W2(0,2,l)), aimag(W2(0,2,l))
! write(1000,*) r_, real(W2(0,3,l)), aimag(W2(0,3,l))
! write(1000,*) r_, real(W2(1,0,l)), aimag(W2(1,0,l))
! write(1000,*) r_, real(W2(1,1,l)), aimag(W2(1,1,l))
! write(1000,*) r_, real(W2(1,2,l)), aimag(W2(1,2,l))
! write(1000,*) r_, real(W2(1,3,l)), aimag(W2(1,3,l))

end subroutine

!------------------------------------------------------------------------------
