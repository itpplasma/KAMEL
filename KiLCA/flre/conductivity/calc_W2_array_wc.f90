!<Calculates array of W2 special functions.

!------------------------------------------------------------------------------

subroutine calc_W2_array ()

use constants, only: dpc, dp, sqrt2p;
use conduct_arrays, only: vE, W2;
use conduct_parameters, only: omega_, omc_, kp_, ks_, nu_, vT_, Vp_, r_, q_;
use flre_sett, only: Nmax;

implicit none;

real(dp), dimension(1:6) :: kp, Vt, nu;
complex(dpc), dimension(1:6) :: dom;
real(dp) :: omega0;

complex(dpc) :: t1, t2, F11m;
complex(8), parameter :: I = cmplx(0.0d0, 1.0d0, 8), one = cmplx(1.0d0, 0.0d0, 8);

integer :: l;

real(dp) :: F_im, F_re;

kp(1) = kp_
Vt(1) = vT_;
nu(1) = nu_;

do l = 2,6
    kp(l) = kp(l-1) * kp(1); 
    Vt(l) = Vt(l-1) * Vt(1);
    nu(l) = nu(l-1) * nu(1);
end do

t1 = cmplx(kp(2)*Vt(2)/nu(2), 0.0d0, 8);

do l = -Nmax,Nmax

    omega0 = dble(l) * omc_ + ks_ * vE + kp_ * Vp_;

    dom(1) = omega_ - omega0;
    dom(2) = dom(1)*dom(1);
    dom(3) = dom(2)*dom(1);
    dom(4) = dom(3)*dom(1);

    t2 = - I*dom(1)/nu(1) + t1;

!      call Hypergeometric1F1_kummer_modified_0_ada (real(one+t2), aimag(one+t2), real(t1), aimag(t1), F_re, F_im);
! 
!      F11m = F_re + I*F_im;

     call hypergeometric1f1_cont_fract_1_modified_0_ada (real(one+t2), aimag(one+t2), real(t1), aimag(t1), F_re, F_im);
! 
     F11m = F_re + I*F_im;

!asymptotic form of expressions: has no fake singularity at kp=0

!complex(dpc), dimension(0:1, 0:3, -Nmax:Nmax) :: W2

W2(0,0,l) =  &
((-I)*sqrt2p*nu(1)*Vt(1)*(-(dom(2)*nu(2)) + 2*nu(4) + 5*kp(2)*nu(2)*Vt(2) -  &
   (3*I)*dom(1)*(nu(3) + kp(2)*nu(1)*Vt(2)) + (3 + F11m)*kp(4)*Vt(4)))/ &
 ((dom(1)*nu(1) + I*kp(2)*Vt(2))*(dom(1)*nu(1) + I*(nu(2) + kp(2)*Vt(2)))* &
  (dom(1)*nu(1) + I*(2*nu(2) + kp(2)*Vt(2))))

W2(0,1,l) =  &
((-I)*sqrt2p*kp(1)*Vt(3)*(dom(1)*(-nu(3) + F11m*kp(2)*nu(1)*Vt(2)) -  &
   I*(2*nu(4) + 3*kp(2)*nu(2)*Vt(2) + kp(4)*Vt(4))))/ &
 ((dom(1)*nu(1) + I*kp(2)*Vt(2))*(dom(1)*nu(1) + I*(nu(2) + kp(2)*Vt(2)))* &
  (dom(1)*nu(1) + I*(2*nu(2) + kp(2)*Vt(2))))

W2(0,2,l) =  &
-((sqrt2p*(dom(1) + I*nu(1))*Vt(3)*(2*nu(4) + 3*kp(2)*nu(2)*Vt(2) -  &
    I*dom(1)*(nu(3) - F11m*kp(2)*nu(1)*Vt(2)) + kp(4)*Vt(4)))/ &
  ((dom(1)*nu(1) + I*kp(2)*Vt(2))*(dom(1)*nu(1) + I*(nu(2) + kp(2)*Vt(2)))* &
   (dom(1)*nu(1) + I*(2*nu(2) + kp(2)*Vt(2)))))

W2(0,3,l) =  &
((-I)*sqrt2p*kp(1)*(F11m*dom(3)*nu(1) - (3 + 2*F11m)*dom(1)*nu(3) +  &
   I*dom(2)*(3*F11m*nu(2) - kp(2)*Vt(2)) -  &
   I*(6*nu(4) + (7 + 2*F11m)*kp(2)*nu(2)*Vt(2) + kp(4)*Vt(4)))*Vt(5))/ &
 ((dom(1)*nu(1) + I*kp(2)*Vt(2))*(dom(1)*nu(1) + I*(nu(2) + kp(2)*Vt(2)))* &
  (dom(1)*nu(1) + I*(2*nu(2) + kp(2)*Vt(2))))

W2(1,0,l) =  &
((-I)*sqrt2p*kp(1)*Vt(3)*(dom(1)*(-nu(3) + F11m*kp(2)*nu(1)*Vt(2)) -  &
   I*(2*nu(4) + 3*kp(2)*nu(2)*Vt(2) + kp(4)*Vt(4))))/ &
 ((dom(1)*nu(1) + I*kp(2)*Vt(2))*(dom(1)*nu(1) + I*(nu(2) + kp(2)*Vt(2)))* &
  (dom(1)*nu(1) + I*(2*nu(2) + kp(2)*Vt(2))))

W2(1,1,l) =  &
((-I)*sqrt2p*dom(1)*Vt(3)*(dom(1)*(-nu(3) + F11m*kp(2)*nu(1)*Vt(2)) -  &
   I*(2*nu(4) + 3*kp(2)*nu(2)*Vt(2) + kp(4)*Vt(4))))/ &
 ((dom(1)*nu(1) + I*kp(2)*Vt(2))*(dom(1)*nu(1) + I*(nu(2) + kp(2)*Vt(2)))* &
  (dom(1)*nu(1) + I*(2*nu(2) + kp(2)*Vt(2))))

W2(1,2,l) =  &
((-I)*sqrt2p*kp(1)*(F11m*dom(3)*nu(1) + I*dom(2)*(F11m*nu(2) - kp(2)*Vt(2)) -  &
   dom(1)*(3*nu(3) + 2*kp(2)*nu(1)*Vt(2)) -  &
   I*(2*nu(4) + 3*kp(2)*nu(2)*Vt(2) + kp(4)*Vt(4)))*Vt(5))/ &
 ((dom(1)*nu(1) + I*kp(2)*Vt(2))*(dom(1)*nu(1) + I*(nu(2) + kp(2)*Vt(2)))* &
  (dom(1)*nu(1) + I*(2*nu(2) + kp(2)*Vt(2))))

W2(1,3,l) =  &
(sqrt2p*dom(1)*((-I)*F11m*dom(3)*nu(1) + I*(3 + 2*F11m)*dom(1)*nu(3) -  &
   6*nu(4) - (7 + 2*F11m)*kp(2)*nu(2)*Vt(2) +  &
   dom(2)*(3*F11m*nu(2) - kp(2)*Vt(2)) - kp(4)*Vt(4))*Vt(5))/ &
 ((dom(1)*nu(1) + I*kp(2)*Vt(2))*(dom(1)*nu(1) + I*(nu(2) + kp(2)*Vt(2)))* &
  (dom(1)*nu(1) + I*(2*nu(2) + kp(2)*Vt(2))))

!usual form of expressions: has fake singularity at kp = 0

! call Hypergeometric1F1_kummer (real(one+t2), aimag(one+t2), real(t1), aimag(t1), F_re, F_im);
!
! F11 = F_re + I*F_im;

!complex(dpc), dimension(0:1, 0:3, -Nmax:Nmax) :: W2

! W2(0,0,l) =  &
! (F11*sqrt2p*Vt(1))/(t2*nu(1))
! 
! W2(0,1,l) =  &
! (sqrt2p*(-I + (F11*dom(1))/(t2*nu(1)))*Vt(1))/kp(1)
! 
! W2(0,2,l) =  &
! (sqrt2p*(F11*dom(2) + I*F11*dom(1)*nu(1) +  &
!    t2*nu(1)*((-I)*(dom(1)) + nu(1)))*Vt(1))/(t2*kp(2)*nu(1))
! 
! W2(0,3,l) =  &
! (sqrt2p*Vt(1)*(F11*dom(3) + I*(3*F11 - t2)*dom(2)*nu(1) +  &
!    (-2*F11 + 3*t2)*dom(1)*nu(2) -  &
!    I*(-2*t2*nu(3) + (2*F11 + t2)*kp(2)*nu(1)*Vt(2))))/(t2*kp(3)*nu(1))
! 
! W2(1,0,l) =  &
! (sqrt2p*(-I + (F11*dom(1))/(t2*nu(1)))*Vt(1))/kp(1)
! 
! W2(1,1,l) =  &
! (sqrt2p*dom(1)*(F11*dom(1) - I*t2*nu(1))*Vt(1))/(t2*kp(2)*nu(1))
! 
! W2(1,2,l) =  &
! (sqrt2p*Vt(1)*((dom(1)*(dom(1) + I*nu(1))*(F11*dom(1) - I*t2*nu(1)))/ &
!     (t2*nu(1)) - I*kp(2)*Vt(2)))/kp(3)
! 
! W2(1,3,l) =  &
! (sqrt2p*dom(1)*Vt(1)*(F11*dom(3) + I*(3*F11 - t2)*dom(2)*nu(1) +  &
!    (-2*F11 + 3*t2)*dom(1)*nu(2) -  &
!    I*(-2*t2*nu(3) + (2*F11 + t2)*kp(2)*nu(1)*Vt(2))))/(t2*kp(4)*nu(1))

!write(999,*) r_, real(F11m), aimag(F11m)

end do

! write(1000,*) r_, real(W2(0,0,:)), aimag(W2(0,0,:))
! write(1001,*) r_, real(W2(0,1,:)), aimag(W2(0,1,:))
! write(1002,*) r_, real(W2(0,2,:)), aimag(W2(0,2,:))
! write(1003,*) r_, real(W2(0,3,:)), aimag(W2(0,3,:))
! write(1011,*) r_, real(W2(1,1,:)), aimag(W2(1,1,:))
! write(1012,*) r_, real(W2(1,2,:)), aimag(W2(1,2,:))
! write(1013,*) r_, real(W2(1,3,:)), aimag(W2(1,3,:))

end subroutine

!------------------------------------------------------------------------------
