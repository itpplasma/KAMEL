!<The module containes basic arrays for evaluation of sigma (k) matrices.

module conduct_arrays

use constants, only: dpc, dp;

implicit none;

real(dp) :: vE, dvE, dvEz, dvEt; !depend on a sort of particles

integer, dimension(1:3) :: dimF;
integer, dimension(1:3,1:6) :: Fs1, Fs2;
complex(dpc), allocatable, dimension(:,:,:) :: F; !r,theta,z; ind; l;
complex(dpc), allocatable, dimension(:,:,:) :: F1; !r,theta,z; ind; l;

real(dp), allocatable, dimension(:) :: L;
real(dp), allocatable, dimension(:,:) :: bico;

real(dp), allocatable, dimension(:,:,:) :: A;
integer, dimension(1:3,1:3) :: s1, s2;

integer, dimension(1:3) :: dimI;
integer, dimension(1:3,1:2) :: Is1;
real(dp), dimension(1:3,1:2) :: I;

real(dp), allocatable, dimension(:) :: factor;

complex(dpc), allocatable, dimension(:,:,:) :: W2; !mu; nu; l;
complex(dpc), allocatable, dimension(:,:,:,:,:,:) :: D; !m1; n1; m2; n2; b; l;

real(dp), allocatable, dimension(:, :, :) :: dEM, dJMI;

end module

!------------------------------------------------------------------------------

subroutine allocate_and_set_conductivity_arrays ()

use flre_sett, only: Nmax, flre_order;
use conduct_arrays, only: dimI, Is1;
use conduct_arrays, only: dimF, Fs1, Fs2, F, F1;
use conduct_arrays, only: L, A, s1, s2, bico;
use conduct_arrays, only: W2, D;
use conduct_arrays, only: factor, dEM, dJMI;
use constants, only: dp;

implicit none;

interface
    real(8) function aco (m, k)
        integer, intent(in) :: m, k;
    end function
end interface

integer :: lc;

allocate (L(-Nmax:Nmax));

do lc = 0,Nmax
    L(lc)  =   lc;
    L(-lc) = - lc;
end do

!dimensions:
dimF(1) = 4; !phi;
dimF(2) = 6; !theta;
dimF(3) = 6; !z;

!1 & 2 indices - powers of w1 and J: for (r(1),theta(2),z(3))
Fs1(:,1) = 0; Fs2(:,1) = 0; !{mu=0, nu=0}
Fs1(:,2) = 1; Fs2(:,2) = 0; !{mu=1, nu=0}
Fs1(:,3) = 0; Fs2(:,3) = 1; !{mu=0, nu=1}
Fs1(:,4) = 2; Fs2(:,4) = 0; !{mu=2, nu=0}
Fs1(:,5) = 1; Fs2(:,5) = 1; !{mu=1, nu=1}
Fs1(:,6) = 3; Fs2(:,6) = 0; !{mu=3, nu=0}

allocate (F(1:3,1:6,-Nmax:Nmax)); !r,theta,z; index, cycl. harm.
allocate (F1(1:3,1:6,-Nmax:Nmax)); !r,theta,z; index, cycl. harm.

dimI(1) = 1; !phi;
dimI(2) = 2; !theta
dimI(3) = 2; !z

Is1(1,1) = 0;
Is1(2,1) = 0; Is1(2,2) = 1;
Is1(3,1) = 0; Is1(3,2) = 1;

allocate (A(0:flre_order,1:3,1:3));

s1 = 0; s1(1,1) = 1;
s2 = 0; s2(2,1) = 1; s2(3,1) = 1;

!factorial:
allocate (factor(0:flre_order));

factor(0) = 1;
do lc = 1,flre_order
    factor(lc) = lc * factor(lc-1);
end do

allocate (D(0:1, 0:flre_order+1, 0:1, 0:flre_order+1, 0:1, -Nmax:Nmax));!m1;n1; m2;n2; b;l;
allocate (W2(0:1, 0:3, -Nmax:Nmax)); !mu; nu; l;

allocate (bico(0:flre_order,0:flre_order));

bico = 0.0d0;

call binomial_coefficients (%val(flre_order), bico);

! print *
! print *, 'check aco and bico:'
! print *
! print *, bico(5,:)
! print *, '34',aco(2,4)

allocate (dEM(1:3, 1:3, 0:flre_order), dJMI(1:3, 1:3, 0:flre_order));

end subroutine

!------------------------------------------------------------------------------

subroutine deallocate_conductivity_arrays ()

use conduct_arrays, only: F, F1, L, bico, A, factor, W2, D, dEM, dJMI;

deallocate (F, F1, L, bico, A, factor, W2, D, dEM, dJMI);

end subroutine

!------------------------------------------------------------------------------

subroutine eval_electric_drift_velocities ()
!only to be called after all apropriate back subs

use constants, only: dp, c;
use conduct_parameters, only: r_, ht_, dht_, hz_, dhz_, B, dB, dPhi0_, ddPhi0_;
use conduct_arrays, only: vE, dvE, dvEz, dvEt;

implicit none;

vE  = c/B*dPhi0_;
dvE = c/B*(ddPhi0_ - dB/B*dPhi0_);

dvEz = - (dht_ * vE + ht_ * dvE);                !vEz = - ht_ * vE
dvEt = (hz_ + r_ * dhz_) * vE + r_ * hz_ * dvE;  !vEt = r * hz_ * vE - covariant

end subroutine

!------------------------------------------------------------------------------

subroutine eval_FGI_arrays () !only to be called after all apropriate back subs

use constants, only: dp, dpc;
use conduct_parameters, only: r_, m_, ht_, dht_, hz_, dhz_;
use conduct_parameters, only: kp_, dkp_, ks_, dks_, omega_;
use conduct_parameters, only: n_, dn_, vT_, dvT_, Vp_, dVp_, omc_, domc_;
use conduct_arrays, only: vE, dvE, dvEz, dvEt;
use conduct_arrays, only: F, F1, I, L;
use flre_sett, only: Nmax;

implicit none;

real(dp) :: v1, v2;
real(dp) :: h_th, dh_th;          !covariant components
real(dp) :: hth, hz, Vth, Vz;     !contravariant components
real(dp) :: dVt, dVz;
real(dp) :: rks, drks, rkp, drkp, drkV, drksvT;

real(dp) :: NaN = 1.0d200, J0, oms;

real(dp), dimension(0:2,0:1) :: P;
real(dp), dimension(1:3,0:1,0:0) :: Q;
real(dp), dimension(1:3,0:2,0:1) :: Z;
complex(dpc), dimension(0:1,0:0) :: R;

complex(dpc), dimension(-Nmax:Nmax) :: P00, R00;

integer :: k;

!useful variables:
v1 = (- dn_ / n_ + 3.0d0 * dvT_ / vT_);
v2 = (domc_ / omc_ - 2.0d0 * dvT_ / vT_);

!covariant:
h_th  = r_ * ht_;
dh_th = ht_ + r_ * dht_;

!contravariant:
hth = ht_ / r_;
hz = hz_;

!contravariant:
Vth = hz / r_ * vE + hth * Vp_;
Vz  = - ht_ * vE + hz * Vp_;

dVt = dh_th * Vp_ + h_th * dVp_ + dvEt;
dVz = dhz_ * Vp_ + hz_ * dVp_ + dvEz;

rks = r_ * ks_;
drks = ks_ + r_ * dks_;
rkp = r_ * kp_;
drkp = kp_ + r_ * dkp_;

drkV = drks * Vp_ + rks * dVp_ - (drkp * vE + rkp * dvE);
drksvT = drks * vT_ + rks * dvT_;

!Q - array:
Q(1,0,0) = omc_;
Q(1,1,0) = 0.0d0;

Q(2,0,0) = Vth;
Q(2,1,0) = hth;

Q(3,0,0) = Vz;
Q(3,1,0) = hz;

!P array: without l*omc term!
P(0,0) =   ks_ * vT_ *vT_ * v1 / omc_;
P(0,1) =   ks_ * v2 / m_;
P(1,0) =   (rkp * omc_ - drkV) / r_ / omc_;
P(2,0) = - (drksvT / vT_) / r_ / omc_;

!R array: without l*omc term!
R(0,0) = ks_ * vE + kp_ * Vp_ - omega_;
R(1,0) = kp_;

!Z array:
!phi:
Z(1,0,0) = omc_;
Z(1,1,0) = 0.0d0;
Z(1,0,1) = 0.0d0;
Z(1,2,0) = 0.0d0;

!theta:
Z(2,0,0) =   hz_ * vT_ * vT_ * v1 / r_ / omc_;
Z(2,1,0) =   (omc_ * ht_ - dVz) / r_ / omc_;
Z(2,0,1) =   hz_ * v2 / r_ / m_;
Z(2,2,0) = - (hz_ * dvT_ + dhz_ * vT_) / vT_ / r_ / omc_;

!z:
Z(3,0,0) = - h_th * vT_ * vT_ * v1 / r_ / omc_;
Z(3,1,0) =   (r_ * omc_ * hz_ + dVt) / r_ / omc_;
Z(3,0,1) = - h_th * v2 / r_ / m_;
Z(3,2,0) =   (h_th * dvT_ + dh_th * vT_) / vT_ / r_ / omc_;

P00 = P(0,0) + L*omc_;
R00 = R(0,0) + L*omc_;

!F coefficients for phi, theta, z
do k = 1,3
    F(k,1,:) = P00*Q(k,0,0) - R00*Z(k,0,0);
    F(k,2,:) = P(1,0)*Q(k,0,0) + P00*Q(k,1,0) - R00*Z(k,1,0) - R(1,0)*Z(k,0,0);
    F(k,3,:) = P(0,1)*Q(k,0,0) - R00*Z(k,0,1);
    F(k,4,:) = P(2,0)*Q(k,0,0) + P(1,0)*Q(k,1,0) - R00*Z(k,2,0) - R(1,0)*Z(k,1,0);
    F(k,5,:) = P(0,1)*Q(k,1,0) - R(1,0)*Z(k,0,1);
    F(k,6,:) = P(2,0)*Q(k,1,0) - R(1,0)*Z(k,2,0);
end do

!redefinition of F for phi component to remove big numbers substraction:
oms = ks_ * (vE - vT_ * vT_ * v1 / omc_) + kp_ * Vp_;
F(1,1,:) =   omc_ * (omega_ - oms);
F(1,2,:) = - drkV / r_;
F(1,3,:) =   omc_ * P(0,1);
F(1,4,:) = - drksvT / vT_ / r_;
F(1,5,:) =   NaN;
F(1,6,:) =   NaN;

!F1 coefficients for phi, theta, z - without adiabatic contribution:
do k = 1,3
    F1(k,1,:) = P00*Q(k,0,0);
    F1(k,2,:) = P(1,0)*Q(k,0,0) + P00*Q(k,1,0);
    F1(k,3,:) = P(0,1)*Q(k,0,0);
    F1(k,4,:) = P(2,0)*Q(k,0,0) + P(1,0)*Q(k,1,0);
    F1(k,5,:) = P(0,1)*Q(k,1,0);
    F1(k,6,:) = P(2,0)*Q(k,1,0);
end do

!I - array:
J0 = m_ * m_ * r_ * omc_;

I(1,1) = J0 * Q(1,0,0);
I(1,2) = NaN; !should not be used!

I(2,1) = J0 * Q(2,0,0);
I(2,2) = J0 * Q(2,1,0);

I(3,1) = J0 * Q(3,0,0);
I(3,2) = J0 * Q(3,1,0);

! write(100,*) r_, real(F), aimag(F)
! write(200,*) r_, I

end subroutine

!------------------------------------------------------------------------------

subroutine eval_A_matrix () !only to be called after all apropriate back subs

use constants, only: dp;
use flre_sett, only: flre_order;
use conduct_parameters, only: r_, ht_, hz_;
use conduct_arrays, only: A;

implicit none;

A(:,1,1) = - 1.0d0;
A(flre_order,1,1) = 0.0d0;
A(:,1,2) =   0.0d0;
A(:,1,3) =   0.0d0;

A(:,2,1) = - hz_ / r_;
A(:,2,2) =   1.0d0;
A(:,2,3) =   0.0d0;

A(:,3,1) =   ht_;
A(:,3,2) =   0.0d0;
A(:,3,3) =   1.0d0;

end subroutine

!------------------------------------------------------------------------------

subroutine print_conduct_arrays ()

use conduct_arrays;
use flre_sett, only: Nmax

implicit none;

print *
print *, 'conductivity arrays:'
print *, 'vE:', vE;
print *, 'dvE:', dvE;
print *, 'dvEz:', dvEz;
print *, 'dvEt:', dvEt;
print *, 'F:', F;
print *, 'I:', I;
print *, 'A:', A;
print *, 'W2:', W2;
print *, 'D:', D;

call flush(101)

end subroutine

!------------------------------------------------------------------------------

real(8) function aco (m, k)

use constants, only: dp;
use conduct_arrays, only: bico;
use flre_sett, only: flre_order;

implicit none;

integer, intent(in) :: m, k;

if (m<-1 .or. m>flre_order-1 .or. k<0 .or. k>m+1) then
    print *, 'warning: illegal index passed to aco function:', m, k;
    aco = 0.0d0;
    return;
end if

if (k==0 .or. m==0 .or. k==m+1) then
    aco = 1.0d0;
    return;
end if

aco = bico(m,k-1) + bico(m,k);

end function

!------------------------------------------------------------------------------

subroutine calc_dEM_dJMI_arrays (r)

use constants, only: dp;
use conduct_arrays, only: dEM, dJMI;
use flre_sett, only: flre_order;

implicit none;

real(dp), intent(in) :: r;

call calc_dEM_dJMI_matrices (r, flre_order, dEM, dJMI);

end subroutine

!------------------------------------------------------------------------------
