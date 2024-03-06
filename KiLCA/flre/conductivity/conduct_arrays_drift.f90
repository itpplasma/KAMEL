!<The module containes basic arrays for the evaluation of sigma (k) matrices.

module conduct_arrays

use constants, only: dp, dpc;

implicit none;

real(dp) :: vE, dvE, dvEz, dvEt; !depend on a sort of particles

integer, parameter :: dimF = 12;
integer, dimension(dimF,2) :: FI;
complex(dpc), allocatable, dimension(:,:,:) :: F; !r,theta,z; ind; l;
complex(dpc), allocatable, dimension(:,:,:) :: F1; !r,theta,z; ind; l;

real(dp), allocatable, dimension(:) :: L;
real(dp), allocatable, dimension(:,:) :: bico;

real(dp), allocatable, dimension(:,:,:) :: A;
integer, dimension(1:3,1:3) :: s1, s2;

real(dp), allocatable, dimension(:) :: factor;

complex(dpc), allocatable, dimension(:,:,:,:)   :: D;     !m; n; b; l;
complex(dpc), allocatable, dimension(:,:,:,:,:) :: Imnkl; !mu; nu; k; l; lc;

real(dp), allocatable, dimension(:, :, :) :: dEM, dJMI;

real(dp),     dimension(0:3,0:2)     :: P;
real(dp),     dimension(1:3,0:2,0:2) :: Q;
real(dp),     dimension(1:3,0:3,0:2) :: Z;
complex(dpc), dimension(0:2,0:2)     :: R;

complex(dpc), allocatable, dimension(:) :: P00, R00;

real(dp), allocatable, dimension(:) :: x;

integer :: mpara, mperp;

complex(dpc), allocatable, dimension(:,:,:,:) :: SImnkl; ! for check

end module

!------------------------------------------------------------------------------

subroutine allocate_and_set_conductivity_arrays ()

use flre_sett, only: Nmax, Nbmax, flre_order;
use conduct_arrays, only: dimF, FI, F, F1;
use conduct_arrays, only: L, A, s1, s2, bico;
use conduct_arrays, only: D, Imnkl;
use conduct_arrays, only: factor, dEM, dJMI;
use conduct_arrays, only: P, Q, Z, R, P00, R00, x;
use conduct_arrays, only: mpara, mperp;
use conduct_arrays, only: SImnkl;
use constants, only: dp;

implicit none;

interface
    real(8) function aco (m, k)
        integer, intent(in) :: m, k;
    end function
end interface

integer :: lc, Nf, lperp;

real(dp) :: NaN = 1.0d100;

allocate (L(-Nmax:Nmax));

do lc = 0,Nmax
    L(lc)  =   lc;
    L(-lc) = - lc;
end do

!1 & 2 indices - power of u and power/2 for lambda:
FI(1,1)  = 0;  FI(1,2)  = 0; !{mu=0, nu=0}: F00
FI(2,1)  = 0;  FI(2,2)  = 1; !{mu=0, nu=1}: F02
FI(3,1)  = 0;  FI(3,2)  = 2; !{mu=0, nu=2}: F04
FI(4,1)  = 1;  FI(4,2)  = 0; !{mu=1, nu=0}: F10
FI(5,1)  = 1;  FI(5,2)  = 1; !{mu=1, nu=1}: F12
FI(6,1)  = 1;  FI(6,2)  = 2; !{mu=1, nu=2}: F14
FI(7,1)  = 2;  FI(7,2)  = 0; !{mu=2, nu=0}: F20
FI(8,1)  = 2;  FI(8,2)  = 1; !{mu=2, nu=1}: F22
FI(9,1)  = 3;  FI(9,2)  = 0; !{mu=3, nu=0}: F30
FI(10,1) = 3;  FI(10,2) = 1; !{mu=3, nu=1}: F32
FI(11,1) = 4;  FI(11,2) = 0; !{mu=4, nu=0}: F40
FI(12,1) = 5;  FI(12,2) = 0; !{mu=5, nu=0}: F50

allocate (F(1:3,1:dimF,-Nmax:Nmax)); !r,theta,z; index, cycl. harm.
allocate (F1(1:3,1:dimF,-Nmax:Nmax)); !r,theta,z; index, cycl. harm.

F = NaN;
F1 = NaN;

allocate (A(0:flre_order,1:3,1:3));

A = NaN;

s1 = 0; s1(1,1) = 1;
s2 = 0; s2(2,1) = 1; s2(3,1) = 1;

!factorial:
Nf = max(flre_order, 2+Nmax+2*Nbmax);
allocate (factor(0:Nf));

factor(0) = 1;
do lc = 1,Nf
    factor(lc) = lc * factor(lc-1);
end do

mpara = 5;

mperp = Nmax + 2*Nbmax + 4;

lperp = 2 * mperp + 1;

allocate (D(0:1, 0:flre_order+1, -Nmax:Nmax, 0:Nbmax)); !m; n; lc; b;

allocate (Imnkl(0:mpara, 0:mpara, 0:mperp, 0:mperp, -Nmax:Nmax)); !mu; nu; k; l; lc;

D = NaN;

Imnkl = NaN;

allocate (bico(0:flre_order,0:flre_order));

bico = 0.0d0;

call binomial_coefficients (%val(flre_order), bico);

allocate (dEM(1:3, 1:3, 0:flre_order), dJMI(1:3, 1:3, 0:flre_order));

dEM = NaN;
dJMI = NaN;

P = NaN;
Q = NaN;
Z = NaN;
R = NaN;

allocate (P00(-Nmax:Nmax), R00(-Nmax:Nmax));

P00 = NaN;
R00 = NaN;

allocate (x(0:Nmax+2*Nbmax));

allocate(SImnkl(0:mpara, 0:mpara, 0:mperp, 0:mperp));

call init_spec_data(mpara, lperp);

end subroutine

!------------------------------------------------------------------------------

subroutine deallocate_conductivity_arrays ()

use conduct_arrays, only: F, F1, L, bico, A, factor, Imnkl, D, dEM, dJMI, P00, R00, x, SImnkl;

deallocate (F, F1, L, bico, A, factor, Imnkl, D, dEM, dJMI, P00, R00, x, SImnkl);

call close_spec_data();

end subroutine

!------------------------------------------------------------------------------

subroutine eval_electric_drift_velocities ()
!only to be called after all apropriate back subs

use constants, only: c;
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

use constants, only: dp;
use conduct_parameters, only: r_, m_, ht_, dht_, hz_, dhz_;
use conduct_parameters, only: kp_, dkp_, ks_, dks_, omega_;
use conduct_parameters, only: n_, dn_, vT_, dvT_, Vp_, dVp_, omc_, domc_;
use conduct_arrays, only: vE, dvE, dvEz, dvEt;
use conduct_arrays, only: F, F1, L;
use conduct_arrays, only: P, Q, Z, R, P00, R00;

implicit none;

real(dp) :: v1, v2;
real(dp) :: h_th, dh_th;                 !covariant components
real(dp) :: hth, hz, dhth, dhz, Vth, Vz; !contravariant components
real(dp) :: dVt, dVz;                    !covariant
real(dp) :: dVth;                        !contravariant
real(dp) :: rks, drks, rkp, drkp, drkV, drksvT;

! aux quantities:
real(dp) :: dr2kpVth2_omc, dhthr2vTkpVth_omc, dhth2r2vT2kp_omc;
real(dp) :: dhthr2Vth2_omc, dhth2r2vTVth_omc, dhthr2hth2vT2_omc;
real(dp) :: dhzr2Vth2_omc, dhzhthr2vTVth_omc, dhzr2hth2vT2_omc;

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
dhth = (dht_ - hth) / r_;
dhz = dhz_;

!contravariant:
Vth = hz / r_ * vE + hth * Vp_;
Vz  = - ht_ * vE + hz * Vp_;

!covariant:
dVt = dh_th * Vp_ + h_th * dVp_ + dvEt;
dVz = dhz_ * Vp_ + hz_ * dVp_ + dvEz;

!contravariant:
dVth = (dVt / r_ - 2.0d0 * Vth) / r_;

rks = r_ * ks_;
drks = ks_ + r_ * dks_;
rkp = r_ * kp_;
drkp = kp_ + r_ * dkp_;

drkV = drks * Vp_ + rks * dVp_ - (drkp * vE + rkp * dvE);
drksvT = drks * vT_ + rks * dvT_;

!(r0^2 kp Vth^2 / omc)':
dr2kpVth2_omc = (2.0d0 * Vth + 2.0d0 * r_ * dVth - r_ * Vth * (domc_ / omc_)) * r_ * Vth / omc_ * kp_ + &
                r_ * r_ * Vth * Vth / omc_ * dkp_;

!(hth r0^2 vT kp Vth / omc)':
dhthr2vTkpVth_omc = r_ / omc_ * (hth * r_ * vT_ * Vth * dkp_ + kp_ * (dhth * r_ * vT_ * Vth + &
                    hth * 2.0d0 * vT_ * Vth + hth * r_ * dvT_ * Vth + hth * r_ * vT_ * dVth - &
                    hth * r_ * vT_ * Vth * (domc_ / omc_)));

!(hth^2 r0^2 vT^2 kp / omc)':
dhth2r2vT2kp_omc = r_ * hth / omc_ * (hth * r_ * vT_ * vT_ * dkp_ + kp_ * (2.0d0 * dhth * r_ * vT_ * vT_ + &
                   hth * 2.0d0 * vT_ * vT_ + hth * r_ * 2.0d0 * vT_ * dvT_ - (domc_ / omc_) * hth * r_ * vT_ * vT_));

!(hth r0^2 Vth^2 /omc)':
dhthr2Vth2_omc = r_ * Vth / omc_ * (dhth * r_ * Vth + hth * 2.0d0 * Vth + hth * r_ * 2.0d0 * dVth - &
                 (domc_ / omc_) * hth * r_ * Vth);

!(hth^2 r0^2 vT Vth / omc)':
dhth2r2vTVth_omc = hth * r_ / omc_ * (2.0d0 * dhth * r_ * vT_ * Vth + hth * 2.0d0 * vT_ * Vth + &
                   hth * r_ * dvT_ * Vth + hth * r_ * vT_ * dVth - (domc_ / omc_) * hth * r_ * vT_ * Vth);

!(hth r0^2 hth^2 vT^2 / omc)':
dhthr2hth2vT2_omc = r_ * hth * hth / omc_ * (2.0d0 * hth * vT_ * vT_ + r_ * 3.0d0 * dhth * vT_ * vT_ + &
                    r_ * hth * 2.0d0 * vT_ * dvT_ - (domc_ / omc_) * r_ * hth * vT_ * vT_);

!(hz r0^2 Vth^2 / omc)':
dhzr2Vth2_omc = r_ * Vth / omc_ * (dhz * r_ * Vth + hz * 2.0d0 * Vth + hz * r_ * 2.0d0 * dVth - &
                (domc_ / omc_) * hz * r_ * Vth);

!(hz hth r0^2 vT Vth / omc)':
dhzhthr2vTVth_omc = r_ / omc_ * (dhz * hth * r_ * vT_ * Vth + hz * dhth * r_ * vT_ * Vth + &
                    hz * hth * 2.0d0 * vT_ * Vth + hz * hth * r_ * dvT_ * Vth + hz * hth * r_ * vT_ * dVth - &
                    (domc_ / omc_) * hz * hth * r_ * vT_ * Vth);

!(hz r0^2 hth^2 vT^2 / omc):
dhzr2hth2vT2_omc = r_ * hth * vT_ / omc_ * (dhz * r_ * hth * vT_ + hz * 2.0d0 * hth * vT_ + &
                   hz * r_ * 2.0d0 * dhth * vT_ + hz * r_ * hth * 2.0d0 * dvT_ - (domc_ / omc_) * hz * r_ * hth * vT_);

!Q - array:
!phi:
Q(1,0,0) = omc_;

Q(1,1,0) = 0.0d0;

Q(1,2,0) = 0.0d0;

Q(1,0,2) = 0.0d0;

!theta:
Q(2,0,0) = Vth - hz / omc_ * Vth * Vth;

Q(2,1,0) = vT_ * (hth - 2.0d0 * hth * hz * Vth / omc_);

Q(2,2,0) = - hz * hth * hth * vT_ * vT_ / omc_;

Q(2,0,2) = hz * vT_ * vT_ * (domc_/ omc_) / 2.0 / omc_ / r_;

!z:
Q(3,0,0) = Vz + hth * r_ * r_ / omc_ * Vth * Vth;

Q(3,1,0) = vT_ * (hz + 2.0d0 * hth * hth * r_ * r_ * Vth / omc_);

Q(3,2,0) = hth * hth * hth * r_ * r_ * vT_ * vT_ / omc_;

Q(3,0,2) = - hth * r_ * vT_ * vT_ * (domc_/ omc_) / 2.0 / omc_;

!P array: without l term!
P(0,0) =   v1 * m_ * r_ * (ks_ + 2.0d0 * hth * r_ * kp_ * Vth / omc_);

P(1,0) =   m_ * (v1 * 2.0d0 * r_ * r_ * vT_ * hth * hth * kp_ / omc_ + omc_ * r_ * kp_ / vT_ - &
           drkV / vT_ - dr2kpVth2_omc / vT_);

P(2,0) = - m_ / vT_ * (drksvT + 2.0d0 * dhthr2vTkpVth_omc);

P(3,0) = - m_ / vT_ * dhth2r2vT2kp_omc;

P(0,2) =   0.5d0 * v2 * m_ * r_ * (ks_ + 2.0d0 * hth * r_ * kp_ * Vth / omc_);

P(1,2) =   v2 * m_ * r_ * r_ * hth * hth * vT_ * kp_ / omc_;

!R array: without l*omc term!
R(0,0) = ks_ * vE + kp_ * Vp_ - ks_ * r_ * Vth * Vth / omc_ - omega_;

R(1,0) = kp_ * vT_ - 2.0d0 * hth * Vth * vT_ * ks_ * r_ / omc_;

R(2,0) = - hth * hth * vT_ * vT_ * ks_ * r_ / omc_;

R(0,2) = ks_ * vT_ * vT_ * (domc_ / omc_) / 2.0d0 / omc_;

!Z array:
!phi:
Z(1,0,0) = m_ * r_ * omc_ * omc_ / vT_ / vT_;

Z(1,1,0) = 0.0d0;

Z(1,2,0) = 0.0d0;

Z(1,3,0) = 0.0d0;

Z(1,0,2) = 0.0d0;

Z(1,1,2) = 0.0d0;

!theta:
Z(2,0,0) =   m_ * v1 * (hz + 2.0d0 * hth * hth * r_ * r_ * Vth / omc_);

Z(2,1,0) =   m_ * (2.0d0 * v1 * r_ * r_ * vT_ * hth * hth * hth / omc_ + (omc_ * r_ * hth - dVz) / vT_ - dhthr2Vth2_omc / vT_);

Z(2,2,0) = - m_ / vT_ * (hz_ * dvT_ + dhz_ * vT_ + 2.0d0 * dhth2r2vTVth_omc);

Z(2,3,0) = - m_ / vT_ * dhthr2hth2vT2_omc;

Z(2,0,2) =   m_ * v2 * (hz / 2.0d0 + hth * hth * r_ * r_ * Vth / omc_);

Z(2,1,2) =   m_ * v2 * hth * hth * hth * r_ * r_ * vT_ / omc_;

!z:
Z(3,0,0) = - m_ * v1 * (h_th - 2.0d0 * hth * hz * r_ * r_ * Vth / omc_);

Z(3,1,0) =   m_ * (2.0d0 * v1 * hz * hth * hth * r_ * r_ * vT_ / omc_ + (omc_ * r_ * hz + dVt) / vT_ - dhzr2Vth2_omc / vT_);

Z(3,2,0) =   m_ / vT_ * (h_th * dvT_ + dh_th * vT_ - 2.0d0 * dhzhthr2vTVth_omc);

Z(3,3,0) = - m_ / vT_ * dhzr2hth2vT2_omc;

Z(3,0,2) = - m_ * v2 * (h_th / 2.0d0 - hz * hth * r_ * r_ * Vth / omc_);

Z(3,1,2) =   m_ * v2 * hz * hth * hth * r_ * r_ * vT_ / omc_;

! now include the principal terms:
P00 = P(0,0) + L * (m_ * r_ * omc_ * omc_ / vT_ / vT_);

R00 = R(0,0) + L * (omc_);

do k = 1,3

    !{mu=0, nu=0}: F00: P00 Q00 - R00 Z00
    F(k,1,:) = P00*Q(k,0,0) - R00*Z(k,0,0);

    !{mu=0, nu=1}: F02: P02 Q00 + P00 Q02 - R02 Z00 - R00 Z02
    F(k,2,:) = P(0,2)*Q(k,0,0) + P00*Q(k,0,2) - R(0,2)*Z(k,0,0) - R00*Z(k,0,2);

    !{mu=0, nu=2}: F04: P02 Q02 - R02 Z02
    F(k,3,:) = P(0,2)*Q(k,0,2) - R(0,2)*Z(k,0,2);

    !{mu=1, nu=0}: F10: P10 Q00 + P00 Q10 - R10 Z00 - R00 Z10
    F(k,4,:) = P(1,0)*Q(k,0,0) + P00*Q(k,1,0) - R(1,0)*Z(k,0,0) - R00*Z(k,1,0);

    !{mu=1, nu=1}: F12: P12 Q00 + P10 Q02 + P02 Q10 - R10 Z02 - R02 Z10 - R00 Z12
    F(k,5,:) = P(1,2)*Q(k,0,0) + P(1,0)*Q(k,0,2) + P(0,2)*Q(k,1,0) - R(1,0)*Z(k,0,2) - R(0,2)*Z(k,1,0) - R00*Z(k,1,2);

    !{mu=1, nu=2}: F14: P12 Q02 - R02 Z12
    F(k,6,:) = P(1,2)*Q(k,0,2) - R(0,2)*Z(k,1,2);

    !{mu=2, nu=0}: F20: P20 Q00 + P10 Q10 + P00 Q20 - R20 Z00 - R10 Z10 - R00 Z20
    F(k,7,:) = P(2,0)*Q(k,0,0) + P(1,0)*Q(k,1,0) + P00*Q(k,2,0) - R(2,0)*Z(k,0,0) - R(1,0)*Z(k,1,0) - R00*Z(k,2,0);

    !{mu=2, nu=1}: F22: P20 Q02 + P12 Q10 + P02 Q20 - R20 Z02 - R10 Z12 - R02 Z20
    F(k,8,:) = P(2,0)*Q(k,0,2) + P(1,2)*Q(k,1,0) + P(0,2)*Q(k,2,0) - R(2,0)*Z(k,0,2) - R(1,0)*Z(k,1,2) - R(0,2)*Z(k,2,0);

    !{mu=3, nu=0}: F30: P30 Q00 + P20 Q10 + P10 Q20 - R20 Z10 - R10 Z20 - R00 Z30
    F(k,9,:) = P(3,0)*Q(k,0,0) + P(2,0)*Q(k,1,0) + P(1,0)*Q(k,2,0) - R(2,0)*Z(k,1,0) - R(1,0)*Z(k,2,0) - R00*Z(k,3,0);

    !{mu=3, nu=1}: F32: P30 Q02 + P12 Q20 - R20 Z12 - R02 Z30
    F(k,10,:) = P(3,0)*Q(k,0,2) + P(1,2)*Q(k,2,0) - R(2,0)*Z(k,1,2) - R(0,2)*Z(k,3,0);

    !{mu=4, nu=0}: F40: P30 Q10 + P20 Q20 - R20 Z20 - R10 Z30
    F(k,11,:) = P(3,0)*Q(k,1,0) + P(2,0)*Q(k,2,0) - R(2,0)*Z(k,2,0) - R(1,0)*Z(k,3,0);

    !{mu=5, nu=0}: F50: P30 Q20 - R20 Z30
    F(k,12,:) = P(3,0)*Q(k,2,0) - R(2,0)*Z(k,3,0);

end do

!redefinition of F for phi component to remove big numbers substraction:
!F(1,1,:) = F(1,1,0);

!F1 coefficients for phi, theta, z - without adiabatic contribution:
do k = 1,3

!    F1(k,1,:) = P00*Q(k,0,0);
!    F1(k,2,:) = P(1,0)*Q(k,0,0) + P00*Q(k,1,0);
!    F1(k,3,:) = P(0,1)*Q(k,0,0);
!    F1(k,4,:) = P(2,0)*Q(k,0,0) + P(1,0)*Q(k,1,0);
!    F1(k,5,:) = P(0,1)*Q(k,1,0);
!    F1(k,6,:) = P(2,0)*Q(k,1,0);

end do

end subroutine

!------------------------------------------------------------------------------

subroutine eval_A_matrix () !only to be called after all apropriate back subs

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

implicit none;

print *
print *, 'conductivity arrays:'
print *, 'vE:', vE;
print *, 'dvE:', dvE;
print *, 'dvEz:', dvEz;
print *, 'dvEt:', dvEt;
print *, 'F:', F;
print *, 'A:', A;
print *, 'Imnl:', Imnkl;
print *, 'D:', D;

call flush(101)

end subroutine

!------------------------------------------------------------------------------

real(8) function aco (m, k)

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
