!<The module used to calculate the genious FLRE correction dreamed to Sergey after 36 hours of continious work and ~30 coffee cups without a meal...

module gal_corr

!for additional current:
real(8), dimension(3) :: comfac
real(8) :: dVeq

end module

!------------------------------------------------------------------------------

subroutine eval_and_set_params_for_additional_current (r, spec, flag_back)

!only allowed to be called at the end of all other calls!..

use constants, only: dp, dpc, c
use background, only: huge_factor
use gal_corr
use core, only: bp_ptr;
use conduct_parameters, only: Vp_, dVp_, ht_, dht_, hz_, dhz_, q_, omc_, m_, B, dPhi0_, r_

implicit none;

real(dp), intent(in) :: r
integer, intent(in) :: spec
character(1), intent(in) :: flag_back

real(dp) :: Vs, Vp, dVp, n0;

!for additional current:
call vs_0_f (r, spec, bp_ptr, Vs);

!the difference between Vp_m and Vp_p is very small:
!vp_mom (Vp);
Vp = Vp_
dVp = dVp_

dVeq = dVp - ht_*ht_/r*Vp - (ht_*hz_/r+dht_*hz_-dhz_*ht_)*Vs;

if (flag_back == 'w') then
    !keep rotational transform and set only derivs to zero
    Vs = c/B*dPhi0_;
    dVeq = 0.0d0;
else if (flag_back == 'h') then
    Vs = 0.0d0;
    Vp = 0.0d0;
    dVeq = 0.0d0;
else if (flag_back == 's') then
    print *, 'warning: case is not implemented!..'
end if

call dens_mom (n0);

comfac(1) = q_*q_/(m_*omc_)/(m_*omc_)/omc_/r_

comfac(2) = q_*q_/(c*m_*omc_)*n0*Vs

comfac(3) = q_*q_/(c*m_*omc_)*n0*Vp

end subroutine

!------------------------------------------------------------------------------

subroutine calc_and_add_galilelian_correction (r, spec, flag_back, ct)

use constants, only: dp, dpc, c, im
use flre_sett, only: flre_order
use gal_corr, only: comfac
use conduct_parameters, only: omega_, ht_, dht_, hz_, dhz_, ks_, r_

implicit none;

real(dp), intent(in) :: r
integer, intent(in) :: spec
character (1), intent(in) :: flag_back
complex(dpc), dimension(1:3,1:3,0:2*flre_order), intent(inout) :: ct
complex(dpc), dimension(3,3)  :: delC0, delC1, delC2 !additional current

complex(dpc) :: cio, Ns, N3, N4
real(dp) :: dhs

integer :: i, j

call eval_and_set_all_conduct_parameters(r, spec, flag_back);

dhs = hz_*dht_ - ht_*dhz_;

cio = c/(im*omega_);
Ns  = c/omega_*ks_;
N3  = cio*(dhs+ht_*hz_/r_);
N4  = cio*hz_*hz_/r_;

call eval_and_set_params_for_additional_current (r, spec, flag_back);

ct(1,3,0) = ct(1,3,0) - comfac(3)*Ns
ct(2,3,0) = ct(2,3,0) + comfac(3)*N4
ct(3,3,0) = ct(3,3,0) + comfac(3)*N3
ct(2,3,1) = ct(2,3,1) + comfac(3)*cio

ct(1,2,0) = ct(1,2,0) - comfac(2)*Ns
ct(2,2,0) = ct(2,2,0) + comfac(2)*N4
ct(3,2,0) = ct(3,2,0) + comfac(2)*N3
ct(2,2,1) = ct(2,2,1) + comfac(2)*cio

call delC_matrices (delC0, delC1, delC2);

do i = 1,3
    do j = 1,3
        ct(j,i,0) = ct(j,i,0) - comfac(1)*delC0(i,j)
        ct(j,i,1) = ct(j,i,1) - comfac(1)*delC1(i,j)
        ct(j,i,2) = ct(j,i,2) - comfac(1)*delC2(i,j)
    end do
end do

end subroutine

!------------------------------------------------------------------------------

subroutine delC_matrices(delC0, delC1, delC2)

use constants, only: dpc

complex(dpc), dimension(3,3), intent(out) :: delC0, delC1, delC2

call delC0_11 (delC0(1,1))
call delC0_12 (delC0(1,2))
call delC0_13 (delC0(1,3))
call delC0_21 (delC0(2,1))
call delC0_22 (delC0(2,2))
call delC0_23 (delC0(2,3))
call delC0_31 (delC0(3,1))
call delC0_32 (delC0(3,2))
call delC0_33 (delC0(3,3))

call delC1_11 (delC1(1,1))
call delC1_12 (delC1(1,2))
call delC1_13 (delC1(1,3))
call delC1_21 (delC1(2,1))
call delC1_22 (delC1(2,2))
call delC1_23 (delC1(2,3))
call delC1_31 (delC1(3,1))
call delC1_32 (delC1(3,2))
call delC1_33 (delC1(3,3))

call delC2_11 (delC2(1,1))
call delC2_12 (delC2(1,2))
call delC2_13 (delC2(1,3))
call delC2_21 (delC2(2,1))
call delC2_22 (delC2(2,2))
call delC2_23 (delC2(2,3))
call delC2_31 (delC2(3,1))
call delC2_32 (delC2(3,2))
call delC2_33 (delC2(3,3))

end subroutine delC_matrices

!------------------------------------------------------------------------------

subroutine calc_correction_matrices_in_cylindrical_system (rdim, rgrd, spec, cct)

use constants, only: dp, dpc, c, im;
use conduct_parameters, only: B, dB, ht_, hz_, dht_, dhz_, q_, omc_, m_, r_, ks_;
use conduct_parameters, only: omega_, n_, dn_, ddn_, vT_, dvT_, ddvT_, Vp_, dVp_, dPhi0_, ddPhi0_;
use background, only: rtor;
use mode_data, only: m, n;

implicit none;

integer, intent(in) :: rdim;
real(dp), dimension(rdim), intent(in) :: rgrd;
integer, intent(in) :: spec;
complex(dpc), dimension(1:3,1:3,0:1,1:rdim), intent(out) :: cct; !(i, j, derivative, grid)

character (1) :: flag_back = 'f';

real(dp), dimension(3,3) :: gt;
real(dp), dimension(3) :: V, dV;
complex(dpc), dimension(3) :: F, T;
complex(dpc) :: Ct, Cz;
real(dp) :: G;
real(dp) :: ht_con, ht_cov;
real(dp) :: nT, dnT, ddnT, Vs, dVs, dnTV;

real(dp) :: kz;

integer :: k;

kz = dble(n) / rtor;

V(1) = 0.0d0;  dV(1) = 0.0d0;
gt(1,1) = 1.0d0;  gt(1,2) = 0.0d0;  gt(1,3) = 0.0d0;  gt(2,1) = 0.0d0;  gt(3,1) = 0.0d0;

do k = 1,rdim

    call eval_and_set_all_conduct_parameters(rgrd(k), spec, flag_back);

    nT = n_ * m_ * vT_ * vT_;

    dnT = m_ * (dn_ * vT_ * vT_ + n_ * 2.0d0 * vT_ * dvT_);

    ddnT = m_ * (ddn_ * vT_ * vT_ + dn_ * 2.0d0 * vT_ * dvT_ + &
                 2.0d0 * (dn_ * vT_ + n_ * dvT_) * dvT_ + 2.0d0 * n_ * vT_ * ddvT_);

    Vs = c / B * (dPhi0_ + dnT / q_ / n_);

    V(2) = hz_ / r_ * Vs + ht_ / r_ * Vp_;  !contravariant V^th

    V(3) = - ht_ * Vs + hz_ * Vp_;          !contravariant V^z

    dVs = - c * dB / B / B  * (dPhi0_ + dnT / q_ / n_) +  c / B * (ddPhi0_ + (ddnT * n_ - dnT * dn_) / q_ / n_ / n_);

    dV(2) = (dhz_ * r_ - hz_) / r_ / r_ * Vs  + hz_ / r_ * dVs + &
            (dht_ * r_ - ht_) / r_ / r_ * Vp_ + ht_ / r_ * dVp_;    !contravariant (V^th)'

    dV(3) = - dht_ * Vs - ht_ * dVs + dhz_ * Vp_ + hz_ * dVp_       !contravariant (V^z)'

    ht_con = ht_ / r_;   !contravariant
    ht_cov = ht_ * r_;   !covariant

    gt(2,2) = 1.0d0 / r_ / r_ - ht_con * ht_con;
    gt(3,3) = 1.0d0 - hz_ * hz_;
    gt(2,3) = - ht_con * hz_;
    gt(3,2) = gt(2,3);

    F = q_ * q_ * n_ / (m_ * omc_ * omega_ * r_) * V;

    G = q_ * q_ / (m_ * omc_) / (m_ * omc_) / omc_ / r_;

    dnTV = ht_cov * (dnT * V(2) + nT * dV(2)) + hz_ * (dnT * V(3) + nT * dV(3));

    Ct = dnT * ht_cov - m  * dnTV / omega_;
    Cz = dnT * hz_    - kz * dnTV / omega_;

    T = im * (m * gt(:,2) + kz * gt(:,3)) * G;

    cct(:,1,0,k) = - r_ * ks_ * F;                       ! E_r

    cct(:,2,0,k) = - T * Cz;                             ! E_theta

    cct(:,3,0,k) =   T * Ct;                             ! E_z

    cct(:,1,1,k) = 0.0d0;                                ! E_r'

    cct(:,2,1,k) = - im * hz_ * F    - gt(:,1) * Cz * G; ! E_theta'

    cct(:,3,1,k) =   im * ht_cov * F + gt(:,1) * Ct * G; ! E_z'

end do

end subroutine

!------------------------------------------------------------------------------
