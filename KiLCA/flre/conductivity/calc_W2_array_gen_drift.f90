!<Calculates array of W2 special functions.

!------------------------------------------------------------------------------

subroutine calc_W2_array (spec)

use constants, only: dpc, dp, sqrt2p;
use conduct_arrays, only: vE, Imnkl, mpara, mperp;
use conduct_parameters, only: omega_, omc_, domc_, nu_, kp_, ks_, ht_, hz_, vT_, Vp_, r_;
use flre_sett, only: Nmax, collmod;

implicit none;

integer, intent(in) :: spec;

integer :: lc, status;

complex(dpc) :: x1, x2, x3, x4;

real(dp) :: Vth; ! contravariant theta component of the particle velocity

logical :: consenergy;

Vth = (hz_ * vE + ht_ * Vp_) / r_;

x1 = kp_ * vT_ / nu_ - 2.0d0 * ks_ * vT_ * ht_ * Vth / nu_ / omc_;

x3 = ks_ * ht_ * ht_ * vT_ * vT_ / nu_ / r_ / omc_;

x4 = ks_ * vT_ * vT_ * (domc_ / omc_) / omc_ / nu_;

!print *, 'r_=', r_, 'spec=', spec;

if (collmod(spec) == 0) then ! conservation of the number of particles

    consenergy = .false.;

else if (collmod(spec) == 1) then ! conservation of the number and energy of particles

    consenergy = .true.;

else if (collmod(spec) == 2) then ! conservation of the number and momentum of particles

    print *, 'not implemented yet!'

else if (collmod(spec) == 3) then ! conservation of the number, energy and momentum of particles

    print *, 'not implemented yet!'

end if

do lc = -Nmax,Nmax ! over cyclotron harmonics

    x2 = (omega_ - dble(lc) * omc_ - ks_ * vE - kp_ * Vp_) / nu_ + ks_ * r_ * Vth * Vth / nu_ / omc_;

    call calc_imnkl_quad(consenergy, mpara, mperp, x1, x2, x3, x4, Imnkl(:,:,:,:,lc), status);

    if(status /= 0) then
        print *, 'error: calc_W2_array: calc_imnkl_quad failed to evaluate the Imnkl:', x1, x2, x3, x4;
    end if

end do

end subroutine

!------------------------------------------------------------------------------
