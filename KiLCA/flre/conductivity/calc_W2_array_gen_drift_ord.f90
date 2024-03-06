!<Calculates array of W2 special functions.

!------------------------------------------------------------------------------

subroutine calc_W2_array (spec)

use constants, only: dpc, dp, sqrt2p;
use conduct_arrays, only: vE, Imnl, SImnl, IImnl, mpara, mperp;
use conduct_parameters, only: omega_, omc_, domc_, nu_, kp_, ks_, ht_, hz_, vT_, Vp_, r_;
use flre_sett, only: Nmax, collmod;

implicit none;

integer, intent(in) :: spec;

integer :: lc, m, n, l, status, lperp;

complex(dpc) :: x1, x2, x3, x4;

real(dp) :: Vth; ! contravariant theta component of the particle velocity

Vth = (hz_ * vE + ht_ * Vp_) / r_;

x1 = kp_ * vT_ / nu_ - 2.0d0 * ks_ * vT_ * ht_ * Vth / nu_ / omc_;

x3 = ks_ * ht_ * ht_ * vT_ * vT_ / nu_ / r_ / omc_;

x4 = ks_ * vT_ * vT_ * (domc_ / omc_) / omc_ / nu_;

lperp = 2 * mperp + 1;

!print *, 'r_=', r_, 'spec=', spec;

do lc = -Nmax,Nmax ! over cyclotron harmonics

    x2 = (omega_ - dble(lc) * omc_ - ks_ * vE - kp_ * Vp_) / nu_ + ks_ * r_ * Vth * Vth / nu_ / omc_;

    if (collmod(spec) == 0) then ! conservation of the number of particles

        !call GreenMom(mpara, lperp, real(x1), real(x2), real(x3), real(x4), SImnl); ! there is a problem for complex omega!!!

        call calc_imnl_quad(mpara, lperp, x1, x2, x3, x4, IImnl, status);

        if(status /= 0) then
            print *, 'error: calc_W2_array: calc_imnl_quad failed to evaluate the Imnl:', x1, x2, x3, x4;
        end if

        do m = 0,mpara
            do n = 0,mpara
                do l = 0,mperp

                    !Imnl(m,n,l,lc) = SImnl(m,n,2*l+1); !Sergey's version
                    Imnl(m,n,l,lc) = IImnl(m,n,2*l+1); !Ivan's version

                end do
            end do
        end do

!          ! compare results:
!         if(maxval(abs((-sqrt2p*SImnl - IImnl)/IImnl)) > 1.0d-3) then
!
!             print *, 'WARNING: a difference between two versions is found:';
!             do m = 0,mpara
!                 do n = 0,mpara
!                     do l = 1,lperp
!                         write(*,'(A1,3(I3),A3,4(E30.16))') '(',m, n, l,') =', IImnl(m,n,l), -sqrt2p * SImnl(m,n,l);
!                     end do
!                 end do
!             end do
!
!         end if

    else if (collmod(spec) == 1) then ! conservation of the number and energy of particles

        print *, 'not implemented yet!'

    else if (collmod(spec) == 2) then ! conservation of the number and momentum of particles

        print *, 'not implemented yet!'

    else if (collmod(spec) == 3) then ! conservation of the number, energy and momentum of particles

        print *, 'not implemented yet!'

    end if

end do

end subroutine

!------------------------------------------------------------------------------
