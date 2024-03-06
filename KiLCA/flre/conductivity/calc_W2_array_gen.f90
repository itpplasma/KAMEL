!<Calculates array of W2 special functions.

!------------------------------------------------------------------------------

subroutine calc_W2_array (spec)

use constants, only: dpc, dp, sqrt2p;
use conduct_arrays, only: vE, W2;
use conduct_parameters, only: omega_, omc_, kp_, ks_, nu_, vT_, Vp_, r_, q_;
use flre_sett, only: Nmax, collmod;

implicit none;

integer, intent(in) :: spec;

real(dp), dimension(1:7) :: Vt;
real(dp) :: omega0;

integer :: l, m, n;

complex(dpc) :: x1, x2;

complex(dpc), dimension(0:3,0:3) :: Imn;

complex(dpc) :: Ip, A, AM, B, BM, DELTA;

Vt(1) = vT_;

do l = 2,7
    Vt(l) = Vt(l-1) * Vt(1);
end do

do l = -Nmax,Nmax

    omega0 = dble(l) * omc_ + ks_ * vE + kp_ * Vp_;

    x1 = cmplx(kp_ * vT_ / nu_, 0.0d0, 8);
    x2 = (omega_ - omega0) / nu_;

    call calc_Imn_array(x1, x2, Imn);

    if (collmod(spec) == 0) then ! conservation of the number of particles

        do m = 0,1
            do n = 0,3
                W2(m,n,l) = sqrt2p / nu_ * Vt(m+n+1) * Imn(m,n);
            end do
        end do

    else if (collmod(spec) == 1) then ! conservation of the number and energy of particles

        A = cmplx(1.0d0, 0.0d0, 8) + (Imn(2,0) - Imn(0,0)) + (Imn(2,0) - Imn(2,2));
        !A = cmplx(1.0d0, 0.0d0, 8) - Imn(0,0) + 2.0d0 * Imn(2,0) - Imn(2,2);
        do m = 0,1
            do n = 0,3

                Ip = Imn(m,n) + ((Imn(m,0) - Imn(m,2)) * (Imn(n,0) - Imn(n,2))) / A;

                W2(m,n,l) = sqrt2p / nu_ * Vt(m+n+1) * Ip;

            end do
        end do

    else if (collmod(spec) == 2) then ! conservation of the number and momentum of particles

        AM = cmplx(1.0d0, 0.0d0, 8) - Imn(1,1);

        do m = 0,1
            do n = 0,3

                Ip = Imn(m,n) + (Imn(m,1) * Imn(n,1)) / AM;

                W2(m,n,l) = sqrt2p / nu_ * Vt(m+n+1) * Ip;

            end do
        end do

    else if (collmod(spec) == 3) then ! conservation of the number, energy and momentum of particles

        A  = cmplx(1.0d0, 0.0d0, 8) + (Imn(2,0) - Imn(0,0)) + (Imn(2,0) - Imn(2,2));
        AM = cmplx(1.0d0, 0.0d0, 8) - Imn(1,1);
        B  = Imn(1,2) - Imn(1,0);
        BM = B;
        DELTA = A * AM - B * BM;

        do m = 0,1
            do n = 0,3

                Ip = Imn(m,n) + (BM * Imn(m,1) * (Imn(n,2) - Imn(n,0)) + &
                                 A  * Imn(m,1) * Imn(n,1) + &
                                 AM * (Imn(m,2) - Imn(m,0)) * (Imn(n,2) - Imn(n,0)) + &
                                 B  * (Imn(m,2) - Imn(m,0)) * Imn(n,1) ) / DELTA;

                W2(m,n,l) = sqrt2p / nu_ * Vt(m+n+1) * Ip;

            end do
        end do

        !check for Onsager symmetry:
        !print *, x1, x2;

        !print *, '12-30: ';
        !print *, W2(1,2,l) / sqrt2p * nu_ / Vt(1+2+1), W2(0,3,l) / sqrt2p * nu_ / Vt(0+3+1);
        !print *, (W2(1,2,l) / sqrt2p * nu_ / Vt(1+2+1) - W2(0,3,l) / sqrt2p * nu_ / Vt(0+3+1)) / &
        !         (W2(1,2,l) / sqrt2p * nu_ / Vt(1+2+1));

        !print *, '11-20: ';
        !print *, W2(1,1,l) / sqrt2p * nu_ / Vt(1+1+1), W2(0,2,l) / sqrt2p * nu_ / Vt(0+2+1);
        !print *, (W2(1,1,l) / sqrt2p * nu_ / Vt(1+1+1) - W2(0,2,l) / sqrt2p * nu_ / Vt(0+2+1)) / &
        !                    (W2(1,1,l) / sqrt2p * nu_ / Vt(1+1+1));

        !31-22 also has been checked: 22 function should be additionally computed

    end if

end do

end subroutine

!------------------------------------------------------------------------------
