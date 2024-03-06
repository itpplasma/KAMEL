!------------------------------------------------------------------------------

subroutine calc_k_matrices (kmat)

use constants, only: dp, dpc, tppoh, Im;
use conduct_parameters, only: n_, m_, q_, vT_;
use conduct_arrays, only: I, Is1, dimI;
use conduct_arrays, only: F, Fs1, Fs2, dimF;
use conduct_arrays, only: A, s1, s2, factor;
use conduct_arrays, only: W2, D;
use flre_sett, only: Nmax, flre_order, rsp;

implicit none;

!fortran: kmat(j,i,n2,n1), C: kmat[n1,n2,i,j]
complex(dpc), dimension(1:3, 1:3, 0:flre_order, 0:flre_order), intent(out) :: kmat;

integer :: n1, n2, k, j, al, be, i1, i2, l;
integer :: s1ka, s2ka, sska, s1jb, s2jb, ssjb;
complex(dpc) :: fac_g, fac_n, fac_ab, fac_mne;
complex(dpc) :: sum_ab, sum_mne, sum_l;
complex(dpc) :: Dv;

fac_g = n_ / (m_ * vT_)**4 / vT_ / tppoh;

do n1 = 0, flre_order !n
    do n2 = 0, flre_order !n'

        fac_n = (-1.0d0)**(n1+1)/factor(n1)/factor(n2)*Im**(n1+n2+1);

        do k = 1,3 !r,theta,z - current density component
            do j = 1,3 !r,theta,z - electric field component

                sum_ab = cmplx(0.0d0, 0.0d0, dpc);

                do al = 1,3 !phi, theta_g, z_g - alpha
                    do be = 1,3 !phi, theta_g, z_g - beta

                        s1ka = s1(k,al);
                        s2ka = s2(k,al);
                        s1jb = s1(j,be);
                        s2jb = s2(j,be);

                        sska = s1ka + s2ka;
                        ssjb = s1jb + s2jb;

                        fac_ab = (-1.0)**(sska)*(Im)**(sska+ssjb)*A(n1,k,al)*A(n2,j,be);

                        sum_mne = cmplx(0.0d0, 0.0d0, dpc);

                        do i2 = 1,dimF(be) !mu,nu
                            do i1 = 1,dimI(al) !eta

                                fac_mne = I(al,i1);

                                sum_l = cmplx(0.0d0, 0.0d0, dpc);

                                do l = -Nmax,Nmax
                                    Dv = D(s1ka, s2ka+n1, s1jb, s2jb+n2, Fs2(be,i2), l);
                                    sum_l = sum_l + F(be,i2,l) * Dv * &
                                            W2(Is1(al,i1), Fs1(be,i2), l);
                                end do

                                sum_mne = sum_mne + fac_mne * sum_l;

                            end do
                        end do

                        sum_ab = sum_ab + fac_ab * sum_mne;

                    end do
                end do

                kmat(j,k,n2,n1) = fac_g * fac_n * sum_ab; !cylindrical system

            end do
        end do
    end do
end do

if (rsp > 0) call transform_k_matrices_to_rsp(kmat); !transformation to rsp system

end subroutine

!------------------------------------------------------------------------------

subroutine transform_k_matrices_to_rsp (kmat)

use constants, only: dp, dpc;
use conduct_arrays, only: bico, dEM, dJMI;
use flre_sett, only: flre_order;

implicit none;

interface
    real(8) function aco (m, k)
        integer, intent(in) :: m, k;
    end function
end interface

!fortran: kmat(j,i,n2,n1), C: kmat[n1,n2,i,j]
complex(dpc), dimension(1:3, 1:3, 0:flre_order, 0:flre_order), intent(out) :: kmat;

complex(dpc), dimension(1:3, 1:3, 0:flre_order, 0:flre_order) :: kcpy;
complex(dpc), dimension(1:3, 1:3) :: mat, summa;
integer :: n1, n2, k, j;

kcpy = kmat;

do k = 0,flre_order
    do j = 0,flre_order

        summa = cmplx(0.0d0, 0.0d0, 8);

        do n1 = 0,flre_order-k
            do n2 = 0,flre_order-j

                !kmat(j,i,n2,n1) - can be optimized by taking one matmul outside!
                mat = matmul(dJMI(:,:,n1), transpose(kcpy(:,:,n2+j,n1+k)));
                summa = summa + aco(n1+k-1,k)*bico(n2+j,j)*matmul(mat, dEM(:,:,n2));

            end do
        end do

        kmat(:,:,j,k) = transpose(summa);

    end do
end do

end subroutine

!------------------------------------------------------------------------------

subroutine calc_dEM_dJMI_matrices (r, N, dEM, dJMI)

use constants, only: dp;
use core, only: bp_ptr;
use background, only: flag_back, huge_factor;

implicit none;

real(dp), intent(in) :: r;
integer, intent(in) :: N;
real(dp), dimension(1:3, 1:3, 0:N), intent(out) :: dEM, dJMI;

real(dp), dimension(0:N,1:2) :: Htz;

real(dp) :: r_, ht, dht, drht, hz, dhz, drhz;
integer :: k;

r_ = r;

call eval_hthz (r, 0, N, bp_ptr, Htz);

if (flag_back == 'w') then

    Htz(1:N,:) = 0.0d0;
    r_ = r*huge_factor;

else if (flag_back == 'h') then

    Htz = 0.0d0;
    Htz(0,2) = 1.0d0 !hz = 1;
    r_ = r*huge_factor;

else if (flag_back == 's') then

    print *, 'warning: this flag_back is not implemented!..'

end if

ht = Htz(0,1);
hz = Htz(0,2);

dEM = 0.0d0;
dJMI = 0.0d0;

dEM(1,1,0) =   1.0d0;
dEM(2,2,0) =   r_ * hz;
dEM(2,3,0) =   r_ * ht;
dEM(3,2,0) = - ht;
dEM(3,3,0) =   hz;

dJMI(1,1,0) =   1.0d0;
dJMI(2,2,0) =   r_ * hz;
dJMI(2,3,0) = - ht;
dJMI(3,2,0) =   r_* ht;
dJMI(3,3,0) =   hz;

do k = 1,N

    dht = Htz(k,1);
    drht = r_ * Htz(k,1) + k * Htz(k-1,1);

    dhz = Htz(k,2);
    drhz = r_ * Htz(k,2) + k * Htz(k-1,2);

    dEM(2,2,k) =   drhz;
    dEM(2,3,k) =   drht;
    dEM(3,2,k) = - dht;
    dEM(3,3,k) =   dhz;

    dJMI(2,2,k) =   drhz;
    dJMI(2,3,k) = - dht;
    dJMI(3,2,k) =   drht;
    dJMI(3,3,k) =   dhz;

end do

! write(100,*) r_, dEM
! write(200,*) r_, dJMI

end subroutine

!------------------------------------------------------------------------------

subroutine calc_dummy_matrices (kmat)

use constants, only: dp, dpc;
use flre_sett, only: flre_order;

implicit none;

!fortran: kmat(j,i,n2,n1), C: kmat[n1,n2,i,j]
complex(dpc), dimension(1:3, 1:3, 0:flre_order, 0:flre_order), intent(out) :: kmat;

kmat = cmplx(123456789.0d123, 987654321.0d213, 8);

end subroutine

!------------------------------------------------------------------------------

subroutine transform_c_matrices_to_rsp (r, cmat)

use constants, only: dp, dpc;
use flre_sett, only: flre_order;

implicit none;

!fortran: cmat(j,i,s), C: cmat[s,i,j]
real(dp), intent(in) :: r;
complex(dpc), dimension(1:3, 1:3, 0:2*flre_order), intent(out) :: cmat;

complex(dpc), dimension(1:3, 1:3, 0:2*flre_order) :: ccpy;
complex(dpc), dimension(1:3, 1:3) :: summa;

real(dp), dimension(1:3, 1:3, 0:2*flre_order) :: dEM, dJMI;
real(dp), dimension(0:2*flre_order, 0:2*flre_order) :: bico;

integer :: s, l;

call calc_dEM_dJMI_matrices (r, 2*flre_order, dEM, dJMI);
call binomial_coefficients (%val(2*flre_order), bico);

ccpy = cmat;

do s = 0,2*flre_order

    summa = cmplx(0.0d0, 0.0d0, 8);

    do l = 0,2*flre_order-s

        summa = summa + bico(l+s,s)*matmul(transpose(ccpy(:,:,l+s)), dEM(:,:,l));

    end do

    cmat(:,:,s) = transpose(matmul(dJMI(:,:,0), summa));

end do

end subroutine

!------------------------------------------------------------------------------

subroutine calc_k1_matrices (kmat) !without adiabatic contribution

use constants, only: dp, dpc, tppoh, Im;
use conduct_parameters, only: n_, m_, q_, vT_;
use conduct_arrays, only: I, Is1, dimI;
use conduct_arrays, only: F1, Fs1, Fs2, dimF;
use conduct_arrays, only: A, s1, s2, factor;
use conduct_arrays, only: W2, D;
use flre_sett, only: Nmax, flre_order, rsp;

implicit none;

!fortran: kmat(j,i,n2,n1), C: kmat[n1,n2,i,j]
complex(dpc), dimension(1:3, 1:3, 0:flre_order, 0:flre_order), intent(out) :: kmat;

integer :: n1, n2, k, j, al, be, i1, i2, l;
integer :: s1ka, s2ka, sska, s1jb, s2jb, ssjb;
complex(dpc) :: fac_g, fac_n, fac_ab, fac_mne;
complex(dpc) :: sum_ab, sum_mne, sum_l;
complex(dpc) :: Dv;

fac_g = n_ / (m_ * vT_)**4 / vT_ / tppoh;

do n1 = 0, flre_order !n
    do n2 = 0, flre_order !n'

        fac_n = (-1.0)**(n1+1)/dble(factor(n1))/dble(factor(n2))*Im**(n1+n2+1);

        do k = 1,3 !r,theta,z - current density component
            do j = 1,3 !r,theta,z - electric field component

                sum_ab = cmplx(0.0d0, 0.0d0, dpc);

                do al = 1,3 !phi, theta_g, z_g - alpha
                    do be = 1,3 !phi, theta_g, z_g - beta

                        s1ka = s1(k,al);
                        s2ka = s2(k,al);
                        s1jb = s1(j,be);
                        s2jb = s2(j,be);

                        sska = s1ka + s2ka;
                        ssjb = s1jb + s2jb;

                        fac_ab = (-1.0)**(sska)*(Im)**(sska+ssjb)*A(n1,k,al)*A(n2,j,be);

                        sum_mne = cmplx(0.0d0, 0.0d0, dpc);

                        do i2 = 1,dimF(be) !mu,nu
                            do i1 = 1,dimI(al) !eta

                                fac_mne = I(al,i1);

                                sum_l = cmplx(0.0d0, 0.0d0, dpc);

                                do l = -Nmax,Nmax
                                    Dv = D(s1ka, s2ka+n1, s1jb, s2jb+n2, Fs2(be,i2), l);
                                    sum_l = sum_l + F1(be,i2,l) * Dv * &
                                            W2(Is1(al,i1), Fs1(be,i2), l);
                                end do

                                sum_mne = sum_mne + fac_mne * sum_l;

                            end do
                        end do

                        sum_ab = sum_ab + fac_ab * sum_mne;

                    end do
                end do

                kmat(j,k,n2,n1) = fac_g * fac_n * sum_ab; !cylindrical system

            end do
        end do
    end do
end do

if (rsp > 0) call transform_k_matrices_to_rsp(kmat); !transformation to rsp system

end subroutine

!------------------------------------------------------------------------------
