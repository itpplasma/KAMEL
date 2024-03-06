!<Calculates array of D special functions.

!------------------------------------------------------------------------------

subroutine calc_D_array ()

use constants, only: dpc, dp, Im;
use conduct_arrays, only: D, x;
use conduct_parameters, only: ks_, vT_, omc_;
use flre_sett, only: Nmax, Nbmax, flre_order;

implicit none;

integer :: l, b, m, n, ind, k;

real(dp) :: x_;

complex(dpc) :: Ifac;

x_ = ks_ * vT_ / omc_;

!if (abs(x_) > 1.0d0) print *, 'calc_D_array: warning: BesselJ argument x_ is greater than 1.0:', x_;

x(0) = 1.0d0;

do l = 1,Nmax+2*Nbmax

    !x(l) = x(l-1) * x_;
    x(l) = (x_)**l;

end do

do l = 0,Nmax

    do b = 0,Nbmax

        k = 2*b+l;

        ind = max(0, k-0); D(0, 0, l, b) = x(ind);

        ind = max(0, k-1); D(0, 1, l, b) = Im*l*x(ind);

        ind = max(0, k-2); D(0, 2, l, b) = (k-l*l)*x(ind);

        ind = max(0, k-1); D(1, 0, l, b) = (k)*x(ind);

        ind = max(0, k-2); D(1, 1, l, b) = Im*l*(-1+k)*x(ind);

        ind = max(0, k-3); D(1, 2, l, b) = (-2+k)*(k-l*l)*x(ind);

        if (flre_order < 3) cycle;

        ind = max(0, k-3); D(0, 3, l, b) = -Im*l*(2-6*b-3*l+l*l)*x(ind);

        ind = max(0, k-4); D(0, 4, l, b) = (12*b*b+(-3+l)*(-2+l)*(-1+l)*l-12*b*(1+(-1+l)*l))*x(ind);

        ind = max(0, k-4); D(1, 3, l, b) = -Im*l*(-3+k)*(2-6*b-3*l+l*l)*x(ind);

        ind = max(0, k-5); D(1, 4, l, b) = (-4+k)*(12*b*b+(-3+l)*(-2+l)*(-1+l)*l-12*b*(1+(-1+l)*l))*x(ind);

        if (flre_order < 5) cycle;

        ind = max(0, k-5); D(0, 5, l, b) = Im*l*(60*b*b+(-4+l)*(-3+l)*(-2+l)*(-1+l)-20*b*(5+(-3+l)*l))*x(ind);

        ind = max(0, k-6); D(0, 6, l, b) = (120*b*b*b-(-5+l)*(-4+l)*(-3+l)*(-2+l)*(-1+l)*l-180*b*b*(2+(-1+l)*l)+&
                                            30*b*(8+(-1+l)*l*(12+(-5+l)*l)))*x(ind);

        ind = max(0, k-6); D(1, 5, l, b) = Im*l*(-5+k)*(60*b*b+(-4+l)*(-3+l)*(-2+l)*(-1+l)-20*b*(5+(-3+l)*l))*x(ind);

        ind = max(0, k-7); D(1, 6, l, b) = (-6+k)*(120*b*b*b-(-5+l)*(-4+l)*(-3+l)*(-2+l)*(-1+l)*l-&
                                           180*b*b*(2+(-1+l)*l)+30*b*(8+(-1+l)*l*(12+(-5+l)*l)))*x(ind);

    end do
end do

do m = 0,1
    do n = 0,flre_order+1

        Ifac = Im**(m + n);

        do b = 0,Nbmax

            l = 0;
            D(m, n, l, b) = Ifac * D(m, n, l, b);

            do l = 1,Nmax

                D(m, n, -l, b) = conjg(D(m, n, l, b));
                D(m, n,  l, b) = Ifac * D(m, n,  l, b);
                D(m, n, -l, b) = Ifac * D(m, n, -l, b);

            end do
        end do

    end do
end do

end subroutine

!------------------------------------------------------------------------------
