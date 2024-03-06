program test

use constants
use conduct_arrays, only: D, factor;
use flre_sett, only: Nmax, Nbmax, flre_order;
implicit none;

external :: calc_D_array;

complex(dpc) :: res, res1, res2;

integer :: m1, n1, m2, n2, l, nu, k, b1, b2;

real(dp) :: x_, rho;

print *, Nmax, Nbmax, flre_order;

call allocate_and_set_conductivity_arrays();

x_ = -0.750d0;
rho = 0.125d0;
m1 = 1;
n1 = 2;
m2 = 0;
n2 = 5;
l = 3;
nu = 2;

do k = 1,1 ! interactive regime

    call calc_D_array(x_);

    ! integral:
    res = cmplx(0.0d0, 0.0d0, dp);

    do b1 = 0,Nbmax
        do b2 = 0,Nbmax

            res = res + 2.0d0**(nu - abs(l) - b1 - b2) * (-1.0d0)**(b1 + b2) * &
                  (factor(nu + abs(l) + b1 + b2) / (factor(b1) * factor(abs(l) + b1) * factor(b2) * factor(abs(l) + b2))) * &
                   (conjg(D(m1, n1, l, b2)) * D(m2, n2, l, b1));

        end do
    enddo

    res = res * rho**(m1 + n1 + m2 + n2);

    print *, 'res =', res;
    print *;

end do

end program
