program test

implicit none;

external :: igamma;

real(8),    parameter :: pi = 3.141592653589793238462643383279502884197
complex(8), parameter :: im = (0.0d0, 1.0d0)

complex(8) :: s, z;
real(8)    :: abs_err = 1.0d-12, rel_err = 1.0d-12;
complex(8) :: G;
real(8)    :: est_err;
integer    :: Nterms;
integer    :: status;

integer :: k, i, j;

real(8) :: dphi;

do k = 1,10 ! interactive regime

    print *, 'input s ='
    read *, s
    print *, 'input z ='
    read *, z

    call incomplete_gamma(s, z, abs_err, rel_err, 0, G, est_err, Nterms, status);

    print *, 's =', s;
    print *, 'z =', z;
    print *, 'G =', G;
    print *, 'err =', est_err;
    print *, 'status =', status;
    print *;

end do

! automatic checks for array of arguments:

s = -2.0 + im;

dphi = 2.0 * pi / 32;

do i = 1,32

    do j = -6,9

        z = cmplx(10,0,8)**j * exp(im * dphi * (i-1));

        call incomplete_gamma(s, z, abs_err, rel_err, G, est_err, Nterms, status);

        write(1000, *) s, z, G, est_err, Nterms, status;

    end do

end do

end program
