program test_det

implicit none;

external :: calc_determinant

complex(8), parameter :: im = (0.0d0, 1.0d0)

complex(8) :: det

integer, parameter :: N = 11;

complex(8), dimension(N,N) :: A;

real(8) :: k, l;

real(8), parameter :: pi    = 3.141592653589793238462643383279502884197
real(8), parameter :: euler = 0.5772156649015328606065120900824024310422

do k = 1,N
    do l = 1,N
        !A(i,j) = k+im*l;
        !A(k,l) = 1e8*(k*l+im*l**2/k**2);
        !A(k,l) = sin(k**l+im*l**k);
        A(k,l) = k**(l+1)-im*l**k;
    end do
end do

call calc_determinant (N, A, det);

print *, det

end program
