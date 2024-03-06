program test

implicit none;

! for mesh:
integer, parameter :: dp = 8;
integer, parameter :: N = 3, deg = 2 * N, nmax = 5;
real(dp), parameter :: delta = 1.0d-3;
real(dp), dimension(-N:N) :: node = (/-3.0d0*delta, -2.0d0*delta, -1.0d0*delta, 0.0d0, 1.0d0*delta, 2.0d0*delta, 3.0d0*delta/);
!real(dp), dimension(-N:N) :: node = (/-4.0d0*delta, -3.0d0*delta, -2.0d0*delta, &
!-1.0d0*delta, 0.0d0, 1.0d0*delta, 2.0d0*delta, 3.0d0*delta, 4.0d0*delta/);
!real(dp), dimension(-N:N) :: node = (/-5.0d0*delta, -4.0d0*delta, -3.0d0*delta, -2.0d0*delta, &
!-1.0d0*delta, 0.0d0, 1.0d0*delta, 2.0d0*delta, 3.0d0*delta, 4.0d0*delta, 5.0d0*delta/);
real(dp), dimension(-N:N,-N:N) :: Qij; ! complex
real(dp), dimension(0:nmax) :: der;

real(dp) :: x, y;

integer :: i, j;
integer :: Dmin = 0, Dmax = nmax;


integer                      :: mpara, mperp;
real(dp), allocatable, dimension(:)       :: dblfac;
real(dp),  allocatable, dimension(:,:)       :: binocoef;


! mesh:
do j = -N,N

    y = node(j);

    do i = j,N

        x = node(i);

        print *, x, y;

        Qij(i,j) = exp(x*x + y*y + x**3 + y**3);
        Qij(j,i) = Qij(i,j);

    end do

end do

x = 0.0d0;

call eval_neville_polynom_ready(node, Qij(:,0), deg, x, Dmin, Dmax, der);

print *, 'nder = ', der;
print *, 'eder = ', 1.0d0;


mpara = 8;

!allocate(MC(0:mpara, 0:mpara, 0:mpara, 0:mpara));

allocate(dblfac(-1:mpara), binocoef(0:mpara, 0:mpara));


! double factorial:
dblfac(-1) = 1.0d0;
dblfac(0)  = 1.0d0;

do i = 1,mpara

    dblfac(i) = dblfac(i-2) * i;

end do

! binomial coefficients C^i_j:
! binomial coefficients:
binocoef = 0.0d0;
do i = 0,mpara

    binocoef(i,0) = 1.0d0;

    do j = 1,i

        binocoef(i,j) = binocoef(i,j-1) * dble(i-j+1) / dble(j);

    end do

end do

print *;
print *, 'dblfac = '
print *, dblfac;

print *;
print *, 'binocoef = '
print *, binocoef(0,:);
print *, binocoef(1,:);
print *, binocoef(2,:);
print *, binocoef(3,:);
print *, binocoef(4,:);
print *, binocoef(5,:);
print *, binocoef(6,:);
print *, binocoef(7,:);
print *, binocoef(8,:);

end program
