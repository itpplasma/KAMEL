program test_eval_ff

implicit none;

integer :: m, n; !poloidal, toridal number
real(8) :: s; !radial variable (as stored in formfactors.dat file). It can be: psi_pol, psi_tor, normalized fluxes, etc...
complex(8) :: ff; !form factors  flux array

integer, parameter :: dimg = 20000;

integer :: k;

m = -2;
n = 1;

!loads data during first call:
call eval_form_factors (m, n, 0.5d0, ff);

do k = 1,dimg
    s = real(k-1)/real(dimg-1);

    call eval_form_factors (m, n, s, ff);

    !write (*,*) s, real(ff), aimag(ff);
    write (10,*) s, real(ff), aimag(ff);
end do

end program

!compilation string:
!gfortran -o test_ff -O2 -fimplicit-none -Wall -Wtabs -O2 -mtune=generic -static test_eval_ff.f90 eval_form_factors.f90
