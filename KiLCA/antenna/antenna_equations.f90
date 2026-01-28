! This module provides information about the antenna. The routines calculate
! the surface current density produced by a single mode in the antenna and
! solve the continuity conditions at the antenna. The current mode numbers
! are found in the module "mode".

!------------------------------------------------------------------------------

subroutine continuity (inner, outer, coeffs)

!This routine solves the continuity conditions at the antenna corresonding
!to the current in the antenna. The mode numbers are imported from the module "mode".

!inner ... matrix containing the solution vectors at the inner side of the antenna;
!dimension 2: number of solution vectors
!dimension 1: components per vector
!
!outer ... matrix containing the solution vectors at the outer side of
!the antenna; dimensions see inner
!
!coeffs ... matrix containing the coefficients for solution superposition
!              found by solving the continuity conditions at the antenna for
!              the given vectors inside and outside;
!
!dimension 2: always 2;
!1 ... coeffs to superpose modes inside;
!2 ... coeffs to superpose modes outside
!
!dimension 1: number of solution vectors to superpose;
!corresponds to SIZE(inner,1) and SIZE(outer,1)

use constants, only: dp, dpc, pi, im, c
use eqs_sett, only: Nwaves, Nfs, hom_sys
use antenna_data, only: ra, flag_debug
use background, only: flag_back
use mode_data, only: omov, det, wd_ptr
use maxwell_equations, only: num_eqs, num_vars, iErsp_state, dim_Ersp_state, iErsp_sys, names_state

implicit none;

complex(dpc), dimension(Nwaves, Nfs), intent(in)  :: inner, outer
complex(dpc), dimension(Nfs, 2), intent(out) :: coeffs

complex(dpc), dimension(Nwaves,2*Nfs) :: mat, tmp
complex(dpc), dimension(Nwaves) :: b, jumps, dev

complex(dpc), dimension(2) :: jsurf, jsurft

integer :: i, k, info
real(dp), parameter :: fpc = 4.0d0*pi/c

complex(dpc) :: delBp, delBs, cio

complex(dpc), dimension(Nwaves,Nwaves) :: Dmat
complex(dpc), dimension(1:Nwaves,0:1) :: Svec

complex(dpc), dimension(num_vars) :: v_sys
complex(dpc), dimension(Nwaves) :: v_state
complex(dpc), dimension(num_eqs) :: rhs

integer, dimension(Nwaves) :: IPIV

integer, dimension(1) :: dimmat;
real(dp) :: err;

!for right side of continuity equations
call current_density (jsurf);
call cyl2rsp (ra, jsurf(1), jsurf(2), jsurft(1), jsurft(2));

!approximated jumps of mfs: used for rhs of the system
delBs = fpc*jsurft(2)
delBp = -fpc*jsurft(1)

! print *, 'delBs ~ ', delBs
! print *, 'delBp ~ ', delBp

cio = c/im/omov

!to find the system: y' = Dy + S0*delta(r-ra) + S1*delta'(r-ra) + ...
!1 stage: D-matrix:

call calc_diff_sys_matrix (ra, flag_back, Dmat);

Svec = cmplx(0.0d0, 0.0d0, dpc)

!2 stage: Svec(:,0) - S0
rhs = cmplx(0.0d0, 0.0d0, dpc)
rhs(5) = -cio*delBp !coeff before delta
rhs(6) =  cio*delBs !coeff before delta

v_state = cmplx(0.0d0, 0.0d0, dpc);

call state2sys (ra, flag_back, v_state, v_sys, rhs);

!set derivatives:
do k = 1,3
    do i=0,dim_Ersp_state(k)-1
        Svec(iErsp_state(k)+i,0) = v_sys(iErsp_sys(k)+i+1)
    end do
end do

if(hom_sys /= 0) then
    if(num_eqs /= 12) then
        print *, 'error: antenna: num_eqs /= 12', num_eqs
    end if

    !3 stage: Svec(:,1) - S1
    rhs = cmplx(0.0d0, 0.0d0, dpc)
    rhs(11) = -cio*delBp  !coeff before delta'

    v_state = cmplx(0.0d0, 0.0d0, dpc)

    call state2sys (ra, flag_back, v_state, v_sys, rhs);

    !set derivatives:
    do k = 1,3
        do i=0,dim_Ersp_state(k)-1
            Svec(iErsp_state(k)+i,1) = v_sys(iErsp_sys(k)+i+1)
        end do
    end do
end if

jumps = Svec(:,0) + matmul(Dmat, Svec(:,1))

if(Nwaves /= 2*Nfs) then
    print *,'warning: antenna: antenna conditions are not implemented yet for this case!..'
end if

dimmat(1) = 2*Nfs*Nwaves;

!build coefficient matrix (put solution vectors in columns)
do i=1,Nfs
    mat(:,i) = -inner(:,i)
    mat(:,i+Nfs) = outer(:,i)
end do

tmp = mat; !for checking solution

b = jumps;

!solve general complex equation:
call zgesv (Nwaves, 1, mat, Nwaves, IPIV, b, Nwaves, info);

if (info /= 0) then
    print *, 'error: antenna: continuity: failed to solve the system: ierr=', info
end if

!check the solution:
if (flag_debug > 0) then
    print *
    print *, 'Jumps of the state vector at the antenna:'

    do i = 1,Nwaves
        print *, 'i = ',i, 'jump of ', names_state(i),' = ', jumps(i)
    end do

    print *
    print *, 'checking the solution of antenna system:'
    do i=1,Nwaves
        print *, 'i = ', i ,' : dev = ', sum(tmp(i,:)*b) - jumps(i)
    end do

    dev = matmul(tmp, b)-jumps
    err = maxval(abs(dev))/maxval(abs(jumps))/maxval(reshape(abs(tmp),dimmat))

    if (err > 1.0d-10) then
        print *
        print *, 'warning: antenna: solution error estimation: ', err
    end if

    print *
    do i=1,Nwaves
        print *, 'solution: i = ', i ,' : b(i) = ', b(i)
    end do

end if

coeffs(:,1) = b(1:Nfs)
coeffs(:,2) = b(Nfs+1:2*Nfs)

call calc_antenna_system_determinant (Nwaves, Nfs, inner, outer, det);

if (flag_debug > 0) then
    print *, 'det=', det
end if

call set_det_in_wd_struct (wd_ptr, real(det), aimag(det));

end subroutine

!------------------------------------------------------------------------------

subroutine calc_antenna_system_determinant (Nwaves, Nfs, inner, outer, det)

use constants, only: dpc

implicit none;

integer, intent(in) :: Nwaves, Nfs
complex(dpc), dimension(Nwaves, Nfs), intent(in)  :: inner, outer
complex(dpc), intent(out) :: det

complex(dpc), dimension(Nwaves,2*Nfs) :: mat
complex(dpc), dimension(Nwaves) :: b
complex(dpc), dimension(Nwaves, Nfs) :: in, out

integer :: i, info

integer, dimension(Nwaves) :: IPIV

!normalize fundamental solutions:
!internal:
mat = cmplx(0.0d0, 0.0d0, 8);
mat(1,1) = inner(1,1);
mat(1,2) = inner(1,2);
mat(2,1) = inner(3,1);
mat(2,2) = inner(3,2);

mat(3,3) = inner(1,1);
mat(3,4) = inner(1,2);
mat(4,3) = inner(3,1);
mat(4,4) = inner(3,2);

b(1) = 1.0d0
b(2) = 0.0d0
b(3) = 0.0d0
b(4) = 1.0d0

!solve general complex equation:
call zgesv (Nwaves, 1, mat, Nwaves, IPIV, b, Nwaves, info);

if (info /= 0) then
    print *, 'error: determinant: failed to solve the system: ierr=', info
end if

in(:,1) = b(1)*inner(:,1) + b(2)*inner(:,2)
in(:,2) = b(3)*inner(:,1) + b(4)*inner(:,2)

! print *, 'in1:'
! print *, in(:,1)
!
! print *, 'in2:'
! print *, in(:,2)

!external:
mat = 0.0d0
mat(1,1) = outer(1,1);
mat(1,2) = outer(1,2);
mat(2,1) = outer(3,1);
mat(2,2) = outer(3,2);

mat(3,3) = outer(1,1);
mat(3,4) = outer(1,2);
mat(4,3) = outer(3,1);
mat(4,4) = outer(3,2);

b(1) = 1.0d0
b(2) = 0.0d0
b(3) = 0.0d0
b(4) = 1.0d0

!solve general complex equation:
call zgesv (Nwaves, 1, mat, Nwaves, IPIV, b, Nwaves, info);

if (info /= 0) then
    print *, 'error: determinant: failed to solve the system: ierr=', info
end if

out(:,1) = b(1)*outer(:,1) + b(2)*outer(:,2)
out(:,2) = b(3)*outer(:,1) + b(4)*outer(:,2)

! print *, 'out1:'
! print *, out(:,1)
!
! print *, 'out2:'
! print *, out(:,2)

!build coefficient matrix (put solution vectors in columns)
do i=1,Nfs
    mat(:,i)     = -in(:,i)
    mat(:,i+Nfs) = out(:,i)
end do

!calculates determinant:
call calc_det_lu (Nwaves, mat, det, IPIV);

end subroutine

!------------------------------------------------------------------------------

subroutine continuity_mhd (Nwaves, Nfs, inner, outer, coeffs)

!This routine solves the continuity conditions at the antenna corresonding
!to the current in the antenna. The mode numbers are imported from the module "mode".

!inner ... matrix containing the solution vectors at the inner side of the antenna;
!dimension 2: number of solution vectors
!dimension 1: components per vector
!
!outer ... matrix containing the solution vectors at the outer side of
!the antenna; dimensions see inner
!
!coeffs ... matrix containing the coefficients for solution superposition
!              found by solving the continuity conditions at the antenna for
!              the given vectors inside and outside;
!
!dimension 2: always 2;
!1 ... coeffs to superpose modes inside;
!2 ... coeffs to superpose modes outside
!
!dimension 1: number of solution vectors to superpose;
!corresponds to SIZE(inner,1) and SIZE(outer,1)

use constants, only: dp, dpc, pi, im, c
use antenna_data, only: ra, flag_debug
use mode_data, only: omov, det, wd_ptr
use maxwell_equations, only: num_eqs, num_vars, iErsp_state, dim_Ersp_state, iErsp_sys, names_state

implicit none;

integer, intent(in) :: Nwaves, Nfs
complex(dpc), dimension(Nwaves, Nfs), intent(in) :: inner, outer
complex(dpc), dimension(Nfs, 2), intent(out) :: coeffs

complex(dpc), dimension(Nwaves,2*Nfs) :: mat, tmp
complex(dpc), dimension(Nwaves) :: b, jumps, dev

complex(dpc), dimension(2) :: jsurf, jsurft

integer :: i, info
real(dp), parameter :: fpc = 4.0d0*pi/c

complex(dpc) :: delBp, delBs, deldBr, cio

integer, dimension(Nwaves) :: IPIV

integer, dimension(1) :: dimmat;
real(dp) :: err;

real(dp) :: kt, kz, ks, kp, k2, kB;
complex(dpc) :: kA, kfac, det2;

dimmat(1) = 2*Nfs*Nwaves;

!for right side of continuity equations
call current_density (jsurf);
call cyl2rsp (ra, jsurf(1), jsurf(2), jsurft(1), jsurft(2));

!approximated jumps of mfs: used for rhs of the system
delBs = fpc*jsurft(2)
delBp = -fpc*jsurft(1)

! print *, 'delBs ~ ', delBs
! print *, 'delBp ~ ', delBp

cio = c/im/omov

!build coefficient matrix (put solution vectors in columns)
do i=1,Nfs
    mat(:,i) = -inner(:,i)
    mat(:,i+Nfs) = outer(:,i)
end do

tmp = mat; !for checking solution

det2 = tmp(1,1)*tmp(2,2) - tmp(1,2)*tmp(2,1);

!ks, kp, kA and other stuff:
call calc_k_vals_sub (ra, kt, kz, ks, kp, k2, kB, kA, kfac);

deldBr = -Im*(ks*delBs+kp*delBp);

! print *, kA, kp
! print *
! print *, deldBr, kfac, k2

!jumps = (/ cmplx(0.d0, 0.d0, 8), deldBr*ra*kfac/k2 /);
jumps = (/ cmplx(0.d0, 0.d0, 8), deldBr /); !another choice of vars in mhd eqn: Br, Br'

b = jumps;

!solve general complex equation:
call zgesv (Nwaves, 1, mat, Nwaves, IPIV, b, Nwaves, info);

if (info /= 0) then
    print *, 'error: antenna: continuity: failed to solve the system: ierr=', info
end if

!check the solution:
if (flag_debug > 0) then
    print *
    print *, 'antenna system matrix:'
    print *, tmp

    print *
    print *, 'jumps at the antenna:'
    do i=1,Nwaves
        print *, 'i = ', i ,' : jumps = ', jumps(i)
    end do

    print *
    print *, 'solution of the antenna system:'
    do i=1,Nwaves
        print *, 'i = ', i ,' : b = ', b(i)
    end do

    print *
    print *, 'checking the solution of the antenna system:'
    do i=1,Nwaves
        print *, 'i = ', i ,' : dev = ', sum(tmp(i,:)*b) - jumps(i)
    end do

    dev = matmul(tmp, b)-jumps
    err = maxval(abs(dev))/maxval(abs(jumps))/maxval(reshape(abs(tmp),dimmat))

    if (err > 1.0d-10) then
        print *
        print *, 'warning: antenna: solution error estimation: ', err
    end if

    print *
    do i=1,Nwaves
        print *, 'solution: i = ', i ,' : b(i) = ', b(i)
    end do

end if

coeffs(:,1) = b(1:Nfs)
coeffs(:,2) = b(Nfs+1:2*Nfs)

!calculate determinant:
call calc_det_lu (Nwaves, tmp, det, IPIV);

if (abs(det-det2) > 1.0d-10) then
    print *, 'dets are different:', det, det2
end if

if (flag_debug > 0) then
    print *, 'det=', det
end if

call set_det_in_wd_struct (wd_ptr, real(det), aimag(det));

end subroutine

!------------------------------------------------------------------------------

subroutine continuity_mhd_zeta (Nwaves, Nfs, inner, outer, coeffs)

!This routine solves the continuity conditions at the antenna corresonding
!to the current in the antenna. The mode numbers are imported from the module "mode".

!inner ... matrix containing the solution vectors at the inner side of the antenna;
!dimension 2: number of solution vectors
!dimension 1: components per vector
!
!outer ... matrix containing the solution vectors at the outer side of
!the antenna; dimensions see inner
!
!coeffs ... matrix containing the coefficients for solution superposition
!              found by solving the continuity conditions at the antenna for
!              the given vectors inside and outside;
!
!dimension 2: always 2;
!1 ... coeffs to superpose modes inside;
!2 ... coeffs to superpose modes outside
!
!dimension 1: number of solution vectors to superpose;
!corresponds to SIZE(inner,1) and SIZE(outer,1)

use constants, only: dp, dpc, pi, im, c
use antenna_data, only: ra, flag_debug
use mode_data, only: omov, det, wd_ptr
use maxwell_equations, only: num_eqs, num_vars, iErsp_state, dim_Ersp_state, iErsp_sys, names_state

implicit none;

integer, intent(in) :: Nwaves, Nfs
complex(dpc), dimension(Nwaves, Nfs), intent(in) :: inner, outer
complex(dpc), dimension(Nfs, 2), intent(out) :: coeffs

complex(dpc), dimension(Nwaves,2*Nfs) :: mat, tmp
complex(dpc), dimension(Nwaves) :: b, jumps, dev

complex(dpc), dimension(2) :: jsurf, jsurft

integer :: i, info
real(dp), parameter :: fpc = 4.0d0*pi/c

complex(dpc) :: delBp, delBs, deldBr, cio

integer, dimension(Nwaves) :: IPIV

integer, dimension(1) :: dimmat;
real(dp) :: err;

real(dp) :: kt, kz, ks, kp, k2, kB;
complex(dpc) :: kA, kfac, det2;

dimmat(1) = 2*Nfs*Nwaves;

!for right side of continuity equations
call current_density (jsurf);
call cyl2rsp (ra, jsurf(1), jsurf(2), jsurft(1), jsurft(2));

!approximated jumps of mfs: used for rhs of the system
delBs = fpc*jsurft(2)
delBp = -fpc*jsurft(1)

! print *, 'delBs ~ ', delBs
! print *, 'delBp ~ ', delBp

cio = c/im/omov

!build coefficient matrix (put solution vectors in columns)
do i=1,Nfs
    mat(:,i) = -inner(:,i)
    mat(:,i+Nfs) = outer(:,i)
end do

tmp = mat; !for checking solution

det2 = tmp(1,1)*tmp(2,2) - tmp(1,2)*tmp(2,1);

!ks, kp, kA and other stuff:
call calc_k_vals_sub (ra, kt, kz, ks, kp, k2, kB, kA, kfac);

deldBr = -Im*(ks*delBs+kp*delBp);

! print *, kA, kp
! print *
! print *, deldBr, kfac, k2, kB, deldBr/(Im*kB)

deldBr = deldBr/(Im*kB); !jump of zeta

!jumps = (/ cmplx(0.d0, 0.d0, 8), deldBr*ra*kfac/k2 /);
jumps = (/ cmplx(0.d0, 0.d0, 8), deldBr /); !another choice of vars in mhd eqn: Br, Br'

b = jumps;

!solve general complex equation:
call zgesv (Nwaves, 1, mat, Nwaves, IPIV, b, Nwaves, info);

if (info /= 0) then
    print *, 'error: antenna: continuity: failed to solve the system: ierr=', info
end if

!check the solution:
if (flag_debug > 0) then
    print *
    print *, 'antenna system matrix:'
    print *, tmp

    print *
    print *, 'jumps at the antenna:'
    do i=1,Nwaves
        print *, 'i = ', i ,' : jumps = ', jumps(i)
    end do

    print *
    print *, 'solution of the antenna system:'
    do i=1,Nwaves
        print *, 'i = ', i ,' : b = ', b(i)
    end do

    print *
    print *, 'checking the solution of the antenna system:'
    do i=1,Nwaves
        print *, 'i = ', i ,' : dev = ', sum(tmp(i,:)*b) - jumps(i)
    end do

    dev = matmul(tmp, b)-jumps
    err = maxval(abs(dev))/maxval(abs(jumps))/maxval(reshape(abs(tmp),dimmat))

    if (err > 1.0d-10) then
        print *
        print *, 'warning: antenna: solution error estimation: ', err
    end if

    print *
    do i=1,Nwaves
        print *, 'solution: i = ', i ,' : b(i) = ', b(i)
    end do

end if

coeffs(:,1) = b(1:Nfs)
coeffs(:,2) = b(Nfs+1:2*Nfs)

!calculate determinant:
call calc_det_lu (Nwaves, tmp, det, IPIV);

if (abs(det-det2) > 1.0d-10) then
    print *, 'dets are different:', det, det2
end if

if (flag_debug > 0) then
    print *, 'det=', det
end if

call set_det_in_wd_struct (wd_ptr, real(det), aimag(det));

end subroutine

!------------------------------------------------------------------------------

subroutine continuity_mhd_hi (Nwaves, Nfs, inner, outer, coeffs)

!This routine solves the continuity conditions at the antenna corresonding
!to the current in the antenna. The mode numbers are imported from the module "mode".

!inner ... matrix containing the solution vectors at the inner side of the antenna;
!dimension 2: number of solution vectors
!dimension 1: components per vector
!
!outer ... matrix containing the solution vectors at the outer side of
!the antenna; dimensions see inner
!
!coeffs ... matrix containing the coefficients for solution superposition
!              found by solving the continuity conditions at the antenna for
!              the given vectors inside and outside;
!
!dimension 2: always 2;
!1 ... coeffs to superpose modes inside;
!2 ... coeffs to superpose modes outside
!
!dimension 1: number of solution vectors to superpose;
!corresponds to SIZE(inner,1) and SIZE(outer,1)

use constants, only: dp, dpc, pi, im, c
use antenna_data, only: ra, flag_debug
use mode_data, only: omov, det, wd_ptr
use maxwell_equations, only: num_eqs, num_vars, iErsp_state, dim_Ersp_state, iErsp_sys, names_state

implicit none;

integer, intent(in) :: Nwaves, Nfs
complex(dpc), dimension(Nwaves, Nfs), intent(in) :: inner, outer
complex(dpc), dimension(Nfs, 2), intent(out) :: coeffs

complex(dpc), dimension(Nwaves,2*Nfs) :: mat, tmp
complex(dpc), dimension(Nwaves) :: b, jumps, dev

complex(dpc), dimension(2) :: jsurf, jsurft

integer :: i, info
real(dp), parameter :: fpc = 4.0d0*pi/c

complex(dpc) :: delBp, delBs, deldBr, cio

integer, dimension(Nwaves) :: IPIV

integer, dimension(1) :: dimmat;
real(dp) :: err;

real(dp) :: kt, kz, ks, kp, k2, kB;
complex(dpc) :: kA, kfac, det2;

dimmat(1) = 2*Nfs*Nwaves;

!for right side of continuity equations
call current_density (jsurf);
call cyl2rsp (ra, jsurf(1), jsurf(2), jsurft(1), jsurft(2));

!approximated jumps of mfs: used for rhs of the system
delBs = fpc*jsurft(2)
delBp = -fpc*jsurft(1)

! print *, 'delBs ~ ', delBs
! print *, 'delBp ~ ', delBp

cio = c/im/omov

!build coefficient matrix (put solution vectors in columns)
do i=1,Nfs
    mat(:,i) = -inner(:,i)
    mat(:,i+Nfs) = outer(:,i)
end do

tmp = mat; !for checking solution

det2 = tmp(1,1)*tmp(2,2) - tmp(1,2)*tmp(2,1);

!ks, kp, kA and other stuff:
call calc_k_vals_sub (ra, kt, kz, ks, kp, k2, kB, kA, kfac);

deldBr = -Im*(ks*delBs+kp*delBp);

! print *, kA, kp
! print *
! print *, deldBr, kfac, k2, kB, deldBr/(Im*kB)

deldBr = ra*deldBr/(Im*kB); !jump of (r*zeta)'

!jumps = (/ cmplx(0.d0, 0.d0, 8), deldBr*ra*kfac/k2 /);
jumps = (/ cmplx(0.d0, 0.d0, 8), deldBr /); !another choice of vars in mhd eqn: Br, Br'

b = jumps;

!solve general complex equation:
call zgesv (Nwaves, 1, mat, Nwaves, IPIV, b, Nwaves, info);

if (info /= 0) then
    print *, 'error: antenna: continuity: failed to solve the system: ierr=', info
end if

!check the solution:
if (flag_debug > 0) then
    print *
    print *, 'antenna system matrix:'
    print *, tmp

    print *
    print *, 'jumps at the antenna:'
    do i=1,Nwaves
        print *, 'i = ', i ,' : jumps = ', jumps(i)
    end do

    print *
    print *, 'solution of the antenna system:'
    do i=1,Nwaves
        print *, 'i = ', i ,' : b = ', b(i)
    end do

    print *
    print *, 'checking the solution of the antenna system:'
    do i=1,Nwaves
        print *, 'i = ', i ,' : dev = ', sum(tmp(i,:)*b) - jumps(i)
    end do

    dev = matmul(tmp, b)-jumps
    err = maxval(abs(dev))/maxval(abs(jumps))/maxval(reshape(abs(tmp),dimmat))

    if (err > 1.0d-10) then
        print *
        print *, 'warning: antenna: solution error estimation: ', err
    end if

    print *
    do i=1,Nwaves
        print *, 'solution: i = ', i ,' : b(i) = ', b(i)
    end do

end if

coeffs(:,1) = b(1:Nfs)
coeffs(:,2) = b(Nfs+1:2*Nfs)

!calculate determinant:
call calc_det_lu (Nwaves, tmp, det, IPIV);

if (abs(det-det2) > 1.0d-10) then
    print *, 'dets are different:', det, det2
end if

if (flag_debug > 0) then
    print *, 'det=', det
end if

call set_det_in_wd_struct (wd_ptr, real(det), aimag(det));

end subroutine

!------------------------------------------------------------------------------

subroutine get_jump_antenna (jumps)

use constants, only: dp, dpc, pi, im, c
use eqs_sett, only: Nwaves, Nfs, hom_sys
use antenna_data, only: ra, flag_debug
use background, only: flag_back
use mode_data, only: omov, det, wd_ptr
use maxwell_equations, only: num_eqs, num_vars, iErsp_state, dim_Ersp_state, iErsp_sys, names_state

implicit none;

complex(dpc), dimension(Nwaves), intent(out) :: jumps

complex(dpc), dimension(2) :: jsurf, jsurft

integer :: i, k
real(dp), parameter :: fpc = 4.0d0*pi/c

complex(dpc) :: delBp, delBs, cio

complex(dpc), dimension(Nwaves,Nwaves) :: Dmat
complex(dpc), dimension(1:Nwaves,0:1) :: Svec

complex(dpc), dimension(num_vars) :: v_sys
complex(dpc), dimension(Nwaves) :: v_state
complex(dpc), dimension(num_eqs) :: rhs

!for right side of continuity equations
call current_density (jsurf);
call cyl2rsp (ra, jsurf(1), jsurf(2), jsurft(1), jsurft(2));

!approximated jumps of mfs: used for rhs of the system
delBs = fpc*jsurft(2)
delBp = -fpc*jsurft(1)

! print *, 'delBs ~ ', delBs
! print *, 'delBp ~ ', delBp

cio = c/im/omov

!to find the system: y' = Dy + S0*delta(r-ra) + S1*delta'(r-ra) + ...
!1 stage: D-matrix:

call calc_diff_sys_matrix (ra, flag_back, Dmat);

Svec = cmplx(0.0d0, 0.0d0, dpc)

!2 stage: Svec(:,0) - S0
rhs = cmplx(0.0d0, 0.0d0, dpc)
rhs(5) = -cio*delBp !coeff before delta
rhs(6) =  cio*delBs !coeff before delta

v_state = cmplx(0.0d0, 0.0d0, dpc);

call state2sys (ra, flag_back, v_state, v_sys, rhs);

!set derivatives:
do k = 1,3
    do i=0,dim_Ersp_state(k)-1
        Svec(iErsp_state(k)+i,0) = v_sys(iErsp_sys(k)+i+1)
    end do
end do

if(hom_sys /= 0) then
    if(num_eqs /= 12) then
        print *, 'error: antenna: num_eqs /= 12', num_eqs
    end if

    !3 stage: Svec(:,1) - S1
    rhs = cmplx(0.0d0, 0.0d0, dpc)
    rhs(11) = -cio*delBp  !coeff before delta'

    v_state = cmplx(0.0d0, 0.0d0, dpc)

    call state2sys (ra, flag_back, v_state, v_sys, rhs);

    !set derivatives:
    do k = 1,3
        do i=0,dim_Ersp_state(k)-1
            Svec(iErsp_state(k)+i,1) = v_sys(iErsp_sys(k)+i+1)
        end do
    end do
end if

jumps = Svec(:,0) + matmul(Dmat, Svec(:,1))

end subroutine

!------------------------------------------------------------------------------
