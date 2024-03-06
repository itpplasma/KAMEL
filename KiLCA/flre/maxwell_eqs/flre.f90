!<The subroutines used to impose stitching conditions at various boundaries for flre case.

!------------------------------------------------------------------------------

subroutine center_equations_flre (Nw, D, code, EB, neq_, nvar_, M, J)

use constants, only: pp;

implicit none;

integer :: neq, nvar, noc;

integer, intent(in) :: Nw, D;
integer(pp), intent(in) :: code;
complex(8), dimension(D,Nw), intent(in) :: EB;

integer, intent(out) :: neq_, nvar_;

complex(8), dimension(Nw/2,Nw), intent(out) :: M;
complex(8), dimension(Nw/2), intent(out) :: J;

integer :: k;

neq = Nw/2;
nvar = Nw;
noc = nvar/2;

!equations for the boundary: M c = J
neq_ = neq;
nvar_ = nvar;

M = cmplx (0.0d0, 0.0d0, 8);
J = cmplx (0.0d0, 0.0d0, 8);

!supress divergent functions: (actually, these basis functions having indices Nwaves/2+1 : Nwaves are not computed)

do k = 1,neq !over equations
    M(k,noc+k) = cmplx (1.0d0, 0.0d0, 8)*1.0d0;
end do

end subroutine

!------------------------------------------------------------------------------

subroutine infinity_equations_flre (Nw, D, code, EB, neq_, nvar_, M, J)

use constants, only: pp;

implicit none;

integer :: neq, nvar, noc;

integer, intent(in) :: Nw, D;
integer(pp), intent(in) :: code;
complex(8), dimension(D,Nw), intent(in) :: EB;

integer, intent(out) :: neq_, nvar_;

complex(8), dimension(Nw/2,Nw), intent(out) :: M;
complex(8), dimension(Nw/2), intent(out) :: J;

integer :: k;

neq = Nw/2;
nvar = Nw;
noc = nvar/2;

!equations for the boundary: M c = J
neq_ = neq;
nvar_ = nvar;

M = cmplx (0.0d0, 0.0d0, 8);
J = cmplx (0.0d0, 0.0d0, 8);

!supress divergent functions: (actually, these basis functions having indices 1 : Nwaves/2 are not computed)

do k = 1,neq !over equations
    M(k,k) = cmplx (1.0d0, 0.0d0, 8)*1.0d0;
end do

end subroutine

!------------------------------------------------------------------------------

subroutine ideal_wall_equations_flre (Nw, D, code, EB, neq_, nvar_, M, J)

use constants, only: pp;

implicit none;

integer :: neq, nvar, noc;

integer, intent(in) :: Nw, D;
integer(pp), intent(in) :: code;
complex(8), dimension(D,Nw), intent(in) :: EB;

integer, intent(out) :: neq_, nvar_;

complex(8), dimension(Nw/2,Nw), intent(out) :: M;
complex(8), dimension(Nw/2), intent(out) :: J;

integer :: k;

neq = Nw/2;
nvar = Nw;
noc = nvar/2;

!equations for the boundary: M c = J
neq_ = neq;
nvar_ = nvar;

M = cmplx (0.0d0, 0.0d0, 8);
J = cmplx (0.0d0, 0.0d0, 8);

!supress divergent functions: (actually, these basis functions having indices Nwaves/2+1 : Nwaves are not computed)

do k = 1,neq !over equations
    M(k,noc+k) = cmplx (1.0d0, 0.0d0, 8)*1.0d0;
end do

end subroutine

!------------------------------------------------------------------------------

subroutine stitching_equations_flre_flre (Nw1, D1, code1, EB1, Nw2, D2, code2, EB2, flg_ant, neq_, nvar_, M, J)

use constants, only: pp;

use maxwell_equations, only: sys_ind;

implicit none;

integer :: neq, nvar, noc;

integer, intent(in) :: Nw1, D1;
integer(pp), intent(in) :: code1;
complex(8), dimension(D1,Nw1), intent(in) :: EB1;

integer, intent(in) :: Nw2, D2;
integer(pp), intent(in) :: code2;
complex(8), dimension(D2,Nw2), intent(in) :: EB2;

integer, intent(in) :: flg_ant;

integer, intent(out) :: neq_, nvar_;

complex(8), dimension((Nw1+Nw2)/2, Nw1+Nw2), intent(out) :: M;
complex(8), dimension((Nw1+Nw2)/2), intent(out) :: J;

!indices of state fields in sys fields array for zone 1:
integer, dimension(Nw1) :: sys_ind_1;

integer :: k, l;

neq = (Nw1+Nw2)/2;
nvar = Nw1+Nw2;
noc = Nw1;

!stitching equations for the boundary: M c = J
neq_ = neq;
nvar_ = nvar;

M = cmplx (0.0d0, 0.0d0, 8);
J = cmplx (0.0d0, 0.0d0, 8);

!below, the same order of flre is assumed in 2 zones:
if (Nw1 /= Nw2) then
    print *, 'error: stitching_equations_flre_flre: flre zones are different:', Nw1, Nw2;
    stop;
end if

!left zone:
call get_sys_ind_array (code1, sys_ind_1);

!equations matrix:
do k = 1,neq !over equations, or state vector components

    do l = 1,noc !over coefficients
        M(k,0  +l) = - EB1(sys_ind_1(k),l);
        M(k,noc+l) =   EB2(sys_ind_1(k),l);
    end do

end do

!rhs vector:
call activate_fortran_modules_for_zone (code1);

if (flg_ant /= 0) then

    call calc_flre_state_jump_at_antenna (J);

end if

call deactivate_fortran_modules_for_zone (code1);

end subroutine

!------------------------------------------------------------------------------

subroutine stitching_equations_flre_hommed (Nw1, D1, code1, EB1, Nw2, D2, code2, EB2, flg_ant, neq_, nvar_, M, J)

use constants, only: pp;

implicit none;

integer, intent(in) :: Nw1, D1;
integer(pp), intent(in) :: code1;
complex(8), dimension(D1,Nw1), intent(in) :: EB1;

integer, intent(in) :: Nw2, D2;
integer(pp), intent(in) :: code2;
complex(8), dimension(D2,Nw2), intent(in) :: EB2;

integer, intent(in) :: flg_ant;

integer, intent(out) :: neq_, nvar_;

complex(8), dimension((Nw1+Nw2)/2, Nw1+Nw2), intent(out) :: M;
complex(8), dimension((Nw1+Nw2)/2), intent(out) :: J;

if (Nw1 == 4) then
    call stitching_equations_flre_N1_hommed (Nw1, D1, code1, EB1, Nw2, D2, code2, EB2, flg_ant, neq_, nvar_, M, J);
else
    call stitching_equations_flre_NN_hommed (Nw1, D1, code1, EB1, Nw2, D2, code2, EB2, flg_ant, neq_, nvar_, M, J);
end if

end subroutine

!------------------------------------------------------------------------------

subroutine stitching_equations_flre_N1_hommed (Nw1, D1, code1, EB1, Nw2, D2, code2, EB2, flg_ant, neq_, nvar_, M, J)

use constants, only: pp;

implicit none;

integer :: neq, nvar, noc;

integer, intent(in) :: Nw1, D1;
integer(pp), intent(in) :: code1;
complex(8), dimension(D1,Nw1), intent(in) :: EB1;

integer, intent(in) :: Nw2, D2;
integer(pp), intent(in) :: code2;
complex(8), dimension(D2,Nw2), intent(in) :: EB2;

integer, intent(in) :: flg_ant;

integer, intent(out) :: neq_, nvar_;

complex(8), dimension((Nw1+Nw2)/2, Nw1+Nw2), intent(out) :: M;
complex(8), dimension((Nw1+Nw2)/2), intent(out) :: J;

integer, dimension(4) :: ind_med = (/2,3,5,6/); !indices of Et, Ez, Bt, Bz in EB hmed arrays
integer, dimension(4) :: ind_flre;              !indices of Et, Ez, Bt, Bz in EB flre arrays

real(8) :: r;                           !boundary position
complex(8), dimension(D1,Nw1) :: EB1cl; !flre basis in the lab cyl frame

integer, dimension(3) :: iErsp_sys, iBrsp_sys;

integer :: k, l;

if (Nw1 /= 4 .or. D1 /= 13 .or. Nw2 /= 4 .or. D2 /= 6) then
    print *, 'error: stitching_equations_flre_N1_hommed: improper input values:', Nw1, D1, Nw2, D2;
    stop;
end if

if (flg_ant == 1) then
    print *, 'error: stitching_equations_flre_N1_hommed: antenna flag must be zero for this boundary:', flg_ant;
    stop;
end if

neq = (Nw1+Nw2)/2;
nvar = Nw1+Nw2;
noc = Nw1;

!stitching equations for the boundary: M c = J (continuity of Et, Ez, Bt, Bz)
neq_ = neq;
nvar_ = nvar;

M = cmplx (0.0d0, 0.0d0, 8);
J = cmplx (0.0d0, 0.0d0, 8);

!calc basis system vectors in lab cyl frame:
call get_right_boundary_of_zone (code1, r);

call calc_flre_basis_in_lab_cyl_frame_with_full_system_vectors (code1, r, EB1, EB1cl);

call get_iErsp_sys_array (code1, iErsp_sys);
call get_iBrsp_sys_array (code1, iBrsp_sys);

ind_flre(1) = iErsp_sys(2);
ind_flre(2) = iErsp_sys(3);
ind_flre(3) = iBrsp_sys(2);
ind_flre(4) = iBrsp_sys(3);

!equations matrix: (N=1 only)
do k = 1,4 !over equations: Et, Ez, Bt, Bz

    do l = 1,Nw1 !over flre coefficients
        M(k,0+l) = - EB1cl(ind_flre(k),l);
    end do

    do l = 1,4 !over hommedium  coefficients
        M(k,Nw1+l) = EB2(ind_med(k),l);
    end do

end do

!no antenna is assumed at this boundary

end subroutine

!------------------------------------------------------------------------------

subroutine stitching_equations_flre_NN_hommed (Nw1, D1, code1, EB1, Nw2, D2, code2, EB2, flg_ant, neq_, nvar_, M, J)

use constants, only: pp;

implicit none;

integer, intent(in) :: Nw1, D1;
integer(pp), intent(in) :: code1;
complex(8), dimension(D1,Nw1), intent(in) :: EB1;

integer, intent(in) :: Nw2, D2;
integer(pp), intent(in) :: code2;
complex(8), dimension(D2,Nw2), intent(in) :: EB2;

integer, intent(in) :: flg_ant;

integer, intent(out) :: neq_, nvar_;

complex(8), dimension((Nw1+Nw2)/2, Nw1+Nw2), intent(out) :: M;
complex(8), dimension((Nw1+Nw2)/2), intent(out) :: J;

print *, 'error: not implemented!..'

end subroutine

!------------------------------------------------------------------------------

subroutine calc_flre_state_jump_at_antenna (jumps)

!everything is in moving frame!

use constants, only: pi, im, c;
use flre_sett, only: Nwaves, hom_sys;
use antenna_data, only: ra;
use background, only: flag_back;
use mode_data, only: omov;
use maxwell_equations, only: num_eqs, num_vars, iErsp_state, dim_Ersp_state, iErsp_sys;

implicit none;

complex(8), dimension(Nwaves), intent(out) :: jumps;

complex(8), dimension(2) :: jsurf, jsurft;

integer :: i, k;
real(8), parameter :: fpc = 4.0d0*pi/c;

complex(8) :: delBp, delBs, cio;

complex(8), dimension(Nwaves,Nwaves) :: Dmat;
complex(8), dimension(1:Nwaves,0:1) :: Svec;

complex(8), dimension(num_vars) :: v_sys;
complex(8), dimension(Nwaves) :: v_state;
complex(8), dimension(num_eqs) :: rhs;

!for right side of continuity equations
call current_density (jsurf);
call cyl2rsp (ra, jsurf(1), jsurf(2), jsurft(1), jsurft(2));

!approximated jumps of mfs: used for rhs of the system
delBs =   fpc*jsurft(2);
delBp = - fpc*jsurft(1);

cio = c/im/omov;

!to find the system: y' = Dy + S0*delta(r-ra) + S1*delta'(r-ra) + ...

!Svec(:,0) - S0
rhs = cmplx(0.0d0, 0.0d0, 8);

rhs(5) = - cio*delBp !coeff before delta
rhs(6) =   cio*delBs !coeff before delta

v_state = cmplx(0.0d0, 0.0d0, 8);

call state2sys (ra, flag_back, v_state, v_sys, rhs);

!set derivatives:
Svec = cmplx(0.0d0, 0.0d0, 8);

do k = 1,3
    do i=0,dim_Ersp_state(k)-1
        Svec(iErsp_state(k)+i,0) = v_sys(iErsp_sys(k)+i+1);
    end do
end do

jumps = Svec(:,0);

if (hom_sys == 0) then
    return;
end if

!homogenious plasma:
if(num_eqs /= 12) then
    print *, 'error: calc_flre_state_jump_at_antenna: num_eqs /= 12', num_eqs;
    stop;
end if

call calc_diff_sys_matrix (ra, flag_back, Dmat);

!Svec(:,1) - S1
rhs = cmplx(0.0d0, 0.0d0, 8);

rhs(11) = -cio*delBp;  !coefficient before delta'

v_state = cmplx(0.0d0, 0.0d0, 8);

call state2sys (ra, flag_back, v_state, v_sys, rhs);

!set derivatives:
do k = 1,3
    do i=0,dim_Ersp_state(k)-1
        Svec(iErsp_state(k)+i,1) = v_sys(iErsp_sys(k)+i+1);
    end do
end do

jumps = jumps + matmul(Dmat, Svec(:,1));

end subroutine

!------------------------------------------------------------------------------

subroutine normalize_flre_basis (D, Nw, dim, basis, iErsp_sys, ind1, ind2, node)

implicit none;

integer, intent(in) :: D, Nw, dim, ind1, ind2, node;
complex(8), dimension(D,Nw,dim) :: basis;
integer, dimension(3) :: iErsp_sys;

complex(8), allocatable, dimension(:,:) :: mat, b, new;
integer, allocatable, dimension(:) :: IPIV;

integer :: i, k, Nfs, info;

if (Nw /= 4) then !only first order of flre is assumed at the moment
    return;
end if

Nfs = Nw/2;

if (ind2-ind1+1 /= Nfs) then
    print *, 'normalize_flre_basis: wrong ind1 or ind2:', ind1, ind2;
    return;
end if

!below everything is hardcoded for N=1 expansion

allocate (mat(Nfs,Nfs), b(Nfs,Nfs), IPIV(Nfs), new(D,Nfs));

mat = cmplx(0.0d0, 0.0d0, 8);

do i = 1,Nfs
    mat(1,i) = basis(iErsp_sys(2)+1,ind1+i,node+1); !Es
    mat(2,i) = basis(iErsp_sys(3)+1,ind1+i,node+1); !Ep
end do

b(1,1) = 1.0d0; b(1,2) = 0.0d0;
b(2,1) = 0.0d0; b(2,2) = 1.0d0;

!solve general complex equation:
call zgesv (Nfs, Nfs, mat, Nfs, IPIV, b, Nfs, info);

if (info /= 0) then
    print *, 'error: normalize_flre_basis: failed to solve the system: ierr=', info;
end if

!basis normalization:
do k = 1,dim
    do i = 1,Nfs
        new(:,i) = b(1,i)*basis(:,ind1+1,k) + b(2,i)*basis(:,ind1+2,k);
    end do

    do i = 1,Nfs
        basis(:,ind1+i,k) = new(:,i);
    end do
end do

deallocate (mat, b, IPIV, new);

end subroutine

!------------------------------------------------------------------------------
