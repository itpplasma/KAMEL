!<The subroutines used to impose stitching conditions at various boundaries for ideal MHD case.

!------------------------------------------------------------------------------

subroutine center_equations_imhd (Nw, D, code, EB, neq_, nvar_, M, J)

use constants, only: pp;

implicit none;

integer, parameter :: neq = 1, nvar = 2;

integer, intent(in) :: Nw, D;
integer(pp), intent(in) :: code;
complex(8), dimension(D,Nw), intent(in) :: EB;

integer, intent(out) :: neq_, nvar_;

complex(8), dimension(neq,nvar), intent(out) :: M;
complex(8), dimension(neq), intent(out) :: J;

if (Nw /= 2 .or. D /= 8 .or. neq /= Nw/2) then
    print *, 'error: center_equations_imhd: improper Nw or D values:', Nw, D;
    stop;
end if

!equations for the boundary: M c = J
neq_ = neq;
nvar_ = nvar;

M = cmplx (0.0d0, 0.0d0, 8);
J = cmplx (0.0d0, 0.0d0, 8);

!supress divergent functions: (actually, the second basis function is not computed)
M(1,2) = cmplx (1.0d0, 0.0d0, 8)*1.0d6;

end subroutine

!------------------------------------------------------------------------------

subroutine infinity_equations_imhd (Nw, D, code, EB, neq_, nvar_, M, J)

use constants, only: pp;

implicit none;

integer, parameter :: neq = 1, nvar = 2;

integer, intent(in) :: Nw, D;
integer(pp), intent(in) :: code;
complex(8), dimension(D,Nw), intent(in) :: EB;

integer, intent(out) :: neq_, nvar_;

complex(8), dimension(neq,nvar), intent(out) :: M;
complex(8), dimension(neq), intent(out) :: J;

if (Nw /= 2 .or. D /= 8 .or. neq /= Nw/2) then
    print *, 'error: infinity_equations_imhd: improper Nw or D values:', Nw, D;
    stop;
end if

!equations for the boundary: M c = J
neq_ = neq;
nvar_ = nvar;

M = cmplx (0.0d0, 0.0d0, 8);
J = cmplx (0.0d0, 0.0d0, 8);

!supress divergent functions: (actually, the first basis function is not computed)
M(1,1) = cmplx (1.0d0, 0.0d0, 8)*1.0d6;

end subroutine

!------------------------------------------------------------------------------

subroutine ideal_wall_equations_imhd (Nw, D, code, EB, neq_, nvar_, M, J)

use constants, only: pp;

implicit none;

integer, parameter :: neq = 1, nvar = 2;

integer, intent(in) :: Nw, D;
integer(pp), intent(in) :: code;
complex(8), dimension(D,Nw), intent(in) :: EB;

integer, intent(out) :: neq_, nvar_;

complex(8), dimension(neq,nvar), intent(out) :: M;
complex(8), dimension(neq), intent(out) :: J;

if (Nw /= 2 .or. D /= 8 .or. neq /= Nw/2) then
    print *, 'error: ideal_wall_equations_imhd: improper Nw or D values:', Nw, D;
    stop;
end if

!equations for the boundary: M c = J
neq_ = neq;
nvar_ = nvar;

M = cmplx (0.0d0, 0.0d0, 8);
J = cmplx (0.0d0, 0.0d0, 8);

!keep only basis with correct BC (Br = 0). Actually, the second basis function is not computed
M(1,2) = cmplx (1.0d0, 0.0d0, 8)*1.0d6;

end subroutine

!------------------------------------------------------------------------------

subroutine stitching_equations_imhd_imhd (Nw1, D1, code1, EB1, Nw2, D2, code2, EB2, flg_ant, neq_, nvar_, M, J)

use constants, only: pp;
use antenna_data, only: ra;

implicit none;

integer, parameter :: neq = 2, nvar = 4, noc = nvar/2;

integer, intent(in) :: Nw1, D1;
integer(pp), intent(in) :: code1;
complex(8), dimension(D1,Nw1), intent(in) :: EB1;

integer, intent(in) :: Nw2, D2;
integer(pp), intent(in) :: code2;
complex(8), dimension(D2,Nw2), intent(in) :: EB2;

integer, intent(in) :: flg_ant;

integer, intent(out) :: neq_, nvar_;

complex(8), dimension(neq,nvar), intent(out) :: M;
complex(8), dimension(neq), intent(out) :: J;

real(8), parameter :: pi  = 3.141592653589793238462;
real(8), parameter :: c   = 29979245800.0;
real(8), parameter :: fpc = 4.0d0*pi/c;
complex(8), parameter :: Im = cmplx(0.0d0, 1.0d0, 8);

!integer, dimension(neq) :: ind = (/5,6/); !indices of Bt, Bz in EB arrays
!integer, dimension(neq) :: ind = (/4,6/); !indices of Br, Bz in EB arrays
integer, dimension(neq) :: ind = (/7,8/);  !indices of r*zeta, (r*zeta)' in EB arrays

complex(8), dimension(2) :: jsurf;

integer :: k, l;

real(8) :: kt, kz, ks, kp, k2, kB;
complex(8) :: kA, kfac, deldBr;

if (Nw1 /= 2 .or. D1 /= 8 .or. Nw2 /= 2 .or. D2 /= 8 .or. neq /= (Nw1+Nw2)/2) then
    print *, 'error: stitching_equations_imhd_imhd: improper input values:', Nw1, D1, Nw2, D2;
    stop;
end if

!stitching equations for the boundary: M c = J
neq_ = neq;
nvar_ = nvar;

M = cmplx (0.0d0, 0.0d0, 8);
J = cmplx (0.0d0, 0.0d0, 8);

!equations matrix:
do k = 1,neq !over equations

    do l = 1,noc !over coefficients
        M(k,0  +l) = - EB1(ind(k),l);
        M(k,noc+l) =   EB2(ind(k),l);
    end do

end do

!rhs vector:
if (flg_ant /= 0) then

    call current_density (jsurf);

    !J(1) = cmplx(0.0d0, 0.0d0, 8); !Br

    J(1) =   fpc*jsurf(2); !delta Bt
    J(2) = - fpc*jsurf(1); !delta Bz

    call calc_k_vals_sub (ra, kt, kz, ks, kp, k2, kB, kA, kfac);

    deldBr = - Im*(kt*J(1) + kz*J(2)); !delta Br'

    J(1) = cmplx(0.0d0, 0.0d0, 8); !delta (r*zeta)
    J(2) = ra*deldBr/(Im*kB);      !delta (r*zeta)'

end if

end subroutine

!------------------------------------------------------------------------------

subroutine stitching_equations_imhd_hommed (Nw1, D1, code1, EB1, Nw2, D2, code2, EB2, flg_ant, neq_, nvar_, M, J)

use constants, only: pp;

implicit none;

integer, parameter :: neq = 3, nvar = 6;

integer, intent(in) :: Nw1, D1;
integer(pp), intent(in) :: code1;
complex(8), dimension(D1,Nw1), intent(in) :: EB1;

integer, intent(in) :: Nw2, D2;
integer(pp), intent(in) :: code2;
complex(8), dimension(D2,Nw2), intent(in) :: EB2;

integer, intent(in) :: flg_ant;

integer, intent(out) :: neq_, nvar_;

complex(8), dimension(neq,nvar), intent(out) :: M;
complex(8), dimension(neq), intent(out) :: J;

real(8), parameter :: pi  = 3.141592653589793238462;
real(8), parameter :: c   = 29979245800.0;
real(8), parameter :: fpc = 4.0d0*pi/c;

integer, dimension(3) :: ind = (/4,5,6/); !indices of Br, Bt, Bz in EB arrays
complex(8), dimension(2) :: jsurf;

integer :: k, l;

if (Nw1 /= 2 .or. D1 /= 8 .or. Nw2 /= 4 .or. D2 /= 6 .or. neq /= (Nw1+Nw2)/2) then
    print *, 'error: stitching_equations_imhd_hommed: improper input values:', Nw1, D1, Nw2, D2;
    stop;
end if

!stitching equations for the boundary: M c = J (continuity of Br, Bt, Bz without antenna)
neq_ = neq;
nvar_ = nvar;

M = cmplx (0.0d0, 0.0d0, 8);
J = cmplx (0.0d0, 0.0d0, 8);

!equations matrix:
do k = 1,neq !over equations

    do l = 1,2 !over imhd coefficients
        M(k,0+l) = - EB1(ind(k),l);
    end do

    do l = 1,4 !over vac coefficients
        M(k,2+l) =   EB2(ind(k),l);
    end do

end do

!rhs vector:
if (flg_ant /= 0) then

    call current_density (jsurf);

    J(1) =   cmplx (0.0, 0.0, 8); !Br
    J(2) =   fpc*jsurf(2);        !Bt
    J(3) = - fpc*jsurf(1);        !Bz

end if

end subroutine

!------------------------------------------------------------------------------

module imhd_data

use constants, only: pp;

integer(pp) :: zone_ptr;

end module

!------------------------------------------------------------------------------

subroutine set_imhd_data_module (zone)

use constants, only: pp;
use imhd_data, only: zone_ptr;

implicit none;

integer(pp), intent(in) :: zone;

zone_ptr = zone;

end subroutine

!------------------------------------------------------------------------------

subroutine calc_k_vals_sub (r, kt, kz, ks, kp, k2, kB, kA, kfac)

use constants, only: dp, dpc, im
use imhd_data, only: zone_ptr;

implicit none;

real(dp), intent(in) :: r;
real(dp), intent(out) :: kt, kz, ks, kp, k2, kB;
complex(dpc), intent(out) :: kA, kfac;

real(dp) :: re_kA, im_kA, re_kfac, im_kfac;

call calc_k_vals (zone_ptr, r, kt, kz, ks, kp, k2, kB, re_kA, im_kA, re_kfac, im_kfac);

kA = re_kA + im*im_kA;

kfac = re_kfac + im*im_kfac;

end subroutine

!------------------------------------------------------------------------------

subroutine normalize_imhd_basis (D, Nw, dim, basis, ind, node)

implicit none;

integer, intent(in) :: D, Nw, dim, ind, node;
complex(8), dimension(D,Nw,dim) :: basis;
complex(8) :: factor;
integer :: k;

factor = basis(4,ind+1,node+1); !Br at the stitching point

do k = 1,dim

    basis(:,ind+1,k) = basis(:,ind+1,k)/factor;

end do

end subroutine

!------------------------------------------------------------------------------
