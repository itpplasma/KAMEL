!<The subroutines used to calculate vacuum fields and impose stitching conditions at various boundaries for vacuum case.

!------------------------------------------------------------------------------

subroutine eval_fields_in_hom_media (m, kz, omega, sigma, r, S, EB)

!evaluates solution of ME in a media with constant conductivity sigma in cyl geometry

implicit none;

integer, intent(in) :: m;  !poloidal wave number
real(8), intent(in) :: kz; !n/R, n - toroidal wave number, R - big radius of the torus
complex(8), intent(in) :: omega, sigma;
real(8), intent(in) :: r;
complex(8), dimension(4), intent(in) :: S; !superposition coefficients

complex(8), dimension(6), intent(out) :: EB;

complex(8), parameter :: I = (0.0d0, 1.0d0);

real(8), parameter :: c = 29979245800.0;

real(8), parameter :: pi = 3.1415926535897932384626;

complex(8) :: agemo, gamma;

complex(8) :: Im, Km, dIm, dKm;

real(8) :: kt;

!Modified Bessel functions:
real(8) :: ZR, ZI, FNU;
integer :: KODE, NZ, IERR;
real(8), dimension(-1:1) :: CYR, CYI;

complex(8) :: gr;

agemo = omega + 4.0d0*pi*I*sigma;

gamma = sqrt(kz*kz - omega*agemo/c/c);

gr = gamma*r;

kt = m/r;

ZR = real(gr);
ZI = aimag(gr);

FNU = abs(m)-1; ! this is wrong for m==0 !!!
KODE = 1;

call zbesi (ZR, ZI, FNU, KODE, 3, CYR, CYI, NZ, IERR);

!if (IERR /= 0) then
!    print *, 'besseli: warning: ierr =', IERR, 'NZ =', NZ, 'arg=', ZR + I*ZI;
!end if

Im = cmplx(CYR(0), CYI(0), 8);
dIm = 0.5*gamma*(cmplx(CYR(-1), CYI(-1), 8) + cmplx(CYR(1), CYI(1), 8));

call zbesk (ZR, ZI, FNU, KODE, 3, CYR, CYI, NZ, IERR);

!if (IERR /= 0) then
!    print *, 'besselk: warning: ierr =', IERR, 'NZ =', NZ, 'arg=', ZR + I*ZI;
!end if

Km = cmplx(CYR(0), CYI(0), 8);

dKm = - 0.5*gamma*(cmplx(CYR(-1), CYI(-1), 8) + cmplx(CYR(1), CYI(1), 8));

!fields:
EB(1) = -I*kz*(S(1)*dIm + S(3)*dKm) + kt*omega/c*(S(2)*Im + S(4)*Km);

EB(2) = kt*kz*(S(1)*Im + S(3)*Km) + I*omega/c*(S(2)*dIm + S(4)*dKm);

EB(3) = gamma*gamma*(S(1)*Im + S(3)*Km);

EB(4) = -kt*agemo/c*(S(1)*Im + S(3)*Km) - I*kz*(S(2)*dIm + S(4)*dKm);

EB(5) = -I*agemo/c*(S(1)*dIm + S(3)*dKm) + kt*kz*(S(2)*Im + S(4)*Km);

EB(6) = gamma*gamma*(S(2)*Im + S(4)*Km);

end subroutine

!------------------------------------------------------------------------------

subroutine eval_basis_in_hom_media (m, kz, omega, sigma, r, EB)

!evaluates basis solution of ME in a media with constant conductivity sigma in cyl geometry
!there are 4 independent solutions with the state vector (E_theta, Ez, B_theta, B_z)
!the function evaluates all (6) field components of all (4) basis solutions
!generally, a moving frame is assumed

implicit none;

integer, intent(in) :: m;  !poloidal wave number
real(8), intent(in) :: kz; !n/R, n - toroidal wave number, R - big radius of the torus
complex(8), intent(in) :: omega, sigma;
real(8), intent(in) :: r;

complex(8), dimension(6,4), intent(out) :: EB;

complex(8), parameter :: I = (0.0d0, 1.0d0);

real(8), parameter :: c = 29979245800.0;

real(8), parameter :: pi = 3.1415926535897932384626;

complex(8) :: agemo, gamma;

complex(8) :: Im, Km, dIm, dKm;

real(8) :: kt;

!Modified Bessel functions:
real(8) :: ZR, ZI, FNU;
integer :: KODE, NZ, IERR, Nnu;
real(8), dimension(-1:1) :: CYR, CYI;

complex(8) :: gr;

agemo = omega + 4.0d0*pi*I*sigma;

gamma = sqrt(kz*kz - omega*agemo/c/c);

gr = gamma*r;

kt = m/r;

ZR = real(gr);
ZI = aimag(gr);

if (m /= 0) then
    FNU = abs(m)-1;
    Nnu = 3;
else
    FNU = 0;
    Nnu = 2;
end if

KODE = 1;

call zbesi (ZR, ZI, FNU, KODE, Nnu, CYR, CYI, NZ, IERR);

!if (IERR /= 0) then
!    print *, 'besseli: warning: ierr =', IERR, 'NZ =', NZ, 'arg=', ZR + I*ZI;
!end if

Im = cmplx(CYR(0), CYI(0), 8);
dIm = 0.5*gamma*(cmplx(CYR(-1), CYI(-1), 8) + cmplx(CYR(1), CYI(1), 8));

if (m == 0) then
    Im = cmplx(CYR(-1), CYI(-1), 8);
    dIm = 0.5*gamma*(2.0*cmplx(CYR(0), CYI(0), 8));
end if

call zbesk (ZR, ZI, FNU, KODE, Nnu, CYR, CYI, NZ, IERR);

!if (IERR /= 0) then
!    print *, 'besselk: warning: ierr =', IERR, 'NZ =', NZ, 'arg=', ZR + I*ZI;
!end if

Km = cmplx(CYR(0), CYI(0), 8);
dKm = - 0.5*gamma*(cmplx(CYR(-1), CYI(-1), 8) + cmplx(CYR(1), CYI(1), 8));

if (m == 0) then
    Km = cmplx(CYR(-1), CYI(-1), 8);
    dKm = - 0.5*gamma*(2.0*cmplx(CYR(0), CYI(0), 8));
end if

!basis fields: take common factors off!

EB(1,1) = - I*kz*dIm;
EB(1,2) =   kt*omega/c*Im;
EB(1,3) = - I*kz*dKm;
EB(1,4) =   kt*omega/c*Km;

EB(2,1) =   kt*kz*Im;
EB(2,2) =   I*omega/c*dIm;
EB(2,3) =   kt*kz*Km;
EB(2,4) =   I*omega/c*dKm;

EB(3,1) =   gamma*gamma*Im;
EB(3,2) =   cmplx(0.0d0, 0.0d0, 8);
EB(3,3) =   gamma*gamma*Km;
EB(3,4) =   cmplx(0.0d0, 0.0d0, 8);

EB(4,1) = - kt*agemo/c*Im;
EB(4,2) = - I*kz*dIm;
EB(4,3) = - kt*agemo/c*Km;
EB(4,4) = - I*kz*dKm;

EB(5,1) = - I*agemo/c*dIm;
EB(5,2) =   kt*kz*Im;
EB(5,3) = - I*agemo/c*dKm;
EB(5,4) =   kt*kz*Km;

EB(6,1) =   cmplx(0.0d0, 0.0d0, 8);
EB(6,2) =   gamma*gamma*Im;
EB(6,3) =   cmplx(0.0d0, 0.0d0, 8);
EB(6,4) =   gamma*gamma*Km;

end subroutine

!------------------------------------------------------------------------------

subroutine center_equations_hommed (Nw, D, code, EB, neq_, nvar_, M, J)

use constants, only: pp;

implicit none;

integer, parameter :: neq = 2, nvar = 4;

integer, intent(in) :: Nw, D;
integer(pp), intent(in) :: code;
complex(8), dimension(D,Nw), intent(in) :: EB;

integer, intent(out) :: neq_, nvar_;

complex(8), dimension(neq,nvar), intent(out) :: M;
complex(8), dimension(neq), intent(out) :: J;

if (Nw /= 4 .or. D /= 6 .or. neq /= Nw/2) then
    print *, 'error: center_equations_hommed: improper Nw or D values:', Nw, D;
    stop;
end if

!equations for the boundary: M c = J
neq_ = neq;
nvar_ = nvar;

M = cmplx (0.0d0, 0.0d0, 8);
J = cmplx (0.0d0, 0.0d0, 8);

!supress K functions:
M(1,3) = cmplx (1.0d0, 0.0d0, 8)*EB(1,3); !normalization: to get S exactly 0.0
M(2,4) = cmplx (1.0d0, 0.0d0, 8)*EB(1,3); !normalization: to get S exactly 0.0

end subroutine

!------------------------------------------------------------------------------

subroutine infinity_equations_hommed (Nw, D, code, EB, neq_, nvar_, M, J)

use constants, only: pp;

implicit none;

integer, parameter :: neq = 2, nvar = 4;

integer, intent(in) :: Nw, D;
integer(pp), intent(in) :: code;
complex(8), dimension(D,Nw), intent(in) :: EB;

integer, intent(out) :: neq_, nvar_;

complex(8), dimension(neq,nvar), intent(out) :: M;
complex(8), dimension(neq), intent(out) :: J;

if (Nw /= 4 .or. D /= 6 .or. neq /= Nw/2) then
    print *, 'error: infinity_equations_hommed: improper Nw or D values:', Nw, D;
    stop;
end if

!equations for the boundary: M c = J
neq_ = neq;
nvar_ = nvar;

M = cmplx (0.0d0, 0.0d0, 8);
J = cmplx (0.0d0, 0.0d0, 8);

!supress I functions:
M(1,1) = cmplx (1.0d0, 0.0d0, 8)*EB(1,1); !normalization: to get S exactly 0.0
M(2,2) = cmplx (1.0d0, 0.0d0, 8)*EB(1,1); !normalization: to get S exactly 0.0

end subroutine

!------------------------------------------------------------------------------

subroutine ideal_wall_equations_hommed (Nw, D, code, EB, neq_, nvar_, M, J)

use constants, only: pp;

implicit none;

integer, parameter :: neq = 2, nvar = 4;

integer, intent(in) :: Nw, D;
integer(pp), intent(in) :: code;
complex(8), dimension(D,Nw), intent(in) :: EB;

integer, intent(out) :: neq_, nvar_;

complex(8), dimension(neq,nvar), intent(out) :: M;
complex(8), dimension(neq), intent(out) :: J;

integer, dimension(2) :: ind = (/2,3/); !indices of Et, Ez, Bt, Bz in EB arrays

integer :: k, l;

if (Nw /= 4 .or. D /= 6 .or. neq /= Nw/2) then
    print *, 'error: ideal_wall_equations_hommed: improper Nw or D values:', Nw, D;
    stop;
end if

!equations for the boundary: M c = J
neq_ = neq;
nvar_ = nvar;

M = cmplx (0.0d0, 0.0d0, 8);
J = cmplx (0.0d0, 0.0d0, 8);

!equations matrix:
do k = 1,neq !over equations

    do l = 1,nvar !over coefficients
        M(k,l) = EB(ind(k),l);
    end do

end do

end subroutine

!------------------------------------------------------------------------------

subroutine stitching_equations_hommed_hommed (Nw1, D1, code1, EB1, Nw2, D2, code2, EB2, flg_ant, neq_, nvar_, M, J)

use constants, only: pp;

implicit none;

integer, parameter :: neq = 4, nvar = 8, noc = nvar/2;

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

integer, dimension(4) :: ind = (/2,3,5,6/); !indices of Et, Ez, Bt, Bz in EB arrays
complex(8), dimension(2) :: jsurf;

integer :: k, l;

if (Nw1 /= 4 .or. D1 /= 6 .or. Nw2 /= 4 .or. D2 /= 6 .or. neq /= (Nw1+Nw2)/2) then
    print *, 'error: stitching_equations_hommed_hommed: improper input values:', Nw1, D1, Nw2, D2;
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
        M(k,0+l) = - EB1(ind(k),l);
        M(k,4+l) =   EB2(ind(k),l);
    end do

end do

!rhs vector:
if (flg_ant /= 0) then

    call current_density (jsurf);

    J(3) =   fpc*jsurf(2); !Bt
    J(4) = - fpc*jsurf(1); !Bz

end if

end subroutine

!------------------------------------------------------------------------------

subroutine calc_coeffs_vacuum_antenna_vacuum (m, kz, omega, r, S)

integer, intent(in) :: m;  !poloidal wave number
real(8), intent(in) :: kz; !n/R, n - toroidal wave number, R - big radius of the torus
complex(8), intent(in) :: omega;
real(8), intent(in) :: r;
complex(8), dimension(8), intent(out) :: S;

complex(8), dimension(2) :: jsurf;

complex(8) :: dBt, dBz;

real(8), parameter :: pi  = 3.141592653589793238462;
real(8), parameter :: c   = 29979245800.0;
real(8), parameter :: fpc = 4.0d0*pi/c;
complex(8), parameter :: I = (0.0d0, 1.0d0);

complex(8) :: agemo, gamma;

complex(8) :: Im, Km, dIm, dKm;

real(8) :: kt;

!Modified Bessel functions:
real(8) :: ZR, ZI, FNU;
integer :: KODE, NZ, IERR;
real(8), dimension(-1:1) :: CYR, CYI;

complex(8) :: gr;

agemo = omega;

gamma = sqrt(kz*kz - omega*agemo/c/c);

gr = gamma*r;

kt = m/r;

ZR = real(gr);
ZI = aimag(gr);

FNU = abs(m)-1; ! problem for m==0 !!!
KODE = 1;

call zbesi (ZR, ZI, FNU, KODE, 3, CYR, CYI, NZ, IERR);

!if (IERR /= 0) then
!    print *, 'besseli: warning: ierr =', IERR, 'NZ =', NZ, 'arg=', ZR + I*ZI;
!end if

Im = cmplx(CYR(0), CYI(0), 8);
dIm = 0.5*gamma*(cmplx(CYR(-1), CYI(-1), 8) + cmplx(CYR(1), CYI(1), 8));

call zbesk (ZR, ZI, FNU, KODE, 3, CYR, CYI, NZ, IERR);

!if (IERR /= 0) then
!    print *, 'besselk: warning: ierr =', IERR, 'NZ =', NZ, 'arg=', ZR + I*ZI;
!end if

Km = cmplx(CYR(0), CYI(0), 8);

dKm = - 0.5*gamma*(cmplx(CYR(-1), CYI(-1), 8) + cmplx(CYR(1), CYI(1), 8));

S = cmplx(0.0d0, 0.0d0, 8);

call current_density (jsurf);

dBt =   fpc*jsurf(2); !Bt
dBz = - fpc*jsurf(1); !Bz

S(2) = dKm/(dIm*Km-dKm*Im)*dBz/gamma/gamma;

S(7) = (gamma*gamma*dBt - kt*kz*dBz)*Im/(dIm*Km-dKm*Im)/gamma/gamma/(I*omega/c);

S(1) = Km/Im*S(7);

S(8) = (dBz +gamma*gamma*Im*S(2))/gamma/gamma/Km;

! //check:
! double *C = new double[2*Nc]; //complex coefficients (solution)
!
! double omega[2] = {real(wd->omov), imag(wd->omov)};   //moving frame
!
! int m = wd->m;
! double kz = (wd->n)/(sd->bs->rtor);
!
! calc_coeffs_vacuum_antenna_vacuum_ (&m, &kz, omega, &(sd->as->ra), C);

end subroutine

!------------------------------------------------------------------------------
