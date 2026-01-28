!------------------------------------------------------------------------------

recursive complex(8) function besselj (nu, zarg, n) result (res)

!Evaluates n-th derivative of a bessel J_nu(zarg) function, nu - integer, zarg - complex!!!
!The subroutines from AMOS are used for.

implicit none;

integer, parameter :: dp = 8, dpc = 8

integer, intent(in) :: nu, n
complex(dpc), intent(in) :: zarg

integer :: k
real(dp) :: rr, ri
integer :: nz, ierr
real(dp), allocatable, dimension(:,:) :: bico

external :: binomial_coefficients

if(n == 0) then
    call zbesj (real(zarg), aimag(zarg), abs(dble(nu)), 1, 1, rr, ri, nz, ierr)

    !if (ierr /= 0) then
    !    print *, 'besselj: warning: ierr =', ierr, 'NZ =', nz, 'arg=', zarg
    !end if

    if (ierr > 0 .and. ierr /= 3) then
        print *, 'besselj: failed to compute besselj: answer is set to unit!'
        rr = 1.0d0
        ri = 0.0d0
    end if

    if (nu < 0) then
        res = cmplx(rr,ri,dp)*(-1.0d0)**nu
    else
        res = cmplx(rr,ri,dp)
    end if

    return
end if

allocate (bico(0:n,0:n));

!computes C^k_n = n!/k!/(n-k)! coefficients for n=0..N, k=0..n
call binomial_coefficients (%val(n), bico);

res = cmplx(0.0d0,0.0d0,dp)
do k=0,n
    res = res + (-1.0)**k*bico(n,k)*besselj(2*k-n+nu,zarg,0)
end do

res = res/(2.0**n)

deallocate (bico);

end function

!------------------------------------------------------------------------------

recursive complex(8) function besseli (nu, zarg, n) result (res)

!evaluates n-th derivative of a modified bessel I_nu(zarg) function, nu - integer, zarg - complex!!!

implicit none;

interface
    recursive complex(8) function besselj (nu, zarg, n)
        integer, intent(in) :: nu, n
        complex(8), intent(in) :: zarg
    end function
end interface

integer, parameter :: dpc = 8

integer, intent(in) :: nu, n
complex(dpc), intent(in) :: zarg

complex(dpc) :: I = (0.0d0, 1.0d0), J_nu

J_nu = besselj (nu, I*zarg, n)

res = I**(n-nu)*J_nu

end function

!------------------------------------------------------------------------------
