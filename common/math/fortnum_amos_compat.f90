! AMOS complex-Bessel compatibility shim backed by the fortnum numerical core.
!
! KiLCA and KIM call the AMOS subroutines zbesj/zbesi/zbesk through the
! historical (ZR, ZI, FNU, KODE, N, CYR, CYI, NZ, IERR) argument list with no
! module use, relying on external linkage. These subroutines keep that exact
! interface so every caller compiles unchanged, but compute the complex Bessel
! sequences with fortnum_special_complex_bessel instead of the bundled AMOS
! sources. order0 = nint(FNU) (all callers pass integer-valued nonnegative
! orders); KODE 2 selects fortnum's scaled variant.
!
! IERR mapping (the callers branch only on IERR > 0, treating 3 as a
! completed-with-warning result for J): fortnum FORTNUM_OK -> 0, a domain
! error -> 1, any other nonzero status -> 4. NZ (AMOS underflow count) is set
! to 0; fortnum does not zero-flush trailing orders.

subroutine zbesj(zr, zi, fnu, kode, n, cyr, cyi, nz, ierr)
    use, intrinsic :: iso_fortran_env, only: dp => real64
    use fortnum_special_complex_bessel, only: bessel_j_complex
    use fortnum_status, only: fortnum_status_t, FORTNUM_OK, FORTNUM_DOMAIN_ERROR
    implicit none
    real(dp), intent(in) :: zr, zi, fnu
    integer, intent(in) :: kode, n
    real(dp), intent(out) :: cyr(n), cyi(n)
    integer, intent(out) :: nz, ierr

    complex(dp) :: z, value
    type(fortnum_status_t) :: status
    integer :: order0, k

    z = cmplx(zr, zi, kind=dp)
    order0 = nint(fnu)
    nz = 0
    ierr = 0

    ! KiLCA only ever requests the single unscaled order (N = 1, KODE = 1).
    do k = 1, n
        call bessel_j_complex(order0 + k - 1, z, value, status)
        cyr(k) = real(value, dp)
        cyi(k) = aimag(value)
        call map_ierr(status, ierr)
    end do
contains
    subroutine map_ierr(st, code)
        type(fortnum_status_t), intent(in) :: st
        integer, intent(inout) :: code
        if (st%code == FORTNUM_OK) return
        if (st%code == FORTNUM_DOMAIN_ERROR) then
            code = 1
        else
            code = 4
        end if
    end subroutine map_ierr
end subroutine zbesj

subroutine zbesi(zr, zi, fnu, kode, n, cyr, cyi, nz, ierr)
    use, intrinsic :: iso_fortran_env, only: dp => real64
    use fortnum_special_complex_bessel, only: bessel_i_complex_array
    use fortnum_status, only: fortnum_status_t, FORTNUM_OK, FORTNUM_DOMAIN_ERROR
    implicit none
    real(dp), intent(in) :: zr, zi, fnu
    integer, intent(in) :: kode, n
    real(dp), intent(out) :: cyr(n), cyi(n)
    integer, intent(out) :: nz, ierr

    complex(dp) :: z, seq(n)
    type(fortnum_status_t) :: status
    integer :: order0, k

    z = cmplx(zr, zi, kind=dp)
    order0 = nint(fnu)
    nz = 0
    ierr = 0

    call bessel_i_complex_array(order0, n, z, kode == 2, seq, status)
    do k = 1, n
        cyr(k) = real(seq(k), dp)
        cyi(k) = aimag(seq(k))
    end do

    if (status%code == FORTNUM_DOMAIN_ERROR) then
        ierr = 1
    else if (status%code /= FORTNUM_OK) then
        ierr = 4
    end if
end subroutine zbesi

subroutine zbesk(zr, zi, fnu, kode, n, cyr, cyi, nz, ierr)
    use, intrinsic :: iso_fortran_env, only: dp => real64
    use fortnum_special_complex_bessel, only: bessel_k_complex_array
    use fortnum_status, only: fortnum_status_t, FORTNUM_OK, FORTNUM_DOMAIN_ERROR
    implicit none
    real(dp), intent(in) :: zr, zi, fnu
    integer, intent(in) :: kode, n
    real(dp), intent(out) :: cyr(n), cyi(n)
    integer, intent(out) :: nz, ierr

    complex(dp) :: z, seq(n)
    type(fortnum_status_t) :: status
    integer :: order0, k

    z = cmplx(zr, zi, kind=dp)
    order0 = nint(fnu)
    nz = 0
    ierr = 0

    call bessel_k_complex_array(order0, n, z, kode == 2, seq, status)
    do k = 1, n
        cyr(k) = real(seq(k), dp)
        cyi(k) = aimag(seq(k))
    end do

    if (status%code == FORTNUM_DOMAIN_ERROR) then
        ierr = 1
    else if (status%code /= FORTNUM_OK) then
        ierr = 4
    end if
end subroutine zbesk
