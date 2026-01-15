!> @file mag_wrapper.f90
!> @brief Magnetic field wrapper module for equilibrium field-line tracing
!>
!> Provides a unified interface for computing magnetic field properties
!> from the equilibrium field (via libneo's field_sub). Computes:
!> - Magnetic field magnitude
!> - Covariant and contravariant components of the unit magnetic field vector
!> - Metric tensor determinant (sqrtg)
!> - Field gradients
!>
!> Uses cylindrical coordinates (R, phi, Z) as input.

module mag_wrapper_m

    use field_sub

    implicit none
    private

    public :: mag_field_components

contains

    !> Compute magnetic field components and geometric quantities
    !>
    !> Given cylindrical coordinates (R, phi, Z), computes the magnetic field
    !> magnitude, covariant/contravariant unit vector components, and their
    !> derivatives. These are needed for field-line tracing.
    !>
    !> @param[in]  x       Coordinates: x(1)=R, x(2)=phi, x(3)=Z
    !> @param[out] bmod    Magnetic field magnitude |B|
    !> @param[out] sqrtg   Square root of metric tensor determinant (= R in cylindrical)
    !> @param[out] bder    Derivatives of log(|B|) w.r.t. coordinates
    !> @param[out] hcovar  Covariant components of unit B-field direction
    !> @param[out] hctrvr  Contravariant components of unit B-field direction
    !> @param[out] hcoder  Derivatives of covariant components: hcoder(i,j) = d(hcovar(j))/dx(i)
    !> @param[out] hctder  Derivatives of contravariant components: hctder(i,j) = d(hctrvr(j))/dx(i)
    subroutine mag_field_components(x, bmod, sqrtg, bder, hcovar, hctrvr, hcoder, hctder)
        double precision, intent(in) :: x(3)
        double precision, intent(out) :: bmod, sqrtg
        double precision, intent(out) :: bder(3), hcovar(3), hctrvr(3)
        double precision, intent(out) :: hcoder(3,3), hctder(3,3)

        double precision :: rbig, ri, fii, zi
        double precision :: br, bf, bz
        double precision :: brr, brf, brz, bfr, bff, bfz, bzr, bzf, bzz
        double precision :: hr, hf, hz

        ! Ensure R > 0 to avoid singularity at axis
        rbig = max(x(1), 1.d-12)
        ri = rbig
        fii = x(2)
        zi = x(3)

        ! Get B-field from equilibrium (via libneo)
        call field_eq(ri, fii, zi, br, bf, bz, &
                      brr, brf, brz, bfr, bff, bfz, bzr, bzf, bzz)

        ! Magnetic field magnitude
        bmod = dsqrt(br**2 + bf**2 + bz**2)

        ! Square root of metric determinant in cylindrical coords
        sqrtg = rbig

        ! Unit vector components
        hr = br / bmod
        hf = bf / bmod
        hz = bz / bmod

        ! Derivatives of log(|B|) = (1/|B|) * d|B|/dx
        bder(1) = (brr*hr + bfr*hf + bzr*hz) / bmod
        bder(2) = (brf*hr + bff*hf + bzf*hz) / bmod
        bder(3) = (brz*hr + bfz*hf + bzz*hz) / bmod

        ! Covariant components: h_i = g_ij * h^j
        ! In cylindrical: g_RR=1, g_phiphi=R^2, g_ZZ=1
        hcovar(1) = hr
        hcovar(2) = hf * rbig  ! = R * h^phi
        hcovar(3) = hz

        ! Contravariant components: h^i = B^i / |B|
        hctrvr(1) = hr
        hctrvr(2) = hf / rbig  ! = B^phi / |B| = (B_phi/R) / |B|
        hctrvr(3) = hz

        ! Derivatives of covariant components
        hcoder(1,1) = brr/bmod - hcovar(1)*bder(1)
        hcoder(2,1) = brf/bmod - hcovar(1)*bder(2)
        hcoder(3,1) = brz/bmod - hcovar(1)*bder(3)
        hcoder(1,2) = (rbig*bfr + bf)/bmod - hcovar(2)*bder(1)
        hcoder(2,2) = rbig*bff/bmod - hcovar(2)*bder(2)
        hcoder(3,2) = rbig*bfz/bmod - hcovar(2)*bder(3)
        hcoder(1,3) = bzr/bmod - hcovar(3)*bder(1)
        hcoder(2,3) = bzf/bmod - hcovar(3)*bder(2)
        hcoder(3,3) = bzz/bmod - hcovar(3)*bder(3)

        ! Derivatives of contravariant components
        hctder(1,1) = brr/bmod - hctrvr(1)*bder(1)
        hctder(2,1) = brf/bmod - hctrvr(1)*bder(2)
        hctder(3,1) = brz/bmod - hctrvr(1)*bder(3)
        hctder(1,2) = (bfr - bf/rbig)/(rbig*bmod) - hctrvr(2)*bder(1)
        hctder(2,2) = bff/(rbig*bmod) - hctrvr(2)*bder(2)
        hctder(3,2) = bfz/(rbig*bmod) - hctrvr(2)*bder(3)
        hctder(1,3) = bzr/bmod - hctrvr(3)*bder(1)
        hctder(2,3) = bzf/bmod - hctrvr(3)*bder(2)
        hctder(3,3) = bzz/bmod - hctrvr(3)*bder(3)

    end subroutine mag_field_components

end module mag_wrapper_m


!> Legacy interface for backward compatibility with existing fouriermodes code
subroutine mag(x, bmod, sqrtg, bder, hcovar, hctrvr, hcoder, hctder)
    use mag_wrapper_m, only: mag_field_components

    implicit none
    double precision, intent(in) :: x(3)
    double precision, intent(out) :: bmod, sqrtg
    double precision, intent(out) :: bder(3), hcovar(3), hctrvr(3)
    double precision, intent(out) :: hcoder(3,3), hctder(3,3)

    call mag_field_components(x, bmod, sqrtg, bder, hcovar, hctrvr, hcoder, hctder)

end subroutine mag
