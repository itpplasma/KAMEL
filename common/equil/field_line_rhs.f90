!> @file field_line_rhs.f90
!> @brief Field-line tracing ODE right-hand side for equilibrium calculations
!>
!> Provides the ODE system for tracing magnetic field lines in toroidal geometry.
!> The state vector contains:
!>   y(1) = R (major radius)
!>   y(2) = Z (vertical position)
!>   y(3) = enclosed cross-sectional area (integrated)
!>   y(4) = toroidal flux or R-accumulator (mode-dependent)
!>   y(5) = volume element or Z-accumulator (mode-dependent)
!>
!> Two modes are supported via the mode flag:
!>   Mode 1 (axis finding): y(4:5) accumulate R and Z for averaging
!>   Mode 2 (profile computation): y(4:5) accumulate flux and volume integrals

module field_line_rhs_m

    implicit none
    private

    public :: set_field_line_mode, get_dz_dphi, field_line_rhs

    !> Module state for field-line integration
    integer, save :: field_line_mode = 1  ! 1 = axis finding, 2 = profile computation
    double precision, save :: dz_dphi_last = 0.d0  ! Last computed dZ/dphi (for Newton iteration)

contains

    !> Set the field-line integration mode
    !> @param[in] mode  1 = axis finding, 2 = profile computation
    subroutine set_field_line_mode(mode)
        integer, intent(in) :: mode
        field_line_mode = mode
    end subroutine set_field_line_mode

    !> Get the last computed dZ/dphi value (used for Newton iteration in axis finding)
    !> @return dz_dphi  The derivative dZ/dphi from the last RHS evaluation
    function get_dz_dphi() result(dz_dphi)
        double precision :: dz_dphi
        dz_dphi = dz_dphi_last
    end function get_dz_dphi

    !> Compute the RHS of the field-line ODE system
    !>
    !> @param[in]  ndim  Dimension of state vector (should be 5)
    !> @param[in]  phi   Toroidal angle (independent variable)
    !> @param[in]  y     State vector
    !> @param[out] dery  Derivatives dy/dphi
    subroutine field_line_rhs(ndim, phi, y, dery)
        implicit none

        integer, intent(in) :: ndim
        double precision, intent(in) :: phi
        double precision, intent(in) :: y(ndim)
        double precision, intent(out) :: dery(ndim)

        double precision :: x(3), bmod, sqrtg
        double precision :: bder(3), hcovar(3), hctrvr(3)
        double precision :: hcoder(3,3), hctder(3,3)

        ! Set up coordinates for mag call: x = (R, phi, Z)
        x(1) = y(1)
        x(2) = phi
        x(3) = y(2)

        ! Get magnetic field components
        call mag(x, bmod, sqrtg, bder, hcovar, hctrvr, hcoder, hctder)

        ! Field-line equations: dR/dphi and dZ/dphi
        ! h^R / h^phi and h^Z / h^phi
        dery(1) = hctrvr(1) / hctrvr(2)  ! dR/dphi
        dery(2) = hctrvr(3) / hctrvr(2)  ! dZ/dphi

        ! Enclosed area differential: dA = R * dZ (cross-sectional area element)
        dery(3) = y(1) * hctrvr(3) / hctrvr(2)

        ! Mode-dependent integrals
        if (field_line_mode == 1) then
            ! Mode 1: Axis finding - accumulate R and Z for averaging
            dery(4) = y(1)   ! Accumulate R (for R-average)
            dery(5) = y(2)   ! Accumulate Z (for Z-average)
        else
            ! Mode 2: Profile computation - accumulate flux and volume
            dery(4) = bmod * y(1) * y(2) * hctrvr(1)  ! Toroidal flux-related
            dery(5) = y(1)**2 * hctrvr(3) / hctrvr(2)  ! Volume element
        end if

        ! Store dZ/dphi for Newton iteration (used in axis finding)
        dz_dphi_last = dery(2)

    end subroutine field_line_rhs

end module field_line_rhs_m


!> Legacy interface for backward compatibility with existing fouriermodes code
!> This module provides the rhs1_mod variables and rhs1 subroutine
module rhs1_mod
    use field_line_rhs_m, only: get_dz_dphi, set_field_line_mode

    implicit none

    !> dZ/dphi from last RHS evaluation (legacy interface)
    double precision :: dz_dphi

    !> Mode switch: 1 = axis finding, 2 = profile computation
    integer :: isw_rhs1 = 1

contains

    !> Update the legacy dz_dphi variable from the module state
    subroutine sync_dz_dphi()
        dz_dphi = get_dz_dphi()
    end subroutine sync_dz_dphi

end module rhs1_mod


!> Legacy subroutine interface for backward compatibility
subroutine rhs1(ndim, phi, y, dery)
    use field_line_rhs_m, only: field_line_rhs, set_field_line_mode
    use rhs1_mod, only: isw_rhs1, dz_dphi, sync_dz_dphi

    implicit none
    integer, intent(in) :: ndim
    double precision, intent(in) :: phi
    double precision, intent(in) :: y(ndim)
    double precision, intent(out) :: dery(ndim)

    ! Sync mode from legacy interface
    call set_field_line_mode(isw_rhs1)

    ! Call the module implementation
    call field_line_rhs(ndim, phi, y, dery)

    ! Sync dz_dphi back to legacy module variable
    call sync_dz_dphi()

end subroutine rhs1
