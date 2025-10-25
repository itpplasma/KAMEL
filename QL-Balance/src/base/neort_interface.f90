!> @brief Module for interfacing KAMEL data with NEO-RT
!> @details Provides data conversion and preparation routines to pass KAMEL
!>          data directly to NEO-RT instead of reading from files
module neort_interface
    use iso_fortran_env, only: dp => real64
    use baseparam_mod, only: EV_TO_ERG => ev

    implicit none

    private

    public :: prepare_plasma_data_for_neort

    real(dp), parameter :: CM3_TO_M3 = 1d-6
    real(dp), parameter :: ERG_TO_EV = 1d0 / EV_TO_ERG

contains

    !> @brief Prepare plasma profile data for NEO-RT from KAMEL arrays
    !> @param[out] plasma_data 2D array (nflux, 6) for NEO-RT plasma input
    !> @param[out] am1 Mass of species 1 [atomic mass units]
    !> @param[out] am2 Mass of species 2 [atomic mass units]
    !> @param[out] Z1 Charge of species 1
    !> @param[out] Z2 Charge of species 2
    subroutine prepare_plasma_data_for_neort(plasma_data, am1, am2, Z1, Z2)
        use baseparam_mod, only: am, Z_i
        use grid_mod, only: npoic, rc, Sc
        use plasma_parameters, only: params

        real(dp), intent(out) :: plasma_data(:, :)
        real(dp), intent(out) :: am1, am2, Z1, Z2

        integer :: ipoi
        real(dp), allocatable :: s(:)

        ! Check dimensions
        if (size(plasma_data, 1) /= npoic) then
            error stop "prepare_neort_plasma_data: plasma_data dimension mismatch"
        end if
        if (size(plasma_data, 2) /= 6) then
            error stop "prepare_neort_plasma_data: plasma_data must have six columns"
        end if

        ! Set species parameters from KAMEL
        am1 = am  ! Main ion mass [u]
        am2 = am  ! Secondary species (same as main for now)
        Z1 = Z_i  ! Main ion charge
        Z2 = Z_i  ! Secondary species charge

        allocate(s(npoic))
        call calculate_s(s, rc)

        ! Fill plasma data array
        ! c.f. e.g.: QL-Balance/src/base/paramscan.f90:257 and plasma.in file in NEO-RT's example base
        do ipoi = 1, npoic
            ! Column 1: Normalized toroidal flux s (0 to 1)
            plasma_data(ipoi, 1) = s(ipoi)

            ! Column 2: Density of species 1 [1/m³]
            ! KAMEL params(1, :) is in [1/cm³], convert to [1/m³]
            plasma_data(ipoi, 2) = params(1, ipoi) / CM3_TO_M3

            ! Column 3: Density of species 2 [1/m³]
            ! Set to zero for single-species case
            plasma_data(ipoi, 3) = 0d0

            ! Column 4: Temperature of species 1 [eV]
            ! KAMEL params(4, :) is ion temperature in [erg], convert to [eV]
            plasma_data(ipoi, 4) = params(4, ipoi) * ERG_TO_EV

            ! Column 5: Temperature of species 2 [eV]
            ! Set to 1 for single-species case
            plasma_data(ipoi, 5) = 1d0

            ! Column 6: Electron temperature [eV]
            ! KAMEL params(3, :) is electron temperature in [erg], convert to [eV]
            plasma_data(ipoi, 6) = params(3, ipoi) * ERG_TO_EV
        end do

    end subroutine prepare_plasma_data_for_neort

    subroutine calculate_s(s, r)
        real(dp), intent(out) :: s(:)
        real(dp), intent(in) :: r(:)

        real(dp) :: rmin, rmax

        rmin = minval(r)
        rmax = maxval(r)
        s = (r - rmin) / (rmax - rmin)
    end subroutine calculate_s

end module neort_interface
