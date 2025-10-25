!> @brief Module for interfacing KAMEL data with NEO-RT
!> @details Provides data conversion and preparation routines to pass KAMEL
!>          data directly to NEO-RT instead of reading from files
module neort_interface
    use iso_fortran_env, only: dp => real64
    use baseparam_mod, only: EV_TO_ERG => ev

    implicit none

    private

    public :: prepare_plasma_data_for_neort
    public :: prepare_profile_data_for_neort

    real(dp), parameter :: CM_TO_M = 1d-2
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

    !> @brief Prepare rotation profile data for NEO-RT from KAMEL arrays
    !> @param[out] profile_data 2D array (nflux, 3) for NEO-RT profile input
    subroutine prepare_profile_data_for_neort(profile_data)
        use grid_mod, only: npoic, rc
        use plasma_parameters, only: params
        use baseparam_mod, only: R0 => rtor, am, p_mass

        real(dp), intent(out) :: profile_data(:, :)

        integer :: ipoi
        real(dp), allocatable :: s(:)
        real(dp) :: omega_tor, v_tor, T_i, m_i, vth, M_t

        ! Check dimensions
        if (size(profile_data, 1) /= npoic) then
            error stop "prepare_neort_profile_data: profile_data dimension mismatch"
        end if
        if (size(profile_data, 2) /= 2) then
            error stop "prepare_neort_profile_data: profile_data must have two columns (s, M_t)"
        end if

        allocate(s(npoic))
        call calculate_s(s, rc)

        ! Fill profile data array
        do ipoi = 1, npoic
            ! Column 1: Normalized toroidal flux s (0 to 1)
            profile_data(ipoi, 1) = s(ipoi)

            ! Column 2: Toroidal Mach number M_t in [m/s]
            ! KAMEL params(2, :) is toroidal angular frequency [rad/s]
            ! KAMEL params(4, :) is ion temperatur [erg]
            omega_tor = params(2, ipoi)
            v_tor = omega_tor * R0  ! toroidal velocity [cm/s]
            T_i = params(4, ipoi)
            m_i = am * p_mass  ! ion mass [g]
            vth = sqrt(2d0 * T_i / m_i)  ! thermal velocity [cm/s]
            M_t = v_tor / vth * CM_TO_M
            profile_data(ipoi, 2) = M_t

            ! Calculate Mach number
            ! TODO: is this necessary?
            ! if (vth > 0.0d0) then
            !     profile_data(ipoi, 2) = v_toroidal / vth
            ! else
            !     profile_data(ipoi, 2) = 0.0d0
            ! end if
        end do

    end subroutine prepare_profile_data_for_neort

    subroutine calculate_s(s, r)
        real(dp), intent(out) :: s(:)
        real(dp), intent(in) :: r(:)

        real(dp) :: rmin, rmax

        rmin = minval(r)
        rmax = maxval(r)
        s = (r - rmin) / (rmax - rmin)
    end subroutine calculate_s

end module neort_interface
