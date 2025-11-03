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

contains

    !> @brief Prepare plasma profile data for NEO-RT from KAMEL arrays
    !> @param[out] plasma_data 2D array (nflux, 6) for NEO-RT plasma input
    subroutine prepare_plasma_data_for_neort(plasma_data, s_tor)
        use plasma_parameters, only: params

        real(dp), dimension(:, :), intent(out) :: plasma_data
        real(dp), dimension(:), intent(in) :: s_tor

        integer :: i, s_size

        s_size = size(s_tor)

        ! Check dimensions
        if (size(plasma_data, 1) /= s_size) then
            error stop "prepare_neort_plasma_data: plasma_data dimension mismatch"
        end if
        if (size(plasma_data, 2) /= 6) then
            error stop "prepare_neort_plasma_data: plasma_data must have six columns"
        end if

        ! Fill plasma data array
        ! c.f. e.g.: QL-Balance/src/base/paramscan.f90:257 and plasma.in file in NEO-RT's example base
        do i = 1, s_size
            ! Column 1: Normalized toroidal flux s (0 to 1)
            plasma_data(i, 1) = s_tor(i)

            ! Column 2: Density of species 1
            plasma_data(i, 2) = params(1, i)

            ! Column 3: Density of species 2
            plasma_data(i, 3) = 0d0

            ! Column 4: Temperature of species 1
            plasma_data(i, 4) = params(4, i)

            ! Column 5: Temperature of species 2
            plasma_data(i, 5) = 1d0

            ! Column 6: Electron temperature
            plasma_data(i, 6) = params(3, i)
        end do

    end subroutine prepare_plasma_data_for_neort

    !> @brief Prepare rotation profile data for NEO-RT from KAMEL arrays
    !> @param[out] profile_data 2D array (nflux, 3) for NEO-RT profile input
    subroutine prepare_profile_data_for_neort(profile_data, s_tor)
        use baseparam_mod, only: R0 => rtor, am, p_mass, btor
        use grid_mod, only: Ercov
        use plasma_parameters, only: params, qsafb

        real(dp), dimension(:, :), intent(out) :: profile_data
        real(dp), dimension(:), intent(in) :: s_tor

        integer :: i, s_size
        real(dp) :: T_i, m_i, vth, M_t, E_r, B_theta, B_mag

        s_size = size(s_tor)

        ! Check dimensions
        if (size(profile_data, 1) /= s_size) then
            error stop "prepare_neort_profile_data: profile_data dimension mismatch"
        end if
        if (size(profile_data, 2) /= 2) then
            error stop "prepare_neort_profile_data: profile_data must have two columns (s, M_t)"
        end if

        ! Fill profile data array
        ! TODO: check this, is is probably wrong
        do i = 1, s_size
            ! Column 1: Normalized toroidal flux s (0 to 1)
            profile_data(i, 1) = s_tor(i)

            ! Calculate thermal velocity and ion mass
            T_i = params(4, i)
            m_i = am * p_mass  ! ion mass
            vth = sqrt(T_i / m_i)  ! thermal velocity

            ! Calculate ExB Mach number: M_t = E_r * B_theta / (B^2 * v_th)
            E_r = Ercov(i)
            B_theta = btor / qsafb(i)
            B_mag = sqrt(btor**2 + B_theta**2)
            M_t = E_r * B_theta / (B_mag**2 * vth)

            ! Column 2: ExB Mach number M_t
            profile_data(i, 2) = M_t
        end do

    end subroutine prepare_profile_data_for_neort

end module neort_interface
