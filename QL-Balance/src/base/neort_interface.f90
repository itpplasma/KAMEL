!> @brief Module for interfacing KAMEL data with NEO-RT
!> @details Provides data conversion and preparation routines to pass KAMEL
!>          data directly to NEO-RT instead of reading from files
module neort_interface
    use iso_fortran_env, only: dp => real64
    use neort_datatypes, only: magfie_data_t, transport_data_t

    implicit none

    private

    public :: prepare_plasma_data_for_neort
    public :: prepare_profile_data_for_neort
    public :: read_equil_file
    public :: calculate_s_tor

contains

    !> @brief Prepare plasma profile data for NEO-RT from KAMEL arrays
    !> @param[out] plasma_data 2D array (nflux, 6) for NEO-RT plasma input
    !> @details Passes KAMEL data to NEO-RT. Data format expected by NEO-RT:
    !>          Column 1: s_tor - normalized toroidal flux (dimensionless)
    !>          Column 2: ni1 - density of species 1 [cm^-3]
    !>          Column 3: ni2 - density of species 2 [cm^-3]
    !>          Column 4: Ti1 - temperature of species 1 [eV]
    !>          Column 5: Ti2 - temperature of species 2 [eV]
    !>          Column 6: Te - electron temperature [eV]
    subroutine prepare_plasma_data_for_neort(plasma_data, s_tor)
        use baseparam_mod, only: EV_TO_ERG => ev
        use plasma_parameters, only: params

        real(dp), dimension(:, :), intent(out) :: plasma_data
        real(dp), dimension(:), intent(in) :: s_tor

        integer :: i, s_size
        real(dp), parameter :: ERG_TO_EV = 1.0_dp / EV_TO_ERG

        s_size = size(s_tor)

        ! Check dimensions
        if (size(plasma_data, 1) /= s_size) then
            error stop "prepare_neort_plasma_data: plasma_data dimension mismatch"
        end if
        if (size(plasma_data, 2) /= 6) then
            error stop "prepare_neort_plasma_data: plasma_data must have six columns"
        end if

        ! Fill plasma data array for NEO-RT. C.f.:
        ! - KAMEL/QL-Balance/src/base/paramscan.f90:257 
        ! - NEO-RT/examples/base/plasma.in
        ! - NEO-RT/doc/running.md
        do i = 1, s_size
            ! Column 1: Normalized toroidal flux s [1], in [0, 1]
            plasma_data(i, 1) = s_tor(i)

            ! Column 2: Density of species 1 [cm^-3]
            plasma_data(i, 2) = params(1, i)

            ! Column 3: Density of species 2 [cm^-3]
            ! Set to 0 for single-species calculations
            plasma_data(i, 3) = 0.0_dp

            ! Column 4: Temperature of species 1 [eV]
            plasma_data(i, 4) = params(4, i) * ERG_TO_EV

            ! Column 5: Temperature of species 2 [eV]
            ! Dummy value for unused species 2
            plasma_data(i, 5) = 1.0_dp

            ! Column 6: Electron temperature [eV]
            plasma_data(i, 6) = params(3, i) * ERG_TO_EV
        end do

    end subroutine prepare_plasma_data_for_neort

    !> @brief Prepare rotation profile data for NEO-RT from KAMEL arrays
    !> @param[out] profile_data 2D array (nflux, 3) for NEO-RT profile input
    subroutine prepare_profile_data_for_neort(profile_data, s_tor)
        use baseparam_mod, only: am, p_mass, btor
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

    !> @brief Read toroidal flux and safety factor from equilibrium file
    !> @details Reads the toroidal flux (phi) and safety factor (q) from the equil_r_q_psi.dat file
    !>          specified in the balance configuration. The file format is:
    !>          - 3 header lines
    !>          - Data columns: r, q, psi_pol, phi_tor, dphi/dpsi, r_geom, V, R_beg, Z_beg, R_min, R_max
    !> @param[out] phi_tor Toroidal flux array (column 4 from equil file)
    !> @param[out] q_prof Safety factor array (column 2 from equil file)
    subroutine read_equil_file(phi_tor, q_prof)
        use control_mod, only: equil_path

        real(dp), dimension(:), intent(out) :: phi_tor
        real(dp), dimension(:), intent(out) :: q_prof

        integer :: iunit, ipoi, ios, phi_size, q_size
        real(dp) :: r_eff, q, psi_pol, phi, dphi_dpsi, r_geom, V, R_beg, Z_beg, R_min, R_max
        integer :: nlines
        character(len=1024) :: errmsg

        phi_size = size(phi_tor)
        q_size = size(q_prof)

        ! Check that both arrays have the same size
        if (phi_size /= q_size) then
            write (errmsg, '(A,I0,A,I0)') &
                "read_toroidal_flux: phi_tor size (", phi_size, ") != q_prof size (", q_size, ")"
            error stop trim(errmsg)
        end if

        ! Open equilibrium file
        open (newunit=iunit, file=trim(equil_path), status='old', action='read', iostat=ios)
        if (ios /= 0) then
            write (errmsg, '(A,A)') "read_toroidal_flux: cannot open file ", trim(equil_path)
            error stop trim(errmsg)
        end if

        ! Skip 3 header lines
        do ipoi = 1, 3
            read (iunit, *, iostat=ios)
            if (ios /= 0) then
                close (iunit)
                error stop "read_toroidal_flux: error reading header lines"
            end if
        end do

        ! Read data lines and extract toroidal flux (column 4) and safety factor (column 2)
        nlines = 0
        do ipoi = 1, phi_size
            read (iunit, *, iostat=ios) r_eff, q, psi_pol, phi, dphi_dpsi, r_geom, V, R_beg, &
                                        Z_beg, R_min, R_max
            if (ios /= 0) then
                close (iunit)
                write (errmsg, '(A,I0,A,I0,A)') "read_toroidal_flux: error reading line ", &
                                                ipoi + 3, " (expected ", phi_size, " data lines)"
                error stop trim(errmsg)
            end if
            phi_tor(ipoi) = phi
            q_prof(ipoi) = q
            nlines = nlines + 1
        end do

        close (iunit)

        ! Verify we read the expected number of lines
        if (nlines /= phi_size) then
            write (errmsg, '(A,I0,A,I0)') &
                "read_toroidal_flux: read ", nlines, " lines, expected ", phi_size
            error stop trim(errmsg)
        end if

    end subroutine read_equil_file

    subroutine calculate_s_tor(s, phi)
        real(dp), dimension(:), intent(out) :: s
        real(dp), dimension(:), intent(in) :: phi

        real(dp) :: phi_min, phi_max

        phi_min = phi(1)
        phi_max = phi(size(phi))
        s = (phi - phi_min) / (phi_max - phi_min)
    end subroutine calculate_s_tor

end module neort_interface
