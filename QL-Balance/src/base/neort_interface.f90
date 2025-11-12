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
        use baseparam_mod, only: am, p_mass, btor, rtor, c
        use grid_mod, only: rc
        use plasma_parameters, only: params, qsaf
        use wave_code_data, only: dPhi0

        real(dp), dimension(:, :), intent(out) :: profile_data
        real(dp), dimension(:), intent(in) :: s_tor

        integer :: i, s_size
        real(dp) :: T_i, m_i, vth, M_t, dPhi_dr, dpsi_pol_dr, Omega_tE

        s_size = size(s_tor)

        ! Check dimensions
        if (size(profile_data, 1) /= s_size) then
            error stop "prepare_neort_profile_data: profile_data dimension mismatch"
        end if
        if (size(profile_data, 2) /= 2) then
            error stop "prepare_neort_profile_data: profile_data must have two columns (s, M_t)"
        end if

        ! Fill profile data array
        do i = 1, s_size
            ! Column 1: Normalized toroidal flux s (0 to 1)
            profile_data(i, 1) = s_tor(i)

            ! Calculate thermal velocity and ion mass
            T_i = params(4, i)
            m_i = am * p_mass  ! ion mass
            vth = sqrt(2 * T_i / m_i)  ! thermal velocity

            ! Calculate toroidal electric precession frequency:
            ! Omega_tE = -c * (dPhi/dr) / (dpsi_pol/dr)
            ! Note: dPhi0 already contains dPhi/dr (electric potential gradient)
            ! and dpsi_pol/dr = B^theta = psi_tor' / q = r * B_tor / q
            ! Psi_tor = r²π * B_tor  =>  Psi_tor' = 2πr * B_tor
            ! our psi_tor == Psi_tor / (2π)
            dPhi_dr = dPhi0(i)
            dpsi_pol_dr = rc(i) * btor / qsaf(i)
            Omega_tE = -c * dPhi_dr / dpsi_pol_dr

            ! Calculate ExB Mach number: M_t = Omega_tE * R_0 / v_th
            M_t = Omega_tE * rtor / vth

            ! Column 2: ExB Mach number M_t
            profile_data(i, 2) = M_t
        end do

    end subroutine prepare_profile_data_for_neort

    !> @brief Read equilibrium data from equilibrium file
    !> @details Reads equilibrium quantities from the equil_r_q_psi.dat file
    !>          specified in the balance configuration. The file format is:
    !>          - 3 header lines
    !>          - Data columns:
    !>            r, q, psi_pol, phi_tor, dphi/dpsi, r_geom, V, R_beg, Z_beg, R_min, R_max
    !> @param[out] r_eff Effective radius array (column 1, optional)
    !> @param[out] q_prof Safety factor array (column 2, optional)
    !> @param[out] psi_pol Poloidal flux array (column 3, optional)
    !> @param[out] phi_tor Toroidal flux array (column 4, optional)
    !> @param[out] dphi_dpsi dphi/dpsi array (column 5, optional)
    !> @param[out] r_geom Geometric radius array (column 6, optional)
    !> @param[out] V Volume array (column 7, optional)
    !> @param[out] R_beg R beginning array (column 8, optional)
    !> @param[out] Z_beg Z beginning array (column 9, optional)
    !> @param[out] R_min R minimum array (column 10, optional)
    !> @param[out] R_max R maximum array (column 11, optional)
    subroutine read_equil_file(r_eff, q_prof, psi_pol, psi_tor, dphi_dpsi, r_geom, V, R_beg, &
                               Z_beg, R_min, R_max)
        use control_mod, only: equil_path

        real(dp), dimension(:), intent(out), optional :: r_eff
        real(dp), dimension(:), intent(out), optional :: q_prof
        real(dp), dimension(:), intent(out), optional :: psi_pol
        real(dp), dimension(:), intent(out), optional :: psi_tor
        real(dp), dimension(:), intent(out), optional :: dphi_dpsi
        real(dp), dimension(:), intent(out), optional :: r_geom
        real(dp), dimension(:), intent(out), optional :: V
        real(dp), dimension(:), intent(out), optional :: R_beg
        real(dp), dimension(:), intent(out), optional :: Z_beg
        real(dp), dimension(:), intent(out), optional :: R_min
        real(dp), dimension(:), intent(out), optional :: R_max

        integer :: iunit, ipoi, ios, npoints
        real(dp) :: r_val, q_val, psi_pol_val, psi_tor_val, dphi_dpsi_val
        real(dp) :: r_geom_val, V_val, R_beg_val, Z_beg_val, R_min_val, R_max_val
        integer :: nlines
        character(len=1024) :: errmsg

        ! Determine array size from first present optional argument
        if (present(r_eff)) then
            npoints = size(r_eff)
        else if (present(q_prof)) then
            npoints = size(q_prof)
        else if (present(psi_pol)) then
            npoints = size(psi_pol)
        else if (present(psi_tor)) then
            npoints = size(psi_tor)
        else if (present(dphi_dpsi)) then
            npoints = size(dphi_dpsi)
        else if (present(r_geom)) then
            npoints = size(r_geom)
        else if (present(V)) then
            npoints = size(V)
        else if (present(R_beg)) then
            npoints = size(R_beg)
        else if (present(Z_beg)) then
            npoints = size(Z_beg)
        else if (present(R_min)) then
            npoints = size(R_min)
        else if (present(R_max)) then
            npoints = size(R_max)
        else
            error stop "read_equil_file: at least one optional argument must be present"
        end if

        ! Check that all present arrays have the same size
        if (present(r_eff) .and. size(r_eff) /= npoints) then
            error stop "read_equil_file: r_eff size mismatch"
        end if
        if (present(q_prof) .and. size(q_prof) /= npoints) then
            error stop "read_equil_file: q_prof size mismatch"
        end if
        if (present(psi_pol) .and. size(psi_pol) /= npoints) then
            error stop "read_equil_file: psi_pol size mismatch"
        end if
        if (present(psi_tor) .and. size(psi_tor) /= npoints) then
            error stop "read_equil_file: psi_tor size mismatch"
        end if
        if (present(dphi_dpsi) .and. size(dphi_dpsi) /= npoints) then
            error stop "read_equil_file: dphi_dpsi size mismatch"
        end if
        if (present(r_geom) .and. size(r_geom) /= npoints) then
            error stop "read_equil_file: r_geom size mismatch"
        end if
        if (present(V) .and. size(V) /= npoints) then
            error stop "read_equil_file: V size mismatch"
        end if
        if (present(R_beg) .and. size(R_beg) /= npoints) then
            error stop "read_equil_file: R_beg size mismatch"
        end if
        if (present(Z_beg) .and. size(Z_beg) /= npoints) then
            error stop "read_equil_file: Z_beg size mismatch"
        end if
        if (present(R_min) .and. size(R_min) /= npoints) then
            error stop "read_equil_file: R_min size mismatch"
        end if
        if (present(R_max) .and. size(R_max) /= npoints) then
            error stop "read_equil_file: R_max size mismatch"
        end if

        ! Open equilibrium file
        open (newunit=iunit, file=trim(equil_path), status='old', action='read', iostat=ios)
        if (ios /= 0) then
            write (errmsg, '(A,A)') "read_equil_file: cannot open file ", trim(equil_path)
            error stop trim(errmsg)
        end if

        ! Skip 3 header lines
        do ipoi = 1, 3
            read (iunit, *, iostat=ios)
            if (ios /= 0) then
                close (iunit)
                error stop "read_equil_file: error reading header lines"
            end if
        end do

        ! Read data lines and store requested quantities
        nlines = 0
        do ipoi = 1, npoints
            read (iunit, *, iostat=ios) r_val, q_val, psi_pol_val, psi_tor_val, dphi_dpsi_val, &
                r_geom_val, V_val, R_beg_val, Z_beg_val, &
                R_min_val, R_max_val
            if (ios /= 0) then
                close (iunit)
                write (errmsg, '(A,I0,A,I0,A)') "read_equil_file: error reading line ", ipoi + 3, &
                    " (expected ", npoints, " data lines)"
                error stop trim(errmsg)
            end if

            ! Store values in arrays if present
            if (present(r_eff)) r_eff(ipoi) = r_val
            if (present(q_prof)) q_prof(ipoi) = q_val
            if (present(psi_pol)) psi_pol(ipoi) = psi_pol_val
            if (present(psi_tor)) psi_tor(ipoi) = psi_tor_val
            if (present(dphi_dpsi)) dphi_dpsi(ipoi) = dphi_dpsi_val
            if (present(r_geom)) r_geom(ipoi) = r_geom_val
            if (present(V)) V(ipoi) = V_val
            if (present(R_beg)) R_beg(ipoi) = R_beg_val
            if (present(Z_beg)) Z_beg(ipoi) = Z_beg_val
            if (present(R_min)) R_min(ipoi) = R_min_val
            if (present(R_max)) R_max(ipoi) = R_max_val

            nlines = nlines + 1
        end do

        close (iunit)

        ! Verify we read the expected number of lines
        if (nlines /= npoints) then
            write (errmsg, '(A,I0,A,I0)') "read_equil_file: read ", nlines, &
                " lines, expected ", npoints
            error stop trim(errmsg)
        end if
    end subroutine read_equil_file

    subroutine calculate_s_tor(s_tor, psi_tor)
        real(dp), dimension(:), intent(out) :: s_tor
        real(dp), dimension(:), intent(in) :: psi_tor

        real(dp) :: psi_tor_min, psi_tor_max

        psi_tor_min = psi_tor(1)
        psi_tor_max = psi_tor(size(psi_tor))
        s_tor = (psi_tor - psi_tor_min) / (psi_tor_max - psi_tor_min)
    end subroutine calculate_s_tor

end module neort_interface
