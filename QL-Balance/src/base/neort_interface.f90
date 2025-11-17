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
    public :: calculate_Omega_tE_splined
    public :: read_equil_file
    public :: calculate_s_tor
    public :: calculate_coarse_s_tor

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
    subroutine prepare_plasma_data_for_neort(plasma_data, r, s_tor)
        use baseparam_mod, only: EV_TO_ERG => ev
        use grid_mod, only: rc
        use plasma_parameters, only: params
        use spline, only: spline_coeff, spline_val

        real(dp), dimension(:, :), intent(out) :: plasma_data
        real(dp), dimension(:), intent(in) :: r
        real(dp), dimension(:), intent(in) :: s_tor

        real(dp), dimension(:, :), allocatable :: ni_of_r_coeffs, Ti_of_r_coeffs, Te_of_r_coeffs
        real(dp), dimension(:, :), allocatable :: ni_splined, Ti_splined, Te_splined
        integer :: i, p_size, c_size, s_size

        real(dp), parameter :: ERG_TO_EV = 1.0_dp / EV_TO_ERG

        p_size = size(params, 2)
        c_size = p_size - 1
        s_size = size(s_tor)

        ! Check dimensions
        if (size(plasma_data, 1) /= s_size .or. size(plasma_data, 2) /= 6) then
            error stop "prepare_neort_plasma_data: plasma_data dimension mismatch"
        end if

        allocate (ni_of_r_coeffs(c_size, 5))
        allocate (Ti_of_r_coeffs(c_size, 5))
        allocate (Te_of_r_coeffs(c_size, 5))
        allocate (ni_splined(s_size, 3))
        allocate (Ti_splined(s_size, 3))
        allocate (Te_splined(s_size, 3))

        ni_of_r_coeffs = spline_coeff(rc, params(1, :))
        ni_splined = spline_val(ni_of_r_coeffs, r)
        Ti_of_r_coeffs = spline_coeff(rc, params(4, :))
        Ti_splined = spline_val(Ti_of_r_coeffs, r)
        Te_of_r_coeffs = spline_coeff(rc, params(3, :))
        Te_splined = spline_val(Te_of_r_coeffs, r)

        ! Fill plasma data array for NEO-RT. C.f.:
        ! - KAMEL/QL-Balance/src/base/paramscan.f90:257
        ! - NEO-RT/examples/base/plasma.in
        ! - NEO-RT/doc/running.md
        do i = 1, s_size
            ! Column 1: Normalized toroidal flux s [1], in [0, 1]
            plasma_data(i, 1) = s_tor(i)

            ! Column 2: Density of species 1 [cm^-3]
            plasma_data(i, 2) = ni_splined(i, 1)

            ! Column 3: Density of species 2 [cm^-3]
            ! Set to 0 for single-species calculations
            plasma_data(i, 3) = 0.0_dp

            ! Column 4: Temperature of species 1 [eV]
            plasma_data(i, 4) = Ti_splined(i, 1) * ERG_TO_EV

            ! Column 5: Temperature of species 2 [eV]
            ! Dummy value for unused species 2
            plasma_data(i, 5) = 1.0_dp

            ! Column 6: Electron temperature [eV]
            plasma_data(i, 6) = Te_splined(i, 1) * ERG_TO_EV
        end do
    end subroutine prepare_plasma_data_for_neort

    !> @brief Prepare rotation profile data for NEO-RT from KAMEL arrays
    !> @param[out] profile_data 2D array (nflux, 3) for NEO-RT profile input
    subroutine prepare_profile_data_for_neort(profile_data, r, s_tor, Omega_tE)
        use baseparam_mod, only: am, p_mass, rtor
        use grid_mod, only: rc
        use plasma_parameters, only: params
        use spline, only: spline_coeff, spline_val

        real(dp), dimension(:, :), intent(out) :: profile_data
        real(dp), dimension(:), intent(in) :: r
        real(dp), dimension(:), intent(in) :: s_tor
        real(dp), dimension(:), intent(in) :: Omega_tE

        real(dp) :: T_i, m_i, vth, M_t
        real(dp), dimension(:, :), allocatable :: Ti_of_r_coeffs, Ti_splined
        integer :: i, p_size, c_size, s_size

        p_size = size(params, 2)
        c_size = p_size - 1
        s_size = size(s_tor)

        allocate (Ti_of_r_coeffs(c_size, 5))
        allocate (Ti_splined(s_size, 3))

        Ti_of_r_coeffs = spline_coeff(rc, params(4, :))
        Ti_splined = spline_val(Ti_of_r_coeffs, r)

        m_i = am * p_mass  ! ion mass

        ! Check dimensions
        if (size(profile_data, 1) /= s_size .or. size(profile_data, 2) /= 2) then
            error stop "prepare_neort_profile_data: profile_data dimension mismatch"
        end if

        ! Fill profile data array
        do i = 1, s_size
            ! Column 1: Normalized toroidal flux s (0 to 1)
            profile_data(i, 1) = s_tor(i)

            ! Calculate thermal velocity
            T_i = Ti_splined(i, 1)  ! ion temperature in erg
            vth = sqrt(2 * T_i / m_i)  ! thermal velocity

            ! Calculate ExB Mach number: M_t = Omega_tE * R_0 / v_th
            M_t = Omega_tE(i) * rtor / vth

            ! Column 2: ExB Mach number M_t
            profile_data(i, 2) = M_t
        end do
    end subroutine prepare_profile_data_for_neort

    subroutine calculate_Omega_tE_splined(Omega_tE, r)
        use baseparam_mod, only: btor, c
        use grid_mod, only: rb, rc
        use plasma_parameters, only: qsaf
        use wave_code_data, only: dPhi0
        use spline, only: spline_coeff, spline_val

        real(dp), dimension(:), intent(out) :: Omega_tE
        real(dp), dimension(:), intent(in) :: r

        real(dp) :: dPhi_dr, dpsi_pol_dr
        real(dp), dimension(:, :), allocatable :: dPhi0_coeffs, q_coeffs
        real(dp), dimension(:, :), allocatable :: dPhi0_splined, q_splined
        integer :: i, rb_size, rc_size, r_size

        rb_size = size(rb)
        rc_size = size(rc)
        r_size = size(r)

        allocate (dPhi0_coeffs(rb_size - 1, 5))
        allocate (q_coeffs(rc_size - 1, 5))

        ! TODO: here, rb, rc and r are intermixed, make sure to handle this correctly!
        dPhi0_coeffs = spline_coeff(rb, dPhi0)
        dPhi0_splined = spline_val(dPhi0_coeffs, r)
        q_coeffs = spline_coeff(rc, qsaf)
        q_splined = spline_val(q_coeffs, r)

        do i = 1, r_size
            ! Calculate toroidal electric precession frequency:
            ! Omega_tE = -c * (dPhi/dr) / (dpsi_pol/dr)
            ! Note: dPhi0 already contains dPhi/dr (electric potential gradient)
            ! and dpsi_pol/dr = B^theta = psi_tor' / q = r * B_tor / q
            ! Psi_tor = r²π * B_tor  =>  Psi_tor' = 2πr * B_tor
            ! our psi_tor == Psi_tor / (2π)
            dPhi_dr = dPhi0_splined(i, 1)
            dpsi_pol_dr = r(i) * btor / q_splined(i, 1)
            Omega_tE(i) = -c * dPhi_dr / dpsi_pol_dr
        end do
    end subroutine calculate_Omega_tE_splined

    !> @brief Read equilibrium data from equilibrium file
    !> @details Reads equilibrium quantities from the equil_r_q_psi.dat file
    !>          specified in the balance configuration. The file format is:
    !>          - Header lines starting with '#' (skipped automatically)
    !>          - Data columns:
    !>            r, q, psi_pol, phi_tor, dphi/dpsi, r_geom, V, R_beg, Z_beg, R_min, R_max
    !>          All arrays are automatically allocated to the correct size based on
    !>          the number of data lines in the file.
    !> @param[out] r_eff Effective radius array (column 1, optional, allocatable)
    !> @param[out] q_prof Safety factor array (column 2, optional, allocatable)
    !> @param[out] psi_pol Poloidal flux array (column 3, optional, allocatable)
    !> @param[out] psi_tor Toroidal flux array (column 4, optional, allocatable)
    !> @param[out] dphi_dpsi dphi/dpsi array (column 5, optional, allocatable)
    !> @param[out] r_geom Geometric radius array (column 6, optional, allocatable)
    !> @param[out] V Volume array (column 7, optional, allocatable)
    !> @param[out] R_beg R beginning array (column 8, optional, allocatable)
    !> @param[out] Z_beg Z beginning array (column 9, optional, allocatable)
    !> @param[out] R_min R minimum array (column 10, optional, allocatable)
    !> @param[out] R_max R maximum array (column 11, optional, allocatable)
    subroutine read_equil_file(r_eff, q_prof, psi_pol, psi_tor, dphi_dpsi, r_geom, V, R_beg, &
                               Z_beg, R_min, R_max)
        use control_mod, only: equil_path

        real(dp), dimension(:), allocatable, intent(out), optional :: r_eff
        real(dp), dimension(:), allocatable, intent(out), optional :: q_prof
        real(dp), dimension(:), allocatable, intent(out), optional :: psi_pol
        real(dp), dimension(:), allocatable, intent(out), optional :: psi_tor
        real(dp), dimension(:), allocatable, intent(out), optional :: dphi_dpsi
        real(dp), dimension(:), allocatable, intent(out), optional :: r_geom
        real(dp), dimension(:), allocatable, intent(out), optional :: V
        real(dp), dimension(:), allocatable, intent(out), optional :: R_beg
        real(dp), dimension(:), allocatable, intent(out), optional :: Z_beg
        real(dp), dimension(:), allocatable, intent(out), optional :: R_min
        real(dp), dimension(:), allocatable, intent(out), optional :: R_max

        integer :: iunit, ipoi, ios, npoints, nheader
        real(dp) :: r_val, q_val, psi_pol_val, psi_tor_val, dphi_dpsi_val
        real(dp) :: r_geom_val, V_val, R_beg_val, Z_beg_val, R_min_val, R_max_val
        character(len=1024) :: errmsg, line

        ! Check that at least one argument is present
        if (.not. (present(r_eff) .or. present(q_prof) .or. present(psi_pol) .or. &
                   present(psi_tor) .or. present(dphi_dpsi) .or. present(r_geom) .or. &
                   present(V) .or. present(R_beg) .or. present(Z_beg) .or. &
                   present(R_min) .or. present(R_max))) then
            error stop "read_equil_file: at least one optional argument must be present"
        end if

        ! Open equilibrium file
        open (newunit=iunit, file=trim(equil_path), status='old', action='read', iostat=ios)
        if (ios /= 0) then
            write (errmsg, '(A,A)') "read_equil_file: cannot open file ", trim(equil_path)
            error stop trim(errmsg)
        end if

        ! Skip all header lines starting with '#' (with optional leading whitespace)
        nheader = 0
        do
            read (iunit, '(A)', iostat=ios) line
            if (ios /= 0) then
                close (iunit)
                error stop "read_equil_file: error reading file (unexpected end during header)"
            end if
            ! Remove leading whitespace and check if line starts with '#'
            line = adjustl(line)
            ! Exit loop when we find a non-empty line that doesn't start with '#'
            if (len_trim(line) > 0 .and. line(1:1) /= '#') then
                ! This is the first data line, backspace to re-read it
                backspace (iunit)
                exit
            end if
            nheader = nheader + 1
        end do

        ! First pass: count the number of data lines
        npoints = 0
        do
            read (iunit, *, iostat=ios) r_val, q_val, psi_pol_val, psi_tor_val, dphi_dpsi_val, &
                r_geom_val, V_val, R_beg_val, Z_beg_val, R_min_val, R_max_val
            if (ios /= 0) exit
            npoints = npoints + 1
        end do

        if (npoints == 0) then
            close (iunit)
            error stop "read_equil_file: no data lines found in file"
        end if

        ! Allocate arrays based on the number of data lines
        if (present(r_eff)) allocate (r_eff(npoints))
        if (present(q_prof)) allocate (q_prof(npoints))
        if (present(psi_pol)) allocate (psi_pol(npoints))
        if (present(psi_tor)) allocate (psi_tor(npoints))
        if (present(dphi_dpsi)) allocate (dphi_dpsi(npoints))
        if (present(r_geom)) allocate (r_geom(npoints))
        if (present(V)) allocate (V(npoints))
        if (present(R_beg)) allocate (R_beg(npoints))
        if (present(Z_beg)) allocate (Z_beg(npoints))
        if (present(R_min)) allocate (R_min(npoints))
        if (present(R_max)) allocate (R_max(npoints))

        ! Rewind to the beginning of the file for second pass
        rewind (iunit)

        ! Skip the known number of header lines
        do ipoi = 1, nheader
            read (iunit, '(A)')
        end do

        ! Second pass: read and store data
        do ipoi = 1, npoints
            read (iunit, *, iostat=ios) r_val, q_val, psi_pol_val, psi_tor_val, dphi_dpsi_val, &
                r_geom_val, V_val, R_beg_val, Z_beg_val, R_min_val, R_max_val
            if (ios /= 0) then
                close (iunit)
                write (errmsg, '(A,I0,A)') "read_equil_file: error reading data line ", ipoi
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
        end do

        close (iunit)
    end subroutine read_equil_file

    subroutine calculate_s_tor(s_tor, psi_tor)
        real(dp), dimension(:), intent(out) :: s_tor
        real(dp), dimension(:), intent(in) :: psi_tor

        real(dp) :: psi_tor_min, psi_tor_max

        psi_tor_min = psi_tor(1)
        psi_tor_max = psi_tor(size(psi_tor))
        s_tor = (psi_tor - psi_tor_min) / (psi_tor_max - psi_tor_min)
    end subroutine calculate_s_tor

    subroutine calculate_coarse_s_tor(s_tor, s_min, s_max, npoints)
        use spline, only: spline_coeff, spline_val

        real(dp), dimension(:), intent(out) :: s_tor
        real(dp), intent(in) :: s_min, s_max
        integer, intent(in) :: npoints

        real(dp) :: ds
        integer :: i

        ds = (s_max - s_min) / npoints

        do i = 1, npoints
            s_tor(i) = s_min + (i - 1) * ds
        end do
    end subroutine calculate_coarse_s_tor

end module neort_interface
