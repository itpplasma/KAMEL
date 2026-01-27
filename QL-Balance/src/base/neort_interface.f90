!> @brief Module for interfacing KAMEL data with NEO-RT
!> @details Provides data conversion and preparation routines to pass KAMEL
!>          data directly to NEO-RT instead of reading from files
module neort_interface
    use iso_fortran_env, only: dp => real64
    use neort_lib, only: config_t, transport_data_t

    implicit none

    private

    type, public :: meta_config_neort_t
        character(len=1024) :: boozer_file
        character(len=1024) :: boozer_pert_file
        integer :: amount_of_s
        type(config_t) :: config
    end type meta_config_neort_t

    public :: apply_ntv_transport
    public :: calculate_Omega_tE
    public :: calculate_coarse_s_tor
    public :: calculate_s_tor
    public :: prepare_plasma_data_for_neort
    public :: prepare_profile_data_for_neort
    public :: read_equil_file
    public :: read_neort_meta_config

contains

    subroutine read_neort_meta_config(config_file, meta_config)
        use wave_code_data, only: m_vals, n_vals

        character(len=*), intent(in) :: config_file
        type(meta_config_neort_t), intent(out) :: meta_config

        character(len=*), parameter :: unset = "__!UNSET!__"

        ! mode numbers
        integer :: poloidal_mode
        integer :: toroidal_mode

        ! I/O error handling
        integer :: iostat
        character(len=512) :: iomsg

        ! namelist content
        ! This is a modified version of the full input format used by NEO-RT
        ! as some values are deliberately unused or fixed plus the following:
        ! explicit paths for the Boozer input files, where the perturbation file is optional
        ! amount of equidistant s to compute on
        ! new meta config values:
        character(len=1024) :: boozer_file
        character(len=1024) :: boozer_pert_file  ! optional
        integer :: amount_of_s
        ! original NEO-RT config values
        real(dp) :: epsmn
        logical :: magdrift
        logical :: nopassing
        logical :: noshear
        logical :: nonlin
        real(dp) :: bfac
        real(dp) :: efac
        integer :: inp_swi
        integer :: vsteps
        integer :: log_level

        namelist /NEORT/ boozer_file, boozer_pert_file, amount_of_s, epsmn, magdrift, &
            nopassing, noshear, nonlin, bfac, efac, inp_swi, vsteps, log_level

        ! some default values
        boozer_pert_file = unset
        amount_of_s = 100
        bfac = 1.0
        efac = 1.0
        inp_swi = 9  ! AUG
        vsteps = 512
        log_level = -1  ! silent

        ! read namelist
        open (22, file=config_file)
        read (22, nml=NEORT, iostat=iostat, iomsg=iomsg)

        if (iostat /= 0) then
            write (*, *) "type_of_run == 'TimeEvolutionNTV' specified, but "// &
                "no NEORT config namelist was found!"
            write (*, *) "Error message: "//trim(iomsg)
            error stop
        end if

        ! mode numbers (from KiLCA input)
        poloidal_mode = abs(m_vals(1))
        toroidal_mode = n_vals(1)

        ! set values from namelist in meta config struct
        meta_config%boozer_file = trim(adjustl(boozer_file))
        meta_config%boozer_pert_file = trim(adjustl(boozer_pert_file))
        meta_config%amount_of_s = amount_of_s
        meta_config%config%epsmn = epsmn
        meta_config%config%m0 = poloidal_mode
        meta_config%config%mph = toroidal_mode
        meta_config%config%comptorque = .true.  ! must be true
        meta_config%config%magdrift = magdrift
        meta_config%config%nopassing = nopassing
        meta_config%config%noshear = noshear
        meta_config%config%pertfile = boozer_pert_file /= unset  ! only if boozer_pert_file was set
        meta_config%config%nonlin = nonlin
        meta_config%config%bfac = bfac
        meta_config%config%efac = efac
        meta_config%config%inp_swi = inp_swi
        meta_config%config%vsteps = vsteps
        meta_config%config%log_level = log_level
    end subroutine read_neort_meta_config

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
        integer :: i, rc_size, s_size

        real(dp), parameter :: ERG_TO_EV = 1.0_dp / EV_TO_ERG

        rc_size = size(rc)
        s_size = size(s_tor)

        ! Check dimensions
        if (size(plasma_data, 1) /= s_size .or. size(plasma_data, 2) /= 6) then
            error stop "prepare_neort_plasma_data: plasma_data dimension mismatch"
        end if

        allocate (ni_of_r_coeffs(rc_size - 1, 5))
        allocate (Ti_of_r_coeffs(rc_size - 1, 5))
        allocate (Te_of_r_coeffs(rc_size - 1, 5))
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
        integer :: i, rc_size, s_size

        rc_size = size(rc)
        s_size = size(s_tor)

        allocate (Ti_of_r_coeffs(rc_size - 1, 5))
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

    subroutine calculate_Omega_tE(Omega_tE, r)
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
    end subroutine calculate_Omega_tE

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
        real(dp), dimension(:), intent(out) :: s_tor
        real(dp), intent(in) :: s_min, s_max
        integer, intent(in) :: npoints

        real(dp) :: ds
        integer :: i

        ds = (s_max - s_min) / (npoints - 1)

        do i = 1, npoints
            s_tor(i) = s_min + (i - 1) * ds
        end do
    end subroutine calculate_coarse_s_tor

    subroutine apply_ntv_transport(r, transport_data)
        use grid_mod, only: rb, torque_ntv
        use PolyLagrangeInterpolation, only: binsrc
        use spline, only: spline_coeff, spline_val

        real(dp), dimension(:), intent(in) :: r
        type(transport_data_t), dimension(:), intent(in) :: transport_data

        real(dp), dimension(:), allocatable :: total_torque  ! in cgs over s
        real(dp), dimension(:), allocatable :: dVds
        real(dp), dimension(:, :), allocatable :: torque_of_r_coeffs
        real(dp), dimension(:, :), allocatable :: torque_splined
        integer :: ntorque, nrb, i, i_separatrix

        ! extract array sizes
        ntorque = size(transport_data)
        nrb = size(rb)

        allocate (total_torque(ntorque))
        allocate (dVds(ntorque))

        do i = 1, ntorque
            total_torque(i) = transport_data(i)%torque%Tco + &
                              transport_data(i)%torque%Tctr + &
                              transport_data(i)%torque%Tt
            dVds(i) = transport_data(i)%torque%dVds
        end do

        allocate (torque_of_r_coeffs(ntorque - 1, 5))
        allocate (torque_splined(nrb, 3))

        ! translate to cgs via dividing by dVds
        torque_of_r_coeffs = spline_coeff(r, total_torque / dVds)

        ! find index/index+1 of last element from r in rb
        call binsrc(rb, 1, nrb, r(ntorque), i_separatrix)

        ! make sure to not hit the separatrix
        i_separatrix = i_separatrix - 1

        ! do not extrapolate
        torque_splined = spline_val(torque_of_r_coeffs, rb(1:i_separatrix))

        ! torque for r > separatrix should be zero
        torque_ntv(1:i_separatrix) = torque_splined(:, 1)
        torque_ntv(i_separatrix+1:nrb) = 0.0_dp
    end subroutine apply_ntv_transport

!================== Writing routines for debugging ==================

    subroutine print_meta_config(meta_config)
        type(meta_config_neort_t), intent(in) :: meta_config

        print *, "NEO-RT meta config:"
        print *, "boozer_file = ", trim(adjustl(meta_config%boozer_file))
        print *, "boozer_pertPfile = ", trim(adjustl(meta_config%boozer_pert_file))
        print *, "amount_of_s = ", meta_config%amount_of_s
        print *, "s = ", meta_config%config%s
        print *, "M_t = ", meta_config%config%M_t
        print *, "qs = ", meta_config%config%qs
        print *, "ms = ", meta_config%config%ms
        print *, "vth = ", meta_config%config%vth
        print *, "epsmn = ", meta_config%config%epsmn
        print *, "m0 = ", meta_config%config%m0
        print *, "mph = ", meta_config%config%mph
        print *, "comptorque = ", meta_config%config%comptorque
        print *, "magdrift = ", meta_config%config%magdrift
        print *, "nopassing = ", meta_config%config%nopassing
        print *, "noshear = ", meta_config%config%noshear
        print *, "pertfile = ", meta_config%config%pertfile
        print *, "nonlin = ", meta_config%config%nonlin
        print *, "bfac = ", meta_config%config%bfac
        print *, "efac = ", meta_config%config%efac
        print *, "inp_swi = ", meta_config%config%inp_swi
        print *, "vsteps = ", meta_config%config%vsteps
        print *, "log_level = ", meta_config%config%log_level
    end subroutine print_meta_config

    !> @brief Write plasma_data array to a file for debugging
    !> @param[in] plasma_data 2D array with plasma data
    !> @param[in] filename Output filename
    !> @param[in] r Radial coordinate array (optional)
    subroutine write_plasma_data_to_file(plasma_data, filename, r)
        real(dp), dimension(:, :), intent(in) :: plasma_data
        character(len=*), intent(in) :: filename
        real(dp), dimension(:), intent(in), optional :: r

        integer :: iunit, i, ios
        character(len=1024) :: errmsg

        open (newunit=iunit, file=filename, status='replace', action='write', iostat=ios)
        if (ios /= 0) then
            write (errmsg, '(A,A)') "write_plasma_data_to_file: cannot open file ", filename
            error stop trim(errmsg)
        end if

        ! Write header
        if (present(r)) then
            write (iunit, '(A)') &
               "# Plasma data: s_tor, ni1 [cm^-3], ni2 [cm^-3], Ti1 [eV], Ti2 [eV], Te [eV], r [cm]"
            write (iunit, '(A,I0,A,I0)') "# npoints = ", size(plasma_data, 1), &
                ", ncols = ", size(plasma_data, 2) + 1
        else
            write (iunit, '(A)') &
                "# Plasma data: s_tor, ni1 [cm^-3], ni2 [cm^-3], Ti1 [eV], Ti2 [eV], Te [eV]"
            write (iunit, '(A,I0,A,I0)') "# npoints = ", size(plasma_data, 1), &
                ", ncols = ", size(plasma_data, 2)
        end if

        ! Write data
        do i = 1, size(plasma_data, 1)
            if (present(r)) then
                write (iunit, '(*(ES23.15E3,2X))') plasma_data(i, :), r(i)
            else
                write (iunit, '(*(ES23.15E3,2X))') plasma_data(i, :)
            end if
        end do

        close (iunit)

        print *, "Plasma data written to: ", trim(filename)
    end subroutine write_plasma_data_to_file

    !> @brief Write profile_data array to a file for debugging
    !> @param[in] profile_data 2D array with profile data
    !> @param[in] filename Output filename
    !> @param[in] r Radial coordinate array (optional)
    subroutine write_profile_data_to_file(profile_data, filename, r)
        real(dp), dimension(:, :), intent(in) :: profile_data
        character(len=*), intent(in) :: filename
        real(dp), dimension(:), intent(in), optional :: r

        integer :: iunit, i, ios
        character(len=1024) :: errmsg

        open (newunit=iunit, file=filename, status='replace', action='write', iostat=ios)
        if (ios /= 0) then
            write (errmsg, '(A,A)') "write_profile_data_to_file: cannot open file ", filename
            error stop trim(errmsg)
        end if

        ! Write header
        if (present(r)) then
            write (iunit, '(A)') "# Profile data: s_tor, M_t, r [cm]"
            write (iunit, '(A,I0,A,I0)') "# npoints = ", size(profile_data, 1), &
                ", ncols = ", size(profile_data, 2) + 1
        else
            write (iunit, '(A)') "# Profile data: s_tor, M_t"
            write (iunit, '(A,I0,A,I0)') "# npoints = ", size(profile_data, 1), &
                ", ncols = ", size(profile_data, 2)
        end if

        ! Write data
        do i = 1, size(profile_data, 1)
            if (present(r)) then
                write (iunit, '(*(ES23.15E3,2X))') profile_data(i, :), r(i)
            else
                write (iunit, '(*(ES23.15E3,2X))') profile_data(i, :)
            end if
        end do

        close (iunit)

        print *, "Profile data written to: ", trim(filename)
    end subroutine write_profile_data_to_file

    !> @brief Write original (unsplined) KAMEL plasma data to a file
    !> @param[in] filename Output filename
    !> @details Writes original KAMEL data from rc grid and params array.
    !>          Format matches plasma_data: ni1, ni2, Ti1, Ti2, Te
    !>          where rc is used as the radial coordinate and params provides
    !>          the plasma quantities (ni from params(1,:), Ti from params(4,:),
    !>          Te from params(3,:)).
    subroutine write_original_plasma_data_to_file(filename)
        use baseparam_mod, only: EV_TO_ERG => ev
        use grid_mod, only: rc
        use plasma_parameters, only: params

        character(len=*), intent(in) :: filename

        integer :: iunit, i, ios, npoints
        character(len=1024) :: errmsg
        real(dp), parameter :: ERG_TO_EV = 1.0_dp / EV_TO_ERG

        npoints = size(rc)

        ! Check dimension consistency
        if (size(params, 2) /= npoints) then
            error stop &
                "write_original_plasma_data_to_file: params second dimension mismatch with rc"
        end if

        open (newunit=iunit, file=filename, status='replace', action='write', iostat=ios)
        if (ios /= 0) then
            write (errmsg, '(A,A)') &
                "write_original_plasma_data_to_file: cannot open file ", filename
            error stop trim(errmsg)
        end if

        ! Write header
        write (iunit, '(A)') &
      "# Original KAMEL plasma data: ni1 [cm^-3], ni2 [cm^-3], Ti1 [eV], Ti2 [eV], Te [eV], rc [cm]"
        write (iunit, '(A,I0,A)') "# npoints = ", npoints, ", ncols = 6"

        ! Write data: same column structure as plasma_data but from original rc/params
        do i = 1, npoints
            write (iunit, '(*(ES23.15E3,2X))') &
                rc(i), &  ! Column 6: rc [cm]
                params(1, i), &  ! Column 1: ni1 [cm^-3]
                0.0_dp, &  ! Column 2: ni2 [cm^-3] (zero for single species)
                params(4, i) * ERG_TO_EV, &  ! Column 3: Ti1 [eV]
                1.0_dp, &  ! Column 4: Ti2 [eV] (dummy value)
                params(3, i) * ERG_TO_EV  ! Column 5: Te [eV]
        end do

        close (iunit)

        print *, "Original plasma data written to: ", trim(filename)
    end subroutine write_original_plasma_data_to_file

    !> @brief Write original (unsplined) KAMEL profile data to a file
    !> @param[in] filename Output filename
    !> @details Writes original KAMEL data from rc grid.
    !>          Format matches profile_data: M_t
    !>          where M_t is calculated from Omega_tE (calculated internally)
    !>          at rc grid points using the same method as calculate_Omega_tE_splined.
    subroutine write_original_profile_data_to_file(filename)
        use baseparam_mod, only: am, p_mass, rtor, btor, c
        use grid_mod, only: rb, rc
        use plasma_parameters, only: params, qsaf
        use wave_code_data, only: dPhi0

        character(len=*), intent(in) :: filename

        integer :: iunit, i, ios, npoints, rb_size
        character(len=1024) :: errmsg
        real(dp) :: T_i, m_i, vth, M_t, dPhi_dr, dpsi_pol_dr, Omega_tE

        npoints = size(rc)
        rb_size = size(rb)

        ! Check dimension consistency
        if (size(params, 2) /= npoints) then
            error stop &
                "write_original_profile_data_to_file: params second dimension mismatch with rc"
        end if

        if (size(qsaf) /= npoints) then
            error stop "write_original_profile_data_to_file: qsaf size mismatch with rc"
        end if

        if (size(dPhi0) /= rb_size) then
            error stop "write_original_profile_data_to_file: dPhi0 size mismatch with rb"
        end if

        m_i = am * p_mass  ! ion mass

        open (newunit=iunit, file=filename, status='replace', action='write', iostat=ios)
        if (ios /= 0) then
            write (errmsg, '(A,A)') "write_original_profile_data_to_file: cannot open file ", &
                filename
            error stop trim(errmsg)
        end if

        ! Write header
        write (iunit, '(A)') "# Original KAMEL profile data: M_t, rc [cm]"
        write (iunit, '(A,I0,A)') "# npoints = ", npoints, ", ncols = 2"

        ! Write data: same column structure as profile_data but from original rc/params
        do i = 1, npoints
          ! Calculate toroidal electric precession frequency (same as in calculate_Omega_tE_splined)
            ! Omega_tE = -c * (dPhi/dr) / (dpsi_pol/dr)
            ! Note: dPhi0 is on rb grid, so we spline it to rc
            dPhi_dr = dPhi0(i)
            dpsi_pol_dr = rc(i) * btor / qsaf(i)
            Omega_tE = -c * dPhi_dr / dpsi_pol_dr

            ! Calculate thermal velocity from original Ti
            T_i = params(4, i)  ! ion temperature in erg
            vth = sqrt(2 * T_i / m_i)  ! thermal velocity

            ! Calculate ExB Mach number: M_t = Omega_tE * R_0 / v_th
            M_t = Omega_tE * rtor / vth

            write (iunit, '(*(ES23.15E3,2X))') &
                rc(i), &  ! Column 1: rc [cm]
                M_t  ! Column 2: M_t (ExB Mach number)
        end do

        close (iunit)

        print *, "Original profile data written to: ", trim(filename)
    end subroutine write_original_profile_data_to_file

    !> @brief Write NEO-RT transport data to a file
    !> @param[in] transport_data Transport data from NEO-RT calculation
    !> @param[in] filename Output filename
    !> @param[in] s_tor Normalized toroidal flux coordinate (optional)
    !> @details Writes comprehensive transport data including summary transport
    !>          coefficients, torque data, and per-harmonic contributions.
    !>          Output format:
    !>          - Summary section: Total transport coefficients and torque
    !>          - Per-harmonic section: Individual resonance contributions
    subroutine write_neort_transport_data_to_file(transport_data, filename, s_tor)
        type(transport_data_t), intent(in) :: transport_data
        character(len=*), intent(in) :: filename
        real(dp), intent(in), optional :: s_tor

        integer :: iunit, i, ios, nharmonics
        character(len=1024) :: errmsg

        open (newunit=iunit, file=filename, status='replace', action='write', iostat=ios)
        if (ios /= 0) then
            write (errmsg, '(A,A)') &
                "write_neort_transport_data_to_file: cannot open file ", filename
            error stop trim(errmsg)
        end if

        ! Write header
        write (iunit, '(A)') "# NEO-RT Transport Data Output"
        write (iunit, '(A)') "#"
        if (present(s_tor)) then
            write (iunit, '(A,ES23.15E3)') "# s_tor = ", s_tor
        end if
        write (iunit, '(A)') "#"

        ! Write summary transport coefficients
        write (iunit, '(A)') "# ========================================="
        write (iunit, '(A)') "# SUMMARY TRANSPORT COEFFICIENTS"
        write (iunit, '(A)') "# ========================================="
        write (iunit, '(A)') "#"
        write (iunit, '(A,ES23.15E3)') "# M_t (Toroidal Mach number) = ", &
            transport_data%summary%M_t
        write (iunit, '(A)') "#"
        write (iunit, '(A)') "# Co-passing particles:"
        write (iunit, '(A,2(ES23.15E3,2X))') "#   Dco(D11, D12) = ", &
            transport_data%summary%Dco(1), transport_data%summary%Dco(2)
        write (iunit, '(A)') "#"
        write (iunit, '(A)') "# Counter-passing particles:"
        write (iunit, '(A,2(ES23.15E3,2X))') "#   Dctr(D11, D12) = ", &
            transport_data%summary%Dctr(1), transport_data%summary%Dctr(2)
        write (iunit, '(A)') "#"
        write (iunit, '(A)') "# Trapped particles:"
        write (iunit, '(A,2(ES23.15E3,2X))') "#   Dt(D11, D12) = ", &
            transport_data%summary%Dt(1), transport_data%summary%Dt(2)
        write (iunit, '(A)') "#"

        ! Write torque summary if available
        if (transport_data%torque%has_torque) then
            write (iunit, '(A)') "# ========================================="
            write (iunit, '(A)') "# TORQUE SUMMARY"
            write (iunit, '(A)') "# ========================================="
            write (iunit, '(A)') "#"
            write (iunit, '(A,ES23.15E3)') "# s = ", transport_data%torque%s
            write (iunit, '(A,ES23.15E3)') "# dVds = ", transport_data%torque%dVds
            write (iunit, '(A,ES23.15E3)') "# M_t = ", transport_data%torque%M_t
            write (iunit, '(A)') "#"
            write (iunit, '(A,ES23.15E3)') "# Tco (Co-passing torque) = ", &
                transport_data%torque%Tco
            write (iunit, '(A,ES23.15E3)') "# Tctr (Counter-passing torque) = ", &
                transport_data%torque%Tctr
            write (iunit, '(A,ES23.15E3)') "# Tt (Trapped torque) = ", &
                transport_data%torque%Tt
            write (iunit, '(A,ES23.15E3)') "# Total torque = ", &
                transport_data%torque%Tco + transport_data%torque%Tctr + transport_data%torque%Tt
            write (iunit, '(A)') "#"
        end if

        ! Write per-harmonic data if available
        if (allocated(transport_data%harmonics)) then
            nharmonics = size(transport_data%harmonics)
            write (iunit, '(A)') "# ========================================="
            write (iunit, '(A)') "# PER-HARMONIC TRANSPORT DATA"
            write (iunit, '(A)') "# ========================================="
            write (iunit, '(A)') "#"
            write (iunit, '(A,I0)') "# Number of harmonics = ", nharmonics
            write (iunit, '(A)') "#"
            write (iunit, '(A)') &
  "# Columns: mth, Dresco(1:2), Dresctr(1:2), Drest(1:2), Tresco, Tresctr, Trest, vminp/vth, &
&vmaxp/vth, vmint/vth, vmaxt/vth"
            write (iunit, '(A)') "#"

            do i = 1, nharmonics
                write (iunit, '(I4,2X,12(ES23.15E3,2X))') &
                    transport_data%harmonics(i)%mth, &
                    transport_data%harmonics(i)%Dresco(1), &
                    transport_data%harmonics(i)%Dresco(2), &
                    transport_data%harmonics(i)%Dresctr(1), &
                    transport_data%harmonics(i)%Dresctr(2), &
                    transport_data%harmonics(i)%Drest(1), &
                    transport_data%harmonics(i)%Drest(2), &
                    transport_data%harmonics(i)%Tresco, &
                    transport_data%harmonics(i)%Tresctr, &
                    transport_data%harmonics(i)%Trest, &
                    transport_data%harmonics(i)%vminp_over_vth, &
                    transport_data%harmonics(i)%vmaxp_over_vth, &
                    transport_data%harmonics(i)%vmint_over_vth, &
                    transport_data%harmonics(i)%vmaxt_over_vth
            end do
        end if

        close (iunit)

        print *, "NEO-RT transport data written to: ", trim(filename)
    end subroutine write_neort_transport_data_to_file

    !> @brief Write torque data to CSV file
    !> @details Writes torque components (Tco, Tctr, Tt) for multiple flux surfaces
    !>          to a CSV file with header row for easy import into analysis tools
    !> @param[in] s_tor Array of normalized toroidal flux coordinates
    !> @param[in] transport_data Array of transport data containing torque information
    !> @param[in] filename Output CSV filename
    subroutine write_torque_csv(s_tor, transport_data, filename)
        use iso_fortran_env, only: dp => real64

        real(dp), dimension(:), intent(in) :: s_tor
        type(transport_data_t), dimension(:), intent(in) :: transport_data
        character(len=*), intent(in) :: filename

        integer :: iunit, ios, i, n
        character(len=1024) :: errmsg

        n = size(s_tor)
        if (size(transport_data) /= n) then
            error stop "write_torque_csv: s_tor and transport_data size mismatch"
        end if

        open (newunit=iunit, file=filename, status='replace', action='write', iostat=ios)
        if (ios /= 0) then
            write (errmsg, '(A,A)') "write_torque_csv: cannot open file ", filename
            error stop trim(errmsg)
        end if

        ! Write header
        write (iunit, '(A)') "s_tor,Tco,Tctr,Tt"

        ! Write data rows
        do i = 1, n
            write (iunit, '(ES23.15E3,A,ES23.15E3,A,ES23.15E3,A,ES23.15E3)') &
                s_tor(i), ",", &
                transport_data(i)%torque%Tco, ",", &
                transport_data(i)%torque%Tctr, ",", &
                transport_data(i)%torque%Tt
        end do

        close (iunit)
    end subroutine write_torque_csv

    subroutine write_plasma_input(path, nplasma, am1, am2, Z1, Z2, plasma)
        character(len=*), intent(in) :: path
        integer, intent(in) :: nplasma
        real(dp), intent(in) :: am1, am2, Z1, Z2
        real(dp), intent(in) :: plasma(:, :)

        integer :: k
        integer, parameter :: fd = 1

        open (fd, file=path, status="replace")
        write (fd, '(A)') " % N am1 am2 Z1 Z2"
        write (fd, '(I0, 4(ES24.16))') nplasma, am1, am2, Z1, Z2
        write (fd, '(A)') " % s ni_1[cm^-3] ni_2[cm^-3] Ti_1[eV] Ti_2[eV] Te[eV]"
        do k = 1, nplasma
            write (fd, '(6(ES24.16))') plasma(k, :)
        end do
        close (fd)
    end subroutine write_plasma_input

    subroutine write_profile_input(path, data)
        character(len=*), intent(in) :: path
        real(8), intent(in) :: data(:, :)

        integer :: k
        integer, parameter :: fd = 1

        open (fd, file=path, status="replace")
        do k = 1, size(data, 1)
            write (fd, '(2(ES24.16))') data(k, 1), data(k, 2)
        end do
        close (fd)
    end subroutine write_profile_input

    subroutine write_transport_data_to_files(data, base_path)
        type(transport_data_t), intent(in) :: data
        character(len=*), intent(in) :: base_path

        integer :: k
        real(8) :: total_D1, total_D2
        integer, parameter :: unit1 = 9
        integer, parameter :: unit2 = 10

        open (unit=unit1, file=trim(adjustl(base_path))//".out", recl=1024)
        write (unit1, *) "# M_t D11co D11ctr D11t D11 D12co D12ctr D12t D12"
        total_D1 = data%summary%Dco(1) + data%summary%Dctr(1) + data%summary%Dt(1)
        total_D2 = data%summary%Dco(2) + data%summary%Dctr(2) + data%summary%Dt(2)
        write (unit1, *) data%summary%M_t, data%summary%Dco(1), data%summary%Dctr(1), &
                        data%summary%Dt(1), total_D1, data%summary%Dco(2), data%summary%Dctr(2), &
                        data%summary%Dt(2), total_D2
        close (unit=unit1)

        if (data%torque%has_torque) then
            open (unit=unit1, file=trim(adjustl(base_path))//"_torque.out", recl=1024)
            write (unit1, *) "# s dVds M_t Tco Tctr Tt"
            write (unit1, *) data%torque%s, data%torque%dVds, data%torque%M_t, data%torque%Tco, &
                             data%torque%Tctr, data%torque%Tt
            close (unit=unit1)
        end if

        open (unit=unit1, file=trim(adjustl(base_path))//"_integral.out", recl=1024)
        open (unit=unit2, file=trim(adjustl(base_path))//"_torque_integral.out", recl=1024)
        do k = 1, size(data%harmonics)
            total_D1 = data%harmonics(k)%Dresco(1) + data%harmonics(k)%Dresctr(1) + &
                       data%harmonics(k)%Drest(1)
            total_D2 = data%harmonics(k)%Dresco(2) + data%harmonics(k)%Dresctr(2) + &
                       data%harmonics(k)%Drest(2)
            write (unit1, *) data%summary%M_t, data%harmonics(k)%mth, data%harmonics(k)%Dresco(1), &
                             data%harmonics(k)%Dresctr(1), data%harmonics(k)%Drest(1), &
                             total_D1, data%harmonics(k)%Dresco(2), data%harmonics(k)%Dresctr(2), &
                             data%harmonics(k)%Drest(2), total_D2, &
                             data%harmonics(k)%vminp_over_vth, data%harmonics(k)%vmaxp_over_vth, &
                             data%harmonics(k)%vmint_over_vth, data%harmonics(k)%vmaxt_over_vth

            write (unit2, *) data%harmonics(k)%mth, data%harmonics(k)%Tresco, &
                             data%harmonics(k)%Tresctr, data%harmonics(k)%Trest
        end do
        close (unit=unit1)
        close (unit=unit2)
    end subroutine write_transport_data_to_files

end module neort_interface
