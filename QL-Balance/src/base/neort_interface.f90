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
        real(dp) :: antenna_factor_exponent
        type(config_t) :: config
    end type meta_config_neort_t

    public :: apply_ntv_transport
    public :: calculate_coarse_s_tor
    public :: calculate_s_tor
    public :: prepare_plasma_data_for_neort
    public :: prepare_profile_data_for_neort
    public :: read_equil_file
    public :: read_neort_meta_config
    public :: meta_config_neort_problem

contains

    subroutine read_neort_meta_config(config_file, meta_config)
        use wave_code_data, only: m_vals, n_vals

        character(len=*), intent(in) :: config_file
        type(meta_config_neort_t), intent(out) :: meta_config

        ! mode numbers
        integer :: poloidal_mode_number
        integer :: toroidal_mode_number

        ! I/O error handling
        integer :: iunit, iostat
        character(len=512) :: iomsg

        ! NEO-RT's own defaults; every key missing from the namelist keeps them
        type(config_t) :: defaults

        ! namelist content
        ! This is a modified version of the full input format used by NEO-RT
        ! as some values are deliberately unused or fixed plus the following:
        ! - explicit paths for the Boozer input files
        ! - amount of equidistant s to compute on
        ! - an exponent for the antenna_factor the resulting torque gets scaled with
        ! new meta config values:
        character(len=1024) :: boozer_file
        character(len=1024) :: boozer_pert_file
        integer :: amount_of_s
        real(dp) :: antenna_factor_exponent
        ! original NEO-RT config values
        real(dp) :: epsmn
        ! m0 is set from the internal poloidal mode
        ! mph is set from the internal toroidal mode
        ! comptorque must be .true.
        logical :: supban
        logical :: magdrift
        integer :: magdrift_passing
        logical :: nopassing
        logical :: noshear
        ! pertfile is .true. if boozer_pert_file is non-empty
        logical :: nonlin
        real(dp) :: bfac
        real(dp) :: efac
        integer :: inp_swi
        integer :: inp_swi_pert
        integer :: vsteps
        integer :: mth_max_abs
        real(dp) :: vmax_over_vth
        integer :: log_level
        ! output_format not needed, nothing is written to disk

        namelist /NEORT/ boozer_file, boozer_pert_file, amount_of_s, antenna_factor_exponent, &
            epsmn, supban, magdrift, magdrift_passing, nopassing, noshear, nonlin, bfac, efac, &
            inp_swi, inp_swi_pert, vsteps, mth_max_abs, vmax_over_vth, log_level

        ! defaults
        boozer_file = ""
        boozer_pert_file = ""
        amount_of_s = 100
        ! antenna_factor is the squared field ratio (compute_antenna_factor_from_Ipar),
        ! so the quasilinear NTV torque ~ (delta B)^2 scales with its first power
        antenna_factor_exponent = 1.0_dp
        epsmn = defaults%epsmn
        supban = defaults%supban
        magdrift = defaults%magdrift
        magdrift_passing = defaults%magdrift_passing
        nopassing = defaults%nopassing
        noshear = defaults%noshear
        nonlin = defaults%nonlin
        bfac = defaults%bfac
        efac = defaults%efac
        inp_swi = defaults%inp_swi
        inp_swi_pert = defaults%inp_swi_pert
        vsteps = defaults%vsteps
        mth_max_abs = defaults%mth_max_abs
        vmax_over_vth = defaults%vmax_over_vth
        log_level = defaults%log_level

        ! read namelist
        open (newunit=iunit, file=config_file, status="old", action="read", iostat=iostat, &
              iomsg=iomsg)
        if (iostat /= 0) then
            write (*, *) "Cannot open "//trim(config_file)//": "//trim(iomsg)
            error stop
        end if
        read (iunit, nml=NEORT, iostat=iostat, iomsg=iomsg)
        close (iunit)

        if (iostat /= 0) then
            write (*, *) "type_of_run == 'TimeEvolutionNTV' specified, but "// &
                "no NEORT config namelist was found!"
            write (*, *) "Error message: "//trim(iomsg)
            error stop
        end if

        ! mode numbers (from KiLCA input)
        poloidal_mode_number = abs(m_vals(1))
        toroidal_mode_number = n_vals(1)

        ! set values from namelist in meta config struct
        meta_config%boozer_file = trim(adjustl(boozer_file))
        meta_config%boozer_pert_file = trim(adjustl(boozer_pert_file))
        meta_config%amount_of_s = amount_of_s
        meta_config%antenna_factor_exponent = antenna_factor_exponent
        meta_config%config%epsmn = epsmn
        meta_config%config%m0 = poloidal_mode_number
        meta_config%config%mph = toroidal_mode_number
        meta_config%config%comptorque = .true.  ! must be true
        meta_config%config%supban = supban
        meta_config%config%magdrift = magdrift
        meta_config%config%magdrift_passing = magdrift_passing
        meta_config%config%nopassing = nopassing
        meta_config%config%noshear = noshear
        meta_config%config%pertfile = len_trim(boozer_pert_file) > 0
        meta_config%config%nonlin = nonlin
        meta_config%config%bfac = bfac
        meta_config%config%efac = efac
        meta_config%config%inp_swi = inp_swi
        meta_config%config%inp_swi_pert = inp_swi_pert
        meta_config%config%vsteps = vsteps
        meta_config%config%mth_max_abs = mth_max_abs
        meta_config%config%vmax_over_vth = vmax_over_vth
        meta_config%config%log_level = log_level
        ! output_format irrelevant, write_output not called by the lib interface
    end subroutine read_neort_meta_config

    !> @brief Describe why a meta config cannot be used; empty if it can.
    function meta_config_neort_problem(meta_config) result(problem)
        type(meta_config_neort_t), intent(in) :: meta_config
        character(len=256) :: problem

        problem = ""
        if (len_trim(meta_config%boozer_file) == 0) then
            problem = "NEORT: boozer_file must be given"
        else if (len_trim(meta_config%boozer_pert_file) == 0) then
            ! Without a perturbation file NEO-RT would use an analytic
            ! perturbation of relative amplitude epsmn with KiLCA's mode numbers,
            ! whose sign convention is not defined for NEO-RT's Boozer angles.
            problem = "NEORT: boozer_pert_file must be given (analytic perturbation unsupported)"
        else if (meta_config%antenna_factor_exponent < 0.0_dp) then
            problem = "NEORT: antenna_factor_exponent must not be negative"
        else if (meta_config%amount_of_s < 2) then
            problem = "NEORT: amount_of_s must be at least 2"
        end if
    end function meta_config_neort_problem

    !> @brief Prepare plasma profile data for NEO-RT from KAMEL arrays
    !> @param[out] plasma_data 2D array (nflux, 6) for NEO-RT plasma input
    !> @details Passes KAMEL data to NEO-RT. Data format expected by NEO-RT:
    !>          Column 1: s_tor - normalized toroidal flux (dimensionless)
    !>          Column 2: ni1   - density of species 1 [cm^-3]
    !>          Column 3: ni2   - density of species 2 [cm^-3]
    !>          Column 4: Ti1   - temperature of species 1 [eV]
    !>          Column 5: Ti2   - temperature of species 2 [eV]
    !>          Column 6: Te    - electron temperature [eV]
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
    !> @param[out] profile_data 2D array (nflux, 2) for NEO-RT profile input:
    !>             Column 1: normalized toroidal flux s [dimensionless]
    !>             Column 2: toroidal Mach number M_t [dimensionless]
    subroutine prepare_profile_data_for_neort(profile_data, r, s_tor)
        use baseparam_mod, only: am, btor, c, p_mass, rtor
        use grid_mod, only: rb, rc
        use plasma_parameters, only: params, qsaf
        use rhs_balance_m, only: compute_unsmoothed_radial_electric_field
        use spline, only: spline_coeff, spline_val

        ! Calculate toroidal electric precession frequency Omega_tE
        ! The prime (') denotes the spatial derivative d/dr:
        !   <Omega_tE>_b = -c / psi'_pol * <dPhi/dr>_b
        !     [Albert thesis (3.103)]
        !   q = psi'_tor / psi'_pol
        !     [Albert thesis (A.14)]
        !   psi_tor = r^2 * B_tor / 2
        !     [Markl2023 (37)]
        !
        ! E_r_boundary is the radial electric field E_r = -dPhi/dr.
        ! Therefore, the implemented electric-field form is:
        !   Omega_tE = c * E_r / psi'_pol
        !
        ! The effective-radius relation gives:
        !   dpsi_pol/dr = psi'_pol = r * B_tor / q
        !   Psi_tor = r^2 * pi * B_tor
        !   Psi_tor' = 2 * pi * r * B_tor
        ! Our psi_tor is Psi_tor / (2*pi), so psi_tor' = r * B_tor.
        !
        ! Source details for the relations above:
        ! Albert defines the bounce-averaged ExB contribution in
        ! [Albert2020, Eq. (3.103), printed p. 71, PDF p. 74].
        ! In Albert's flux-coordinate convention:
        !   psi'_pol = sqrt(g) * B^theta
        !   q = psi'_tor / psi'_pol
        ! [Albert2020, Eqs. (A.12) and (A.14), printed p. 122, PDF p. 125].
        ! The shorthand dpsi_pol/dr = B^theta omits the Jacobian factor;
        ! the implemented denominator is the flux derivative psi'_pol.
        !
        ! Markl defines the effective radius through:
        !   r = sqrt(2 * psi_tor / B_ref)
        ! [Markl2023 (37), PDF p. 10].
        ! This assumes a circular toroidal cross-section and fixed B_ref.
        ! Therefore, psi'_tor = r * B_ref and psi'_pol = r * B_ref / q.
        ! In this convention, KAMEL's btor is the signed reference field.
        ! Thus, B_tor in the shorthand above means this fixed B_ref,
        ! not a radius-dependent local toroidal field.
        ! The effective-radius conversion is implemented in
        ! $KAMEL/common/equil/equil_profiles.f90:421-423 and
        ! $KAMEL/KIM/src/util/IO_collection.f90:153-156.
        !
        ! Albert defines lowercase fluxes per radian by dividing full fluxes
        ! by 2*pi [Albert2020, Eqs. (A.8)-(A.9), printed p. 120, PDF p. 123].
        ! The accepted-state E_r profile is splined from the QL-Balance boundary
        ! grid rb to the NEO-RT radial grid r before applying this formula.

        real(dp), dimension(:, :), intent(out) :: profile_data
        real(dp), dimension(:), intent(in) :: r
        real(dp), dimension(:), intent(in) :: s_tor

        real(dp) :: E_r, T_i, m_i, vth, dpsi_pol_dr, Omega_tE, M_t
        real(dp), dimension(:), allocatable :: E_r_boundary
        real(dp), dimension(:, :), allocatable :: E_r_coeffs, q_coeffs
        real(dp), dimension(:, :), allocatable :: Ti_of_r_coeffs
        real(dp), dimension(:, :), allocatable :: E_r_splined
        real(dp), dimension(:, :), allocatable :: q_splined, Ti_splined
        integer :: i, rb_size, rc_size, s_size

        rb_size = size(rb)
        rc_size = size(rc)
        s_size = size(s_tor)

        ! Check dimensions
        if (size(profile_data, 1) /= s_size .or. size(profile_data, 2) /= 2) then
            error stop "prepare_neort_profile_data: profile_data dimension mismatch"
        end if

        allocate (E_r_boundary(rb_size))
        allocate (E_r_coeffs(rb_size - 1, 5))
        allocate (q_coeffs(rc_size - 1, 5))
        allocate (Ti_of_r_coeffs(rc_size - 1, 5))
        allocate (E_r_splined(s_size, 3))
        allocate (q_splined(s_size, 3))
        allocate (Ti_splined(s_size, 3))

        call compute_unsmoothed_radial_electric_field(params, E_r_boundary)
        E_r_coeffs = spline_coeff(rb, E_r_boundary)
        E_r_splined = spline_val(E_r_coeffs, r)
        q_coeffs = spline_coeff(rc, qsaf)
        q_splined = spline_val(q_coeffs, r)
        Ti_of_r_coeffs = spline_coeff(rc, params(4, :))
        Ti_splined = spline_val(Ti_of_r_coeffs, r)

        m_i = am * p_mass  ! ion mass

        ! Fill profile data array
        do i = 1, s_size
            ! Column 1: Normalized toroidal flux s (0 to 1)
            profile_data(i, 1) = s_tor(i)

            ! Calculate electric precession frequency from the accepted-state field.
            E_r = E_r_splined(i, 1)
            dpsi_pol_dr = r(i) * btor / q_splined(i, 1)
            Omega_tE = c * E_r / dpsi_pol_dr

            ! Calculate thermal velocity [Albert2016, above (51)]
            ! Source detail: Albert2016, PDF p. 14, text before Eq. (51).
            T_i = Ti_splined(i, 1)  ! ion temperature in erg
            vth = sqrt(2 * T_i / m_i)  ! thermal velocity

            ! Calculate toroidal Mach number (ExB Mach number) [Albert2016 below (57)]
            ! Source detail: Albert2016, PDF p. 16, text following Eq. (57).
            ! rtor is the major radius (R)
            M_t = Omega_tE * rtor / vth

            ! Column 2: toroidal Mach number M_t
            profile_data(i, 2) = M_t
        end do
    end subroutine prepare_profile_data_for_neort

    !> @brief Read effective radius and toroidal flux from the equilibrium file
    !> @details Reads the equil_r_q_psi.dat file named by the balance configuration:
    !>          header lines start with '#', data lines begin with the columns
    !>          r [cm], q, psi_pol, psi_tor [G cm^2 per radian]. Further columns
    !>          (7 in equil_profiles.f90, 11 in the profile preprocessor) are ignored.
    !> @param[out] r_eff Effective radius (column 1), allocated here
    !> @param[out] psi_tor Toroidal flux (column 4), allocated here
    subroutine read_equil_file(r_eff, psi_tor)
        use control_mod, only: equil_path

        real(dp), dimension(:), allocatable, intent(out) :: r_eff
        real(dp), dimension(:), allocatable, intent(out) :: psi_tor

        integer :: iunit, ipoi, ios, npoints, nheader
        real(dp) :: r_val, q_val, psi_pol_val, psi_tor_val
        character(len=1024) :: errmsg, line

        open (newunit=iunit, file=trim(equil_path), status='old', action='read', iostat=ios)
        if (ios /= 0) then
            write (errmsg, '(A,A)') "read_equil_file: cannot open file ", trim(equil_path)
            error stop trim(errmsg)
        end if

        ! Skip header lines starting with '#' (with optional leading whitespace)
        nheader = 0
        do
            read (iunit, '(A)', iostat=ios) line
            if (ios /= 0) then
                close (iunit)
                error stop "read_equil_file: error reading file (unexpected end during header)"
            end if
            line = adjustl(line)
            if (len_trim(line) > 0 .and. line(1:1) /= '#') then
                backspace (iunit)
                exit
            end if
            nheader = nheader + 1
        end do

        ! First pass: count the number of data lines
        npoints = 0
        do
            read (iunit, *, iostat=ios) r_val, q_val, psi_pol_val, psi_tor_val
            if (ios /= 0) exit
            npoints = npoints + 1
        end do

        if (npoints == 0) then
            close (iunit)
            error stop "read_equil_file: no data lines found in file"
        end if

        allocate (r_eff(npoints), psi_tor(npoints))

        ! Second pass: read and store data
        rewind (iunit)
        do ipoi = 1, nheader
            read (iunit, '(A)')
        end do
        do ipoi = 1, npoints
            read (iunit, *, iostat=ios) r_val, q_val, psi_pol_val, psi_tor_val
            if (ios /= 0) then
                close (iunit)
                write (errmsg, '(A,I0,A)') "read_equil_file: error reading data line ", ipoi
                error stop trim(errmsg)
            end if
            r_eff(ipoi) = r_val
            psi_tor(ipoi) = psi_tor_val
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

    subroutine apply_ntv_transport(r, transport_data, antenna_factor_exponent)
        use PolyLagrangeInterpolation, only: binsrc
        use control_mod, only: data_verbosity
        use grid_mod, only: rb, torque_density_ntv
        use writeData_m, only: write_NTV_torque_diagnostic
        use spline, only: spline_coeff, spline_val
        use wave_code_data, only: antenna_factor
        use time_evolution, only: antenna_factor_max, time_ind

        real(dp), dimension(:), intent(in) :: r
        type(transport_data_t), dimension(:), intent(in) :: transport_data
        real(dp), intent(in) :: antenna_factor_exponent

        real(dp), dimension(:), allocatable :: total_torque_per_s  ! dT_phi/ds [dyn cm]
        real(dp), dimension(:), allocatable :: dVds  ! dV/ds [cm³]
        real(dp), dimension(:), allocatable :: total_torque_density  ! dT_phi/dV [dyn cm-²]
        real(dp), dimension(:, :), allocatable :: torque_of_r_coeffs
        real(dp), dimension(:, :), allocatable :: torque_splined
        integer :: ntorque, nrb, i, i_separatrix
        real(dp) :: antenna_factor_scaling

        ! extract array sizes
        ntorque = size(transport_data)
        nrb = size(rb)

        allocate (total_torque_per_s(ntorque))
        allocate (dVds(ntorque))

        do i = 1, ntorque
            total_torque_per_s(i) = transport_data(i)%torque%Tco + &
                                    transport_data(i)%torque%Tctr + &
                                    transport_data(i)%torque%Tt
            dVds(i) = transport_data(i)%torque%dVds
        end do

        allocate (total_torque_density(ntorque))
        allocate (torque_of_r_coeffs(ntorque - 1, 5))
        allocate (torque_splined(nrb, 3))

        ! Convert flux-surface torque derivative dT_phi/ds [dyn cm] to
        ! volume torque density dT_phi/dV [dyn cm-²] using dV/ds [cm³].
        total_torque_density = total_torque_per_s / dVds

        if (data_verbosity >= 2) then
            call write_NTV_torque_diagnostic(r, total_torque_density, time_ind)
        end if

        ! spline
        torque_of_r_coeffs = spline_coeff(r, total_torque_density)

        ! find index/index+1 of last element from r in rb
        call binsrc(rb, 1, nrb, r(ntorque), i_separatrix)

        ! make sure to not hit the separatrix
        i_separatrix = i_separatrix - 1

        ! do not extrapolate
        torque_splined = spline_val(torque_of_r_coeffs, rb(1:i_separatrix))

        ! torque for r > separatrix should be zero
        torque_density_ntv(1:i_separatrix) = torque_splined(:, 1)
        torque_density_ntv(i_separatrix + 1:nrb) = 0.0_dp

        ! rescale ntv by antenna factor proportion and
        ! exponent (for different collisionality regimes)
        antenna_factor_scaling = (antenna_factor / antenna_factor_max)**antenna_factor_exponent
        torque_density_ntv = torque_density_ntv * antenna_factor_scaling
    end subroutine apply_ntv_transport

end module neort_interface
