module time_evolution_ntv
    use iso_fortran_env, only: dp => real64
    use time_evolution, only: TimeEvolution_t
    use neort_datatypes, only: magfie_data_t, transport_data_t

    implicit none

    type, public, extends(TimeEvolution_t) :: TimeEvolutionNTV_t
    contains
        procedure :: init_balance => initTimeEvolutionNTV
        procedure :: run_balance => runTimeEvolutionNTV
    end type

    private :: initTimeEvolutionNTV
    private :: runTimeEvolutionNTV
    private :: doStep

    ! for NEO-RT
    real(dp), dimension(:), allocatable :: s_tor
    real(dp), dimension(:, :), allocatable :: plasma_data
    real(dp), dimension(:, :), allocatable :: profile_data
    real(dp) :: am1, am2, Z1, Z2, efac, bfac

    ! from NEO-RT
    type(magfie_data_t) :: magfie_data
    type(transport_data_t), dimension(:), allocatable :: transport_data

contains

    subroutine initTimeEvolutionNTV(this)
        use grid_mod, only: npoic
        use logger, only: set_log_level

        class(TimeEvolutionNTV_t), intent(inout) :: this
        real(dp), dimension(:), allocatable :: phi_tor

        call this%TimeEvolution_t%init_balance
        this%runType = "TimeEvolutionNTV"

        ! NEO-RT
        ! allocate (phi_tor(npoic))
        ! allocate (s_tor(npoic))
        ! allocate (plasma_data(npoic, 6))
        ! allocate (profile_data(npoic, 2))
        ! allocate (transport_data(npoic))

        ! call read_toroidal_flux(phi_tor)
        ! call calculate_s_tor(s_tor, phi_tor)
        ! ONLY FOR NOW, DEBUG!
        allocate (phi_tor(dim_for_now))
        allocate (s_tor(dim_for_now))
        allocate (plasma_data(dim_for_now, 6))
        allocate (profile_data(dim_for_now, 2))
        allocate (transport_data(dim_for_now))
        s_tor(1) = 0.3d0
        s_tor(2) = 0.4d0
        s_tor(3) = 0.5d0
        s_tor(4) = 0.6d0
        s_tor(5) = 0.7d0

        ! ! set parameters
        ! s = 0.5d0
        ! efac = 1d0
        ! bfac = 1d0
        ! ! pass same mass and charge for species 1 and 2
        ! am1 = am
        ! am2 = am
        ! Z1 = Z_i
        ! Z2 = Z_i

        call set_log_level(4)  ! for development purposes

        ! ! TODO: loop over

        ! call read_and_set_control("neo-rt/driftorbit")  ! NEO-RT config, TODO: can be done without file
        ! ! TODO: make sure that NEO-RT does not use globals that are not set without this call!
        ! call do_magfie_init("neo-rt/in_file")  ! Boozer field file !! CAVEAT: sets other stuff as well (like R0) , must be kept with file, no other way around
        ! ! call do_magfie_pert_init("neo-rt/in_file_pert") ! Boozer perturbed field file ! Maybe TODO:
        ! call init_profiles(R0)

        ! call prepare_plasma_data_for_neort(plasma_data)
        ! call prepare_profile_data_for_neort(profile_data)

        ! ! calculates only on a fixed s!! => loop over grid
        ! call init_plasma_input(s, npoic, am1, am2, Z1, Z2, plasma_data)
        ! call init_profile_input(s, R0, efac, bfac, profile_data)

        ! call init()
    end subroutine initTimeEvolutionNTV

    subroutine runTimeEvolutionNTV(this)
        use time_evolution, only: time_ind, Nstorage

        class(TimeEvolutionNTV_t), intent(inout) :: this

        do time_ind = 1, Nstorage
            call doStep(this)
        end do
    end subroutine runTimeEvolutionNTV

    subroutine doStep(this)
        use baseparam_mod, only: R0 => rtor, btor
        use do_magfie_mod, only: do_magfie_init
        use do_magfie_pert_mod, only: do_magfie_pert_init
        use neort, only: read_and_set_control, init, compute_transport
        use neort_interface, only: prepare_plasma_data_for_neort, prepare_profile_data_for_neort
        use neort_profiles, only: init_profiles, init_plasma_input, init_profile_input
        use time_evolution, only: time_ind, doStepBase => doStep

        class(TimeEvolutionNTV_t), intent(inout) :: this

        integer :: s_idx
        integer :: s_size

        call doStepBase(this%TimeEvolution_t)

        ! NEO-RT
        s_size = size(s_tor)
        do s_idx = 1, s_size
            call read_and_set_control("neo-rt/driftorbit")  ! NEO-RT config, TODO: can be done without file
            ! TODO: make sure that NEO-RT does not use globals that are not set without this call!
            call do_magfie_init("neo-rt/in_file")  ! Boozer field file !! CAVEAT: sets other stuff as well (like R0) , must be kept with file, no other way around
            ! call do_magfie_pert_init("neo-rt/in_file_pert") ! Boozer perturbed field file ! Maybe TODO:
            call init_profiles(R0)

            call prepare_plasma_data_for_neort(plasma_data, s_tor)
            call prepare_profile_data_for_neort(profile_data, s_tor)

            call init_plasma_input(s_tor(s_idx), s_size, am1, am2, Z1, Z2, plasma_data)
            call init_profile_input(s_tor(s_idx), R0, efac, bfac, profile_data)

            call init()

            ! TODO: omp loop over all s, do init and compute_transport for each s here
            ! call init_thermodynamic_forces(psi_pr, q)

            call compute_transport(transport_data(time_ind))

            ! TODO: Apply NEO-RT transport coefficients back to KAMEL
            ! This would involve updating the transport coefficient arrays in grid_mod
            ! For example:
            ! call apply_ntv_transport(D11_ntv, D12_ntv, torque_ntv)
        end do
    end subroutine doStep

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

end module
