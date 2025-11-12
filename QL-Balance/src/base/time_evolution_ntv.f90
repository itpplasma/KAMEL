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

    ! from NEO-RT
    type(magfie_data_t) :: magfie_data
    type(transport_data_t), dimension(:), allocatable :: transport_data

contains

    subroutine initTimeEvolutionNTV(this)
        use do_magfie_mod, only: do_magfie_init
        use do_magfie_pert_mod, only: do_magfie_pert_init
        use driftorbit, only: pertfile
        use grid_mod, only: npoic
        use logger, only: set_log_level
        use neort, only: read_and_set_control
        use neort_interface, only: read_equil_file, calculate_s_tor, &
                                   prepare_plasma_data_for_neort, prepare_profile_data_for_neort

        class(TimeEvolutionNTV_t), intent(inout) :: this
        real(dp), dimension(:), allocatable :: psi_tor

        integer :: npoic_save

        call this%TimeEvolution_t%init_balance
        this%runType = "TimeEvolutionNTV"

        ! TODO: don't init NEO-RT at all here, only in omp loop during time evolution
        ! NEO-RT
        npoic_save = npoic
        npoic = 10000 ! for now
        allocate (psi_tor(npoic))
        allocate (s_tor(npoic))
        allocate (plasma_data(npoic, 6))
        allocate (profile_data(npoic, 2))
        allocate (transport_data(npoic))

        call read_equil_file(psi_tor=psi_tor)
        call calculate_s_tor(s_tor, psi_tor)
        npoic = npoic_save
        ! ONLY FOR NOW, DEBUG!
        ! s_tor(1) = 0.3d0
        ! s_tor(2) = 0.4d0
        ! s_tor(3) = 0.5d0
        ! s_tor(4) = 0.6d0
        ! s_tor(5) = 0.7d0

        call read_and_set_control("neo-rt/driftorbit")  ! NEO-RT config, TODO: can be done without file
        ! TODO: make sure that NEO-RT does not use globals that are not set without this call!
        call do_magfie_init("neo-rt/in_file")  ! Boozer field file !! CAVEAT: sets other stuff as well (like R0) , must be kept with file, no other way around
        if (pertfile) call do_magfie_pert_init("neo-rt/in_file_pert") ! Boozer perturbed field file ! Maybe TODO:
        call set_log_level(5)  ! for development purposes
        ! call set_log_level(-1)
        call prepare_plasma_data_for_neort(plasma_data, s_tor)
        call prepare_profile_data_for_neort(profile_data, s_tor)

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
        use time_evolution, only: doStepBase => doStep

        class(TimeEvolutionNTV_t), intent(inout) :: this

        integer :: s_size, s_idx

        call doStepBase(this%TimeEvolution_t)

        ! NEO-RT
        s_size = size(s_tor)
        ! TODO: omp loop over all s, do init and compute_transport for each s here
        do s_idx = 1, s_size
            call neo_rt(s_tor(s_idx), s_size, transport_data(s_idx))
            ! TODO: Apply NEO-RT transport coefficients back to KAMEL
            ! This would involve updating the transport coefficient arrays in grid_mod
            ! For example:
            ! call apply_ntv_transport(D11_ntv, D12_ntv, torque_ntv)
        end do
    end subroutine doStep

    subroutine neo_rt(s, s_size, transport_data_)
        use baseparam_mod, only: am, Z_i
        use do_magfie_mod, only: R0, bfac, psi_pr, q, do_magfie
        use do_magfie_pert_mod, only: do_magfie_pert_amp
        use driftorbit, only: pertfile, efac, nopassing, sign_vpar, etamin, etamax
        use neort, only: set_s, compute_transport, set_to_trapped_region
        use neort_freq, only: init_canon_freq_trapped_spline, init_canon_freq_passing_spline
        use neort_magfie, only: init_flux_surface_average
        use neort_profiles, only: init_plasma_input, init_profile_input, init_thermodynamic_forces
        use time_evolution, only: doStepBase => doStep

        real(8), intent(in) :: s
        integer, intent(in) :: s_size
        type(transport_data_t), intent(out) :: transport_data_

        ! for do_magfie, all dummies
        integer, parameter :: dim = 3
        real(8) :: x(dim)
        real(8) :: bmod
        real(8) :: sqrtg
        real(8), dimension(dim) :: bder
        real(8), dimension(dim) :: hcovar
        real(8), dimension(dim) :: hctrvr
        real(8), dimension(dim) :: hcurl
        complex(8) :: dummy

        ! for init_plasma_input
        real(dp) :: am1, am2, Z1, Z2

        call set_s(s)  ! global s

        x(1) = s
        x(2) = 0.0
        x(3) = 0.0
        ! sets most of the globals in module do_magfie_mod
        call do_magfie(x, bmod, sqrtg, bder, hcovar, hctrvr, hcurl)
        if (pertfile) call do_magfie_pert_amp(x, dummy)
        ! note: init_profiles is useless here, as it is done again in init_profile_input

        ! pass same mass and charge for species 1 and 2
        am1 = am
        am2 = am
        Z1 = Z_i
        Z2 = Z_i
        call init_plasma_input(s, s_size, am1, am2, Z1, Z2, plasma_data)
        call init_profile_input(s, R0, efac, bfac, profile_data)

        ! subroutine init
        call init_flux_surface_average(s)
        call init_canon_freq_trapped_spline()  ! sets etamin and etamax
        if (.not. nopassing) call init_canon_freq_passing_spline()
        sign_vpar = 1
        call set_to_trapped_region(etamin, etamax)  ! sets etamin and etamax again
        ! psi_pr is the torodial flux at plasma boundary, fixed for all s
        call init_thermodynamic_forces(psi_pr, q)  ! psi_pr=const., q is set in do_magfie

        call compute_transport(transport_data_)
    end subroutine neo_rt

end module
