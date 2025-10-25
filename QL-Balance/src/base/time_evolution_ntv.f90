module time_evolution_ntv
    use iso_fortran_env, only: dp => real64
    use time_evolution, only: TimeEvolution_t

    implicit none

    type, public, extends(TimeEvolution_t) :: TimeEvolutionNTV_t
        ! for NEO-RT
        real(dp), allocatable :: plasma_data(:, :)
        real(dp) :: am1, am2, Z1, Z2
    contains
        procedure :: init_balance => initTimeEvolutionNTV
        procedure :: run_balance => runTimeEvolutionNTV
    end type

    private :: initTimeEvolutionNTV
    private :: runTimeEvolutionNTV

contains

    subroutine initTimeEvolutionNTV(this)
        use do_magfie_mod, only: R0, s, bfac, do_magfie_init
        use do_magfie_pert_mod, only: do_magfie_pert_init
        use neort_profiles, only: init_profiles
        use logger, only: set_log_level
        use neort, only: read_and_set_control, read_and_init_plasma_input, &
                         read_and_init_profile_input, init, check_magfie
        use driftorbit, only: efac
        use neort_interface, only: prepare_plasma_data_for_neort

        class(TimeEvolutionNTV_t), intent(inout) :: this

        call this % TimeEvolution_t % init_balance
        this % runType = "TimeEvolutionNTV"

        ! NEO-RT
        allocate(this % plasma_data(npoic, 6))
        call set_log_level(4)  ! for development purposes

        call read_and_set_control("neo-rt/driftorbit") ! NEO-RT config
        call do_magfie_init("neo-rt/in_file") ! Boozer field file
        ! call do_magfie_pert_init("neo-rt/in_file_pert") ! Boozer perturbed field file
        call init_profiles(R0) ! minor stuff
        call read_and_init_profile_input("neo-rt/profile.in", s, R0, efac, bfac)
        call prepare_plasma_data_for_neort(this % plasma_data, this % am1, this % am2, this % Z1, &
                                           this % Z2)

        call init
        call check_magfie("neo-rt/qlb")
    end subroutine initTimeEvolutionNTV

    subroutine runTimeEvolutionNTV(this)
        use time_evolution, only: time_ind, Nstorage, doStep
        use neort_profiles, only: init_thermodynamic_forces

        class(TimeEvolutionNTV_t), intent(inout) :: this

        do time_ind = 1, Nstorage
            call doStep(this % TimeEvolution_t)
            ! TODO: prepare data for NEO-RT
            ! TODO: call NEO-RT

            ! update NEO-RT inputs
            ! call init_thermodynamic_forces(psi_pr, q) ! from subroutine init
            ! call init_plasma_input(s, nplasma, am1, am2, Z1, Z2, plasma)
            ! call init_profile_input(s, R0, efac, bfac, data)

            ! call NEO-RT to compute transport
            ! call compute_transport("neo-rt/qlb")
        end do

    end subroutine runTimeEvolutionNTV

end module
