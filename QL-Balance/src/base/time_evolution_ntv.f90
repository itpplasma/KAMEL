module time_evolution_ntv
    use time_evolution, only: TimeEvolution_t

    implicit none

    type, public, extends(TimeEvolution_t) :: TimeEvolutionNTV_t
        contains
            procedure :: init_balance => initTimeEvolutionNTV
            procedure :: run_balance => runTimeEvolutionNTV
    end type

    private :: initTimeEvolutionNTV
    private :: runTimeEvolutionNTV

contains

    subroutine initTimeEvolutionNTV(this)
        use do_magfie_mod, only: R0, s, bfac, do_magfie_init
        use neort_profiles, only: init_profiles
        use logger, only: set_log_level
        use neort, only: read_and_set_control, read_and_init_plasma_input, &
                         read_and_init_profile_input, init
        use driftorbit, only: efac

        class(TimeEvolutionNTV_t), intent(inout) :: this

        call this%TimeEvolution_t%init_balance
        this%runType = "TimeEvolutionNTV"

        ! NEO-RT
        call read_and_set_control("neo-rt/driftorbit")
        call do_magfie_init("neo-rt/in_file")
        call init_profiles(R0)
        call read_and_init_plasma_input("neo-rt/plasma.in", s)
        call read_and_init_profile_input("neo-rt/profile.in", s, R0, efac, bfac)

        call init
    end subroutine initTimeEvolutionNTV

    subroutine runTimeEvolutionNTV(this)
        class(TimeEvolutionNTV_t), intent(inout) :: this

        call this%TimeEvolution_t%run_balance()
    end subroutine runTimeEvolutionNTV
end module
