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
        class(TimeEvolutionNTV_t), intent(inout) :: this

        call this%TimeEvolution_t%init_balance
        this%runType = "TimeEvolutionNTV"
    end subroutine initTimeEvolutionNTV

    subroutine runTimeEvolutionNTV(this)
        class(TimeEvolutionNTV_t), intent(inout) :: this

        call this%TimeEvolution_t%run_balance()
    end subroutine runTimeEvolutionNTV
end module
