module balanceFactory

    implicit none

    public :: balance_t, fromBalanceFactoryGetBalance

    type, abstract :: balance_t
        character(100) :: runType = ""
        contains
            procedure(init), deferred :: initBalance
            procedure(run), deferred :: runBalance
    end type balance
    
    abstract interface
        subroutine init(this)
            import :: balance_t
            class(balance_t), intent(inout) :: this
        end subroutine run
    end interface
    
    abstract interface
        subroutine run(this)
            import :: balance_t
            class(balance_t), intent(inout) :: this
        end subroutine run
    end interface


    type, extends(balance_t) :: SingleStep
        contains
            procedure :: initBalance => initSingleStep
            procedure :: runBalance => runSingleStep
    end type

    type, extends(balance_t) :: TimeEvolution
        contains
            procedure :: initBalance => initTimeEvolution
            procedure :: runBalance => runTimeEvolution
    end type 

    type, extends(balance_t) :: ParameterScan
        contains
            procedure :: initBalance => initParameterScan
            procedure :: runBalance => runParameterScan
    end type balanceParameterScan

    contains

    subroutine initSingleStep(this)
        class(SingleStep), intent(inout) :: this
        this%runType = "SingleStep"
    end subroutine

    subroutine runSingleStep(this)
        class(SingleStep), intent(inout) :: this
        write(*,*) "Running SingleStep"
    end subroutine

    subroutine initTimeEvolution(this)
        class(TimeEvolution), intent(inout) :: this
        this%runType = "TimeEvolution"
    end subroutine

    subroutine runTimeEvolution(this)
        class(TimeEvolution), intent(inout) :: this
        write(*,*) "Running TimeEvolution"
    end subroutine

    subroutine initParameterScan(this)
        class(ParameterScan), intent(inout) :: this
        this%runType = "ParameterScan"
    end subroutine

    subroutine runParameterScan(this)
        class(ParameterScan), intent(inout) :: this
        write(*,*) "Running ParameterScan"
    end subroutine

    subroutine fromBalanceFactoryGetBalance(type, balanceInstance)

        character(100), intent(in) :: type
        class(balance_t), allocatable, intent(out) :: balanceInstance

        select case(type)
            case("SingleStep")
                allocate(balanceInstance, source=SingleStep())
            case("TimeEvolution")
                allocate(balanceInstance, source=TimeEvolution())
            case("ParameterScan")
                allocate(balanceInstance, source=ParameterScan())
            case default
                write(*,*) "Invalid balance type"
        end select

    end subroutine


end module