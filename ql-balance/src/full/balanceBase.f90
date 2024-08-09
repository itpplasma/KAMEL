module balanceBase

    implicit none

    public :: balance_t

    type, abstract :: balance_t
        character(100) :: runType = ""
        contains
            procedure(init), deferred :: initBalance
            procedure(run), deferred :: runBalance
    end type balance_t
    
    abstract interface
        subroutine init(this)
            import :: balance_t
            class(balance_t), intent(inout) :: this
        end subroutine
    end interface
    
    abstract interface
        subroutine run(this)
            import :: balance_t
            class(balance_t), intent(inout) :: this
        end subroutine
    end interface


end module