module kim_base_m

    implicit none

    public :: kim_t

    type, abstract :: kim_t
        character(100) :: run_type = ""
        contains
            procedure(init), deferred :: init
            procedure(run), deferred :: run
    end type kim_t

    abstract interface
        subroutine init(this)
            import :: kim_t
            class(kim_t), intent(inout) :: this
        end subroutine
    end interface

    abstract interface
        subroutine run(this)
            import :: kim_t
            class(kim_t), intent(inout) :: this
        end subroutine
    end interface

end module