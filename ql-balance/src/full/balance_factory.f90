module balance_factory

    implicit none

    type, abstract :: balance_t
    contains
        procedure :: run => balance_run
    end type balance

    abstract interface
        function run(this) result(s)
            import :: balance_t
            class(balance_t), intent(inout) :: this
        end function run
    end interface

end module