module rt_poisson

    use kim_base, only: kim_t

    implicit none

    type, extends(kim_t) :: poisson_t
        contains
            procedure :: init => init_poisson
            procedure :: run => run_poisson
    end type poisson_t

    contains

    subroutine init_poisson(this)

        implicit none

        class(poisson_t), intent(inout) :: this

        this%run_type = "poisson"
        print *, " ____  __.___   _____  "
        print *, "|    |/ _|   | /     \  "
        print *, "|      < |   |/  \ /  \ "
        print *, "|    |  \|   /    Y    \"
        print *, "|____|__ \___\____|__  /"
        print *, "        \/           \/ "

        call generate_grids

    end subroutine

    subroutine run_poisson(this)

        use kernels, only: fill_rho_kernels

        implicit none
        class(poisson_t), intent(inout) :: this

        ! calculate kernels
        ! write kernels
        call fill_rho_kernels

    
    end subroutine

end module
