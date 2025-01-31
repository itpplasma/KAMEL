module WKB_dispersion

    use kim_base, only: kim_t

    implicit none

    type, extends(kim_t) :: WKB_dispersion_t
        contains
            procedure :: init => init_wkb_dispersion
            procedure :: run => run_wkb_dispersion
    end type WKB_dispersion_t

    contains

    subroutine init_wkb_dispersion(this)

        implicit none

        class(WKB_dispersion_t), intent(inout) :: this

        this%run_type = "WKB_dispersion"

        call generate_grids

    end subroutine

    subroutine run_wkb_dispersion(this)

        implicit none

        class(WKB_dispersion_t), intent(inout) :: this

        print *, "Running WKB_dispersion"

    end subroutine

end module