module kim_mod

    implicit none

    contains

    subroutine from_kim_factory_get_kim(type_of_run, kim_instance)

        use kim_base, only: kim_t
        use WKB_dispersion, only: WKB_dispersion_t

        implicit none

        character(100), intent(in) :: type_of_run
        class(kim_t), allocatable, intent(out) :: kim_instance

        select case(trim(type_of_run))
            case("WKB_dispersion")
                allocate(kim_instance, source=WKB_dispersion_t())
            case default
                print *, "Invalid kim type of run " // trim(type_of_run)
                print *, "Options are: WKB_dispersion"
                stop "Due to invalid type of run"
        end select

    end subroutine

end module