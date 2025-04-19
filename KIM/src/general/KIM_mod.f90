module kim_mod

    implicit none

    contains

    subroutine from_kim_factory_get_kim(type_of_run, kim_instance)

        use kim_base, only: kim_t
        use rt_WKB_dispersion, only: WKB_dispersion_t
        use rt_poisson, only: poisson_t
        use rt_reduced, only: reduced_t

        implicit none

        character(100), intent(in) :: type_of_run
        class(kim_t), allocatable, intent(out) :: kim_instance

        select case(trim(type_of_run))
            case("standard")
                allocate(kim_instance, source=poisson_t())
            case("reduced")
                allocate(kim_instance, source=reduced_t())
            case("WKB_dispersion")
                allocate(kim_instance, source=WKB_dispersion_t())
            case default
                print *, "Invalid kim type of run " // trim(type_of_run)
                print *, "Options are: WKB_dispersion"
                stop "Due to invalid type of run"
        end select

    end subroutine

end module