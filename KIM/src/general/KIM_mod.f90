module kim_mod_m

    implicit none

    contains

    subroutine from_kim_factory_get_kim(type_of_run, kim_instance)

        use kim_base_m, only: kim_t
        use rt_WKB_dispersion_m, only: WKB_dispersion_t
        use rt_electrostatic_m, only: electrostatic_t
        use rt_electromagnetic_m, only: electromagnetic_t
        use rt_flr2_benchmark_m, only: flr2_benchmark_t

        implicit none

        character(100), intent(in) :: type_of_run
        class(kim_t), allocatable, intent(out) :: kim_instance

        select case(trim(type_of_run))
            case("electrostatic")
                allocate(kim_instance, source=electrostatic_t())
            case("flr2_benchmark")
                allocate(kim_instance, source=flr2_benchmark_t())
            case("WKB_dispersion")
                allocate(kim_instance, source=WKB_dispersion_t())
            case("electromagnetic")
                allocate(kim_instance, source=electromagnetic_t())
            case default
                print *, "Invalid kim type of run " // trim(type_of_run)
                print *, "Options are: electrostatic, electromagnetic, flr2_benchmark, WKB_dispersion"
                stop "Due to invalid type of run"
        end select

    end subroutine

end module
