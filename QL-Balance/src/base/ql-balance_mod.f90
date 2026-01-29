
module balance_mod

    implicit none

    contains

    subroutine from_balance_factory_get_balance(type_of_run, balance_instance)

        use balance_base, only: balance_t
        use singleStep, only: SingleStep_t
        use time_evolution, only: TimeEvolution_t
        use time_evolution_stellarator, only: time_evolution_stellarator_t
        use paramscan_mod, only: ParameterScan_t

        implicit none

        character(100), intent(in) :: type_of_run
        class(balance_t), allocatable, intent(out) :: balance_instance

        select case(trim(type_of_run))
            case("SingleStep")
                allocate(balance_instance, source=SingleStep_t())
            case("TimeEvolution")
                allocate(balance_instance, source=TimeEvolution_t())
            case("ParameterScan")
                allocate(balance_instance, source=ParameterScan_t())
            case("TimeEvolutionStellarator")
                allocate(balance_instance, source=time_evolution_stellarator_t())
            !case("TimeEvolutionParameterScan") ! TODO
            !case("ConstantPsi") ! TODO
            case default
                print *, "Invalid balance type of run " // trim(type_of_run)
                print *, "Options are: SingleStep, TimeEvolution, ParameterScan"
                stop "Due to invalid type of run"
        end select

    end subroutine


end module
