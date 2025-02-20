program test_D_one_over_nu

    use balanceBase, only: balance_t
    use balance_mod
    use parallelTools

    implicit none
    class(balance_t), allocatable :: balance_instance
    character(100) :: typ_of_run = "SingleStep"

    call read_config
    call from_balance_factory_get_balance(typ_of_run, balance_instance)
    call initMPI
    call balance_instance%init_balance()

    call get_D_one_over_nu

    call MPI_finalize(ierror)

end program