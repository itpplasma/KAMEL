!> @file
!> This is the main program file.

!> @details This program runs the balance code that solves the balance equations described in Heyn et al. NF2014.
!> The solution of the balance code includes envoking KiLCA. Note that KiLCA will read the profiles still from
!> ASCII files, in contrast to the balance code, which can read them from a HDF5 file.
program ql_balance

    use balance_base
    use plasma_parameters
    use balance_mod

    implicit none

    class(balance_t), allocatable :: balance_instance

    write(*,*) ' ________  .____             __________        .__                              '
    write(*,*) ' \_____  \ |    |            \______   \_____  |  | _____    ____   ____  ____  '
    write(*,*) '  /  / \  \|    |      ______ |    |  _/\__  \ |  | \__  \  /    \_/ ___\/ __ \ '
    write(*,*) ' /   \_/.  \    |___  /_____/ |    |   \ / __ \|  |__/ __ \|   |  \  \__\  ___/ '
    write(*,*) ' \_____\ \_/_______ \         |______  /(____  /____(____  /___|  /\___  >___  >'
    write(*,*) '        \__>       \/                \/      \/          \/     \/     \/    \/ '

    call read_config
    call from_balance_factory_get_balance(type_of_run, balance_instance)

    call balance_instance%init_balance()
    call balance_instance%run_balance()

end program ql_balance
