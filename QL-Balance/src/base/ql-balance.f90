!> @file
!> This is the main program file.

!> @details This program runs the balance code that solves the balance equations described in Heyn et al. NF2014.
!> The solution of the balance code includes envoking KiLCA. Note that KiLCA will read the profiles still from
!> ASCII files, in contrast to the balance code, which can read them from a HDF5 file.
program ql_balance

    use balance_base
    use plasma_parameters
    use balance_mod
    use control_mod, only: ihdf5IO
    use run_info_m, only: write_run_info
    use h5mod, only: prepare_output_file

    implicit none

    class(balance_t), allocatable :: balance_instance

    write(*,*) ' ________  .____             __________        .__                              '
    write(*,*) ' \_____  \ |    |            \______   \_____  |  | _____    ____   ____  ____  '
    write(*,*) '  /  / \  \|    |      ______ |    |  _/\__  \ |  | \__  \  /    \_/ ___\/ __ \ '
    write(*,*) ' /   \_/.  \    |___  /_____/ |    |   \ / __ \|  |__/ __ \|   |  \  \__\  ___/ '
    write(*,*) ' \_____\ \_/_______ \         |______  /(____  /____(____  /___|  /\___  >___  >'
    write(*,*) '        \__>       \/                \/      \/          \/     \/     \/    \/ '

    call read_config
    if (ihdf5IO .eq. 1) then
        ! Streaming output (fresh file) unless a restart reads from path2out.
        call prepare_output_file
        ! Start time, git commit and namelist of this run, see run_info.f90.
        call write_run_info
    end if
    call from_balance_factory_get_balance(type_of_run, balance_instance)

    call balance_instance%init_balance()
    call balance_instance%run_balance()

end program ql_balance
