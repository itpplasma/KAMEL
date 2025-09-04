program kim_main

    use kim_base_m, only: kim_t
    use kim_mod_m, only: from_kim_factory_get_kim
    use config_m, only: type_of_run
    use omp_lib, only: omp_get_wtime
    use KIM_kinds_m, only: dp
    use config_display_m, only: display_kim_banner

    implicit none

    real(dp) :: t_start, t_finish
    class(kim_t), allocatable :: kim_instance

    ! Display the beautiful KIM banner at startup
    call display_kim_banner()
    
    call kim_init
    print *, "KIM initialized."
    call from_kim_factory_get_kim(type_of_run, kim_instance)

    t_start = omp_get_wtime()

    call kim_instance%init()
    call kim_instance%run()
    
    !call solve_debye_in_kr_space

    t_finish =  omp_get_wtime()

    write(*,*) ' Time: ', (t_finish - t_start), ' s'

end program