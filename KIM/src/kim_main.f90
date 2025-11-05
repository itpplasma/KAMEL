program kim_main

    use kim_base_m, only: kim_t
    use kim_mod_m, only: from_kim_factory_get_kim
    use config_m, only: type_of_run, hdf5_output
    use omp_lib, only: omp_get_wtime, omp_get_max_threads
    use KIM_kinds_m, only: dp
    use config_display_m, only: display_kim_banner
    use IO_collection_m, only: deinitialize_hdf5_output, h5id
    use KAMEL_hdf5_tools, only: h5_add

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

    if (hdf5_output) then 
        call h5_add(h5id, 'runtime', t_finish - t_start, 'KIM runtime', 's')
        call h5_add(h5id, 'omp_num_threads', omp_get_max_threads(), 'Number of OpenMP threads', '1')
        call deinitialize_hdf5_output()
    end if

end program