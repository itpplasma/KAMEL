program kim_main

    use kim_base, only: kim_t
    use kim_mod, only: from_kim_factory_get_kim
    use config, only: type_of_run
    use omp_lib, only: omp_get_wtime

    implicit none

    real :: t_start, t_finish
    integer :: ierr = 0
    class(kim_t), allocatable :: kim_instance

    
    call kim_init
    call from_kim_factory_get_kim(type_of_run, kim_instance)

    t_start = omp_get_wtime()

!    call generate_k_space_grid(100, .true.)  
    
    call generate_grids

    call solve_debye_in_kr_space

    !if (artificial_debye_case .eqv. .true.) then
    !    write(*,*) ' === Artificial Debye case ==='
    !    call fill_spline_kernel_debye(.true.)
    !else
    !    call basis_transformation_of_kernels(.true.)
    !end if

    !call test_sparse_solver(ierr)
    !if (ierr /= 0) then
    !    write(*,*) 'Unit tests failed with sum of errors: ', ierr
    !    stop
    !end if

    ! solve poisson's equation with spline solver
    !call solve_poisson

    t_finish =  omp_get_wtime()

    write(*,*) ' Time: ', (t_finish - t_start), ' s'

end program