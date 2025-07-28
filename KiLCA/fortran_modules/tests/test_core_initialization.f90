program test_core_initialization
    use iso_fortran_env, only: real64, int32, int64
    use iso_c_binding
    use kilca_types_m
    use kilca_core_m
    use kilca_settings_m
    implicit none
    
    ! Test variables
    integer :: test_status = 0
    character(len=256) :: test_path = "/test/project/path/"
    
    ! Test 1: Basic initialization
    call test_basic_initialization()
    
    ! Test 2: Mode independent data initialization
    call test_mode_independent_initialization()
    
    ! Test 3: Mode dependent antenna initialization
    call test_mode_dependent_antenna_initialization()
    
    ! Test 4: Mode dependent eigenmode initialization
    call test_mode_dependent_eigmode_initialization()
    
    ! Test 5: Mode dependent interface initialization
    call test_mode_dependent_interface_initialization()
    
    ! Test 6: Pointer precision check
    call test_pointer_precision_check()
    
    if (test_status == 0) then
        print *, "All initialization tests PASSED!"
    else
        print *, "Some tests FAILED. Status:", test_status
        stop 1
    end if
    
contains

    subroutine test_basic_initialization()
        type(core_data_t), pointer :: cd
        integer :: ierr
        
        print *, "Testing basic initialization..."
        
        ! Create core data structure
        call core_data_create(cd, test_path, ierr)
        if (ierr /= KILCA_SUCCESS) then
            print *, "FAIL: Could not create core_data"
            test_status = test_status + 1
            return
        end if
        
        ! Verify path was set correctly
        if (cd%path2project /= test_path) then
            print *, "FAIL: Path not set correctly"
            test_status = test_status + 1
        end if
        
        ! Clean up
        call core_data_destroy(cd, ierr)
        
        print *, "test_basic_initialization completed"
    end subroutine test_basic_initialization
    
    subroutine test_mode_independent_initialization()
        type(core_data_t), pointer :: cd
        integer :: ierr
        character(len=256) :: vacuum_path, flre_path
        
        print *, "Testing mode independent initialization..."
        
        ! Test vacuum case
        vacuum_path = "/test/vacuum/project/"
        call core_data_create(cd, vacuum_path, ierr)
        if (ierr /= KILCA_SUCCESS) then
            print *, "FAIL: Could not create vacuum core_data"
            test_status = test_status + 1
            return
        end if
        
        ! Initialize mode independent data
        call calc_and_set_mode_independent_core_data(cd, ierr)
        
        ! For this test, we expect failure since files don't exist
        if (ierr == KILCA_SUCCESS) then
            print *, "WARN: Unexpected success - files shouldn't exist"
        end if
        
        call core_data_destroy(cd, ierr)
        
        ! Test FLRE case
        flre_path = "/test/flre/project/"
        call core_data_create(cd, flre_path, ierr)
        if (ierr /= KILCA_SUCCESS) then
            print *, "FAIL: Could not create flre core_data"
            test_status = test_status + 1
            return
        end if
        
        ! Initialize mode independent data
        call calc_and_set_mode_independent_core_data(cd, ierr)
        
        ! For this test, we expect failure since files don't exist
        if (ierr == KILCA_SUCCESS) then
            print *, "WARN: Unexpected success - files shouldn't exist"
        end if
        
        call core_data_destroy(cd, ierr)
        
        print *, "test_mode_independent_initialization completed"
    end subroutine test_mode_independent_initialization
    
    subroutine test_mode_dependent_antenna_initialization()
        type(core_data_t), pointer :: cd
        type(settings_t), pointer :: sd
        type(antenna_t), pointer :: ant
        integer :: ierr
        integer, allocatable :: modes(:)
        
        print *, "Testing mode dependent antenna initialization..."
        
        ! Create and set up core data
        call core_data_create(cd, test_path, ierr)
        if (ierr /= KILCA_SUCCESS) return
        
        ! Create settings
        call settings_create(sd, test_path, ierr)
        if (ierr /= KILCA_SUCCESS) then
            call core_data_destroy(cd, ierr)
            return
        end if
        
        ! Set up antenna settings
        call settings_get_antenna(sd, ant, ierr)
        ant%dma = 2  ! 2 modes
        allocate(modes(4))
        modes = [1, 2, 3, 4]  ! m1=1, n1=2, m2=3, n2=4
        call antenna_set_parameters(ant, 10.0_dp, 1.0_dp, 1000.0_dp, &
                                   cmplx(50.0e6_dp, 0.0_dp, dp), 2, modes, ierr)
        
        ! Associate settings with core data
        cd%sd => sd
        
        ! Initialize mode dependent data
        call calc_and_set_mode_dependent_core_data_antenna(cd, ierr)
        
        ! Verify mode array was allocated
        if (.not. allocated(cd%mda)) then
            print *, "FAIL: Mode array not allocated"
            test_status = test_status + 1
        else
            if (size(cd%mda) /= 2) then
                print *, "FAIL: Mode array size incorrect"
                test_status = test_status + 1
            end if
        end if
        
        ! Clean up
        deallocate(modes)
        call settings_destroy(sd, ierr)
        call core_data_destroy(cd, ierr)
        
        print *, "test_mode_dependent_antenna_initialization completed"
    end subroutine test_mode_dependent_antenna_initialization
    
    subroutine test_mode_dependent_eigmode_initialization()
        type(core_data_t), pointer :: cd
        type(settings_t), pointer :: sd
        type(antenna_t), pointer :: ant
        type(eigmode_sett_t), pointer :: es
        integer :: ierr
        integer, allocatable :: modes(:)
        
        print *, "Testing mode dependent eigenmode initialization..."
        
        ! Create and set up core data
        call core_data_create(cd, test_path, ierr)
        if (ierr /= KILCA_SUCCESS) return
        
        ! Create settings
        call settings_create(sd, test_path, ierr)
        if (ierr /= KILCA_SUCCESS) then
            call core_data_destroy(cd, ierr)
            return
        end if
        
        ! Set up antenna settings
        call settings_get_antenna(sd, ant, ierr)
        ant%dma = 1  ! 1 mode
        allocate(modes(2))
        modes = [1, 1]  ! m=1, n=1
        call antenna_set_parameters(ant, 10.0_dp, 1.0_dp, 1000.0_dp, &
                                   cmplx(50.0e6_dp, 0.0_dp, dp), 1, modes, ierr)
        
        ! Set up eigenmode settings
        call settings_get_eigmode(sd, es, ierr)
        es%search_flag = 1  ! Frequency scan
        
        ! Associate settings with core data
        cd%sd => sd
        
        ! Initialize mode dependent data
        call calc_and_set_mode_dependent_core_data_eigmode(cd, ierr)
        
        ! Verify mode array was allocated
        if (.not. allocated(cd%mda)) then
            print *, "FAIL: Mode array not allocated for eigmode"
            test_status = test_status + 1
        end if
        
        ! Clean up
        deallocate(modes)
        call settings_destroy(sd, ierr)
        call core_data_destroy(cd, ierr)
        
        print *, "test_mode_dependent_eigmode_initialization completed"
    end subroutine test_mode_dependent_eigmode_initialization
    
    subroutine test_mode_dependent_interface_initialization()
        type(core_data_t), pointer :: cd
        type(settings_t), pointer :: sd
        type(antenna_t), pointer :: ant
        integer :: ierr
        integer, allocatable :: modes(:)
        
        print *, "Testing mode dependent interface initialization..."
        
        ! Create and set up core data
        call core_data_create(cd, test_path, ierr)
        if (ierr /= KILCA_SUCCESS) return
        
        ! Create settings
        call settings_create(sd, test_path, ierr)
        if (ierr /= KILCA_SUCCESS) then
            call core_data_destroy(cd, ierr)
            return
        end if
        
        ! Set up antenna settings
        call settings_get_antenna(sd, ant, ierr)
        ant%dma = 2  ! 2 modes
        allocate(modes(4))
        modes = [1, 2, 3, 4]  ! m1=1, n1=2, m2=3, n2=4
        call antenna_set_parameters(ant, 10.0_dp, 1.0_dp, 1000.0_dp, &
                                   cmplx(50.0e6_dp, 0.0_dp, dp), 2, modes, ierr)
        
        ! Associate settings with core data
        cd%sd => sd
        
        ! Initialize mode dependent data for interface
        call calc_and_set_mode_dependent_core_data_antenna_interface(cd, ierr)
        
        ! Verify mode array was allocated
        if (.not. allocated(cd%mda)) then
            print *, "FAIL: Mode array not allocated for interface"
            test_status = test_status + 1
        end if
        
        ! Test single mode interface
        call core_data_destroy(cd, ierr)
        call core_data_create(cd, test_path, ierr)
        cd%sd => sd
        
        ! Initialize with specific m,n values
        call calc_and_set_mode_dependent_core_data_antenna_interface_mn(cd, 2, 3, 0, ierr)
        
        ! Verify single mode was allocated
        if (.not. allocated(cd%mda)) then
            print *, "FAIL: Mode array not allocated for single mode"
            test_status = test_status + 1
        else
            if (cd%dim /= 1) then
                print *, "FAIL: Single mode dimension incorrect"
                test_status = test_status + 1
            end if
        end if
        
        ! Clean up
        deallocate(modes)
        call settings_destroy(sd, ierr)
        call core_data_destroy(cd, ierr)
        
        print *, "test_mode_dependent_interface_initialization completed"
    end subroutine test_mode_dependent_interface_initialization
    
    subroutine test_pointer_precision_check()
        integer :: pp_fortran, pp_c
        integer :: ierr
        
        print *, "Testing pointer precision check..."
        
        ! Get Fortran pointer precision
        call get_pointer_precision(pp_fortran, ierr)
        
        ! Get C pointer size (in bytes)
        pp_c = c_sizeof(c_null_ptr)
        
        ! Verify they match
        if (pp_fortran /= pp_c) then
            print *, "WARN: Pointer size mismatch - Fortran:", pp_fortran, "C:", pp_c
            ! This is a warning, not a failure, as it depends on system
        else
            print *, "Pointer sizes match:", pp_fortran, "bytes"
        end if
        
        print *, "test_pointer_precision_check completed"
    end subroutine test_pointer_precision_check

end program test_core_initialization