program test_kilca_zone
    use iso_fortran_env, only: real64, int32, error_unit
    use kilca_zone_m
    use kilca_types_m
    use kilca_constants_m
    use kilca_settings_m
    use kilca_background_m
    use kilca_mode_m, only: wave_data_t, wave_data_create, wave_data_destroy
    implicit none
    
    ! Test counters
    integer :: total_tests, passed_tests, failed_tests
    logical :: test_passed
    
    ! Initialize test counters
    total_tests = 0
    passed_tests = 0
    failed_tests = 0
    
    print *, "========================================================"
    print *, "KiLCA Zone Management System - Unit Tests"
    print *, "========================================================"
    print *, ""
    
    ! Test 1: Zone type creation and initialization
    call test_zone_creation()
    
    ! Test 2: Zone boundary management
    call test_zone_boundaries()
    
    ! Test 3: Zone type and model settings
    call test_zone_types()
    
    ! Test 4: Zone radial grid management
    call test_zone_radial_grid()
    
    ! Test 5: Zone basis field allocation
    call test_zone_basis_allocation()
    
    ! Test 6: Zone indexing functions
    call test_zone_indexing()
    
    ! Test 7: Zone file I/O operations
    call test_zone_file_operations()
    
    ! Test 8: Zone validation and error handling
    call test_zone_validation()
    
    ! Final summary
    print *, ""
    print *, "========================================================"
    print *, "TEST SUMMARY"
    print *, "========================================================"
    print *, "Total tests run:    ", total_tests
    print *, "Tests passed:       ", passed_tests
    print *, "Tests failed:       ", failed_tests
    print *, "Success rate:       ", real(passed_tests)/real(total_tests)*100.0, "%"
    print *, ""
    
    if (failed_tests == 0) then
        print *, "*** ALL TESTS PASSED! ***"
    else
        print *, "*** SOME TESTS FAILED! ***"
        stop 1
    end if
    
contains

    !---------------------------------------------------------------------------
    ! Test 1: Zone type creation and initialization
    !---------------------------------------------------------------------------
    subroutine test_zone_creation()
        type(zone_extended_t) :: zone
        type(settings_t), target :: settings
        type(background_t), target :: background
        type(wave_data_t), target :: wave_data
        integer :: ierr
        
        call start_test("Zone creation and initialization")
        
        ! Create dependencies
        call settings_create(settings, "dummy", ierr)
        call background_create(background, ierr)
        call wave_data_create(wave_data, 1, 1, cmplx(1.0_real64, 0.0_real64), cmplx(1.0_real64, 0.0_real64))
        
        ! Test zone creation
        call zone_create(zone, settings, background, wave_data, "test_path", 1, ierr)
        
        test_passed = (ierr == 0)
        if (test_passed) then
            test_passed = (zone%index == 1)
            test_passed = test_passed .and. (associated(zone%settings))
            test_passed = test_passed .and. (associated(zone%background))
            test_passed = test_passed .and. (associated(zone%wave_data))
            test_passed = test_passed .and. (len_trim(zone%path) > 0)
        end if
        
        ! Clean up
        call zone_destroy(zone)
        call wave_data_destroy(wave_data)
        call background_destroy(background)
        call settings_destroy(settings)
        
        call end_test(test_passed)
    end subroutine test_zone_creation
    
    !---------------------------------------------------------------------------
    ! Test 2: Zone boundary management
    !---------------------------------------------------------------------------
    subroutine test_zone_boundaries()
        type(zone_extended_t) :: zone
        type(settings_t), target :: settings
        type(background_t), target :: background
        type(wave_data_t), target :: wave_data
        real(real64) :: r1_test, r2_test
        integer :: bc1_test, bc2_test, ierr
        
        call start_test("Zone boundary management")
        
        ! Create test zone
        call settings_create(settings, "dummy", ierr)
        call background_create(background, ierr)
        call wave_data_create(wave_data, 1, 1, cmplx(1.0_real64, 0.0_real64), cmplx(1.0_real64, 0.0_real64))
        call zone_create(zone, settings, background, wave_data, "test_path", 1, ierr)
        
        ! Test boundary setting
        r1_test = 0.1_real64
        r2_test = 0.9_real64
        bc1_test = BOUNDARY_CENTER
        bc2_test = BOUNDARY_INTERFACE
        
        call zone_set_boundaries(zone, r1_test, r2_test, bc1_test, bc2_test, ierr)
        
        test_passed = (ierr == 0)
        if (test_passed) then
            test_passed = (abs(zone%r1 - r1_test) < 1.0e-14_real64)
            test_passed = test_passed .and. (abs(zone%r2 - r2_test) < 1.0e-14_real64)
            test_passed = test_passed .and. (zone%bc1 == bc1_test)
            test_passed = test_passed .and. (zone%bc2 == bc2_test)
        end if
        
        ! Test validation of invalid boundaries
        call zone_set_boundaries(zone, r2_test, r1_test, bc1_test, bc2_test, ierr)
        test_passed = test_passed .and. (ierr /= 0)  ! Should fail for r1 > r2
        
        ! Clean up
        call zone_destroy(zone)
        call wave_data_destroy(wave_data)
        call background_destroy(background)
        call settings_destroy(settings)
        
        call end_test(test_passed)
    end subroutine test_zone_boundaries
    
    !---------------------------------------------------------------------------
    ! Test 3: Zone type and model settings
    !---------------------------------------------------------------------------
    subroutine test_zone_types()
        type(zone_extended_t) :: zone
        type(settings_t), target :: settings
        type(background_t), target :: background
        type(wave_data_t), target :: wave_data
        integer :: ierr
        
        call start_test("Zone type and model settings")
        
        ! Create test zone
        call settings_create(settings, "dummy", ierr)
        call background_create(background, ierr)
        call wave_data_create(wave_data, 1, 1, cmplx(1.0_real64, 0.0_real64), cmplx(1.0_real64, 0.0_real64))
        call zone_create(zone, settings, background, wave_data, "test_path", 1, ierr)
        
        ! Test plasma model setting
        call zone_set_plasma_model(zone, PLASMA_MODEL_FLRE, 1, ierr)
        
        test_passed = (ierr == 0)
        if (test_passed) then
            test_passed = (zone%medium == PLASMA_MODEL_FLRE)
            test_passed = test_passed .and. (zone%version == 1)
        end if
        
        ! Test invalid model
        call zone_set_plasma_model(zone, -1, 1, ierr)
        test_passed = test_passed .and. (ierr /= 0)  ! Should fail for invalid model
        
        ! Clean up
        call zone_destroy(zone)
        call wave_data_destroy(wave_data)
        call background_destroy(background)
        call settings_destroy(settings)
        
        call end_test(test_passed)
    end subroutine test_zone_types
    
    !---------------------------------------------------------------------------
    ! Test 4: Zone radial grid management
    !---------------------------------------------------------------------------
    subroutine test_zone_radial_grid()
        type(zone_extended_t) :: zone
        type(settings_t), target :: settings
        type(background_t), target :: background
        type(wave_data_t), target :: wave_data
        real(real64), allocatable :: r_grid(:)
        integer :: grid_dim, ierr, i
        
        call start_test("Zone radial grid management")
        
        ! Create test zone
        call settings_create(settings, "dummy", ierr)
        call background_create(background, ierr)
        call wave_data_create(wave_data, 1, 1, cmplx(1.0_real64, 0.0_real64), cmplx(1.0_real64, 0.0_real64))
        call zone_create(zone, settings, background, wave_data, "test_path", 1, ierr)
        
        ! Create test radial grid
        grid_dim = 100
        allocate(r_grid(grid_dim))
        do i = 1, grid_dim
            r_grid(i) = real(i-1, real64) / real(grid_dim-1, real64)
        end do
        
        ! Set radial grid
        call zone_set_radial_grid(zone, grid_dim, r_grid, ierr)
        
        test_passed = (ierr == 0)
        if (test_passed) then
            test_passed = (zone%dim == grid_dim)
            test_passed = test_passed .and. allocated(zone%r)
            if (associated(zone%r)) then
                test_passed = test_passed .and. (size(zone%r) == grid_dim)
                test_passed = test_passed .and. (abs(zone%r(1) - r_grid(1)) < 1.0e-14_real64)
                test_passed = test_passed .and. (abs(zone%r(grid_dim) - r_grid(grid_dim)) < 1.0e-14_real64)
            end if
        end if
        
        ! Clean up
        deallocate(r_grid)
        call zone_destroy(zone)
        call wave_data_destroy(wave_data)
        call background_destroy(background)
        call settings_destroy(settings)
        
        call end_test(test_passed)
    end subroutine test_zone_radial_grid
    
    !---------------------------------------------------------------------------
    ! Test 5: Zone basis field allocation
    !---------------------------------------------------------------------------
    subroutine test_zone_basis_allocation()
        type(zone_extended_t) :: zone
        type(settings_t), target :: settings
        type(background_t), target :: background
        type(wave_data_t), target :: wave_data
        integer :: nwaves, ncomps, dim, ierr
        
        call start_test("Zone basis field allocation")
        
        ! Create test zone
        call settings_create(settings, "dummy", ierr)
        call background_create(background, ierr)
        call wave_data_create(wave_data, 1, 1, cmplx(1.0_real64, 0.0_real64), cmplx(1.0_real64, 0.0_real64))
        call zone_create(zone, settings, background, wave_data, "test_path", 1, ierr)
        
        ! Set zone dimensions
        nwaves = 2
        ncomps = 6
        dim = 50
        zone%Nwaves = nwaves
        zone%Ncomps = ncomps
        zone%dim = dim
        
        ! Allocate basis fields
        call zone_allocate_basis_fields(zone, ierr)
        
        test_passed = (ierr == 0)
        if (test_passed) then
            test_passed = allocated(zone%basis)
            test_passed = test_passed .and. allocated(zone%EB)
            test_passed = test_passed .and. allocated(zone%S)
            
            if (allocated(zone%basis)) then
                test_passed = test_passed .and. (size(zone%basis) == 2*ncomps*nwaves*dim)
            end if
            if (allocated(zone%EB)) then
                test_passed = test_passed .and. (size(zone%EB) == 2*ncomps*dim)
            end if
            if (allocated(zone%S)) then
                test_passed = test_passed .and. (size(zone%S) == nwaves)
            end if
        end if
        
        ! Clean up
        call zone_destroy(zone)
        call wave_data_destroy(wave_data)
        call background_destroy(background)
        call settings_destroy(settings)
        
        call end_test(test_passed)
    end subroutine test_zone_basis_allocation
    
    !---------------------------------------------------------------------------
    ! Test 6: Zone indexing functions
    !---------------------------------------------------------------------------
    subroutine test_zone_indexing()
        type(zone_extended_t) :: zone
        integer :: idx1, idx2, expected_idx
        
        call start_test("Zone indexing functions")
        
        ! Set up test zone dimensions
        zone%Ncomps = 6
        zone%Nwaves = 2
        zone%dim = 10
        
        ! Test basis indexing function
        ! ib(node, sol, comp, part) where node=0:dim-1, sol=0:Nwaves-1, comp=0:Ncomps-1, part=0:1
        idx1 = zone_ib(zone, 0, 0, 0, 0)  ! First element
        expected_idx = 1  ! Fortran 1-based
        test_passed = (idx1 == expected_idx)
        
        idx2 = zone_ib(zone, 1, 1, 2, 1)  ! Complex element
        expected_idx = 2 + 2*(2 + 6*(1 + 2*1))  ! Following the C++ formula but 1-based
        test_passed = test_passed .and. (idx2 == expected_idx)
        
        ! Test EB field indexing function
        idx1 = zone_iEB(zone, 0, 0, 0)  ! First element
        expected_idx = 1  ! Fortran 1-based
        test_passed = test_passed .and. (idx1 == expected_idx)
        
        idx2 = zone_iEB(zone, 2, 3, 1)  ! Complex element
        expected_idx = 2 + 2*(3 + 6*2)  ! Following the C++ formula but 1-based
        test_passed = test_passed .and. (idx2 == expected_idx)
        
        call end_test(test_passed)
    end subroutine test_zone_indexing
    
    !---------------------------------------------------------------------------
    ! Test 7: Zone file I/O operations
    !---------------------------------------------------------------------------
    subroutine test_zone_file_operations()
        type(zone_extended_t) :: zone
        type(settings_t), target :: settings
        type(background_t), target :: background
        type(wave_data_t), target :: wave_data
        character(len=256) :: test_filename
        integer :: ierr
        
        call start_test("Zone file I/O operations")
        
        ! Create test zone
        call settings_create(settings, "dummy", ierr)
        call background_create(background, ierr)
        call wave_data_create(wave_data, 1, 1, cmplx(1.0_real64, 0.0_real64), cmplx(1.0_real64, 0.0_real64))
        call zone_create(zone, settings, background, wave_data, "test_path", 1, ierr)
        
        ! Test zone filename generation
        call zone_get_filename(1, test_filename, ierr)
        test_passed = (ierr == 0)
        test_passed = test_passed .and. (len_trim(test_filename) > 0)
        test_passed = test_passed .and. (index(test_filename, "zone_1") > 0)
        
        ! Test zone type determination (stub implementation)
        call zone_determine_type("dummy_file.in", zone%medium, ierr)
        test_passed = test_passed .and. (ierr == 0)  ! Should succeed with stub
        
        ! Clean up
        call zone_destroy(zone)
        call wave_data_destroy(wave_data)
        call background_destroy(background)
        call settings_destroy(settings)
        
        call end_test(test_passed)
    end subroutine test_zone_file_operations
    
    !---------------------------------------------------------------------------
    ! Test 8: Zone validation and error handling
    !---------------------------------------------------------------------------
    subroutine test_zone_validation()
        type(zone_extended_t) :: zone
        integer :: ierr
        
        call start_test("Zone validation and error handling")
        
        ! Test zone validation with uninitialized zone
        call zone_validate(zone, ierr)
        test_passed = (ierr /= 0)  ! Should fail for uninitialized zone
        
        ! Test zone print without initialization (should not crash)
        call zone_print(zone)
        test_passed = test_passed .and. .true.  ! If we get here, it didn't crash
        
        call end_test(test_passed)
    end subroutine test_zone_validation
    
    !---------------------------------------------------------------------------
    ! Test utilities
    !---------------------------------------------------------------------------
    subroutine start_test(name)
        character(len=*), intent(in) :: name
        total_tests = total_tests + 1
        write(*,'(A,A)', advance='no') "Testing ", name
        write(*,'(A)', advance='no') " ... "
    end subroutine start_test
    
    subroutine end_test(passed)
        logical, intent(in) :: passed
        if (passed) then
            print *, "PASSED"
            passed_tests = passed_tests + 1
        else
            print *, "FAILED"
            failed_tests = failed_tests + 1
        end if
    end subroutine end_test

end program test_kilca_zone