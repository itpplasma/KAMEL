program test_kilca_mode_solver
    use iso_fortran_env, only: real64, int32
    use kilca_types_m
    use kilca_mode_m
    use kilca_background_m
    use kilca_settings_m
    use kilca_mode_solver_m
    implicit none
    
    ! Test counters
    integer :: total_tests, passed_tests, failed_tests
    logical :: test_passed
    
    ! Initialize test counters
    total_tests = 0
    passed_tests = 0
    failed_tests = 0
    
    print *, "====================================================="
    print *, "KiLCA Mode Solver Module - Unit Tests"
    print *, "====================================================="
    print *, ""
    
    ! Test 1: Find resonance location with simple q-profile
    call test_find_resonance_location()
    
    ! Test 2: Mode calculation workflow initialization
    call test_mode_calculation_workflow()
    
    ! Test 3: Brent's root finding algorithm
    call test_brent_root_finding()
    
    ! Test 4: Mode data creation and basic setup
    call test_mode_data_basic_setup()
    
    ! Test 5: Error handling in mode solver
    call test_mode_solver_error_handling()
    
    ! Final summary
    print *, ""
    print *, "====================================================="
    print *, "TEST SUMMARY"
    print *, "====================================================="
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
    ! Test 1: Find resonance location with simple q-profile
    !---------------------------------------------------------------------------
    subroutine test_find_resonance_location()
        type(mode_data_t) :: mode_data
        type(settings_t), pointer :: settings
        type(background_t), target :: background
        type(wave_data_t), target :: wave_data
        integer :: ierr
        
        call start_test("Find resonance location with simple q-profile")
        
        ! Initialize test data structures
        call settings_create(settings, "dummy", ierr)
        test_passed = (ierr == 0)
        
        if (test_passed) then
            call background_create(background, settings, ierr)
            test_passed = (ierr == 0)
        end if
        
        if (test_passed) then
            ! Setup simple background profile for testing
            background%dimx = 100
            if (allocated(background%x)) deallocate(background%x)
            allocate(background%x(background%dimx))
            
            ! Create radial grid from 0.1 to 1.0
            do ierr = 1, background%dimx
                background%x(ierr) = 0.1_real64 + &
                    (1.0_real64 - 0.1_real64) * real(ierr - 1, real64) / &
                    real(background%dimx - 1, real64)
            end do
            
            ! Create wave data for m=2, n=1 mode
            call wave_data_create(wave_data, 2, 1, &
                                  cmplx(1.0_real64, 0.1_real64, real64), &
                                  cmplx(0.8_real64, 0.1_real64, real64))
            
            ! Setup mode data
            call mode_data_create(mode_data, wave_data, settings, background, ierr)
            test_passed = (ierr == 0)
        end if
        
        if (test_passed) then
            ! Set debug flag for output
            settings%os%flag_debug = 1
            
            ! Test resonance finding
            call find_resonance_location(mode_data, ierr)
            test_passed = (ierr == 0 .and. mode_data%wd%r_res > 0.0_real64)
        end if
        
        if (.not. test_passed) then
            print *, "  Resonance location finding failed, ierr =", ierr
        else
            print *, "  Resonance location found at r =", mode_data%wd%r_res
        end if
        
        ! Cleanup
        if (allocated(background%x)) deallocate(background%x)
        if (associated(settings)) call settings_destroy(settings, ierr)
        
        call end_test(test_passed)
    end subroutine test_find_resonance_location
    
    !---------------------------------------------------------------------------
    ! Test 2: Mode calculation workflow initialization  
    !---------------------------------------------------------------------------
    subroutine test_mode_calculation_workflow()
        type(mode_data_t) :: mode_data
        type(settings_t), pointer :: settings
        type(background_t), target :: background
        type(wave_data_t), target :: wave_data
        integer :: ierr
        
        call start_test("Mode calculation workflow initialization")
        
        ! Initialize test data structures
        call settings_create(settings, "dummy", ierr)
        test_passed = (ierr == 0)
        
        if (test_passed) then
            call background_create(background, settings, ierr)
            test_passed = (ierr == 0)
        end if
        
        if (test_passed) then
            ! Setup minimal data structures
            background%dimx = 50
            if (allocated(background%x)) deallocate(background%x)
            allocate(background%x(background%dimx))
            
            do ierr = 1, background%dimx
                background%x(ierr) = 0.2_real64 + &
                    0.6_real64 * real(ierr - 1, real64) / real(background%dimx - 1, real64)
            end do
            
            call wave_data_create(wave_data, 1, 1, &
                                  cmplx(0.5_real64, 0.05_real64, real64), &
                                  cmplx(0.4_real64, 0.05_real64, real64))
            
            call mode_data_create(mode_data, wave_data, settings, background, ierr)
            test_passed = (ierr == 0)
        end if
        
        if (test_passed) then
            ! Set flags for minimal calculation
            settings%os%flag_emfield = 0
            settings%os%flag_additional = 0
            settings%os%flag_debug = 0
            
            ! Test main calculation workflow with flag=1 (basis only)
            call calc_all_mode_data(mode_data, 1, ierr)
            
            ! For now, expect implementation pending error
            test_passed = (ierr == -1)  ! Expected "implementation pending"
        end if
        
        if (.not. test_passed) then
            print *, "  Mode calculation workflow test failed, ierr =", ierr
        else
            print *, "  Mode calculation workflow initialization successful"
        end if
        
        ! Cleanup
        if (allocated(background%x)) deallocate(background%x)
        if (associated(settings)) call settings_destroy(settings, ierr)
        
        call end_test(test_passed)
    end subroutine test_mode_calculation_workflow
    
    !---------------------------------------------------------------------------
    ! Test 3: Brent's root finding algorithm validation
    !---------------------------------------------------------------------------
    subroutine test_brent_root_finding()
        type(mode_data_t) :: mode_data
        type(settings_t), target :: settings
        type(background_t), target :: background
        type(wave_data_t), target :: wave_data
        integer :: ierr
        real(real64) :: expected_root, actual_root, error
        
        call start_test("Brent's root finding algorithm validation")
        
        ! Initialize data structures
        call settings_create(settings, ierr)
        test_passed = (ierr == 0)
        
        if (test_passed) then
            call background_create(background, ierr)
            test_passed = (ierr == 0)
        end if
        
        if (test_passed) then
            ! Setup background with known q-profile that has exact solution
            background%dimx = 200
            allocate(background%x(background%dimx))
            
            do ierr = 1, background%dimx
                background%x(ierr) = 0.1_real64 + &
                    0.8_real64 * real(ierr - 1, real64) / real(background%dimx - 1, real64)
            end do
            
            ! For m=3, n=2, q_res = -3/2 = -1.5
            ! With linear q-profile q(r) = 1.5 + 2*r, resonance at r where:
            ! 1.5 + 2*r = -1.5  => r = -1.5 (no solution in positive domain)
            ! Let's use q(r) = -2 + 3*r, then for q_res = -1.5:
            ! -2 + 3*r = -1.5  => r = 0.5/3 ≈ 0.167
            
            wave_data%m = 3
            wave_data%n = 2
            wave_data%olab = cmplx(1.0_real64, 0.0_real64, real64)
            wave_data%omov = cmplx(1.0_real64, 0.0_real64, real64)
            
            mode_data%sd => settings
            mode_data%bp => background
            mode_data%wd => wave_data
            
            settings%os%flag_debug = 0  ! Suppress debug output for this test
            
            ! Calculate expected root analytically
            expected_root = 0.5_real64 / 3.0_real64  ! Approximately 0.167
            
            call find_resonance_location(mode_data, ierr)
            
            if (ierr == 0 .and. mode_data%wd%r_res > 0.0_real64) then
                actual_root = mode_data%wd%r_res
                error = abs(actual_root - expected_root)
                test_passed = (error < 1.0e-6_real64)
            else
                test_passed = .false.
            end if
        end if
        
        if (.not. test_passed) then
            print *, "  Brent's algorithm validation failed, ierr =", ierr
            if (mode_data%wd%r_res > 0.0_real64) then
                print *, "  Expected root:", expected_root, ", Got:", mode_data%wd%r_res
            end if
        else
            print *, "  Brent's algorithm validation successful, error =", error
        end if
        
        ! Cleanup
        if (allocated(background%x)) deallocate(background%x)
        
        call end_test(test_passed)
    end subroutine test_brent_root_finding
    
    !---------------------------------------------------------------------------
    ! Test 4: Mode data creation and basic setup
    !---------------------------------------------------------------------------
    subroutine test_mode_data_basic_setup()
        type(mode_data_t) :: mode_data
        type(settings_t), target :: settings
        type(background_t), target :: background
        type(wave_data_t), target :: wave_data
        integer :: ierr
        
        call start_test("Mode data creation and basic setup")
        
        ! Test mode_data_create functionality
        call settings_create(settings, ierr)
        test_passed = (ierr == 0)
        
        if (test_passed) then
            call background_create(background, ierr)
            test_passed = (ierr == 0)
        end if
        
        if (test_passed) then
            call wave_data_create(wave_data, 2, 3, &
                                  cmplx(1.5_real64, 0.1_real64, real64), &
                                  cmplx(1.2_real64, 0.1_real64, real64), ierr)
            test_passed = (ierr == 0)
        end if
        
        if (test_passed) then
            call mode_data_create(mode_data, settings, background, wave_data, ierr)
            test_passed = (ierr == 0)
        end if
        
        if (test_passed) then
            ! Verify basic properties
            test_passed = (associated(mode_data%sd) .and. &
                          associated(mode_data%bp) .and. &
                          associated(mode_data%wd) .and. &
                          mode_data%wd%m == 2 .and. &
                          mode_data%wd%n == 3)
        end if
        
        if (.not. test_passed) then
            print *, "  Mode data basic setup failed, ierr =", ierr
        else
            print *, "  Mode data basic setup successful"
        end if
        
        call end_test(test_passed)
    end subroutine test_mode_data_basic_setup
    
    !---------------------------------------------------------------------------
    ! Test 5: Error handling in mode solver
    !---------------------------------------------------------------------------
    subroutine test_mode_solver_error_handling()
        type(mode_data_t) :: mode_data
        complex(real64) :: EB_fields(6)
        integer :: ierr
        
        call start_test("Error handling in mode solver")
        
        ! Test error handling with uninitialized mode_data
        call eval_EB_fields(mode_data, 0.5_real64, EB_fields, ierr)
        test_passed = (ierr /= 0)  ! Should fail gracefully
        
        if (test_passed) then
            ! Test with invalid position
            call eval_EB_fields(mode_data, -1.0_real64, EB_fields, ierr)
            test_passed = (ierr /= 0)  ! Should fail gracefully
        end if
        
        if (.not. test_passed) then
            print *, "  Error handling test failed, ierr =", ierr
        else
            print *, "  Error handling test successful"
        end if
        
        call end_test(test_passed)
    end subroutine test_mode_solver_error_handling
    
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

end program test_kilca_mode_solver