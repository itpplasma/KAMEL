program test_kilca_mode
    use iso_fortran_env, only: real64, int32
    use kilca_types_m
    use kilca_mode_m
    use kilca_settings_m
    use kilca_background_m
    implicit none
    
    ! Test counters
    integer :: total_tests, passed_tests, failed_tests
    logical :: test_passed
    
    ! Test data
    type(wave_data_t) :: wave
    type(mode_data_t) :: mode
    type(zone_t), allocatable :: test_zone
    type(settings_t), pointer :: settings => null()
    type(background_t), pointer :: background => null()
    real(real64) :: tol
    integer :: ierr, m, n
    complex(real64) :: omega_lab, omega_mov
    
    ! Initialize test counters
    total_tests = 0
    passed_tests = 0
    failed_tests = 0
    tol = epsilon(1.0_real64) * 100.0_real64
    
    print *, "===================================================="
    print *, "KiLCA Mode Module - Unit Tests"
    print *, "===================================================="
    print *, ""
    
    ! Test 1: Wave data creation and destruction
    call start_test("Wave data creation and destruction")
    
    m = 1
    n = 10
    omega_lab = cmplx(1.0e4_real64, 100.0_real64, real64)
    omega_mov = cmplx(0.9e4_real64, 90.0_real64, real64)
    
    call wave_data_create(wave, m, n, omega_lab, omega_mov)
    test_passed = (wave%m == m) .and. (wave%n == n) .and. &
                  (abs(wave%omega_lab - omega_lab) < tol) .and. &
                  (abs(wave%omega_mov - omega_mov) < tol)
    
    call wave_data_destroy(wave)
    call end_test(test_passed)
    
    ! Test 2: Zone creation with boundary conditions
    call start_test("Zone creation with boundary conditions")
    
    allocate(test_zone)
    call zone_create(test_zone, 0.0_real64, 1.0_real64, BOUNDARY_CENTER, BOUNDARY_INTERFACE, &
                     PLASMA_MODEL_VACUUM, ierr)
    
    test_passed = (ierr == 0) .and. &
                  (abs(test_zone%r1 - 0.0_real64) < tol) .and. &
                  (abs(test_zone%r2 - 1.0_real64) < tol) .and. &
                  (test_zone%bc1 == BOUNDARY_CENTER) .and. &
                  (test_zone%bc2 == BOUNDARY_INTERFACE) .and. &
                  (test_zone%medium == PLASMA_MODEL_VACUUM)
    
    call zone_destroy(test_zone)
    deallocate(test_zone)
    call end_test(test_passed)
    
    ! Test 3: Mode data initialization
    call start_test("Mode data initialization")
    
    ! Initialize settings and background first
    call settings_create(settings, "./", ierr)
    allocate(background)
    call background_create(background, settings, ierr)
    
    ! Create wave data
    call wave_data_create(wave, m, n, omega_lab, omega_mov)
    
    ! Initialize mode
    call mode_data_create(mode, wave, settings, background, ierr)
    
    test_passed = (ierr == 0) .and. &
                  associated(mode%sd) .and. &
                  associated(mode%bp) .and. &
                  associated(mode%wd)
    
    call mode_data_destroy(mode)
    call wave_data_destroy(wave)
    call background_destroy(background)
    deallocate(background)
    call settings_destroy(settings, ierr)
    call end_test(test_passed)
    
    ! Test 4: Zone allocation in mode
    call start_test("Zone allocation in mode structure")
    
    call settings_create(settings, "./", ierr)
    allocate(background)
    call background_create(background, settings, ierr)
    call wave_data_create(wave, m, n, omega_lab, omega_mov)
    call mode_data_create(mode, wave, settings, background, ierr)
    
    ! Allocate zones
    mode%Nzones = 3
    call mode_allocate_zones(mode, ierr)
    
    test_passed = (ierr == 0) .and. &
                  allocated(mode%zones) .and. &
                  (size(mode%zones) == 3)
    
    call mode_data_destroy(mode)
    call wave_data_destroy(wave)
    call background_destroy(background)
    deallocate(background)
    call settings_destroy(settings, ierr)
    call end_test(test_passed)
    
    ! Test 5: Field indexing functions
    call start_test("Field indexing functions")
    
    ! Test basis indexing
    ! ib(node, sol, comp, part) for complex basis array
    test_passed = (ib_index(0, 0, 0, 0) == 1) .and. &  ! First real part
                  (ib_index(0, 0, 0, 1) == 2) .and. &  ! First imag part
                  (ib_index(1, 0, 0, 0) > 2)           ! Next node
    
    ! Test EB field indexing
    test_passed = test_passed .and. &
                  (iEB_index(0, 0, 0) == 1) .and. &
                  (iEB_index(0, 0, 1) == 2)
    
    call end_test(test_passed)
    
    ! Test 6: Zone boundary type strings
    call start_test("Zone boundary type string conversion")
    
    test_passed = (trim(get_boundary_string(BOUNDARY_CENTER)) == "center") .and. &
                  (trim(get_boundary_string(BOUNDARY_INFINITY)) == "infinity") .and. &
                  (trim(get_boundary_string(BOUNDARY_IDEALWALL)) == "idealwall") .and. &
                  (trim(get_boundary_string(BOUNDARY_INTERFACE)) == "interface") .and. &
                  (trim(get_boundary_string(BOUNDARY_ANTENNA)) == "antenna")
    
    call end_test(test_passed)
    
    ! Test 7: Plasma model type strings
    call start_test("Plasma model type string conversion")
    
    test_passed = (trim(get_medium_string(PLASMA_MODEL_VACUUM)) == "vacuum") .and. &
                  (trim(get_medium_string(PLASMA_MODEL_MEDIUM)) == "medium") .and. &
                  (trim(get_medium_string(PLASMA_MODEL_IMHD)) == "imhd") .and. &
                  (trim(get_medium_string(PLASMA_MODEL_RMHD)) == "rmhd") .and. &
                  (trim(get_medium_string(4)) == "flre")  ! Use explicit 4 for FLRE to match C++
    
    call end_test(test_passed)
    
    ! Test 8: Mode path generation
    call start_test("Mode path generation")
    
    call settings_create(settings, "./", ierr)
    allocate(background)
    call background_create(background, settings, ierr)
    call wave_data_create(wave, 1, 10, omega_lab, omega_mov)
    call mode_data_create(mode, wave, settings, background, ierr)
    
    call mode_set_paths(mode, "./test_project/", ierr)
    
    test_passed = (ierr == 0) .and. &
                  (len_trim(mode%path2linear) > 0) .and. &
                  (len_trim(mode%path2dispersion) > 0) .and. &
                  (len_trim(mode%path2poincare) > 0)
    
    if (.not. test_passed) then
        print *, "  Linear path: ", trim(mode%path2linear)
        print *, "  Dispersion path: ", trim(mode%path2dispersion)
        print *, "  Poincare path: ", trim(mode%path2poincare)
    end if
    
    call mode_data_destroy(mode)
    call wave_data_destroy(wave)
    call background_destroy(background)
    deallocate(background)
    call settings_destroy(settings, ierr)
    call end_test(test_passed)
    
    ! Final summary
    print *, ""
    print *, "===================================================="
    print *, "TEST SUMMARY"
    print *, "===================================================="
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

end program test_kilca_mode