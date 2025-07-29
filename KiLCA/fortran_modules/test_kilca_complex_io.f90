program test_kilca_complex_io
    use iso_fortran_env, only: real64, int32, output_unit, error_unit
    use kilca_complex_m
    implicit none
    
    logical :: test_passed
    complex(real64) :: z1, z2, z3
    character(len=200) :: str, filename
    integer :: io_unit, ios
    real(real64) :: re, im
    
    test_passed = .true.
    
    print *, "Testing kilca_complex_m I/O functionality..."
    print *, ""
    
    ! Test 1: Basic string formatting
    print *, "Testing string formatting..."
    
    z1 = cmplx_make(3.14159_real64, -2.71828_real64)
    str = cmplx_to_string(z1)
    print *, "Default format: ", trim(str)
    
    ! Test scientific format
    str = cmplx_to_string(z1, '(ES15.8,SP,ES15.8,"i")')
    print *, "Scientific format: ", trim(str)
    
    ! Test engineering format
    str = cmplx_to_string(z1, '(EN15.6,SP,EN15.6,"i")')
    print *, "Engineering format: ", trim(str)
    
    ! Test compact format
    str = cmplx_to_string(z1, '(F8.4,SP,F8.4,"i")')
    print *, "Compact format: ", trim(str)
    
    ! Test 2: Print functionality
    print *, ""
    print *, "Testing print functionality..."
    
    call cmplx_print(z1)
    call cmplx_print(z1, "z1")
    
    z2 = cmplx_I
    call cmplx_print(z2, "imaginary unit")
    
    z3 = cmplx_polar(2.0_real64, 0.785398_real64)  ! 45 degrees
    call cmplx_print(z3, "polar(2, pi/4)")
    
    ! Test 3: File I/O
    print *, ""
    print *, "Testing file I/O..."
    
    filename = "test_complex_io.dat"
    
    ! Write complex numbers to file
    open(newunit=io_unit, file=filename, status='replace', action='write', iostat=ios)
    if (ios /= 0) then
        print *, "FAIL: Could not open file for writing"
        test_passed = .false.
    else
        ! Write in different formats
        write(io_unit, '(A)') "# Complex number I/O test file"
        
        ! Format 1: Real and imaginary parts separated
        z1 = cmplx_make(1.23456789_real64, -9.87654321_real64)
        write(io_unit, '(A,ES20.12,1X,ES20.12)') "Format1: ", real(z1, real64), aimag(z1)
        
        ! Format 2: Parentheses format (like C++)
        write(io_unit, '(A,"(",ES20.12,", ",ES20.12,")")') "Format2: ", real(z1, real64), aimag(z1)
        
        ! Format 3: Using cmplx_to_string
        write(io_unit, '(A,A)') "Format3: ", trim(cmplx_to_string(z1))
        
        ! Write array of complex numbers
        write(io_unit, '(A)') "# Array of complex numbers"
        do ios = 1, 5
            z2 = cmplx_polar(real(ios, real64), real(ios, real64) * 0.1_real64)
            write(io_unit, '(I3,2X,ES20.12,2X,ES20.12)') ios, real(z2, real64), aimag(z2)
        end do
        
        close(io_unit)
        print *, "PASS: Successfully wrote complex numbers to file"
    end if
    
    ! Read complex numbers from file
    open(newunit=io_unit, file=filename, status='old', action='read', iostat=ios)
    if (ios /= 0) then
        print *, "FAIL: Could not open file for reading"
        test_passed = .false.
    else
        ! Skip header
        read(io_unit, '(A)') str
        
        ! Read format 1
        read(io_unit, '(A9,ES20.12,1X,ES20.12)') str, re, im
        z2 = cmplx_make(re, im)
        if (abs(z2 - z1) > epsilon(1.0_real64) * 10.0_real64) then
            print *, "FAIL: Format 1 read mismatch"
            test_passed = .false.
        else
            print *, "PASS: Format 1 read correctly"
        end if
        
        ! Read format 2 (with parentheses)
        read(io_unit, '(A9,1X,ES20.12,2X,ES20.12)') str, re, im
        z2 = cmplx_make(re, im)
        if (abs(z2 - z1) > epsilon(1.0_real64) * 10.0_real64) then
            print *, "FAIL: Format 2 read mismatch"
            test_passed = .false.
        else
            print *, "PASS: Format 2 read correctly"
        end if
        
        close(io_unit)
    end if
    
    ! Clean up test file
    open(newunit=io_unit, file=filename, status='old', iostat=ios)
    if (ios == 0) then
        close(io_unit, status='delete')
    end if
    
    ! Test 4: Error handling
    print *, ""
    print *, "Testing error handling..."
    
    ! Test with special values
    z1 = cmplx_make(huge(1.0_real64), 0.0_real64)
    call cmplx_print(z1, "Infinity real part")
    
    z2 = cmplx_make(0.0_real64, -huge(1.0_real64))
    call cmplx_print(z2, "Negative infinity imag part")
    
    z3 = cmplx_make(tiny(1.0_real64), tiny(1.0_real64))
    call cmplx_print(z3, "Tiny values")
    
    print *, ""
    if (test_passed) then
        print *, "All I/O tests PASSED!"
    else
        print *, "Some I/O tests FAILED!"
        stop 1
    end if
    
end program test_kilca_complex_io