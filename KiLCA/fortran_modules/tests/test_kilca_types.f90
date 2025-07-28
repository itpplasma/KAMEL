program test_kilca_types
    use iso_fortran_env, only: real32, real64, int32, int64, int8
    use iso_c_binding
    use kilca_types_m
    implicit none
    
    ! Test variables
    integer :: test_status = 0
    
    ! Test 1: Verify basic type definitions exist
    call test_basic_types()
    
    ! Test 2: Verify numeric kinds
    call test_numeric_kinds()
    
    ! Test 3: Verify string length parameters
    call test_string_parameters()
    
    ! Test 4: Verify error codes
    call test_error_codes()
    
    ! Test 5: Verify type compatibility with C
    call test_c_compatibility()
    
    ! Test 6: Verify complex number types
    call test_complex_types()
    
    ! Test 7: Verify physical constants
    call test_physical_constants()
    
    ! Test 8: Verify mathematical constants
    call test_mathematical_constants()
    
    ! Test 9: Verify array dimension constants
    call test_array_constants()
    
    ! Test 10: Verify debug and control flags
    call test_debug_flags()
    
    if (test_status == 0) then
        print *, "All tests PASSED!"
    else
        print *, "Some tests FAILED. Status:", test_status
        stop 1
    end if
    
contains

    subroutine test_basic_types()
        ! Test that uchar and schar equivalents exist
        type(uchar_t) :: uc
        type(schar_t) :: sc
        
        ! Test assignment and range
        uc = uchar_from_int(255)
        sc = schar_from_int(-128)
        
        if (int_from_uchar(uc) /= 255) then
            print *, "FAIL: uchar_t assignment"
            test_status = test_status + 1
        end if
        
        if (int_from_schar(sc) /= -128) then
            print *, "FAIL: schar_t assignment"
            test_status = test_status + 1
        end if
        
        print *, "test_basic_types completed"
    end subroutine test_basic_types
    
    subroutine test_numeric_kinds()
        ! Test that all numeric kinds are properly defined
        real(sp) :: real_single
        real(dp) :: real_double
        real(qp) :: real_quad
        integer(i32) :: int_32
        integer(i64) :: int_64
        
        ! Verify precision
        if (precision(real_single) < 6) then
            print *, "FAIL: Single precision too low"
            test_status = test_status + 1
        end if
        
        if (precision(real_double) < 15) then
            print *, "FAIL: Double precision too low"
            test_status = test_status + 1
        end if
        
        ! Verify integer ranges
        if (huge(int_32) < 2147483647) then
            print *, "FAIL: 32-bit integer range"
            test_status = test_status + 1
        end if
        
        print *, "test_numeric_kinds completed"
    end subroutine test_numeric_kinds
    
    subroutine test_string_parameters()
        ! Test string length parameters
        character(len=MAX_PATH_LEN) :: path
        character(len=MAX_NAME_LEN) :: name
        character(len=MAX_LINE_LEN) :: line
        
        ! Verify lengths are reasonable
        if (MAX_PATH_LEN < 256) then
            print *, "FAIL: MAX_PATH_LEN too small"
            test_status = test_status + 1
        end if
        
        if (MAX_NAME_LEN < 64) then
            print *, "FAIL: MAX_NAME_LEN too small"
            test_status = test_status + 1
        end if
        
        if (MAX_LINE_LEN < 132) then
            print *, "FAIL: MAX_LINE_LEN too small"
            test_status = test_status + 1
        end if
        
        print *, "test_string_parameters completed"
    end subroutine test_string_parameters
    
    subroutine test_error_codes()
        ! Test error code definitions
        if (KILCA_SUCCESS /= 0) then
            print *, "FAIL: KILCA_SUCCESS should be 0"
            test_status = test_status + 1
        end if
        
        if (KILCA_ERROR_MEMORY == KILCA_SUCCESS) then
            print *, "FAIL: Error codes not unique"
            test_status = test_status + 1
        end if
        
        if (KILCA_ERROR_FILE == KILCA_ERROR_MEMORY) then
            print *, "FAIL: Error codes not unique"
            test_status = test_status + 1
        end if
        
        print *, "test_error_codes completed"
    end subroutine test_error_codes
    
    subroutine test_c_compatibility()
        ! Test C interoperability
        real(c_double) :: c_real
        complex(c_double_complex) :: c_cmplx
        integer(c_int) :: c_integer
        
        ! Verify compatibility with our types
        if (dp /= c_double) then
            print *, "FAIL: dp not compatible with c_double"
            test_status = test_status + 1
        end if
        
        if (i32 /= c_int) then
            print *, "FAIL: i32 not compatible with c_int"
            test_status = test_status + 1
        end if
        
        print *, "test_c_compatibility completed"
    end subroutine test_c_compatibility
    
    subroutine test_complex_types()
        ! Test complex number types
        complex(sp) :: c_single
        complex(dp) :: c_double
        complex(qp) :: c_quad
        
        ! Test complex constants
        complex(dp) :: test_zero, test_one, test_i
        
        test_zero = cmplx_zero
        test_one = cmplx_one
        test_i = cmplx_i
        
        if (abs(test_zero) > epsilon(1.0_dp)) then
            print *, "FAIL: Complex zero not zero"
            test_status = test_status + 1
        end if
        
        if (abs(test_one - (1.0_dp, 0.0_dp)) > epsilon(1.0_dp)) then
            print *, "FAIL: Complex one incorrect"
            test_status = test_status + 1
        end if
        
        if (abs(test_i - (0.0_dp, 1.0_dp)) > epsilon(1.0_dp)) then
            print *, "FAIL: Complex i incorrect"
            test_status = test_status + 1
        end if
        
        print *, "test_complex_types completed"
    end subroutine test_complex_types
    
    subroutine test_physical_constants()
        ! Test physical constants
        real(dp) :: test_val
        
        ! Test speed of light
        test_val = c_light
        if (abs(test_val - 29979245800.0_dp) > 1.0_dp) then
            print *, "FAIL: Speed of light incorrect"
            test_status = test_status + 1
        end if
        
        ! Test Boltzmann constant
        test_val = k_boltz
        if (abs(test_val - 1.60216428e-12_dp) > 1.0e-20_dp) then
            print *, "FAIL: Boltzmann constant incorrect"
            test_status = test_status + 1
        end if
        
        ! Test masses
        if (abs(m_proton - 1.67262158e-24_dp) > 1.0e-32_dp) then
            print *, "FAIL: Proton mass incorrect"
            test_status = test_status + 1
        end if
        
        if (abs(m_electron - 9.10938185917485e-28_dp) > 1.0e-40_dp) then
            print *, "FAIL: Electron mass incorrect"
            test_status = test_status + 1
        end if
        
        ! Test charge
        if (abs(e_charge - 4.8032e-10_dp) > 1.0e-14_dp) then
            print *, "FAIL: Elementary charge incorrect"
            test_status = test_status + 1
        end if
        
        ! Test adiabatic constant
        if (abs(gamma_adiabatic - 5.0_dp/3.0_dp) > epsilon(1.0_dp)) then
            print *, "FAIL: Adiabatic constant incorrect"
            test_status = test_status + 1
        end if
        
        print *, "test_physical_constants completed"
    end subroutine test_physical_constants
    
    subroutine test_mathematical_constants()
        ! Test mathematical constants
        real(dp) :: test_val
        
        ! Test pi
        test_val = pi
        if (abs(test_val - 3.141592653589793238462643383279502884197_dp) > 1.0e-15_dp) then
            print *, "FAIL: Pi value incorrect"
            test_status = test_status + 1
        end if
        
        ! Test Euler's constant
        test_val = euler
        if (abs(test_val - 0.5772156649015328606065120900824024310422_dp) > 1.0e-15_dp) then
            print *, "FAIL: Euler constant incorrect"
            test_status = test_status + 1
        end if
        
        ! Test sqrt(2*pi)
        test_val = sqrt_2pi
        if (abs(test_val - sqrt(2.0_dp * pi)) > epsilon(1.0_dp)) then
            print *, "FAIL: sqrt(2*pi) incorrect"
            test_status = test_status + 1
        end if
        
        print *, "test_mathematical_constants completed"
    end subroutine test_mathematical_constants
    
    subroutine test_array_constants()
        ! Test array dimension and indexing constants
        if (MAX_MODES < 1) then
            print *, "FAIL: MAX_MODES too small"
            test_status = test_status + 1
        end if
        
        if (MAX_ZONES < 1) then
            print *, "FAIL: MAX_ZONES too small"
            test_status = test_status + 1
        end if
        
        if (MAX_GRID_POINTS < 100) then
            print *, "FAIL: MAX_GRID_POINTS too small"
            test_status = test_status + 1
        end if
        
        print *, "test_array_constants completed"
    end subroutine test_array_constants
    
    subroutine test_debug_flags()
        ! Test debug and control flags
        if (DEBUG_FLAG /= 0 .and. DEBUG_FLAG /= 1) then
            print *, "FAIL: DEBUG_FLAG not binary"
            test_status = test_status + 1
        end if
        
        if (CALC_EIGEN_DECOMPOSITION /= 0 .and. CALC_EIGEN_DECOMPOSITION /= 1) then
            print *, "FAIL: CALC_EIGEN_DECOMPOSITION not binary"
            test_status = test_status + 1
        end if
        
        if (SORT_DISPERSION_PROFILES /= 0 .and. SORT_DISPERSION_PROFILES /= 1) then
            print *, "FAIL: SORT_DISPERSION_PROFILES not binary"
            test_status = test_status + 1
        end if
        
        if (USE_JACOBIAN_IN_ODE_SOLVER /= 0 .and. USE_JACOBIAN_IN_ODE_SOLVER /= 1) then
            print *, "FAIL: USE_JACOBIAN_IN_ODE_SOLVER not binary"
            test_status = test_status + 1
        end if
        
        if (USE_SPLINES_IN_RHS_EVALUATION /= 0 .and. USE_SPLINES_IN_RHS_EVALUATION /= 1) then
            print *, "FAIL: USE_SPLINES_IN_RHS_EVALUATION not binary"
            test_status = test_status + 1
        end if
        
        ! Test plasma model constants
        if (PLASMA_MODEL_VACUUM /= 0) then
            print *, "FAIL: PLASMA_MODEL_VACUUM should be 0"
            test_status = test_status + 1
        end if
        
        if (PLASMA_MODEL_MEDIUM /= 1) then
            print *, "FAIL: PLASMA_MODEL_MEDIUM should be 1"
            test_status = test_status + 1
        end if
        
        if (PLASMA_MODEL_IMHD /= 2) then
            print *, "FAIL: PLASMA_MODEL_IMHD should be 2"
            test_status = test_status + 1
        end if
        
        print *, "test_debug_flags completed"
    end subroutine test_debug_flags

end program test_kilca_types