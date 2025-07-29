program test_kilca_constants
    use iso_fortran_env, only: real64
    use kilca_constants_m
    implicit none
    
    real(real64) :: test_value
    complex(real64) :: test_complex
    logical :: test_passed
    
    test_passed = .true.
    
    ! Test mathematical constants
    print *, "Testing mathematical constants..."
    
    ! Test pi
    test_value = pi
    if (abs(test_value - 3.141592653589793238462643383279502884197_real64) > epsilon(1.0_real64)) then
        print *, "FAIL: pi value incorrect"
        test_passed = .false.
    else
        print *, "PASS: pi = ", pi
    end if
    
    ! Test Euler's constant
    test_value = eul
    if (abs(test_value - 0.5772156649015328606065120900824024310422_real64) > epsilon(1.0_real64)) then
        print *, "FAIL: Euler constant value incorrect"
        test_passed = .false.
    else
        print *, "PASS: eul = ", eul
    end if
    
    ! Test sqrt(2*pi)
    test_value = sqrt2pi
    if (abs(test_value - sqrt(2.0_real64 * pi)) > epsilon(1.0_real64)) then
        print *, "FAIL: sqrt2pi value incorrect"
        test_passed = .false.
    else
        print *, "PASS: sqrt2pi = ", sqrt2pi
    end if
    
    print *, ""
    print *, "Testing physical constants..."
    
    ! Test Boltzmann constant
    test_value = boltz
    if (abs(test_value - 1.60216428e-12_real64) > epsilon(test_value)) then
        print *, "FAIL: Boltzmann constant value incorrect"
        test_passed = .false.
    else
        print *, "PASS: boltz = ", boltz, " erg/eV"
    end if
    
    ! Test speed of light
    test_value = c_light
    if (abs(test_value - 29979245800.0_real64) > epsilon(test_value)) then
        print *, "FAIL: Speed of light value incorrect"
        test_passed = .false.
    else
        print *, "PASS: c_light = ", c_light, " cm/s"
    end if
    
    ! Test proton mass
    test_value = m_p
    if (abs(test_value - 1.67262158e-24_real64) > epsilon(test_value)) then
        print *, "FAIL: Proton mass value incorrect"
        test_passed = .false.
    else
        print *, "PASS: m_p = ", m_p, " g"
    end if
    
    ! Test electron mass
    test_value = m_e
    if (abs(test_value - 9.10938185917485e-28_real64) > epsilon(test_value)) then
        print *, "FAIL: Electron mass value incorrect"
        test_passed = .false.
    else
        print *, "PASS: m_e = ", m_e, " g"
    end if
    
    ! Test elementary charge
    test_value = e_charge
    if (abs(test_value - 4.8032e-10_real64) > epsilon(test_value)) then
        print *, "FAIL: Elementary charge value incorrect"
        test_passed = .false.
    else
        print *, "PASS: e_charge = ", e_charge, " esu"
    end if
    
    ! Test adiabatic constant
    test_value = gamma_adiabatic
    if (abs(test_value - 5.0_real64/3.0_real64) > epsilon(test_value)) then
        print *, "FAIL: Adiabatic constant value incorrect"
        test_passed = .false.
    else
        print *, "PASS: gamma_adiabatic = ", gamma_adiabatic
    end if
    
    print *, ""
    print *, "Testing complex constants..."
    
    ! Test complex zero
    test_complex = cmplx_zero
    if (abs(test_complex) > epsilon(1.0_real64)) then
        print *, "FAIL: Complex zero value incorrect"
        test_passed = .false.
    else
        print *, "PASS: cmplx_zero = ", cmplx_zero
    end if
    
    ! Test complex one
    test_complex = cmplx_one
    if (abs(test_complex - cmplx(1.0_real64, 0.0_real64, real64)) > epsilon(1.0_real64)) then
        print *, "FAIL: Complex one value incorrect"
        test_passed = .false.
    else
        print *, "PASS: cmplx_one = ", cmplx_one
    end if
    
    ! Test complex i
    test_complex = cmplx_i
    if (abs(test_complex - cmplx(0.0_real64, 1.0_real64, real64)) > epsilon(1.0_real64)) then
        print *, "FAIL: Complex i value incorrect"
        test_passed = .false.
    else
        print *, "PASS: cmplx_i = ", cmplx_i
    end if
    
    print *, ""
    print *, "Testing code settings constants..."
    
    ! Test code settings
    if (DEBUG_FLAG /= 0) then
        print *, "FAIL: DEBUG_FLAG value incorrect"
        test_passed = .false.
    else
        print *, "PASS: DEBUG_FLAG = ", DEBUG_FLAG
    end if
    
    if (USE_JACOBIAN_IN_ODE_SOLVER /= 1) then
        print *, "FAIL: USE_JACOBIAN_IN_ODE_SOLVER value incorrect"
        test_passed = .false.
    else
        print *, "PASS: USE_JACOBIAN_IN_ODE_SOLVER = ", USE_JACOBIAN_IN_ODE_SOLVER
    end if
    
    print *, ""
    print *, "Testing solver method constants..."
    
    ! Test solver methods
    if (NORMAL_EXACT /= 0 .or. EIG_EXP_EXACT /= 1 .or. EIG_EXP_PHYS /= 2) then
        print *, "FAIL: Solver method constants incorrect"
        test_passed = .false.
    else
        print *, "PASS: Solver method constants correct"
    end if
    
    print *, ""
    print *, "Testing plasma model constants..."
    
    ! Test plasma models
    if (PLASMA_MODEL_VACUUM /= 0 .or. PLASMA_MODEL_FLRE /= 4) then
        print *, "FAIL: Plasma model constants incorrect"
        test_passed = .false.
    else
        print *, "PASS: Plasma model constants correct"
    end if
    
    print *, ""
    print *, "Testing boundary condition constants..."
    
    ! Test boundary conditions
    if (BOUNDARY_CENTER /= 0 .or. BOUNDARY_ANTENNA /= 4) then
        print *, "FAIL: Boundary condition constants incorrect"
        test_passed = .false.
    else
        print *, "PASS: Boundary condition constants correct"
    end if
    
    print *, ""
    print *, "Testing string length constants..."
    
    ! Test string lengths
    if (MAX_PATH_LENGTH /= 1024 .or. MAX_NAME_LENGTH /= 64) then
        print *, "FAIL: String length constants incorrect"
        test_passed = .false.
    else
        print *, "PASS: String length constants correct"
    end if
    
    print *, ""
    print *, "Testing numerical precision constants..."
    
    ! Test numerical tolerances
    if (abs(DEFAULT_REL_TOL - 1.0e-12_real64) > epsilon(1.0_real64)) then
        print *, "FAIL: DEFAULT_REL_TOL incorrect"
        test_passed = .false.
    else
        print *, "PASS: DEFAULT_REL_TOL = ", DEFAULT_REL_TOL
    end if
    
    if (abs(INTEGRATION_TOL - 1.0e-3_real64) > epsilon(1.0_real64)) then
        print *, "FAIL: INTEGRATION_TOL incorrect"
        test_passed = .false.
    else
        print *, "PASS: INTEGRATION_TOL = ", INTEGRATION_TOL
    end if
    
    if (abs(COLLISION_FREQ_FACTOR - 1.4e-7_real64) > epsilon(COLLISION_FREQ_FACTOR)) then
        print *, "FAIL: COLLISION_FREQ_FACTOR incorrect"
        test_passed = .false.
    else
        print *, "PASS: COLLISION_FREQ_FACTOR = ", COLLISION_FREQ_FACTOR
    end if
    
    print *, ""
    print *, "Testing coordinate system constants..."
    
    ! Test coordinate systems
    if (COORD_CYLINDRICAL /= 0 .or. COORD_CARTESIAN /= 1 .or. COORD_SPHERICAL /= 2) then
        print *, "FAIL: Coordinate system constants incorrect"
        test_passed = .false.
    else
        print *, "PASS: Coordinate system constants correct"
    end if
    
    if (CYL_COMPONENTS /= "rtz") then
        print *, "FAIL: CYL_COMPONENTS incorrect"
        test_passed = .false.
    else
        print *, "PASS: CYL_COMPONENTS = ", CYL_COMPONENTS
    end if
    
    print *, ""
    if (test_passed) then
        print *, "All tests PASSED!"
    else
        print *, "Some tests FAILED!"
        stop 1
    end if
    
end program test_kilca_constants