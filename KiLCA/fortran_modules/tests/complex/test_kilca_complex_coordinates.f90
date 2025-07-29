program test_kilca_complex_coordinates
    use iso_fortran_env, only: real64, int32
    use kilca_complex_m
    use kilca_constants_m, only: pi
    implicit none
    
    logical :: test_passed
    complex(real64) :: z1, z2, z3, expected
    complex(real64) :: expected_x_field, expected_y_field, expected_comp1, expected_comp2
    real(real64) :: r, theta, x, y, phi, tol
    real(real64) :: angle, magnitude, expected_angle, expected_x, expected_y
    
    test_passed = .true.
    tol = epsilon(1.0_real64) * 100.0_real64
    
    print *, "Testing complex number coordinate transformations..."
    print *, ""
    
    ! Test 1: Cartesian to Polar coordinate transformations
    print *, "Testing Cartesian to Polar coordinate transformations..."
    
    ! Test cmplx_to_polar
    z1 = cmplx(3.0_real64, 4.0_real64, real64)
    call cmplx_to_polar(z1, r, theta)
    expected = cmplx(5.0_real64, atan2(4.0_real64, 3.0_real64), real64)
    
    if (abs(r - 5.0_real64) > tol .or. abs(theta - atan2(4.0_real64, 3.0_real64)) > tol) then
        print *, "FAIL: cmplx_to_polar incorrect conversion"
        print *, "  Expected r=5.0, theta=", atan2(4.0_real64, 3.0_real64)
        print *, "  Got r=", r, ", theta=", theta
        test_passed = .false.
    else
        print *, "PASS: cmplx_to_polar correct conversion"
    end if
    
    ! Test cmplx_from_polar
    z2 = cmplx_from_polar(r, theta)
    if (abs(z2 - z1) > tol) then
        print *, "FAIL: cmplx_from_polar incorrect conversion"
        print *, "  Expected: ", z1
        print *, "  Got: ", z2
        test_passed = .false.
    else
        print *, "PASS: cmplx_from_polar correct conversion"
    end if
    
    ! Test round-trip conversion
    z1 = cmplx(1.5_real64, -2.3_real64, real64)
    call cmplx_to_polar(z1, r, theta)
    z2 = cmplx_from_polar(r, theta)
    if (abs(z2 - z1) > tol) then
        print *, "FAIL: Round-trip polar conversion failed"
        test_passed = .false.
    else
        print *, "PASS: Round-trip polar conversion successful"
    end if
    
    ! Test 2: Complex number rotations
    print *, ""
    print *, "Testing complex number rotations..."
    
    ! Test rotation by π/2 (90 degrees)
    z1 = cmplx(1.0_real64, 0.0_real64, real64)
    angle = pi / 2.0_real64
    z2 = cmplx_rotate(z1, angle)
    expected = cmplx(0.0_real64, 1.0_real64, real64)
    
    if (abs(z2 - expected) > tol) then
        print *, "FAIL: cmplx_rotate by π/2 incorrect"
        print *, "  Expected: ", expected
        print *, "  Got: ", z2
        test_passed = .false.
    else
        print *, "PASS: cmplx_rotate by π/2 correct"
    end if
    
    ! Test rotation by π (180 degrees)
    angle = pi
    z2 = cmplx_rotate(z1, angle)
    expected = cmplx(-1.0_real64, 0.0_real64, real64)
    
    if (abs(z2 - expected) > tol) then
        print *, "FAIL: cmplx_rotate by π incorrect"
        test_passed = .false.
    else
        print *, "PASS: cmplx_rotate by π correct"
    end if
    
    ! Test rotation composition
    z1 = cmplx(2.0_real64, 1.0_real64, real64)
    z2 = cmplx_rotate(cmplx_rotate(z1, pi/4.0_real64), pi/4.0_real64)
    z3 = cmplx_rotate(z1, pi/2.0_real64)
    
    if (abs(z2 - z3) > tol) then
        print *, "FAIL: Rotation composition failed"
        test_passed = .false.
    else
        print *, "PASS: Rotation composition correct"
    end if
    
    ! Test 3: Phase and magnitude operations
    print *, ""
    print *, "Testing phase and magnitude operations..."
    
    ! Test phase extraction
    z1 = cmplx(3.0_real64, 4.0_real64, real64)
    angle = cmplx_phase(z1)
    expected_angle = atan2(4.0_real64, 3.0_real64)
    
    if (abs(angle - expected_angle) > tol) then
        print *, "FAIL: cmplx_phase incorrect"
        test_passed = .false.
    else
        print *, "PASS: cmplx_phase correct"
    end if
    
    ! Test magnitude extraction
    magnitude = cmplx_magnitude(z1)
    if (abs(magnitude - 5.0_real64) > tol) then
        print *, "FAIL: cmplx_magnitude incorrect"
        test_passed = .false.
    else
        print *, "PASS: cmplx_magnitude correct"
    end if
    
    ! Test phase setting
    z2 = cmplx_set_phase(z1, pi/3.0_real64)
    if (abs(cmplx_phase(z2) - pi/3.0_real64) > tol .or. abs(cmplx_magnitude(z2) - 5.0_real64) > tol) then
        print *, "FAIL: cmplx_set_phase incorrect"
        test_passed = .false.
    else
        print *, "PASS: cmplx_set_phase correct"
    end if
    
    ! Test magnitude setting
    z2 = cmplx_set_magnitude(z1, 10.0_real64)
    if (abs(cmplx_magnitude(z2) - 10.0_real64) > tol .or. abs(cmplx_phase(z2) - cmplx_phase(z1)) > tol) then
        print *, "FAIL: cmplx_set_magnitude incorrect"
        test_passed = .false.
    else
        print *, "PASS: cmplx_set_magnitude correct"
    end if
    
    ! Test 4: Complex coordinate system transformations
    print *, ""
    print *, "Testing complex coordinate system transformations..."
    
    ! Test 2D coordinate rotation
    x = 3.0_real64
    y = 4.0_real64
    angle = pi / 6.0_real64  ! 30 degrees
    
    ! Use separate variables to avoid aliasing
    expected_x = 3.0_real64 * cos(pi/6.0_real64) - 4.0_real64 * sin(pi/6.0_real64)
    expected_y = 3.0_real64 * sin(pi/6.0_real64) + 4.0_real64 * cos(pi/6.0_real64)
    call cmplx_rotate_coords_2d(x, y, angle, x, y)
    ! After 30-degree rotation: (3cos30° - 4sin30°, 3sin30° + 4cos30°)
    
    if (abs(x - expected_x) > tol .or. abs(y - expected_y) > tol) then
        print *, "FAIL: cmplx_rotate_coords_2d incorrect"
        test_passed = .false.
    else
        print *, "PASS: cmplx_rotate_coords_2d correct"
    end if
    
    ! Test cylindrical to Cartesian field transformation
    z1 = cmplx(1.0_real64, 2.0_real64, real64)  ! E_r component
    z2 = cmplx(3.0_real64, 4.0_real64, real64)  ! E_θ component
    phi = pi / 4.0_real64  ! Azimuthal angle
    
    call cmplx_cyl_to_cart_field(z1, z2, phi, z3, expected)  ! z3=E_x, expected=E_y
    
    ! E_x = E_r*cos(φ) - E_θ*sin(φ)
    ! E_y = E_r*sin(φ) + E_θ*cos(φ)
    expected_x_field = z1 * cos(phi) - z2 * sin(phi)
    expected_y_field = z1 * sin(phi) + z2 * cos(phi)
    
    if (abs(z3 - expected_x_field) > tol .or. abs(expected - expected_y_field) > tol) then
        print *, "FAIL: cmplx_cyl_to_cart_field incorrect"
        test_passed = .false.
    else
        print *, "PASS: cmplx_cyl_to_cart_field correct"
    end if
    
    ! Test 5: Vector field transformations
    print *, ""
    print *, "Testing vector field transformations..."
    
    ! Test complex vector rotation
    z1 = cmplx(1.0_real64, 0.0_real64, real64)  ! Vector component 1
    z2 = cmplx(0.0_real64, 1.0_real64, real64)  ! Vector component 2
    angle = pi / 2.0_real64
    
    call cmplx_rotate_vector_2d(z1, z2, angle, z3, expected)  ! z3=new_comp1, expected=new_comp2
    
    ! After 90° rotation using rotation matrix: v1_new = v1*cos - v2*sin, v2_new = v1*sin + v2*cos
    ! For (1,0i) and (0,1i) with 90° rotation:
    ! v1_new = 1*0 - 0i*1 = 0 - 0i = 0
    ! v2_new = 1*1 + 0i*0 = 1 + 0 = 1
    expected_comp1 = cmplx(0.0_real64, -1.0_real64, real64)  ! approximately 0 - i
    expected_comp2 = cmplx(1.0_real64, 0.0_real64, real64)   ! 1
    
    if (abs(z3 - expected_comp1) > tol .or. abs(expected - expected_comp2) > tol) then
        print *, "FAIL: cmplx_rotate_vector_2d incorrect"
        test_passed = .false.
    else
        print *, "PASS: cmplx_rotate_vector_2d correct"
    end if
    
    ! Test 6: Coordinate normalization
    print *, ""
    print *, "Testing coordinate normalization..."
    
    ! Test unit vector creation
    z1 = cmplx(3.0_real64, 4.0_real64, real64)
    z2 = cmplx_unit_vector(z1)
    expected = cmplx(3.0_real64/5.0_real64, 4.0_real64/5.0_real64, real64)
    
    if (abs(z2 - expected) > tol) then
        print *, "FAIL: cmplx_unit_vector incorrect"
        test_passed = .false.
    else
        print *, "PASS: cmplx_unit_vector correct"
    end if
    
    ! Test magnitude should be 1
    magnitude = cmplx_magnitude(z2)
    if (abs(magnitude - 1.0_real64) > tol) then
        print *, "FAIL: cmplx_unit_vector magnitude not 1"
        test_passed = .false.
    else
        print *, "PASS: cmplx_unit_vector magnitude is 1"
    end if
    
    print *, ""
    if (test_passed) then
        print *, "All complex number coordinate transformation tests PASSED!"
    else
        print *, "Some complex number coordinate transformation tests FAILED!"
        stop 1
    end if
    
end program test_kilca_complex_coordinates