program test_velocity_integration
    use iso_fortran_env, only: real64
    use kilca_plasma_physics_m
    implicit none
    
    integer :: test_status = 0
    real(real64), parameter :: tolerance = 1.0e-6_real64
    
    print *, "==========================================="
    print *, "Testing Velocity Space Integration [RED PHASE - SHOULD FAIL]"
    print *, "==========================================="
    
    call test_gaussian_quadrature()
    call test_maxwell_integration()
    call test_landau_damping()
    call test_adaptive_integration()
    call test_multi_dimensional_integration()
    
    if (test_status == 0) then
        print *, ""
        print *, "Velocity integration tests PASSED (unexpected in RED phase)"
        stop 0
    else
        print *, ""
        print *, "Velocity integration tests FAILED:", test_status, "test(s) failed (expected in RED phase)"
        stop 1
    end if

contains

    !> Test Gaussian quadrature for velocity integration
    subroutine test_gaussian_quadrature()
        real(real64) :: integral, expected
        real(real64) :: v_min, v_max, v_thermal
        integer :: n_points
        
        print *, ""
        print *, "Testing Gaussian quadrature integration..."
        print *, ""
        
        ! Test 1: Integrate constant function
        v_min = -5.0_real64
        v_max = 5.0_real64
        v_thermal = 1.0_real64
        n_points = 32
        
        integral = velocity_space_integrate(constant_func, v_min, v_max, n_points)
        expected = 10.0_real64  ! (v_max - v_min) * 1.0
        
        if (abs(integral - expected) > tolerance) then
            print *, "FAIL: Constant function integration. Got:", integral, "Expected:", expected
            test_status = test_status + 1
        else
            print *, "PASS: Constant function integration"
        end if
        
        ! Test 2: Integrate linear function
        integral = velocity_space_integrate(linear_func, v_min, v_max, n_points)
        expected = 0.0_real64  ! Symmetric integral of linear function
        
        if (abs(integral - expected) > tolerance) then
            print *, "FAIL: Linear function integration. Got:", integral, "Expected:", expected
            test_status = test_status + 1
        else
            print *, "PASS: Linear function integration"
        end if
        
        ! Test 3: Integrate quadratic function
        integral = velocity_space_integrate(quadratic_func, v_min, v_max, n_points)
        expected = 83.33333333_real64  ! ∫(v²)dv from -5 to 5 = 2*5³/3 ≈ 83.33
        
        if (abs(integral - expected)/abs(expected) > 0.01_real64) then
            print *, "FAIL: Quadratic function integration. Got:", integral, "Expected:", expected
            test_status = test_status + 1
        else
            print *, "PASS: Quadratic function integration"
        end if
        
    end subroutine test_gaussian_quadrature
    
    !> Test Maxwell velocity distribution integration
    subroutine test_maxwell_integration()
        real(real64) :: integral, expected
        real(real64) :: v_thermal, normalization
        real(real64), parameter :: pi = 3.14159265358979323846_real64
        
        print *, ""
        print *, "Testing Maxwell distribution integration..."
        print *, ""
        
        ! Test 1: Normalization of Maxwell distribution
        v_thermal = 1.0_real64
        
        ! Integrate Maxwell distribution - should be normalized to 1
        integral = integrate_maxwell_1d(-10.0_real64, 10.0_real64, v_thermal, 64)
        expected = 1.0_real64
        
        if (abs(integral - expected) > 0.01_real64) then
            print *, "FAIL: Maxwell normalization. Got:", integral, "Expected:", expected
            test_status = test_status + 1
        else
            print *, "PASS: Maxwell distribution normalized"
        end if
        
        ! Test 2: Mean velocity should be zero
        integral = integrate_maxwell_velocity(-10.0_real64, 10.0_real64, v_thermal, 64)
        expected = 0.0_real64
        
        if (abs(integral) > 0.01_real64) then
            print *, "FAIL: Maxwell mean velocity. Got:", integral, "Expected:", expected
            test_status = test_status + 1
        else
            print *, "PASS: Maxwell mean velocity is zero"
        end if
        
        ! Test 3: RMS velocity squared should be v_thermal²/2 for 1D Maxwell
        integral = integrate_maxwell_velocity_squared(-10.0_real64, 10.0_real64, v_thermal, 64)
        expected = v_thermal**2 / 2.0_real64  ! For 1D: <v²> = kT/m = v_th²/2
        
        if (abs(integral - expected)/expected > 0.01_real64) then
            print *, "FAIL: Maxwell RMS velocity. Got:", integral, "Expected:", expected
            test_status = test_status + 1
        else
            print *, "PASS: Maxwell RMS velocity correct"
        end if
        
    end subroutine test_maxwell_integration
    
    !> Test Landau damping rate calculation
    subroutine test_landau_damping()
        real(real64) :: gamma, expected
        real(real64) :: omega, k, v_thermal
        real(real64) :: phase_velocity
        
        print *, ""
        print *, "Testing Landau damping calculation..."
        print *, ""
        
        ! Test 1: Weak damping regime (v_phase >> v_thermal)
        v_thermal = 1.0_real64
        k = 1.0_real64
        omega = 5.0_real64  ! Phase velocity = 5 >> v_thermal
        
        gamma = landau_damping_rate(omega, k, v_thermal)
        phase_velocity = omega / k
        
        ! For v_phase >> v_th, γ ∝ exp(-(v_phase/v_th)²)
        expected = -omega * exp(-(phase_velocity/v_thermal)**2)
        
        if (abs(gamma - expected)/abs(expected) > 0.1_real64) then
            print *, "FAIL: Weak Landau damping. Got:", gamma, "Expected:", expected
            test_status = test_status + 1
        else
            print *, "PASS: Weak Landau damping rate"
        end if
        
        ! Test 2: Strong damping regime (v_phase ~ v_thermal)
        omega = 1.0_real64  ! Phase velocity = 1 ~ v_thermal
        
        gamma = landau_damping_rate(omega, k, v_thermal)
        
        ! Should be significant damping
        if (abs(gamma) < 0.1_real64) then
            print *, "FAIL: Strong Landau damping too weak. Got:", gamma
            test_status = test_status + 1
        else
            print *, "PASS: Strong Landau damping significant"
        end if
        
    end subroutine test_landau_damping
    
    !> Test adaptive integration methods
    subroutine test_adaptive_integration()
        real(real64) :: integral, expected
        real(real64) :: error_estimate
        
        print *, ""
        print *, "Testing adaptive integration..."
        print *, ""
        
        ! Test 1: Adaptive integration of oscillatory function
        integral = adaptive_integrate_oscillatory(-10.0_real64, 10.0_real64, 1.0_real64, &
                                                  1.0e-6_real64, error_estimate)
        
        ! Oscillatory function should integrate to near zero
        if (abs(integral) > 0.1_real64) then
            print *, "FAIL: Adaptive oscillatory integration. Got:", integral
            test_status = test_status + 1
        else
            print *, "PASS: Adaptive oscillatory integration"
        end if
        
        ! Test 2: Error estimate should be small
        if (error_estimate > 1.0e-5_real64) then
            print *, "FAIL: Adaptive error estimate too large:", error_estimate
            test_status = test_status + 1
        else
            print *, "PASS: Adaptive error estimate acceptable"
        end if
        
    end subroutine test_adaptive_integration
    
    !> Test multi-dimensional velocity integration
    subroutine test_multi_dimensional_integration()
        real(real64) :: integral_3d, expected
        real(real64) :: v_thermal
        
        print *, ""
        print *, "Testing 3D velocity space integration..."
        print *, ""
        
        v_thermal = 1.0_real64
        
        ! Test: 3D Maxwell distribution normalization
        integral_3d = integrate_maxwell_3d(v_thermal, 32)
        expected = 1.0_real64  ! Should be normalized
        
        if (abs(integral_3d - expected) > 0.01_real64) then
            print *, "FAIL: 3D Maxwell normalization. Got:", integral_3d, "Expected:", expected
            test_status = test_status + 1
        else
            print *, "PASS: 3D Maxwell distribution normalized"
        end if
        
    end subroutine test_multi_dimensional_integration
    
    !---------------------------------------------------------------------------
    ! Test functions for integration
    !---------------------------------------------------------------------------
    
    function constant_func(v) result(val)
        real(real64), intent(in) :: v
        real(real64) :: val
        val = 1.0_real64
    end function constant_func
    
    function linear_func(v) result(val)
        real(real64), intent(in) :: v
        real(real64) :: val
        val = v
    end function linear_func
    
    function quadratic_func(v) result(val)
        real(real64), intent(in) :: v
        real(real64) :: val
        val = v * v
    end function quadratic_func
    
    !---------------------------------------------------------------------------
    ! Helper functions for Maxwell integration
    !---------------------------------------------------------------------------
    
    function integrate_maxwell_1d(v_min, v_max, v_thermal, n_points) result(integral)
        real(real64), intent(in) :: v_min, v_max, v_thermal
        integer, intent(in) :: n_points
        real(real64) :: integral
        real(real64) :: dv, v, sum_val
        integer :: i
        
        ! Simple trapezoidal integration of Maxwell distribution
        dv = (v_max - v_min) / real(n_points - 1, real64)
        sum_val = 0.0_real64
        
        do i = 1, n_points
            v = v_min + real(i-1, real64) * dv
            if (i == 1 .or. i == n_points) then
                sum_val = sum_val + 0.5_real64 * maxwell_velocity_distribution(v, v_thermal)
            else
                sum_val = sum_val + maxwell_velocity_distribution(v, v_thermal)
            end if
        end do
        
        integral = sum_val * dv
    end function integrate_maxwell_1d
    
    function integrate_maxwell_velocity(v_min, v_max, v_thermal, n_points) result(integral)
        real(real64), intent(in) :: v_min, v_max, v_thermal
        integer, intent(in) :: n_points
        real(real64) :: integral
        real(real64) :: dv, v, sum_val
        integer :: i
        
        ! Trapezoidal integration of v * Maxwell distribution
        dv = (v_max - v_min) / real(n_points - 1, real64)
        sum_val = 0.0_real64
        
        do i = 1, n_points
            v = v_min + real(i-1, real64) * dv
            if (i == 1 .or. i == n_points) then
                sum_val = sum_val + 0.5_real64 * v * maxwell_velocity_distribution(v, v_thermal)
            else
                sum_val = sum_val + v * maxwell_velocity_distribution(v, v_thermal)
            end if
        end do
        
        integral = sum_val * dv
    end function integrate_maxwell_velocity
    
    function integrate_maxwell_velocity_squared(v_min, v_max, v_thermal, n_points) result(integral)
        real(real64), intent(in) :: v_min, v_max, v_thermal
        integer, intent(in) :: n_points
        real(real64) :: integral
        real(real64) :: dv, v, sum_val
        integer :: i
        
        ! Trapezoidal integration of v² * Maxwell distribution
        dv = (v_max - v_min) / real(n_points - 1, real64)
        sum_val = 0.0_real64
        
        do i = 1, n_points
            v = v_min + real(i-1, real64) * dv
            if (i == 1 .or. i == n_points) then
                sum_val = sum_val + 0.5_real64 * v * v * maxwell_velocity_distribution(v, v_thermal)
            else
                sum_val = sum_val + v * v * maxwell_velocity_distribution(v, v_thermal)
            end if
        end do
        
        integral = sum_val * dv
    end function integrate_maxwell_velocity_squared
    
    function adaptive_integrate_oscillatory(v_min, v_max, v_thermal, tol, error_est) result(integral)
        real(real64), intent(in) :: v_min, v_max, v_thermal, tol
        real(real64), intent(out) :: error_est
        real(real64) :: integral
        real(real64) :: integral_n, integral_2n
        integer :: n
        
        ! Simple adaptive integration using Richardson extrapolation
        n = 32
        integral_n = integrate_oscillatory_helper(v_min, v_max, v_thermal, n)
        integral_2n = integrate_oscillatory_helper(v_min, v_max, v_thermal, 2*n)
        
        ! Richardson extrapolation
        integral = (4.0_real64 * integral_2n - integral_n) / 3.0_real64
        error_est = abs(integral_2n - integral_n) / 3.0_real64
    end function adaptive_integrate_oscillatory
    
    function integrate_oscillatory_helper(v_min, v_max, v_thermal, n_points) result(integral)
        real(real64), intent(in) :: v_min, v_max, v_thermal
        integer, intent(in) :: n_points
        real(real64) :: integral
        
        integral = velocity_space_integrate(oscillatory_test_func, v_min, v_max, n_points)
    end function integrate_oscillatory_helper
    
    function oscillatory_test_func(v) result(f)
        real(real64), intent(in) :: v
        real(real64) :: f
        f = sin(10.0_real64 * v) * exp(-(v/1.0_real64)**2)  ! Using v_thermal=1 as constant
    end function oscillatory_test_func
    
    function integrate_maxwell_3d(v_thermal, n_points) result(integral)
        real(real64), intent(in) :: v_thermal
        integer, intent(in) :: n_points
        real(real64) :: integral
        real(real64), parameter :: pi = 3.14159265358979323846_real64
        real(real64) :: v_max, dv, v, sum_val
        integer :: i
        
        ! For 3D Maxwell in spherical coordinates
        v_max = 6.0_real64 * v_thermal
        
        ! Trapezoidal integration
        dv = v_max / real(n_points - 1, real64)
        sum_val = 0.0_real64
        
        do i = 1, n_points
            v = real(i-1, real64) * dv
            ! 3D Maxwell-Boltzmann distribution in spherical coords
            if (i == 1 .or. i == n_points) then
                sum_val = sum_val + 0.5_real64 * 4.0_real64 * pi * v * v * &
                    (1.0_real64/(v_thermal * sqrt(2.0_real64 * pi)))**3 * &
                    exp(-(v/v_thermal)**2 / 2.0_real64)
            else
                sum_val = sum_val + 4.0_real64 * pi * v * v * &
                    (1.0_real64/(v_thermal * sqrt(2.0_real64 * pi)))**3 * &
                    exp(-(v/v_thermal)**2 / 2.0_real64)
            end if
        end do
        
        integral = sum_val * dv
    end function integrate_maxwell_3d

end program test_velocity_integration