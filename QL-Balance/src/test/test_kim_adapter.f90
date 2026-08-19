program test_kim_adapter
    !
    ! Unit tests for the KIM wave code adapter module.
    !
    ! These tests verify:
    !   1. Module compilation and linkage (the adapter + KIM_lib resolve)
    !   2. interp_complex_profile helper correctness
    !   3. wave_code selector variable default
    !   4. Rescaled QL-Balance profiles change the subsequent KIM solution
    !
    use kim_wave_code_adapter_m, only: interp_complex_profile, kim_initialize, &
        kim_update_profiles, kim_run_for_all_modes, kim_Br_modes, kim_vac_Br
    use kim_wave_code_adapter_m, only: kim_Bparallel_modes, kim_get_wave_fields
    use control_mod, only: wave_code, kim_config_path, kim_profiles_from_balance, &
        type_of_run, kim_run_type
    use wave_code_data, only: dim_mn, m_vals, n_vals, r, n, Te, Ti, q, &
        Vth, Vz, dPhi0, Bp
    use plasma_parameters, only: params_b
    use grid_mod, only: Ercov
    use baseparam_mod, only: ev

    implicit none

    integer :: num_passed, num_failed
    external :: allocate_wave_code_data

    num_passed = 0
    num_failed = 0

    print *, "========================================"
    print *, "  KIM Adapter Module Tests"
    print *, "========================================"
    print *, ""

    call test_wave_code_default()
    call test_interp_complex_constant()
    call test_interp_complex_linear()
    call test_rescaled_profiles_change_kim_solution()

    print *, ""
    print *, "========================================"
    print *, "  Test Summary"
    print *, "========================================"
    print '(A,I3,A)', "  Passed: ", num_passed, " tests"
    print '(A,I3,A)', "  Failed: ", num_failed, " tests"
    print *, "========================================"

    if (num_failed > 0) then
        print *, "TESTS FAILED"
        stop 1
    else
        print *, "ALL TESTS PASSED"
    end if

contains

    subroutine assert_equal_real(actual, expected, tol, label)
        real(8), intent(in) :: actual, expected, tol
        character(len=*), intent(in) :: label

        if (abs(actual - expected) > tol) then
            print '(A,A,A,ES15.8,A,ES15.8)', "  FAIL: ", label, &
                " got ", actual, " expected ", expected
            num_failed = num_failed + 1
        else
            print '(A,A)', "  PASS: ", label
            num_passed = num_passed + 1
        end if
    end subroutine

    subroutine assert_equal_complex(actual, expected, tol, label)
        complex(8), intent(in) :: actual, expected
        real(8), intent(in) :: tol
        character(len=*), intent(in) :: label

        if (abs(actual - expected) > tol) then
            print '(A,A,A,2ES15.8,A,2ES15.8)', "  FAIL: ", label, &
                " got (", real(actual), aimag(actual), &
                ") expected (", real(expected), aimag(expected), ")"
            num_failed = num_failed + 1
        else
            print '(A,A)', "  PASS: ", label
            num_passed = num_passed + 1
        end if
    end subroutine

    ! ------------------------------------------------------------------
    ! Test: wave_code defaults to 'KiLCA'
    ! ------------------------------------------------------------------
    subroutine test_wave_code_default()
        print *, "--- test_wave_code_default ---"

        if (trim(wave_code) == 'KiLCA') then
            print '(A)', "  PASS: wave_code default is 'KiLCA'"
            num_passed = num_passed + 1
        else
            print '(A,A,A)', "  FAIL: wave_code default is '", &
                trim(wave_code), "', expected 'KiLCA'"
            num_failed = num_failed + 1
        end if
    end subroutine

    ! ------------------------------------------------------------------
    ! Test: interp_complex_profile with constant function
    ! A constant complex function should interpolate exactly.
    ! ------------------------------------------------------------------
    subroutine test_interp_complex_constant()
        integer, parameter :: n_old = 10, n_new = 3
        real(8) :: r_old(n_old), r_new(n_new)
        complex(8) :: z_old(n_old), z_new(n_new)
        complex(8), parameter :: val = (2.5d0, -1.3d0)
        real(8), parameter :: tol = 1.0d-10
        integer :: i

        print *, "--- test_interp_complex_constant ---"

        ! Old grid: [10, 20, 30, ..., 100]
        do i = 1, n_old
            r_old(i) = 10.0d0 * dble(i)
            z_old(i) = val
        end do

        ! New grid: [15, 25, 35] (within old grid bounds)
        r_new(1) = 15.0d0
        r_new(2) = 25.0d0
        r_new(3) = 35.0d0

        call interp_complex_profile(n_old, r_old, z_old, n_new, r_new, z_new)

        do i = 1, n_new
            call assert_equal_complex(z_new(i), val, tol, &
                "constant interp at new grid point")
        end do
    end subroutine

    ! ------------------------------------------------------------------
    ! Test: interp_complex_profile with linear function
    ! z(r) = (1+2i)*r should interpolate exactly with polynomial interp.
    ! ------------------------------------------------------------------
    subroutine test_interp_complex_linear()
        integer, parameter :: n_old = 10, n_new = 5
        real(8) :: r_old(n_old), r_new(n_new)
        complex(8) :: z_old(n_old), z_new(n_new)
        complex(8) :: expected
        real(8), parameter :: tol = 1.0d-8
        integer :: i

        print *, "--- test_interp_complex_linear ---"

        ! Old grid: uniform [10, 20, 30, ..., 100]
        do i = 1, n_old
            r_old(i) = 10.0d0 * dble(i)
            z_old(i) = cmplx(1.0d0, 2.0d0, kind=8) * r_old(i)
        end do

        ! New grid: non-coincident points
        r_new(1) = 15.0d0
        r_new(2) = 25.0d0
        r_new(3) = 55.0d0
        r_new(4) = 75.0d0
        r_new(5) = 95.0d0

        call interp_complex_profile(n_old, r_old, z_old, n_new, r_new, z_new)

        do i = 1, n_new
            expected = cmplx(1.0d0, 2.0d0, kind=8) * r_new(i)
            call assert_equal_complex(z_new(i), expected, tol, &
                "linear interp at new grid point")
        end do
    end subroutine

    ! ------------------------------------------------------------------
    ! Regression for #136: a parameter-scan profile rescaling must cross
    ! the adapter boundary and change KIM's next wave solution.
    ! ------------------------------------------------------------------
    subroutine test_rescaled_profiles_change_kim_solution()
        integer, parameter :: npts = 40
        integer, parameter :: m_mode = -6, n_mode = 2
        real(8) :: r_grid(npts), frac, response_change, response_scale
        real(8), allocatable :: base_n(:)
        complex(8), allocatable :: br_rescaled(:)
        integer :: i

        print *, "--- test_rescaled_profiles_change_kim_solution ---"

        dim_mn = 1
        allocate(m_vals(dim_mn), n_vals(dim_mn))
        m_vals = m_mode
        n_vals = n_mode

        do i = 1, npts
            r_grid(i) = 3.0d0 + 64.0d0 * dble(i - 1) / dble(npts - 1)
        end do
        call allocate_wave_code_data(npts, r_grid)

        do i = 1, npts
            frac = (r(i) - 3.0d0) / 64.0d0
            n(i) = 5.0d13 * (1.0d0 - 0.9d0 * frac)
            Te(i) = 100.0d0 + 1900.0d0 * (1.0d0 - frac)
            Ti(i) = 0.9d0 * Te(i)
            q(i) = 1.0d0 + 3.0d0 * frac
        end do
        Vth = 0.0d0
        Vz = 0.0d0
        dPhi0 = 0.0d0

        kim_config_path = "KIM_config_em_small.nml"
        kim_profiles_from_balance = .true.
        kim_run_type = "electromagnetic"
        type_of_run = "ParameterScan"

        call kim_initialize(npts, r_grid)
        kim_vac_Br = (0.0d0, 0.0d0)

        allocate(base_n(npts))
        base_n = n
        allocate(params_b(4, npts), Ercov(npts))
        params_b(1, :) = 1.6d0 * base_n
        params_b(2, :) = 0.0d0
        params_b(3, :) = Te * ev
        params_b(4, :) = Ti * ev
        Ercov = 0.0d0

        call kim_update_profiles()
        call kim_run_for_all_modes()

        ! B_parallel is a first-class KIM output and must reach the
        ! wave-code contract as Bp (the RSP parallel component).
        call kim_get_wave_fields(1)
        if (maxval(abs(Bp - kim_Bparallel_modes(:, 1))) <= 1.0d-12) then
            print '(A)', "  PASS: KIM Bparallel reaches wave_code_data Bp"
            num_passed = num_passed + 1
        else
            print '(A)', "  FAIL: KIM Bparallel was not copied to Bp"
            num_failed = num_failed + 1
        end if

        allocate(br_rescaled(npts))
        br_rescaled = kim_Br_modes(:, 1)

        params_b(1, :) = base_n
        call kim_update_profiles()
        call kim_run_for_all_modes()

        response_change = maxval(abs(kim_Br_modes(:, 1) - br_rescaled))
        response_scale = max(maxval(abs(br_rescaled)), 1.0d-30)
        if (response_change > 1.0d-6 * response_scale) then
            print '(A,ES15.8)', "  PASS: rescaled density changed KIM Br by ", &
                response_change
            num_passed = num_passed + 1
        else
            print '(A,ES15.8)', "  FAIL: KIM Br remained frozen; max change = ", &
                response_change
            num_failed = num_failed + 1
        end if
    end subroutine

end program test_kim_adapter
