program test_rhs_balance
    !
    ! Unit tests for rhs_balance module helper functions
    !
    ! Test strategy: Each helper function is tested with known inputs
    ! against analytically computed expected outputs.
    !
    ! Tolerance: 1e-12 (relative) for floating-point comparisons
    !
    use QLBalance_kinds, only: dp
    use rhs_balance_m, only: thermodynamic_forces_t, &
                             compute_thermodynamic_forces, &
                             compute_particle_fluxes, &
                             compute_total_heat_fluxes, &
                             species_fluxes_t, &
                             transport_fluxes_t
    use baseparam_mod, only: e_charge

    implicit none

    integer :: num_passed, num_failed

    num_passed = 0
    num_failed = 0

    print *, "========================================"
    print *, "  RHS Balance Module Tests"
    print *, "========================================"
    print *, ""

    call test_thermodynamic_forces()
    call test_particle_fluxes()
    call test_heat_fluxes()

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

    subroutine test_thermodynamic_forces()
        !
        ! Test compute_thermodynamic_forces with known inputs
        !
        ! Physics:
        !   A_noE_1e = d(ln n)/dr - 1.5 * d(ln Te)/dr
        !   A_noE_2e = d(ln Te)/dr
        !   A_noE_1i = d(ln n)/dr - 1.5 * d(ln Ti)/dr
        !   A_noE_2i = d(ln Ti)/dr
        !   A_1e = A_noE_1e + e*Ercov/Te
        !   A_1i = A_noE_1i - Z_i*e*Ercov/Ti
        !
        type(thermodynamic_forces_t) :: forces
        real(dp) :: ddr_n, ddr_Te, ddr_Ti
        real(dp) :: n_b, Te_b, Ti_b
        real(dp) :: Ercov_val, Z_i_val
        real(dp) :: expected_A_noE_1e, expected_A_noE_2e
        real(dp) :: expected_A_noE_1i, expected_A_noE_2i
        real(dp) :: expected_A_1e, expected_A_1i
        character(len=*), parameter :: test_name = "test_thermodynamic_forces"

        print *, "Running: ", test_name

        ! Test case 1: Simple values
        ! Density gradient: 1e12 cm^-3 / cm = 1e12 cm^-4
        ! Temperature gradients: Te = 100 eV, Ti = 80 eV (in erg)
        ! Using CGS units as in the original code

        n_b = 1.0e13_dp         ! density at boundary [cm^-3]
        Te_b = 1.6022e-10_dp    ! 100 eV in erg
        Ti_b = 1.2818e-10_dp    ! 80 eV in erg

        ddr_n = 1.0e12_dp       ! dn/dr [cm^-4]
        ddr_Te = 1.6022e-11_dp  ! dTe/dr [erg/cm]
        ddr_Ti = 1.2818e-11_dp  ! dTi/dr [erg/cm]

        Z_i_val = 1.0_dp              ! hydrogen
        Ercov_val = 1.0e3_dp          ! radial E-field [statV/cm]

        ! Expected values (analytical calculation):
        ! A_noE_1e = ddr_n/n_b - 1.5*ddr_Te/Te_b
        !          = 1e12/1e13 - 1.5*1.6022e-11/1.6022e-10
        !          = 0.1 - 1.5*0.1 = 0.1 - 0.15 = -0.05
        expected_A_noE_1e = ddr_n/n_b - 1.5_dp*ddr_Te/Te_b

        ! A_noE_2e = ddr_Te/Te_b = 1.6022e-11/1.6022e-10 = 0.1
        expected_A_noE_2e = ddr_Te/Te_b

        ! A_noE_1i = ddr_n/n_b - 1.5*ddr_Ti/Ti_b
        expected_A_noE_1i = ddr_n/n_b - 1.5_dp*ddr_Ti/Ti_b

        ! A_noE_2i = ddr_Ti/Ti_b
        expected_A_noE_2i = ddr_Ti/Ti_b

        ! A_1e = A_noE_1e + e*Ercov/Te
        expected_A_1e = expected_A_noE_1e + Ercov_val * e_charge / Te_b

        ! A_1i = A_noE_1i - Z_i*e*Ercov/Ti
        expected_A_1i = expected_A_noE_1i - Ercov_val * e_charge * Z_i_val / Ti_b

        ! Call the function under test
        call compute_thermodynamic_forces( &
            ddr_n, ddr_Te, ddr_Ti, &
            n_b, Te_b, Ti_b, &
            Ercov_val, Z_i_val, &
            forces)

        ! Verify results
        call assert_equal(forces%e%A1_noE, expected_A_noE_1e, "A_noE_1e")
        call assert_equal(forces%e%A2, expected_A_noE_2e, "A_noE_2e")
        call assert_equal(forces%i%A1_noE, expected_A_noE_1i, "A_noE_1i")
        call assert_equal(forces%i%A2, expected_A_noE_2i, "A_noE_2i")
        call assert_equal(forces%e%A1, expected_A_1e, "A_1e")
        call assert_equal(forces%i%A1, expected_A_1i, "A_1i")

        ! Test case 2: Zero electric field
        Ercov_val = 0.0_dp
        expected_A_1e = expected_A_noE_1e  ! Should equal A_noE when E=0
        expected_A_1i = expected_A_noE_1i

        call compute_thermodynamic_forces( &
            ddr_n, ddr_Te, ddr_Ti, &
            n_b, Te_b, Ti_b, &
            Ercov_val, Z_i_val, &
            forces)

        call assert_equal(forces%e%A1, expected_A_1e, "A_1e (E=0)")
        call assert_equal(forces%i%A1, expected_A_1i, "A_1i (E=0)")

        print *, "  PASSED: ", test_name

    end subroutine test_thermodynamic_forces


    subroutine test_particle_fluxes()
        !
        ! Test compute_particle_fluxes with known inputs
        !
        ! Physics:
        !   gamma_a_e = -(D11_a*A_noE_1e + D12_a*A_noE_2e) * n
        !   gamma_ql_e = -(D11_ql*A_1e + D12_ql*A_noE_2e) * n
        !   gamma_e = gamma_a_e + gamma_ql_e
        !
        type(thermodynamic_forces_t) :: forces
        type(transport_fluxes_t) :: fluxes
        real(dp) :: n_b, Z_i_val
        real(dp) :: D11_a_e, D12_a_e, D11_ql_e, D12_ql_e
        real(dp) :: D11_a_i, D12_a_i, D11_ql_i, D12_ql_i
        real(dp) :: expected_gamma_a_e, expected_gamma_ql_e, expected_gamma_e
        real(dp) :: expected_gamma_a_i, expected_gamma_ql_i, expected_gamma_i
        character(len=*), parameter :: test_name = "test_particle_fluxes"

        print *, "Running: ", test_name

        ! Set up test inputs
        n_b = 1.0e13_dp
        Z_i_val = 1.0_dp

        ! Thermodynamic forces (arbitrary test values)
        forces%e%A1_noE = -0.05_dp
        forces%e%A2 = 0.1_dp
        forces%i%A1_noE = -0.04_dp
        forces%i%A2 = 0.1_dp
        forces%e%A1 = -0.03_dp
        forces%i%A1 = -0.06_dp

        ! Diffusion coefficients (individual values)
        D11_a_e = 1.0e4_dp;  D12_a_e = 0.5e4_dp
        D11_ql_e = 2.0e4_dp; D12_ql_e = 1.0e4_dp
        D11_a_i = 0.8e4_dp;  D12_a_i = 0.4e4_dp
        D11_ql_i = 1.5e4_dp; D12_ql_i = 0.8e4_dp

        ! Expected values (analytical calculation)
        expected_gamma_a_e = -(D11_a_e*forces%e%A1_noE + D12_a_e*forces%e%A2) * n_b
        expected_gamma_ql_e = -(D11_ql_e*forces%e%A1 + D12_ql_e*forces%e%A2) * n_b
        expected_gamma_e = expected_gamma_a_e + expected_gamma_ql_e

        expected_gamma_a_i = -(D11_a_i*forces%i%A1_noE + D12_a_i*forces%i%A2) * n_b / Z_i_val
        expected_gamma_ql_i = -(D11_ql_i*forces%i%A1 + D12_ql_i*forces%i%A2) * n_b / Z_i_val
        expected_gamma_i = expected_gamma_a_i + expected_gamma_ql_i

        ! Call the function under test
        call compute_particle_fluxes(forces, n_b, Z_i_val, &
            D11_a_e, D12_a_e, D11_ql_e, D12_ql_e, &
            D11_a_i, D12_a_i, D11_ql_i, D12_ql_i, &
            fluxes)

        ! Verify results
        call assert_equal(fluxes%e%Gamma_a, expected_gamma_a_e, "gamma_a_e")
        call assert_equal(fluxes%e%Gamma_ql, expected_gamma_ql_e, "gamma_ql_e")
        call assert_equal(fluxes%e%Gamma_tot, expected_gamma_e, "gamma_e")
        call assert_equal(fluxes%i%Gamma_a, expected_gamma_a_i, "gamma_a_i")
        call assert_equal(fluxes%i%Gamma_ql, expected_gamma_ql_i, "gamma_ql_i")
        call assert_equal(fluxes%i%Gamma_tot, expected_gamma_i, "gamma_i")

        print *, "  PASSED: ", test_name

    end subroutine test_particle_fluxes


    subroutine test_heat_fluxes()
        !
        ! Test compute_total_heat_fluxes with known inputs
        !
        ! Physics:
        !   Q_e = -(D12_a*A_noE_1e + D21_ql*A_1e + D22*A_noE_2e) * n * Te
        !   Q_i = -(D12_a*A_noE_1i + D21_ql*A_1i + D22*A_noE_2i) * n * Ti / Z_i
        !
        ! Note: D12 (anomalous) couples to A_noE_1 (no E-field)
        !       D21 (quasi-linear) couples to A_1 (with E-field)
        !       D22 couples to temperature gradient force A_noE_2
        !
        type(thermodynamic_forces_t) :: forces
        real(dp) :: n_b, Te_b, Ti_b, Z_i_val
        real(dp) :: D12_a_e, D21_ql_e, D22_a_e, D22_ql_e
        real(dp) :: D12_a_i, D21_ql_i, D22_a_i, D22_ql_i, D22_nc
        real(dp) :: Q_e, Q_i
        real(dp) :: expected_Q_e, expected_Q_i
        real(dp) :: D22_e, D22_i
        character(len=*), parameter :: test_name = "test_heat_fluxes"

        print *, "Running: ", test_name

        ! Set up test inputs
        n_b = 1.0e13_dp
        Te_b = 1.6022e-10_dp  ! 100 eV
        Ti_b = 1.2818e-10_dp  ! 80 eV
        Z_i_val = 1.0_dp

        ! Thermodynamic forces
        forces%e%A1_noE = -0.05_dp
        forces%e%A2 = 0.1_dp
        forces%i%A1_noE = -0.04_dp
        forces%i%A2 = 0.1_dp
        forces%e%A1 = -0.03_dp
        forces%i%A1 = -0.06_dp

        ! Diffusion coefficients (individual values)
        D12_a_e = 0.5e4_dp;  D22_a_e = 1.5e4_dp
        D21_ql_e = 1.0e4_dp; D22_ql_e = 1.5e4_dp
        D12_a_i = 0.4e4_dp;  D22_a_i = 1.0e4_dp
        D21_ql_i = 0.8e4_dp; D22_ql_i = 1.0e4_dp
        D22_nc = 0.5e4_dp

        ! Combined D22 coefficients for expected calculation
        D22_e = D22_a_e + D22_ql_e
        D22_i = D22_a_i + D22_nc + D22_ql_i

        ! Expected values: Q = -(D12_a*A_noE_1 + D21_ql*A_1 + D22*A_noE_2) * n * T
        expected_Q_e = -(D12_a_e*forces%e%A1_noE + D21_ql_e*forces%e%A1 + D22_e*forces%e%A2) * n_b * Te_b
        expected_Q_i = -(D12_a_i*forces%i%A1_noE + D21_ql_i*forces%i%A1 + D22_i*forces%i%A2) * n_b / Z_i_val * Ti_b

        ! Call the function under test
        call compute_total_heat_fluxes(forces, n_b, Te_b, Ti_b, Z_i_val, &
            D12_a_e, D21_ql_e, D22_a_e, D22_ql_e, &
            D12_a_i, D21_ql_i, D22_a_i, D22_nc, D22_ql_i, &
            Q_e, Q_i)

        ! Verify results
        call assert_equal(Q_e, expected_Q_e, "Q_e")
        call assert_equal(Q_i, expected_Q_i, "Q_i")

        print *, "  PASSED: ", test_name

    end subroutine test_heat_fluxes


    subroutine assert_equal(actual, expected, name)
        !
        ! Assert that two values are equal within tolerance
        !
        real(dp), intent(in) :: actual, expected
        character(len=*), intent(in) :: name
        real(dp), parameter :: tol = 1.0d-12
        real(dp) :: rel_diff

        if (abs(expected) > 1.0e-30_dp) then
            rel_diff = abs((actual - expected) / expected)
        else
            rel_diff = abs(actual - expected)
        end if

        if (rel_diff > tol) then
            print '(A,A,A)', "    FAIL: ", name, " mismatch"
            print '(A,E22.15)', "      Expected: ", expected
            print '(A,E22.15)', "      Actual:   ", actual
            print '(A,E22.15)', "      Rel diff: ", rel_diff
            num_failed = num_failed + 1
        else
            num_passed = num_passed + 1
        end if

    end subroutine assert_equal

end program test_rhs_balance
