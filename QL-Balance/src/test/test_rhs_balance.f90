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
                             compute_heat_fluxes, &
                             compute_total_fluxes_at_point, &
                             compute_time_derivatives, &
                             apply_boundary_conditions, &
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
    call test_total_fluxes_at_point()
    call test_time_derivatives()
    call test_boundary_conditions()

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
        call assert_equal(forces%A_noE_1e, expected_A_noE_1e, "A_noE_1e")
        call assert_equal(forces%A_noE_2e, expected_A_noE_2e, "A_noE_2e")
        call assert_equal(forces%A_noE_1i, expected_A_noE_1i, "A_noE_1i")
        call assert_equal(forces%A_noE_2i, expected_A_noE_2i, "A_noE_2i")
        call assert_equal(forces%A_1e, expected_A_1e, "A_1e")
        call assert_equal(forces%A_1i, expected_A_1i, "A_1i")

        ! Test case 2: Zero electric field
        Ercov_val = 0.0_dp
        expected_A_1e = expected_A_noE_1e  ! Should equal A_noE when E=0
        expected_A_1i = expected_A_noE_1i

        call compute_thermodynamic_forces( &
            ddr_n, ddr_Te, ddr_Ti, &
            n_b, Te_b, Ti_b, &
            Ercov_val, Z_i_val, &
            forces)

        call assert_equal(forces%A_1e, expected_A_1e, "A_1e (E=0)")
        call assert_equal(forces%A_1i, expected_A_1i, "A_1i (E=0)")

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
        forces%A_noE_1e = -0.05_dp
        forces%A_noE_2e = 0.1_dp
        forces%A_noE_1i = -0.04_dp
        forces%A_noE_2i = 0.1_dp
        forces%A_1e = -0.03_dp
        forces%A_1i = -0.06_dp

        ! Diffusion coefficients (individual values)
        D11_a_e = 1.0e4_dp;  D12_a_e = 0.5e4_dp
        D11_ql_e = 2.0e4_dp; D12_ql_e = 1.0e4_dp
        D11_a_i = 0.8e4_dp;  D12_a_i = 0.4e4_dp
        D11_ql_i = 1.5e4_dp; D12_ql_i = 0.8e4_dp

        ! Expected values (analytical calculation)
        expected_gamma_a_e = -(D11_a_e*forces%A_noE_1e + D12_a_e*forces%A_noE_2e) * n_b
        expected_gamma_ql_e = -(D11_ql_e*forces%A_1e + D12_ql_e*forces%A_noE_2e) * n_b
        expected_gamma_e = expected_gamma_a_e + expected_gamma_ql_e

        expected_gamma_a_i = -(D11_a_i*forces%A_noE_1i + D12_a_i*forces%A_noE_2i) * n_b / Z_i_val
        expected_gamma_ql_i = -(D11_ql_i*forces%A_1i + D12_ql_i*forces%A_noE_2i) * n_b / Z_i_val
        expected_gamma_i = expected_gamma_a_i + expected_gamma_ql_i

        ! Call the function under test
        call compute_particle_fluxes(forces, n_b, Z_i_val, &
            D11_a_e, D12_a_e, D11_ql_e, D12_ql_e, &
            D11_a_i, D12_a_i, D11_ql_i, D12_ql_i, &
            fluxes)

        ! Verify results
        call assert_equal(fluxes%e%gamma_a, expected_gamma_a_e, "gamma_a_e")
        call assert_equal(fluxes%e%gamma_ql, expected_gamma_ql_e, "gamma_ql_e")
        call assert_equal(fluxes%e%gamma, expected_gamma_e, "gamma_e")
        call assert_equal(fluxes%i%gamma_a, expected_gamma_a_i, "gamma_a_i")
        call assert_equal(fluxes%i%gamma_ql, expected_gamma_ql_i, "gamma_ql_i")
        call assert_equal(fluxes%i%gamma, expected_gamma_i, "gamma_i")

        print *, "  PASSED: ", test_name

    end subroutine test_particle_fluxes


    subroutine test_heat_fluxes()
        !
        ! Test compute_heat_fluxes with known inputs
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
        forces%A_noE_1e = -0.05_dp
        forces%A_noE_2e = 0.1_dp
        forces%A_noE_1i = -0.04_dp
        forces%A_noE_2i = 0.1_dp
        forces%A_1e = -0.03_dp
        forces%A_1i = -0.06_dp

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
        expected_Q_e = -(D12_a_e*forces%A_noE_1e + D21_ql_e*forces%A_1e + D22_e*forces%A_noE_2e) * n_b * Te_b
        expected_Q_i = -(D12_a_i*forces%A_noE_1i + D21_ql_i*forces%A_1i + D22_i*forces%A_noE_2i) * n_b / Z_i_val * Ti_b

        ! Call the function under test
        call compute_heat_fluxes(forces, n_b, Te_b, Ti_b, Z_i_val, &
            D12_a_e, D21_ql_e, D22_e, &
            D12_a_i, D21_ql_i, D22_i, &
            Q_e, Q_i)

        ! Verify results
        call assert_equal(Q_e, expected_Q_e, "Q_e")
        call assert_equal(Q_i, expected_Q_i, "Q_i")

        print *, "  PASSED: ", test_name

    end subroutine test_heat_fluxes


    subroutine test_total_fluxes_at_point()
        !
        ! Test compute_total_fluxes_at_point with known inputs
        !
        ! Physics: Computes diffusive, convective, and nonlinear fluxes
        ! for all 4 balance equations (n, vphi, Te, Ti)
        !
        real(dp) :: gamma_e, gamma_e_nl, gamma_i
        real(dp) :: Q_e, Q_e_nl, Q_i, Q_i_nl
        real(dp) :: ddr_n, ddr_n_nl, ddr_Te, ddr_Te_nl, ddr_Ti, ddr_Ti_nl, ddr_vphi
        real(dp) :: n_b, Te_b, Ti_b, Z_i_val, Sb_val
        real(dp) :: dae11_val, dqle11_val, dae22_val, dqle22_val
        real(dp) :: dai22_val, dni22_val, dqli22_val, dqli21_val
        real(dp) :: visca_val, gpp_av_val
        real(dp) :: flux_dif(4), flux_con(4), flux_con_nl(4)
        real(dp) :: expected_flux_dif(4), expected_flux_con(4), expected_flux_con_nl(4)
        real(dp) :: dfluxvphi
        character(len=*), parameter :: test_name = "test_total_fluxes_at_point"

        print *, "Running: ", test_name

        ! Set up test inputs (representative values)
        n_b = 1.0e13_dp
        Te_b = 1.6022e-10_dp  ! 100 eV
        Ti_b = 1.2818e-10_dp  ! 80 eV
        Z_i_val = 1.0_dp
        Sb_val = 1.0e3_dp     ! Surface element

        gamma_e = 1.0e15_dp
        gamma_e_nl = 1.2e15_dp
        gamma_i = 0.9e15_dp
        Q_e = 5.0e3_dp
        Q_e_nl = 5.5e3_dp
        Q_i = 4.0e3_dp
        Q_i_nl = 4.5e3_dp

        ddr_n = 1.0e12_dp
        ddr_n_nl = 1.1e12_dp
        ddr_Te = 1.6022e-11_dp
        ddr_Te_nl = 1.7e-11_dp
        ddr_Ti = 1.2818e-11_dp
        ddr_Ti_nl = 1.4e-11_dp
        ddr_vphi = 1.0e4_dp

        dae11_val = 1.0e4_dp
        dqle11_val = 2.0e4_dp
        dae22_val = 3.0e4_dp
        dqle22_val = 4.0e4_dp
        dai22_val = 2.5e4_dp
        dni22_val = 0.5e4_dp
        dqli22_val = 1.5e4_dp
        dqli21_val = 1.0e4_dp
        visca_val = 1.0e5_dp
        gpp_av_val = 1.0_dp

        ! Expected values (analytical calculation)
        ! Equation 1: Particle flux
        expected_flux_dif(1) = -Sb_val*ddr_n*(dae11_val &
                               + dqle11_val*(1._dp + Ti_b/Te_b/Z_i_val))
        expected_flux_con(1) = (Sb_val*gamma_e - expected_flux_dif(1))/n_b
        expected_flux_con_nl(1) = (Sb_val*gamma_e_nl - &
                                  (-Sb_val*ddr_n_nl*(dae11_val &
                                  + dqle11_val*(1._dp + Ti_b/Te_b/Z_i_val))))/n_b

        ! Equation 2: Toroidal momentum flux
        dfluxvphi = -visca_val*ddr_vphi*n_b/Z_i_val*gpp_av_val
        expected_flux_dif(2) = Sb_val*dfluxvphi
        expected_flux_con(2) = 0._dp
        expected_flux_con_nl(2) = 0._dp

        ! Equation 3: Electron heat flux
        expected_flux_dif(3) = -Sb_val*(dae22_val + dqle22_val)*n_b*ddr_Te
        expected_flux_con(3) = (Sb_val*Q_e - expected_flux_dif(3))/Te_b
        expected_flux_con_nl(3) = (Sb_val*Q_e_nl - &
                                  (-Sb_val*(dae22_val + dqle22_val)*n_b*ddr_Te_nl))/Te_b

        ! Equation 4: Ion heat flux
        expected_flux_dif(4) = -Sb_val*(dai22_val + dni22_val + dqli22_val &
                               - 2.5_dp*dqli21_val)*n_b/Z_i_val*ddr_Ti
        expected_flux_con(4) = (Sb_val*Q_i - expected_flux_dif(4))/Ti_b
        expected_flux_con_nl(4) = (Sb_val*Q_i_nl - &
                                  (-Sb_val*(dai22_val + dni22_val + dqli22_val - 2.5_dp*dqli21_val) &
                                  *n_b/Z_i_val*ddr_Ti_nl))/Ti_b

        ! Call the function under test
        call compute_total_fluxes_at_point( &
            gamma_e, gamma_e_nl, gamma_i, Q_e, Q_e_nl, Q_i, Q_i_nl, &
            ddr_n, ddr_n_nl, ddr_Te, ddr_Te_nl, ddr_Ti, ddr_Ti_nl, ddr_vphi, &
            n_b, Te_b, Ti_b, Z_i_val, Sb_val, &
            dae11_val, dqle11_val, dae22_val, dqle22_val, &
            dai22_val, dni22_val, dqli22_val, dqli21_val, &
            visca_val, gpp_av_val, &
            flux_dif, flux_con, flux_con_nl)

        ! Verify results
        call assert_equal(flux_dif(1), expected_flux_dif(1), "flux_dif(1)")
        call assert_equal(flux_con(1), expected_flux_con(1), "flux_con(1)")
        call assert_equal(flux_con_nl(1), expected_flux_con_nl(1), "flux_con_nl(1)")
        call assert_equal(flux_dif(2), expected_flux_dif(2), "flux_dif(2)")
        call assert_equal(flux_con(2), expected_flux_con(2), "flux_con(2)")
        call assert_equal(flux_con_nl(2), expected_flux_con_nl(2), "flux_con_nl(2)")
        call assert_equal(flux_dif(3), expected_flux_dif(3), "flux_dif(3)")
        call assert_equal(flux_con(3), expected_flux_con(3), "flux_con(3)")
        call assert_equal(flux_con_nl(3), expected_flux_con_nl(3), "flux_con_nl(3)")
        call assert_equal(flux_dif(4), expected_flux_dif(4), "flux_dif(4)")
        call assert_equal(flux_con(4), expected_flux_con(4), "flux_con(4)")
        call assert_equal(flux_con_nl(4), expected_flux_con_nl(4), "flux_con_nl(4)")

        print *, "  PASSED: ", test_name

    end subroutine test_total_fluxes_at_point


    subroutine test_time_derivatives()
        !
        ! Test compute_time_derivatives with known inputs
        !
        ! Physics:
        !   dot_params(2) -> d(omega)/dt = dot_params(2)*Z_i/n*2/gpp_av
        !   dot_params(3) -> d(Te)/dt = (-Te*dot_n + d(nTe)/dt/1.5)/n
        !   dot_params(4) -> d(Ti)/dt = (-Ti*dot_n + d(nTi)/dt/1.5)/n
        !
        real(dp) :: dot_params_in(4), dot_params_out(4), expected_out(4)
        real(dp) :: n_val, Te_val, Ti_val, Z_i_val, gpp_av_avg
        character(len=*), parameter :: test_name = "test_time_derivatives"

        print *, "Running: ", test_name

        ! Set up test inputs
        n_val = 1.0e13_dp
        Te_val = 1.6022e-10_dp
        Ti_val = 1.2818e-10_dp
        Z_i_val = 1.0_dp
        gpp_av_avg = 2.0_dp

        ! Raw time derivatives (arbitrary test values)
        dot_params_in(1) = 1.0e10_dp     ! d(n)/dt
        dot_params_in(2) = 5.0e8_dp      ! d(L)/dt (momentum)
        dot_params_in(3) = 2.0e-2_dp     ! d(nTe)/dt
        dot_params_in(4) = 1.5e-2_dp     ! d(nTi)/dt

        ! Expected values (analytical calculation)
        ! Density - no conversion
        expected_out(1) = dot_params_in(1)

        ! Momentum -> rotation frequency
        expected_out(2) = dot_params_in(2)*Z_i_val/n_val*2._dp/gpp_av_avg

        ! d(nTe)/dt -> d(Te)/dt
        expected_out(3) = (-Te_val*dot_params_in(1) + dot_params_in(3)/1.5_dp)/n_val

        ! d(nTi)/dt -> d(Ti)/dt
        expected_out(4) = (-Ti_val*dot_params_in(1) + dot_params_in(4)/1.5_dp)/n_val

        ! Call the function under test
        call compute_time_derivatives( &
            dot_params_in, n_val, Te_val, Ti_val, Z_i_val, gpp_av_avg, &
            dot_params_out)

        ! Verify results
        call assert_equal(dot_params_out(1), expected_out(1), "dot_n")
        call assert_equal(dot_params_out(2), expected_out(2), "dot_omega")
        call assert_equal(dot_params_out(3), expected_out(3), "dot_Te")
        call assert_equal(dot_params_out(4), expected_out(4), "dot_Ti")

        print *, "  PASSED: ", test_name

    end subroutine test_time_derivatives


    subroutine test_boundary_conditions()
        !
        ! Test that boundary conditions are correctly applied
        ! (zero flux at inner boundary)
        !
        character(len=*), parameter :: test_name = "test_boundary_conditions"

        integer, parameter :: nbaleqs = 4
        integer, parameter :: npoints = 3
        integer :: ieq, ipoi
        real(dp) :: fluxes_dif(nbaleqs, npoints)
        real(dp) :: fluxes_con(nbaleqs, npoints)
        real(dp) :: fluxes_con_nl(nbaleqs, npoints)
        real(dp) :: expected_dif(nbaleqs, npoints)
        real(dp) :: expected_con(nbaleqs, npoints)
        real(dp) :: expected_con_nl(nbaleqs, npoints)

        print *, "Running: ", test_name

        do ieq = 1, nbaleqs
            do ipoi = 1, npoints
                fluxes_dif(ieq, ipoi) = real(ieq * ipoi, dp)
                fluxes_con(ieq, ipoi) = -real(ieq * ipoi, dp)
                fluxes_con_nl(ieq, ipoi) = real(ieq + ipoi, dp)
            end do
        end do

        expected_dif = fluxes_dif
        expected_con = fluxes_con
        expected_con_nl = fluxes_con_nl
        expected_dif(:, 1) = 0.0_dp
        expected_con(:, 1) = 0.0_dp
        expected_con_nl(:, 1) = 0.0_dp

        call apply_boundary_conditions(fluxes_dif, fluxes_con, fluxes_con_nl, nbaleqs)

        do ieq = 1, nbaleqs
            do ipoi = 1, npoints
                call assert_equal(fluxes_dif(ieq, ipoi), expected_dif(ieq, ipoi), "fluxes_dif")
                call assert_equal(fluxes_con(ieq, ipoi), expected_con(ieq, ipoi), "fluxes_con")
                call assert_equal(fluxes_con_nl(ieq, ipoi), expected_con_nl(ieq, ipoi), "fluxes_con_nl")
            end do
        end do

        print *, "  PASSED: ", test_name

    end subroutine test_boundary_conditions


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
