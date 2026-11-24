program test_neort_grid
    use iso_fortran_env, only: dp => real64
    use, intrinsic :: ieee_arithmetic, only: ieee_value, ieee_quiet_nan
    use neort_interface, only: apply_ntv_transport, boozer_s_range_problem, &
                               calculate_coarse_s_tor, calculate_s_tor, inner_s_of_grid, &
                               outer_s_of_grid, poloidal_flux_problem, read_equil_file, &
                               s_grid_problem
    use control_mod, only: equil_path

    implicit none

    real(dp), parameter :: tolerance = 1.0e-12_dp
    integer :: num_failed

    num_failed = 0
    call test_s_is_flux_over_boundary_flux()
    call test_equil_file_with_seven_columns()
    call test_poloidal_flux_consistency()
    call test_boozer_s_range()
    call test_s_max_ends_grid()
    call test_s_grid_range()
    call test_s_min_starts_grid()
    call test_torque_zero_outside_grid()

    if (num_failed > 0) then
        print '(A,I0)', "FAILED assertions: ", num_failed
        stop 1
    end if
    print *, "ALL TESTS PASSED"

contains

    subroutine test_s_is_flux_over_boundary_flux()
        ! The first surface of an equilibrium file is not the axis and its last
        ! surface is not the Boozer boundary: s must be psi_tor / psi_boundary.
        real(dp), dimension(3) :: psi_tor, s_tor

        psi_tor = [0.5_dp, 2.0_dp, 3.0_dp]
        call calculate_s_tor(s_tor, psi_tor, psi_tor_boundary=-4.0_dp)

        call assert_equal(s_tor(1), 0.125_dp, "s of first surface")
        call assert_equal(s_tor(2), 0.5_dp, "s of middle surface")
        call assert_equal(s_tor(3), 0.75_dp, "s of last surface")
    end subroutine test_s_is_flux_over_boundary_flux

    subroutine test_equil_file_with_seven_columns()
        ! Only r and psi_tor are needed; both the 7-column writer in
        ! equil_profiles.f90 and the 11-column preprocessor output must work.
        real(dp), dimension(:), allocatable :: r_eff, psi_tor
        integer :: iunit

        equil_path = "test_neort_grid_equil.dat"
        open (newunit=iunit, file=equil_path, status="replace", action="write")
        write (iunit, '(A)') " # nrad = 3"
        write (iunit, '(A)') " # r, q, psi_pol, psi_tor, dphi/dpsi, r_geom, volume"
        write (iunit, '(A)') " 1.0 0.9 -0.1 10.0 -9.0 1.1 5.0"
        write (iunit, '(A)') " 2.0 1.1 -0.4 40.0 -9.0 2.1 20.0"
        write (iunit, '(A)') " 3.0 1.5 -0.9 90.0 -9.0 3.1 45.0"
        close (iunit)

        call read_equil_file(r_eff=r_eff, psi_tor=psi_tor)

        call assert_true(size(r_eff) == 3, "three surfaces read")
        if (size(r_eff) == 3) then
            call assert_equal(r_eff(3), 3.0_dp, "last radius")
            call assert_equal(psi_tor(2), 40.0_dp, "middle toroidal flux")
        end if

        open (newunit=iunit, file=equil_path, status="old")
        close (iunit, status="delete")
    end subroutine test_equil_file_with_seven_columns

    subroutine test_poloidal_flux_consistency()
        ! dpsi_pol/ds is signed: the Boozer file and KAMEL's equilibrium must
        ! share the toroidal orientation and agree in magnitude.
        real(dp), dimension(3) :: neort, kamel

        neort = [-1.0e6_dp, -2.0e6_dp, -4.0e6_dp]
        kamel = neort * 1.02_dp
        call assert_true(len_trim(poloidal_flux_problem(neort, kamel)) == 0, &
                         "2 % deviation accepted")

        kamel(2) = -2.5e6_dp
        call assert_true(len_trim(poloidal_flux_problem(neort, kamel)) > 0, &
                         "25 % deviation rejected")

        kamel = -neort
        call assert_true(len_trim(poloidal_flux_problem(neort, kamel)) > 0, &
                         "opposite orientation rejected")

        ! A unit error (T m^2 against G cm^2) gives a mismatch of ~1e8, which
        ! the message must still print.
        kamel = neort * 1.0e-8_dp
        call assert_true(index(poloidal_flux_problem(neort, kamel), "1.000E+08") > 0, &
                         "large mismatch printed")

        kamel = neort
        kamel(2) = ieee_value(1.0_dp, ieee_quiet_nan)
        call assert_true(len_trim(poloidal_flux_problem(neort, kamel)) > 0, &
                         "NaN in KAMEL's derivative rejected")
        kamel = neort
        neort(3) = ieee_value(1.0_dp, ieee_quiet_nan)
        call assert_true(len_trim(poloidal_flux_problem(neort, kamel)) > 0, &
                         "NaN in NEO-RT's derivative rejected")
    end subroutine test_poloidal_flux_consistency

    subroutine test_boozer_s_range()
        ! NEO-RT limits s to the Boozer file's surfaces and reports the s it used.
        real(dp), dimension(3) :: s_requested, s_evaluated

        s_requested = [0.1_dp, 0.5_dp, 0.99_dp]
        s_evaluated = s_requested
        call assert_true(len_trim(boozer_s_range_problem(s_requested, s_evaluated)) == 0, &
                         "s inside the Boozer file accepted")

        s_evaluated(3) = 0.99975_dp
        s_requested(3) = 1.0_dp
        call assert_true(len_trim(boozer_s_range_problem(s_requested, s_evaluated)) > 0, &
                         "s beyond the last Boozer surface rejected")

        s_requested = [1.0e-4_dp, 0.5_dp, 0.99_dp]
        s_evaluated = [2.5e-4_dp, 0.5_dp, 0.99_dp]
        call assert_true(len_trim(boozer_s_range_problem(s_requested, s_evaluated)) > 0, &
                         "s inside the first Boozer surface rejected")
    end subroutine test_boozer_s_range

    subroutine test_s_max_ends_grid()
        ! The default s_max = 1 keeps the grid ending at s(rsepar), a smaller
        ! s_max moves its outer end inward with the same number of points.
        real(dp), dimension(5) :: s_tor

        call calculate_coarse_s_tor(s_tor, 0.2_dp, outer_s_of_grid(1.0_dp, 1.0_dp), 5)
        call assert_equal(s_tor(5), 1.0_dp, "default grid ends at s(rsepar) = 1")
        call calculate_coarse_s_tor(s_tor, 0.2_dp, outer_s_of_grid(0.98_dp, 1.0_dp), 5)
        call assert_equal(s_tor(5), 0.98_dp, "default grid ends at s(rsepar) < 1")

        call calculate_coarse_s_tor(s_tor, 0.2_dp, outer_s_of_grid(1.0_dp, 0.6_dp), 5)
        call assert_equal(s_tor(1), 0.2_dp, "s_max keeps the inner end")
        call assert_equal(s_tor(2), 0.3_dp, "s_max shrinks the spacing")
        call assert_equal(s_tor(5), 0.6_dp, "s_max ends the grid")
    end subroutine test_s_max_ends_grid

    subroutine test_s_grid_range()
        ! Only max(s(rmin), s_min) .. min(s(rsepar), s_max) has to lie in 0 < s <= 1.
        call assert_true(len_trim(s_grid_problem(0.0022_dp, 1.002_dp, 0.0_dp, 0.99_dp)) == 0, &
                         "s_max cuts off an rsepar beyond s = 1")
        call assert_true(len_trim(s_grid_problem(0.0022_dp, 1.002_dp, 0.0_dp, 1.0_dp)) == 0, &
                         "default s_max ends the grid at s = 1")
        call assert_true(len_trim(s_grid_problem(0.0022_dp, 0.9973_dp, 0.0_dp, 1.0_dp)) == 0, &
                         "rsepar inside the boundary")
        call assert_true(len_trim(s_grid_problem(0.0_dp, 0.9973_dp, 0.0_dp, 1.0_dp)) > 0, &
                         "s(rmin) = 0 rejected")
        call assert_true(len_trim(s_grid_problem(0.95_dp, 0.9973_dp, 0.0_dp, 0.9_dp)) > 0, &
                         "s_max below s(rmin) rejected")
        call assert_true(len_trim(s_grid_problem(0.5_dp, 0.4_dp, 0.0_dp, 1.0_dp)) > 0, &
                         "s(rsepar) below s(rmin) rejected")
    end subroutine test_s_grid_range

    subroutine test_s_min_starts_grid()
        ! The default s_min = 0 keeps the grid starting at s(rmin). A larger s_min
        ! moves the inner end outward, e.g. to the first surface of a half-grid
        ! Boozer file, and also rescues s(rmin) = 0.
        call assert_equal(inner_s_of_grid(0.0036_dp, 0.0_dp), 0.0036_dp, &
                          "default grid starts at s(rmin)")
        call assert_equal(inner_s_of_grid(0.0036_dp, 0.005_dp), 0.005_dp, &
                          "s_min starts the grid")
        call assert_equal(inner_s_of_grid(0.0036_dp, 0.001_dp), 0.0036_dp, &
                          "s_min below s(rmin) has no effect")
        call assert_true(len_trim(s_grid_problem(0.0036_dp, 0.9973_dp, 0.005_dp, 1.0_dp)) == 0, &
                         "s_min above s(rmin) accepted")
        call assert_true(len_trim(s_grid_problem(0.0_dp, 0.9973_dp, 0.005_dp, 1.0_dp)) == 0, &
                         "s_min rescues s(rmin) = 0")
        call assert_true(len_trim(s_grid_problem(0.0022_dp, 0.9973_dp, 0.95_dp, 0.9_dp)) > 0, &
                         "s_min above the outer end rejected")
    end subroutine test_s_min_starts_grid

    subroutine test_torque_zero_outside_grid()
        ! NTV torque on the balance grid: spline of the NEO-RT torque density
        ! between the first and the last NEO-RT surface, zero outside.
        use neort_lib, only: transport_data_t
        use grid_mod, only: rb, torque_density_ntv, torque_density_ntv_unscaled
        use control_mod, only: data_verbosity
        use time_evolution, only: antenna_factor_max
        use wave_code_data, only: antenna_factor

        type(transport_data_t) :: transport_data(5)
        integer :: i

        data_verbosity = 0
        antenna_factor_max = 1.0_dp
        antenna_factor = 1.0_dp
        rb = [(real(i, dp), i = 1, 8)]
        allocate (torque_density_ntv(8))
        do i = 1, 5
            transport_data(i)%torque%dVds = 10.0_dp
            transport_data(i)%torque%Tco = 20.0_dp
            transport_data(i)%torque%Tctr = 0.0_dp
            transport_data(i)%torque%Tt = 0.0_dp
        end do

        ! grid starting inside the balance grid, as with s_min > s(rmin)
        call apply_ntv_transport([2.5_dp, 3.5_dp, 4.5_dp, 5.5_dp, 6.5_dp], transport_data, &
                                 1.0_dp)
        call assert_true(all(torque_density_ntv_unscaled(1:2) == 0.0_dp), &
                         "zero torque inside the first NEO-RT surface")
        call assert_true(all(abs(torque_density_ntv_unscaled(3:6) - 2.0_dp) < tolerance), &
                         "torque density between the NEO-RT surfaces")
        call assert_true(all(torque_density_ntv_unscaled(7:8) == 0.0_dp), &
                         "zero torque beyond the last NEO-RT surface")

        ! grid starting at rmin = rb(1), the default
        call apply_ntv_transport([1.0_dp, 2.5_dp, 4.0_dp, 5.5_dp, 7.0_dp], transport_data, &
                                 1.0_dp)
        call assert_true(abs(torque_density_ntv_unscaled(1) - 2.0_dp) < tolerance, &
                         "torque kept at rmin")
    end subroutine test_torque_zero_outside_grid

    subroutine assert_true(condition, label)
        logical, intent(in) :: condition
        character(len=*), intent(in) :: label

        if (.not. condition) then
            print '(A,A)', "FAIL: ", label
            num_failed = num_failed + 1
        end if
    end subroutine assert_true

    subroutine assert_equal(actual, expected, label)
        real(dp), intent(in) :: actual, expected
        character(len=*), intent(in) :: label

        if (abs(actual - expected) > tolerance*max(abs(expected), 1.0_dp)) then
            print '(A,A)', "FAIL: ", label
            print '(A,ES22.14)', "  expected: ", expected
            print '(A,ES22.14)', "  actual:   ", actual
            num_failed = num_failed + 1
        end if
    end subroutine assert_equal

end program test_neort_grid
