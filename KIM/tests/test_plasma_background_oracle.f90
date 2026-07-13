program test_plasma_background_oracle
    use KIM_kinds_m, only: dp

    implicit none

    integer, parameter :: max_rows = 32
    integer :: row_kind(max_rows), row_case(max_rows), row_species(max_rows)
    real(dp) :: row_value(20, max_rows)
    integer :: nrows, failures

    failures = 0
    call read_oracle(row_kind, row_case, row_species, row_value, nrows)
    call test_fixture_inventory(row_kind, row_case, row_species, row_value, &
                                nrows, failures)
    call test_derived_background(row_kind, row_case, row_species, row_value, nrows, failures)
    call test_array_profile_quasineutrality(failures)
    call test_cell_center_recomputation(row_kind, row_case, row_species, &
                                        row_value, nrows, failures)

    if (nrows /= 17) then
        print *, 'FAIL: expected 17 plasma-background oracle rows, got ', nrows
        failures = failures + 1
    end if

    if (failures /= 0) then
        print *, 'FAIL: plasma-background checks failed: ', failures
        error stop
    end if

    print *, 'PASS: ', nrows, ' independent plasma-background oracle rows'
    print *, 'All tests PASSED'

contains

    subroutine read_oracle(kinds, cases, species, values, count)
        integer, intent(out) :: kinds(max_rows), cases(max_rows), species(max_rows)
        real(dp), intent(out) :: values(20, max_rows)
        integer, intent(out) :: count
        integer :: unit, ios, reserved
        character(len=2048) :: line

        count = 0
        open(newunit=unit, file='plasma_background.dat', status='old', &
             action='read', iostat=ios)
        if (ios /= 0) error stop 'cannot open plasma_background.dat'
        do
            read(unit, '(A)', iostat=ios) line
            if (ios < 0) exit
            if (ios > 0) error stop 'cannot read plasma-background oracle line'
            if (len_trim(line) == 0 .or. line(1:1) == '#') cycle
            count = count + 1
            if (count > max_rows) error stop 'too many plasma-background rows'
            read(line, *, iostat=ios) kinds(count), cases(count), species(count), &
                reserved, values(:, count)
            if (ios /= 0) error stop 'malformed plasma-background oracle row'
        end do
        close(unit)
    end subroutine read_oracle

    subroutine test_fixture_inventory(kinds, cases, species, values, count, failed)
        integer, intent(in) :: kinds(max_rows), cases(max_rows), species(max_rows), count
        real(dp), intent(in) :: values(20, max_rows)
        integer, intent(inout) :: failed
        integer :: expected_count(6), actual_count(6), irow, left_row, right_row

        expected_count = [3, 5, 5, 1, 2, 1]
        actual_count = 0
        left_row = 0
        right_row = 0
        do irow = 1, count
            if (kinds(irow) >= 1 .and. kinds(irow) <= 6) &
                actual_count(kinds(irow)) = actual_count(kinds(irow)) + 1
            if (kinds(irow) == 5 .and. species(irow) == 0) left_row = irow
            if (kinds(irow) == 5 .and. species(irow) == 1) right_row = irow
        end do
        if (any(actual_count /= expected_count)) then
            print *, 'FAIL: plasma-background fixture kind counts = ', actual_count
            failed = failed + 1
        end if

        if (left_row == 0 .or. right_row == 0) then
            print *, 'FAIL: periodic seam fixture endpoints are missing'
            failed = failed + 1
        else
            call expect_close('periodic seam separation', &
                values(1, right_row) - values(1, left_row), 18.0_dp, &
                5.0e-12_dp, failed)
            do irow = 2, 13
                call expect_close('periodic seam primitive/derivative', &
                    values(irow, left_row), values(irow, right_row), &
                    5.0e-12_dp, failed)
            end do
        end if

        do irow = 1, count
            if (kinds(irow) == 1 .and. cases(irow) == 2) &
                call expect_close('resonance fixture k_parallel', values(14, irow), &
                                  0.0_dp, 5.0e-12_dp, failed)
            if (kinds(irow) == 1 .and. cases(irow) == 3) &
                call expect_close('zero-rotation fixture', values(15, irow), &
                                  0.0_dp, 5.0e-12_dp, failed)
            if (kinds(irow) == 6) &
                call expect_close('equilibrium force residual', values(9, irow), &
                                  0.0_dp, 5.0e-12_dp, failed)
        end do
    end subroutine test_fixture_inventory

    subroutine test_derived_background(kinds, cases, species, values, count, failed)
        use config_m, only: number_of_ion_species, rescale_density, ion_flr_scale_factor
        use setup_m, only: collisions_off, omega, mphi_max
        use species_m, only: plasma, calculate_plasma_backs, &
                             calculate_thermodynamic_forces_and_susc
        use grid_m, only: rg_grid

        integer, intent(in) :: kinds(max_rows), cases(max_rows), species(max_rows), count
        real(dp), intent(in) :: values(20, max_rows)
        integer, intent(inout) :: failed
        integer :: icase, irow, sp, ns, mphi
        real(dp), parameter :: relative_tolerance = 5.0e-12_dp

        rescale_density = .false.
        ion_flr_scale_factor = 1.0_dp
        collisions_off = .false.
        omega = 0.0_dp
        mphi_max = 1
        rg_grid%npts_b = 1
        rg_grid%npts_c = 0

        do icase = 1, 2
            ns = 0
            do irow = 1, count
                if (kinds(irow) == 2 .and. cases(irow) == icase) &
                    ns = max(ns, species(irow) + 1)
            end do
            if (ns < 2) error stop 'oracle case has no ion species'

            call reset_global_plasma()
            number_of_ion_species = ns - 1
            plasma%n_species = ns
            plasma%grid_size = 1
            allocate(plasma%spec(0:ns - 1))
            allocate(plasma%B0(1), plasma%kp(1), plasma%om_E(1), plasma%Er(1))
            plasma%B0 = -18000.0_dp
            plasma%Er = -0.5_dp

            do irow = 1, count
                if (kinds(irow) /= 2 .or. cases(irow) /= icase) cycle
                sp = species(irow)
                allocate(plasma%spec(sp)%n(1), plasma%spec(sp)%T(1))
                allocate(plasma%spec(sp)%dndr(1), plasma%spec(sp)%dTdr(1))
                allocate(plasma%spec(sp)%Vpar(1))
                plasma%spec(sp)%n = values(1, irow)
                plasma%spec(sp)%T = values(2, irow)
                plasma%spec(sp)%Zspec = nint(values(3, irow))
                plasma%spec(sp)%Aspec = nint(values(4, irow))
                plasma%spec(sp)%mass = values(5, irow)
                plasma%spec(sp)%dndr = values(15, irow)
                plasma%spec(sp)%dTdr = values(16, irow)
                plasma%spec(sp)%Vpar = values(17, irow)
                plasma%kp = values(18, irow)
                plasma%om_E = values(19, irow)
            end do

            call calculate_plasma_backs(plasma)
            call calculate_thermodynamic_forces_and_susc(plasma)

            do irow = 1, count
                if (kinds(irow) /= 2 .or. cases(irow) /= icase) cycle
                sp = species(irow)
                mphi = nint(values(20, irow))
                call expect_close('thermal velocity', plasma%spec(sp)%vT(1), &
                    values(6, irow), relative_tolerance, failed)
                call expect_close('cyclotron frequency', plasma%spec(sp)%omega_c(1), &
                    values(7, irow), relative_tolerance, failed)
                call expect_close('Larmor radius', plasma%spec(sp)%rho_L(1), &
                    values(8, irow), relative_tolerance, failed)
                call expect_close('Debye length', plasma%spec(sp)%lambda_D(1), &
                    values(9, irow), relative_tolerance, failed)
                call expect_close('collision frequency', plasma%spec(sp)%nu(1), &
                    values(10, irow), relative_tolerance, failed)
                call expect_close('A1', plasma%spec(sp)%A1(1), &
                    values(11, irow), relative_tolerance, failed)
                call expect_close('A2', plasma%spec(sp)%A2(1), &
                    values(12, irow), relative_tolerance, failed)
                call expect_close('x1', plasma%spec(sp)%x1(1), &
                    values(13, irow), relative_tolerance, failed)
                call expect_close('x2', plasma%spec(sp)%x2(1, mphi), &
                    values(14, irow), relative_tolerance, failed)
            end do

            do irow = 1, count
                if (kinds(irow) /= 3 .or. cases(irow) /= icase) cycle
                sp = species(irow)
                call expect_close('collision summary', plasma%spec(sp)%nu(1), &
                    values(1, irow), relative_tolerance, failed)
            end do
        end do

        call reset_global_plasma()
    end subroutine test_derived_background

    subroutine test_array_profile_quasineutrality(failed)
        use config_m, only: number_of_ion_species
        use species_m, only: plasma, set_profiles_from_arrays

        integer, intent(inout) :: failed
        integer, parameter :: npts = 4
        real(dp) :: r(npts), ne(npts), te(npts), ti(npts), q(npts), er(npts)
        real(dp) :: ion_density(npts, 2), ion_temperature(npts, 2)
        real(dp) :: residual

        call reset_global_plasma()
        number_of_ion_species = 2
        plasma%n_species = 3
        allocate(plasma%spec(0:2))
        plasma%spec(0)%Zspec = -1
        plasma%spec(1)%Zspec = 1
        plasma%spec(2)%Zspec = 6
        r = [1.0_dp, 2.0_dp, 3.0_dp, 4.0_dp]
        ne = [2.0e13_dp, 2.1e13_dp, 2.2e13_dp, 2.3e13_dp]
        te = 1200.0_dp
        ti = 800.0_dp
        q = 3.0_dp
        er = 0.0_dp

        call set_profiles_from_arrays(r, ne, te, ti, q, er, npts)
        residual = maxval(abs(plasma%spec(1)%n + 6.0_dp*plasma%spec(2)%n - ne))
        if (residual > 50.0_dp*epsilon(1.0_dp)*maxval(ne)) then
            print *, 'FAIL: array profile reconstruction violates quasineutrality: ', residual
            failed = failed + 1
        end if

        ion_density(:, 2) = [1.0e12_dp, 1.1e12_dp, 1.2e12_dp, 1.3e12_dp]
        ion_density(:, 1) = ne - 6.0_dp*ion_density(:, 2)
        ion_temperature(:, 1) = [800.0_dp, 810.0_dp, 820.0_dp, 830.0_dp]
        ion_temperature(:, 2) = [500.0_dp, 510.0_dp, 520.0_dp, 530.0_dp]
        call set_profiles_from_arrays(r, ne, te, ti, q, er, npts, &
            ion_density_in=ion_density, ion_temperature_in=ion_temperature)
        residual = maxval(abs(plasma%spec(1)%n + 6.0_dp*plasma%spec(2)%n - ne))
        if (residual > 50.0_dp*epsilon(1.0_dp)*maxval(ne)) then
            print *, 'FAIL: supplied multispecies profiles violate quasineutrality: ', residual
            failed = failed + 1
        end if
        call expect_close('supplied impurity density', plasma%spec(2)%n(3), &
            ion_density(3, 2), 0.0_dp, failed)
        call expect_close('distinct main-ion temperature', plasma%spec(1)%T(2), &
            ion_temperature(2, 1), 0.0_dp, failed)
        call expect_close('distinct impurity temperature', plasma%spec(2)%T(2), &
            ion_temperature(2, 2), 0.0_dp, failed)
        call reset_global_plasma()
    end subroutine test_array_profile_quasineutrality

    subroutine test_cell_center_recomputation(kinds, cases, species, values, count, failed)
        use config_m, only: number_of_ion_species, rescale_density, ion_flr_scale_factor
        use setup_m, only: collisions_off, omega, mphi_max
        use species_m, only: plasma, calculate_plasma_backs, compute_rg_cell_centers, &
                             calculate_thermodynamic_forces_and_susc
        use grid_m, only: rg_grid

        integer, intent(in) :: kinds(max_rows), cases(max_rows), species(max_rows), count
        real(dp), intent(in) :: values(20, max_rows)
        integer, intent(inout) :: failed
        integer :: sp, irow
        real(dp) :: expected(12, 0:1)

        call reset_global_plasma()
        number_of_ion_species = 1
        rescale_density = .false.
        ion_flr_scale_factor = 1.0_dp
        collisions_off = .false.
        omega = 0.0_dp
        mphi_max = 1
        plasma%n_species = 2
        plasma%grid_size = 2
        allocate(plasma%spec(0:1))
        allocate(plasma%ks(2), plasma%kp(2), plasma%Er(2), plasma%om_E(2), plasma%B0(2))
        plasma%ks = 0.1_dp
        plasma%kp = 0.03_dp
        plasma%Er = -0.5_dp
        plasma%om_E = 20000.0_dp
        plasma%B0 = [-16000.0_dp, -20000.0_dp]
        rg_grid%npts_b = 2
        rg_grid%npts_c = 1

        expected = 0.0_dp
        do irow = 1, count
            if (kinds(irow) /= 2 .or. cases(irow) /= 1) cycle
            sp = species(irow)
            expected(:, sp) = values(6:17, irow)
            plasma%spec(sp)%mass = values(5, irow)
            plasma%spec(sp)%Zspec = nint(values(3, irow))
            plasma%spec(sp)%Aspec = nint(values(4, irow))
        end do

        do sp = 0, 1
            allocate(plasma%spec(sp)%n(2), plasma%spec(sp)%dndr(2))
            allocate(plasma%spec(sp)%T(2), plasma%spec(sp)%dTdr(2))
            allocate(plasma%spec(sp)%Vpar(2))
            plasma%spec(sp)%n = 2.0e13_dp
            plasma%spec(sp)%dndr = expected(10, sp)
            if (sp == 0) then
                plasma%spec(sp)%T = [200.0_dp, 1800.0_dp]
            else
                plasma%spec(sp)%T = [400.0_dp, 1200.0_dp]
            end if
            plasma%spec(sp)%dTdr = expected(11, sp)
            plasma%spec(sp)%Vpar = expected(12, sp)
        end do

        call calculate_plasma_backs(plasma)
        call compute_rg_cell_centers(plasma)
        call calculate_thermodynamic_forces_and_susc(plasma)
        do sp = 0, 1
            call expect_close('cell-center thermal speed', plasma%spec(sp)%vT_cc(1), &
                expected(1, sp), 5.0e-12_dp, failed)
            call expect_close('cell-center cyclotron frequency', plasma%spec(sp)%omega_c_cc(1), &
                expected(2, sp), 5.0e-12_dp, failed)
            call expect_close('cell-center Larmor radius', plasma%spec(sp)%rho_L_cc(1), &
                expected(3, sp), 5.0e-12_dp, failed)
            call expect_close('cell-center Debye length', plasma%spec(sp)%lambda_D_cc(1), &
                expected(4, sp), 5.0e-12_dp, failed)
            call expect_close('cell-center collision frequency', plasma%spec(sp)%nu_cc(1), &
                expected(5, sp), 5.0e-12_dp, failed)
            call expect_close('cell-center x1', plasma%spec(sp)%x1_cc(1), &
                expected(8, sp), 5.0e-12_dp, failed)
            call expect_close('cell-center x2', plasma%spec(sp)%x2_cc(1, 1), &
                expected(9, sp), 5.0e-12_dp, failed)
        end do
        call reset_global_plasma()
    end subroutine test_cell_center_recomputation

    subroutine reset_global_plasma()
        use species_m, only: plasma

        if (allocated(plasma%spec)) deallocate(plasma%spec)
        if (allocated(plasma%om_E)) deallocate(plasma%om_E)
        if (allocated(plasma%ks)) deallocate(plasma%ks)
        if (allocated(plasma%kp)) deallocate(plasma%kp)
        if (allocated(plasma%B0)) deallocate(plasma%B0)
        if (allocated(plasma%q)) deallocate(plasma%q)
        if (allocated(plasma%dqdr)) deallocate(plasma%dqdr)
        if (allocated(plasma%Er)) deallocate(plasma%Er)
        if (allocated(plasma%r_grid)) deallocate(plasma%r_grid)
        if (allocated(plasma%ks_cc)) deallocate(plasma%ks_cc)
        if (allocated(plasma%Er_cc)) deallocate(plasma%Er_cc)
        if (allocated(plasma%om_E_cc)) deallocate(plasma%om_E_cc)
        plasma%n_species = 0
        plasma%grid_size = 0
    end subroutine reset_global_plasma

    subroutine expect_close(label, actual, expected, tolerance, failed)
        character(len=*), intent(in) :: label
        real(dp), intent(in) :: actual, expected, tolerance
        integer, intent(inout) :: failed
        real(dp) :: error

        error = abs(actual - expected)/max(1.0_dp, abs(expected))
        if (error > tolerance) then
            print *, 'FAIL: ', trim(label)
            print *, ' actual/expected/scaled error = ', actual, expected, error
            failed = failed + 1
        end if
    end subroutine expect_close

end program test_plasma_background_oracle
