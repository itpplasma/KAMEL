program test_kim_solver_flr2
    use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
    use KIM_kinds_m, only: dp
    use config_m, only: flr2_ion_flr, output_path, write_diagnostics_dat
    use fields_m, only: EBdat
    use kim_diagnostics_m, only: integrate_Ipar_window
    use kim_solver_m, only: kim_profiles_t, kim_results_t, kim_solver_t, KIM_OK
    implicit none

    integer, parameter :: n_profile = 48
    ! Deliberately differ from the configured m=-6 to verify that the first
    ! in-process solve rebuilds its mode-dependent KIM background.
    integer, parameter :: m_mode = -5, n_mode = 2
    type(kim_solver_t) :: kim
    type(kim_results_t) :: result
    type(kim_profiles_t) :: profiles
    integer :: i, stat
    real(dp) :: fraction
    complex(dp) :: expected_phi_boundary, factor_of_phi
    logical :: passed
    character(len=512) :: argument

    call get_environment_variable('KIM_EARLY_FAILURE', argument)
    if (trim(argument) == '1') then
        call run_flr2_early_failure_child()
        error stop 'FLR2 early-failure child unexpectedly returned'
    end if

    passed = .true.
    allocate(profiles%r(n_profile), profiles%n(n_profile))
    allocate(profiles%Te(n_profile), profiles%Ti(n_profile))
    allocate(profiles%q(n_profile), profiles%Er(n_profile))
    do i = 1, n_profile
        profiles%r(i) = 3.0_dp + 64.0_dp*real(i - 1, dp)/real(n_profile - 1, dp)
        fraction = (profiles%r(i) - 3.0_dp)/64.0_dp
        profiles%n(i) = 5.0e13_dp*(1.0_dp - 0.8_dp*fraction)
        profiles%Te(i) = 250.0_dp + 1250.0_dp*(1.0_dp - fraction)
        profiles%Ti(i) = 0.8_dp*profiles%Te(i)
        profiles%q(i) = 1.2_dp + 3.4_dp*fraction
        profiles%Er(i) = 25.0_dp + 5.0_dp*fraction
    end do

    call kim%init('KIM_config_flr2_small.nml', run_type='flr2', &
                  profiles=profiles, stat=stat)
    call check('init returns KIM_OK', stat == KIM_OK, passed)
    call check('optional KIM_FLR2 group was read', .not. flr2_ion_flr, passed)
    if (stat /= KIM_OK) call done(passed)

    write_diagnostics_dat = .true.

    call kim%solve(m_mode, n_mode, stat)
    call check('solve returns KIM_OK', stat == KIM_OK, passed)
    if (stat /= KIM_OK) call done(passed)
    result = kim%results()

    call check('field grid populated', allocated(result%r_field), passed)
    if (allocated(result%r_field)) then
        call check('field grid matches KIM plasma grid', &
                   allocated(result%r_plasma) &
                   .and. size(result%r_field) == size(result%r_plasma) &
                   .and. maxval(abs(result%r_field - result%r_plasma)) < 1.0e-12_dp, &
                   passed)
        call check('field grid strictly increasing', &
                   all(result%r_field(2:) > result%r_field(:size(result%r_field) - 1)), &
                   passed)
    end if

    call check_complex_field('Br', result%Br, result%r_field, passed)
    call check_complex_field('Phi', result%Phi, result%r_field, passed)
    call check_complex_field('Es', result%Es, result%r_field, passed)
    call check_complex_field('jpar', result%jpar, result%r_field, passed)
    call check_complex_field('jpar_e', result%jpar_e, result%r_field, passed)
    call check_complex_field('jpar_i', result%jpar_i, result%r_field, passed)
    if (allocated(result%jpar) .and. allocated(result%jpar_e) .and. allocated(result%jpar_i)) then
        call check('FLR2 species currents sum to total', &
                   maxval(abs(result%jpar-result%jpar_e-result%jpar_i)) &
                   /max(maxval(abs(result%jpar)), tiny(1.0_dp)) < 1.0e-10_dp, passed)
    end if
    if (allocated(result%Br)) then
        call check('constant Br drive preserved', &
                   maxval(abs(result%Br - cmplx(1.0_dp, 0.0_dp, dp))) < 1.0e-12_dp, &
                   passed)
    end if
    if (allocated(result%Phi)) then
        call check('FLR2 potential is non-trivial', maxval(abs(result%Phi)) > 0.0_dp, passed)
    end if
    call check('signed B0z background populated', allocated(result%B0z), passed)
    if (allocated(result%Phi) .and. allocated(result%B0z)) then
        factor_of_phi = cmplx(0.0_dp, 1.0_dp, dp) &
                        *real(m_mode + n_mode*profiles%q(1), dp) &
                        /(profiles%q(1)*(-profiles%Er(1)))
        expected_phi_boundary = -(cmplx(1.0_dp, 0.0_dp, dp)*165.0_dp/result%B0z(1)) &
                                /factor_of_phi
        call check('left boundary preserves signed B0z mapping', &
                   abs(result%Phi(1) - expected_phi_boundary) &
                   < 1.0e-10_dp*max(1.0_dp, abs(expected_phi_boundary)), passed)
    end if
    if (allocated(result%jpar)) then
        call check('FLR2 current is non-trivial', maxval(abs(result%jpar)) > 0.0_dp, passed)
    end if
    call check('FLR2 diagnostics file exists', &
               file_exists(trim(output_path)//'kim_diagnostics.dat'), passed)
    call check('FLR2 principal coefficient file exists', &
               file_exists(trim(output_path)//'principal_coefficient.dat'), passed)
    if (allocated(EBdat%r_grid) .and. allocated(EBdat%jpar_i)) then
        call check_flr2_output_files(trim(output_path), EBdat, passed)
    end if
    call check('results retain requested mode', &
               result%m == m_mode .and. result%n == n_mode, passed)
    if (allocated(result%kp) .and. allocated(result%B0) &
        .and. allocated(result%B0z) .and. allocated(result%B0th)) then
        call check('first solve background uses requested mode', &
                   abs(result%kp(1) &
                       - (real(m_mode, dp)/result%r_plasma(1) &
                          *result%B0th(1)/result%B0(1) &
                          + real(n_mode, dp)/165.0_dp &
                            *result%B0z(1)/result%B0(1))) < 1.0e-12_dp, passed)
    end if

    call assert_flr2_early_failure_clears_file()

    call kim%solve(m_mode - 1, n_mode, stat)
    call check('second mode solve returns KIM_OK', stat == KIM_OK, passed)
    if (stat == KIM_OK) then
        result = kim%results()
        call check('second solve background uses requested mode', &
                   allocated(result%kp) .and. allocated(result%B0) &
                   .and. allocated(result%B0z) .and. allocated(result%B0th) &
                   .and. abs(result%kp(1) &
                       - (real(m_mode - 1, dp)/result%r_plasma(1) &
                          *result%B0th(1)/result%B0(1) &
                          + real(n_mode, dp)/165.0_dp &
                            *result%B0z(1)/result%B0(1))) < 1.0e-12_dp, passed)
    end if

    call kim%finalize()
    call done(passed)

contains

    subroutine assert_flr2_early_failure_clears_file()
        character(len=512) :: executable, command, out_file, principal_file
        integer :: child_status, iunit
        logical :: exists, principal_exists

        out_file = trim(output_path)//'kim_diagnostics.dat'
        principal_file = trim(output_path)//'principal_coefficient.dat'
        open(newunit=iunit, file=trim(out_file), status='replace', action='write')
        write(iunit, '(A)') '# seeded successful diagnostics'
        write(iunit, '(A)') '1 2 3 4 5 6'
        close(iunit)
        open(newunit=iunit, file=trim(principal_file), status='replace', action='write')
        write(iunit, '(A)') '# seeded successful principal coefficient'
        write(iunit, '(A)') '1 2 3'
        close(iunit)

        call get_command_argument(0, executable)
        command = 'KIM_EARLY_FAILURE=1 "'//trim(executable)//'"'
        call execute_command_line(trim(command), wait=.true., exitstat=child_status)
        if (child_status == 0) error stop 'FLR2 early-failure child unexpectedly succeeded'

        inquire(file=trim(out_file), exist=exists)
        if (exists) error stop 'FLR2 early failure left stale diagnostics file'
        inquire(file=trim(principal_file), exist=principal_exists)
        if (principal_exists) error stop 'FLR2 early failure left stale principal coefficient file'
        print *, 'PASS: FLR2 early failure clears stale diagnostics and principal coefficient'
    end subroutine assert_flr2_early_failure_clears_file

    subroutine run_flr2_early_failure_child()
        use setup_m, only: type_br_field

        type(kim_solver_t) :: child_kim
        type(kim_profiles_t) :: child_profiles
        integer :: i, child_stat
        real(dp) :: fraction

        allocate(child_profiles%r(n_profile), child_profiles%n(n_profile))
        allocate(child_profiles%Te(n_profile), child_profiles%Ti(n_profile))
        allocate(child_profiles%q(n_profile), child_profiles%Er(n_profile))
        do i = 1, n_profile
            child_profiles%r(i) = 3.0_dp + 64.0_dp*real(i - 1, dp)/real(n_profile - 1, dp)
            fraction = (child_profiles%r(i) - 3.0_dp)/64.0_dp
            child_profiles%n(i) = 5.0e13_dp*(1.0_dp - 0.8_dp*fraction)
            child_profiles%Te(i) = 250.0_dp + 1250.0_dp*(1.0_dp - fraction)
            child_profiles%Ti(i) = 0.8_dp*child_profiles%Te(i)
            child_profiles%q(i) = 1.2_dp + 3.4_dp*fraction
            child_profiles%Er(i) = 25.0_dp + 5.0_dp*fraction
        end do
        call child_kim%init('KIM_config_flr2_small.nml', run_type='flr2', &
                            profiles=child_profiles, stat=child_stat)
        if (child_stat /= KIM_OK) error stop 'FLR2 early-failure child init failed'
        write_diagnostics_dat = .false.
        type_br_field = 99
        call child_kim%solve(m_mode, n_mode, child_stat)
        error stop 'FLR2 early-failure child reached solve return'
    end subroutine run_flr2_early_failure_child

    logical function file_exists(path)
        character(*), intent(in) :: path
        inquire(file=trim(path), exist=file_exists)
    end function file_exists

    subroutine check_flr2_output_files(path, fields, all_passed)
        use fields_m, only: EBdat_t

        character(*), intent(in) :: path
        type(EBdat_t), intent(in) :: fields
        logical, intent(inout) :: all_passed

        character(len=512) :: header
        real(dp) :: row(3), diagnostics(6), extra, residual
        complex(dp) :: Ipar, Ipar_e, Ipar_i
        integer :: iunit, ios, data_ios, extra_ios, i
        character(len=512) :: data_line
        logical :: valid

        valid = .true.
        open(newunit=iunit, file=trim(path)//'principal_coefficient.dat', &
             status='old', action='read', iostat=ios)
        if (ios == 0) then
            read(iunit, '(A)', iostat=ios) header
            do i = 1, size(fields%r_grid)
                read(iunit, *, iostat=ios) row
                if (ios /= 0 .or. abs(row(1)-fields%r_grid(i)) > 1.0e-12_dp .or. &
                    .not. ieee_is_finite(row(1)) .or. .not. ieee_is_finite(row(2)) .or. &
                    .not. ieee_is_finite(row(3))) then
                    valid = .false.
                    exit
                end if
            end do
            read(iunit, '(A)', iostat=ios) data_line
            if (ios == 0) valid = .false.
            close(iunit)
        else
            valid = .false.
        end if
        call check('principal coefficient follows the Phi grid and is finite', valid, all_passed)

        valid = .false.
        open(newunit=iunit, file=trim(path)//'kim_diagnostics.dat', &
             status='old', action='read', iostat=ios)
        if (ios == 0) then
            read(iunit, '(A)', iostat=ios) header
            read(iunit, '(A)', iostat=ios) data_line
            valid = ios == 0
            if (ios == 0) then
                read(data_line, *, iostat=data_ios) diagnostics
                valid = data_ios == 0 .and. all(ieee_is_finite(diagnostics))
            end if
            if (valid) then
                read(data_line, *, iostat=data_ios) diagnostics, extra
                valid = data_ios /= 0
            end if
            read(iunit, '(A)', iostat=extra_ios) data_line
            valid = valid .and. extra_ios /= 0
            close(iunit)
        end if
        if (valid) then
            Ipar = cmplx(diagnostics(2), diagnostics(3), dp)
            Ipar_e = cmplx(diagnostics(4), diagnostics(5), dp)
            Ipar_i = integrate_Ipar_window(size(fields%r_grid), fields%r_grid, fields%jpar_i, &
                                            fields%r_resonance, fields%dx_asis)
            call check('FLR2 dqle22 is positive', diagnostics(1) > 0.0_dp, all_passed)
            call check('FLR2 br_abs_res is positive', diagnostics(6) > 0.0_dp, all_passed)
            call check('FLR2 integrated total current is nonzero', abs(Ipar) > 0.0_dp, all_passed)
            call check('FLR2 integrated electron current is nonzero', &
                       abs(Ipar_e) > 0.0_dp, all_passed)
            if (abs(Ipar) > 0.0_dp) then
                residual = abs(Ipar - Ipar_e - Ipar_i) / abs(Ipar)
            else
                residual = huge(1.0_dp)
            end if
            call check('FLR2 diagnostics currents satisfy total=e+i', &
                       residual < 1.0e-10_dp, all_passed)
        end if
        call check('FLR2 diagnostics file has six finite values', valid, all_passed)
    end subroutine check_flr2_output_files

    subroutine check(label, condition, all_passed)
        character(*), intent(in) :: label
        logical, intent(in) :: condition
        logical, intent(inout) :: all_passed

        if (condition) then
            write(*,*) 'PASS: ', label
        else
            write(*,*) 'FAIL: ', label
            all_passed = .false.
        end if
    end subroutine check

    subroutine check_complex_field(label, values, grid, all_passed)
        character(*), intent(in) :: label
        complex(dp), allocatable, intent(in) :: values(:)
        real(dp), allocatable, intent(in) :: grid(:)
        logical, intent(inout) :: all_passed
        logical :: valid

        valid = allocated(values) .and. allocated(grid)
        if (valid) valid = size(values) == size(grid)
        if (valid) then
            valid = all(ieee_is_finite(real(values))) &
                    .and. all(ieee_is_finite(aimag(values)))
        end if
        call check(trim(label)//' populated, sized, and finite', valid, all_passed)
    end subroutine check_complex_field

    subroutine done(all_passed)
        logical, intent(in) :: all_passed

        if (all_passed) then
            write(*,*) 'KIM FLR2 end-to-end test PASSED'
            stop 0
        else
            write(*,*) 'KIM FLR2 end-to-end test FAILED'
            stop 1
        end if
    end subroutine done

end program test_kim_solver_flr2
