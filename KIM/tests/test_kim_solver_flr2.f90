program test_kim_solver_flr2
    use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
    use KIM_kinds_m, only: dp
    use config_m, only: flr2_ion_flr
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
    call check_complex_field('jpar', result%jpar, result%r_field, passed)
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
