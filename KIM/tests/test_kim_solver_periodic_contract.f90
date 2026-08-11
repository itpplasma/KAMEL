program test_kim_solver_periodic_contract
    ! End-to-end contract test for the periodic KIM response returned through
    ! kim_solver_t.  This is intentionally structural: it proves that the
    ! local periodic result contains the same physical field families needed
    ! by the QL-Balance adapter, without freezing a platform-dependent golden.
    use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
    use KIM_kinds_m, only: dp
    use kim_solver_m, only: kim_solver_t, kim_results_t, kim_profiles_t, KIM_OK

    integer, parameter :: npts = 40
    integer, parameter :: m_mode = -6, n_mode = 2
    type(kim_solver_t) :: kim
    type(kim_results_t) :: res
    type(kim_profiles_t) :: prof
    integer :: ierr, i
    real(dp) :: frac
    logical :: all_passed

    all_passed = .true.

    allocate(prof%r(npts), prof%n(npts), prof%Te(npts), prof%Ti(npts), &
             prof%q(npts), prof%Er(npts))
    do i = 1, npts
        prof%r(i) = 2.0_dp + 58.0_dp * real(i - 1, dp) / real(npts - 1, dp)
        frac = (prof%r(i) - 2.0_dp) / 58.0_dp
        prof%n(i) = 2.0e13_dp * (1.1_dp - prof%r(i) / 100.0_dp)
        prof%Te(i) = 1.0e3_dp * (1.2_dp - prof%r(i) / 100.0_dp)
        prof%Ti(i) = prof%Te(i)
        prof%q(i) = 1.5_dp + 2.5_dp * frac
        prof%Er(i) = -0.5_dp
    end do

    call kim%init('KIM_config_em_small.nml', run_type='electrostatic_periodic', &
                  profiles=prof, stat=ierr)
    call check('periodic init returns KIM_OK', ierr == KIM_OK, all_passed)
    if (ierr == KIM_OK) then
        call kim%solve(m=m_mode, n=n_mode, stat=ierr)
        call check('periodic solve returns KIM_OK', ierr == KIM_OK, all_passed)
    end if

    if (ierr == KIM_OK) then
    res = kim%results()
        call check('local field grid populated', allocated(res%r_field), all_passed)
        if (allocated(res%r_field)) then
            call check('local field grid has useful size', size(res%r_field) > 1, all_passed)
        end if
        call check_complex_field('Phi', res%Phi, res%r_field, all_passed)
        call check_complex_field('Br', res%Br, res%r_field, all_passed)
        call check_complex_field('Bparallel', res%Bparallel, res%r_field, all_passed)
        call check_complex_field('Er', res%Er, res%r_field, all_passed)
        call check_complex_field('Etheta', res%Etheta, res%r_field, all_passed)
        call check_complex_field('Ez', res%Ez, res%r_field, all_passed)
        call check_complex_field('Es', res%Es, res%r_field, all_passed)
        call check_complex_field('Ep', res%Ep, res%r_field, all_passed)
        call check_complex_field('jpar', res%jpar, res%r_field, all_passed)
        call check_complex_field('jpar_e', res%jpar_e, res%r_field, all_passed)
        call check_complex_field('jpar_i', res%jpar_i, res%r_field, all_passed)
        call check('result mode matches request', res%m == m_mode .and. res%n == n_mode, all_passed)
    end if

    call kim%finalize()
    if (all_passed) then
        print *, 'All periodic KIM result-contract checks PASSED'
        stop 0
    else
        print *, 'Periodic KIM result-contract checks FAILED'
        stop 1
    end if

contains

    subroutine check(name, ok, passed)
        character(*), intent(in) :: name
        logical, intent(in) :: ok
        logical, intent(inout) :: passed
        if (ok) then
            print *, 'PASS: ', name
        else
            print *, 'FAIL: ', name
            passed = .false.
        end if
    end subroutine check

    subroutine check_complex_field(name, field, grid, passed)
        character(*), intent(in) :: name
        complex(dp), allocatable, intent(in) :: field(:)
        real(dp), allocatable, intent(in) :: grid(:)
        logical, intent(inout) :: passed
        logical :: ok
        ok = allocated(field) .and. allocated(grid)
        if (ok) ok = size(field) == size(grid)
        if (ok) ok = all(ieee_is_finite(real(field, dp))) .and. &
                     all(ieee_is_finite(aimag(field)))
        call check('periodic '//trim(name)//' populated, sized, finite', ok, passed)
    end subroutine check_complex_field

end program test_kim_solver_periodic_contract
