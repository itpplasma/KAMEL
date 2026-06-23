program test_kim_solver_em
    !> End-to-end keystone test for the kim_solver_m seam: drive a real
    !> electromagnetic solve through init -> solve(m,n) -> results on a small
    !> bench config with in-memory profiles, and assert the run path produces
    !> sane output (populated, finite, monotonic field grid; non-trivial Br).
    !>
    !> This is the deferred Phase-1 happy-path test. It exercises the EM run
    !> path (kernel.f90, fields_mod, electromagnetic_solver -- 0% covered) and
    !> is the safety net the QL-Balance adapter migration (Phase 2) routes
    !> through. Assertions are structural, not a hardcoded golden value: a
    !> tight numeric golden for a physics solve is platform-fragile.
    use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
    use KIM_kinds_m, only: dp
    use kim_solver_m, only: kim_solver_t, kim_results_t, kim_profiles_t, KIM_OK

    implicit none

    integer, parameter :: npts = 40
    integer, parameter :: m_mode = -6, n_mode = 2

    type(kim_solver_t) :: kim
    type(kim_results_t) :: res
    type(kim_profiles_t) :: prof
    integer :: ierr, i
    integer(8) :: t0, t1, rate
    real(dp) :: frac
    logical :: all_passed

    all_passed = .true.

    ! In-memory profiles over the configured domain [r_min, r_plas] = [3, 67] cm,
    ! in CGS (n in 1/cm^3, T in eV). q rises 1 -> 4, crossing the m/n = 3
    ! resonance mid-domain so it sits away from the boundaries.
    allocate(prof%r(npts), prof%n(npts), prof%Te(npts), &
             prof%Ti(npts), prof%q(npts), prof%Er(npts))
    do i = 1, npts
        prof%r(i)  = 3.0_dp + 64.0_dp * real(i - 1, dp) / real(npts - 1, dp)
        frac       = (prof%r(i) - 3.0_dp) / 64.0_dp
        prof%n(i)  = 5.0e13_dp * (1.0_dp - 0.9_dp * frac)
        prof%Te(i) = 100.0_dp + 1900.0_dp * (1.0_dp - frac)
        prof%Ti(i) = 0.9_dp * prof%Te(i)
        prof%q(i)  = 1.0_dp + 3.0_dp * frac
        prof%Er(i) = 0.0_dp
    end do

    call kim%init('KIM_config_em_small.nml', run_type='electromagnetic', &
                  profiles=prof, stat=ierr)
    call check('init returns KIM_OK', ierr == KIM_OK, all_passed)
    if (ierr /= KIM_OK) call done(all_passed)

    call system_clock(t0, rate)
    call kim%solve(m=m_mode, n=n_mode, stat=ierr)
    call system_clock(t1)
    print '(A,F8.2,A)', ' solve wall time: ', real(t1 - t0, dp) / real(rate, dp), ' s'
    call check('solve returns KIM_OK', ierr == KIM_OK, all_passed)
    if (ierr /= KIM_OK) call done(all_passed)

    res = kim%results()

    call check('field grid populated', &
               allocated(res%r_field) .and. size(res%r_field) > 1, all_passed)
    if (allocated(res%r_field)) then
        call check('field grid strictly increasing', &
                   all(res%r_field(2:) > res%r_field(:size(res%r_field) - 1)), all_passed)
    end if

    call check('Br populated, matches field grid', &
               allocated(res%Br) .and. allocated(res%r_field) .and. &
               size(res%Br) == size(res%r_field), all_passed)
    if (allocated(res%Br)) then
        call check('Br all finite', &
                   all(ieee_is_finite(real(res%Br))) .and. &
                   all(ieee_is_finite(aimag(res%Br))), all_passed)
        ! A degenerate solve would return only the boundary value everywhere;
        ! a real interior response varies across the grid.
        call check('Br varies across grid (non-degenerate interior)', &
                   maxval(abs(res%Br)) > minval(abs(res%Br)), all_passed)
        print '(A,ES12.4,A,ES12.4)', ' |Br| range: min = ', &
            minval(abs(res%Br)), '  max = ', maxval(abs(res%Br))
    end if

    call check('results mode == requested (m,n)', &
               res%m == m_mode .and. res%n == n_mode, all_passed)

    ! Background quantities the QL-Balance adapter reads from results() --
    ! guard them so the Phase-2 migration's background reads stay covered.
    call check('plasma grid populated', &
               allocated(res%r_plasma) .and. size(res%r_plasma) > 1, all_passed)
    if (allocated(res%r_plasma)) then
        call check_bg('kp',   res%kp,   size(res%r_plasma), all_passed)
        call check_bg('ks',   res%ks,   size(res%r_plasma), all_passed)
        call check_bg('om_E', res%om_E, size(res%r_plasma), all_passed)
        call check_bg('nu_e', res%nu_e, size(res%r_plasma), all_passed)
        call check_bg('nu_i', res%nu_i, size(res%r_plasma), all_passed)
        call check_bg('B0',   res%B0,   size(res%r_plasma), all_passed)
        call check_bg('B0z',  res%B0z,  size(res%r_plasma), all_passed)
        call check_bg('B0th', res%B0th, size(res%r_plasma), all_passed)
    end if

    call kim%finalize()
    call done(all_passed)

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

    subroutine check_bg(name, arr, n_plasma, passed)
        !> A background quantity must be allocated on the plasma grid and finite.
        character(*), intent(in) :: name
        real(dp), allocatable, intent(in) :: arr(:)
        integer, intent(in) :: n_plasma
        logical, intent(inout) :: passed
        logical :: ok
        ok = allocated(arr)
        if (ok) ok = size(arr) == n_plasma
        if (ok) ok = all(ieee_is_finite(arr))
        call check('background '//name//' populated, sized, finite', ok, passed)
    end subroutine check_bg

    subroutine done(passed)
        logical, intent(in) :: passed
        if (passed) then
            print *, 'All kim_solver EM keystone checks PASSED'
            stop 0
        else
            print *, 'kim_solver EM keystone checks FAILED'
            stop 1
        end if
    end subroutine done

end program test_kim_solver_em
