program test_kim_solver
    !> Lifecycle/contract tests for the kim_solver_m API seam.
    !>
    !> These exercise the handle's own logic (status codes, ordering guards,
    !> no-halt-on-bad-input, safe teardown) without standing up a full physics
    !> run -- they are fixture-free and fast. A happy-path solve test (init ->
    !> solve -> results on a real bench config) is deferred until a runnable
    !> fixture + golden value exist (see the Phase 1 notes in the design doc).
    use kim_solver_m, only: kim_solver_t, kim_results_t, &
                            KIM_OK, KIM_NOT_SETUP, KIM_BAD_CONFIG
    implicit none

    logical :: all_passed
    all_passed = .true.

    call test_solve_before_init(all_passed)
    call test_init_rejects_missing_config(all_passed)
    call test_finalize_safe_when_unset(all_passed)
    call test_results_default_empty(all_passed)

    if (all_passed) then
        print *, 'All kim_solver tests PASSED'
        stop 0
    else
        print *, 'Some kim_solver tests FAILED'
        stop 1
    end if

contains

    subroutine test_solve_before_init(passed)
        !> solve() on a fresh handle reports KIM_NOT_SETUP rather than running.
        logical, intent(inout) :: passed
        type(kim_solver_t) :: kim
        integer :: ierr

        call kim%solve(m=2, n=1, stat=ierr)
        call report('solve before init -> KIM_NOT_SETUP', ierr == KIM_NOT_SETUP, passed)
    end subroutine test_solve_before_init

    subroutine test_init_rejects_missing_config(passed)
        !> A missing config path returns KIM_BAD_CONFIG -- it must not abort
        !> the process via the delegated namelist open().
        logical, intent(inout) :: passed
        type(kim_solver_t) :: kim
        integer :: ierr

        call kim%init('/nonexistent/path/__kim_no_such_config__.nml', stat=ierr)
        call report('init missing config -> KIM_BAD_CONFIG (no crash)', &
                    ierr == KIM_BAD_CONFIG, passed)
    end subroutine test_init_rejects_missing_config

    subroutine test_finalize_safe_when_unset(passed)
        !> finalize() on a never-initialised handle is a safe no-op.
        logical, intent(inout) :: passed
        type(kim_solver_t) :: kim

        call kim%finalize()
        call report('finalize safe on un-setup handle', &
                    kim%status_code() == KIM_OK, passed)
    end subroutine test_finalize_safe_when_unset

    subroutine test_results_default_empty(passed)
        !> A handle that has not solved returns empty results (no field grid,
        !> mode 0) rather than stale or uninitialised data.
        logical, intent(inout) :: passed
        type(kim_solver_t) :: kim
        type(kim_results_t) :: res

        res = kim%results()
        call report('default results are empty', &
                    (.not. allocated(res%r_field)) .and. res%m == 0, passed)
    end subroutine test_results_default_empty

    subroutine report(name, ok, passed)
        character(*), intent(in) :: name
        logical, intent(in) :: ok
        logical, intent(inout) :: passed
        if (ok) then
            print *, 'PASS: ', name
        else
            print *, 'FAIL: ', name
            passed = .false.
        end if
    end subroutine report

end program test_kim_solver
