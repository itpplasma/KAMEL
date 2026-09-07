!> Unit test for kilca_sysmat_profiles_m, focused on the argsort-based grid
!> rearrangement (translated from the C++ sort_index_doubles + the i/j
!> flat-index M-rearrangement loop). The fakes (test_sysmat_profiles_fakes)
!> add two out-of-sorted-order points during the (faked) adaptive-grid
!> refinement and verify, via the fake spline_calc, that every grid point's
!> M row landed at the correct post-sort offset.
program test_sysmat_profiles
    use kilca_sysmat_profiles_m, only: sysmat_profiles_create
    use kilca_sysmat_profiles_m, only: sysmat_profiles_destroy
    use, intrinsic :: iso_c_binding, only: c_int, c_intptr_t
    use sysmat_test_state, only: get_sysmat_test_failures
    implicit none

    integer(c_intptr_t) :: handle
    integer(c_int) :: failures

    ! r1=1.0, r2=3.0, rm=0 (forces the midpoint xa(2)=2.0): the 3 boundary
    ! samples are r={1.0, 2.0, 3.0}; the fake adaptive-grid step adds r=2.5
    ! and r=1.5, giving 5 unsorted points whose argsort the fake spline_calc
    ! verifies.
    handle = sysmat_profiles_create(2_c_int, 'f', '.', 3_c_int, &
                                    16_c_int, 1.0d-3, 0_c_int, 1.0d0, 3.0d0, 0.0d0)
    if (handle == 0_c_intptr_t) then
        write (*, '(a)') "FAIL: sysmat_profiles_create returned null handle"
        stop 1
    end if

    call sysmat_profiles_destroy(handle)

    failures = get_sysmat_test_failures()
    if (failures == 0) then
    write (*, '(a)') "PASS: kilca_sysmat_profiles_m grid sort/rearrangement matches expected layout"
    else
        write (*, '(a,i0)') "FAILED: ", failures
        stop 1
    end if
end program test_sysmat_profiles
