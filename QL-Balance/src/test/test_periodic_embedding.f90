program test_periodic_embedding
    use QLBalance_kinds, only: dp
    use periodic_embedding_m, only: compact_transition, embed_complex_profile
    implicit none

    real(dp), parameter :: lo = 2.0_dp, hi = 4.0_dp, width = 0.5_dp
    real(dp), parameter :: tol = 1.0e-14_dp
    real(dp) :: r, w_left, w_right, prev
    integer :: i
    real(dp), parameter :: local_r(6) = [1.5_dp, 2.0_dp, 2.5_dp, 3.0_dp, 3.5_dp, 4.0_dp]
    real(dp), parameter :: global_r(5) = [1.0_dp, 2.0_dp, 3.0_dp, 4.0_dp, 5.0_dp]
    complex(dp) :: local_z(6), global_z(5)
    real(dp) :: weights(5)

    if (compact_transition(lo, lo, hi, width) /= 1.0_dp) error stop 'core lower endpoint'
    if (compact_transition(hi, lo, hi, width) /= 1.0_dp) error stop 'core upper endpoint'
    if (compact_transition(lo-width, lo, hi, width) /= 0.0_dp) error stop 'left support endpoint'
    if (compact_transition(hi+width, lo, hi, width) /= 0.0_dp) error stop 'right support endpoint'
    if (compact_transition(lo-2.0_dp*width, lo, hi, width) /= 0.0_dp) error stop 'left outside support'
    if (compact_transition(hi+2.0_dp*width, lo, hi, width) /= 0.0_dp) error stop 'right outside support'

    w_left = compact_transition(lo-width/2.0_dp, lo, hi, width)
    w_right = compact_transition(hi+width/2.0_dp, lo, hi, width)
    if (abs(w_left - w_right) > tol) error stop 'transition symmetry'
    if (.not. (w_left > 0.0_dp .and. w_left < 1.0_dp)) error stop 'transition range'

    prev = 0.0_dp
    do i = 0, 100
        r = lo - width + width * real(i, dp) / 100.0_dp
        if (compact_transition(r, lo, hi, width) + tol < prev) error stop 'left transition monotonicity'
        prev = compact_transition(r, lo, hi, width)
    end do
    prev = 1.0_dp
    do i = 0, 100
        r = hi + width * real(i, dp) / 100.0_dp
        if (compact_transition(r, lo, hi, width) > prev + tol) error stop 'right transition monotonicity'
        prev = compact_transition(r, lo, hi, width)
    end do

    ! KIM's periodic Fourier grid excludes the duplicated upper endpoint.
    ! It is sufficient that every global point with nonzero compact weight is
    ! covered by the local grid; zero-weight support endpoints need no value.
    local_z = cmplx(2.0_dp, -1.0_dp, dp)
    call embed_complex_profile(local_r, local_z, global_r, lo, hi, width, global_z, weights)
    if (abs(global_z(2)-local_z(1)) > tol .or. abs(global_z(3)-local_z(1)) > tol &
            .or. abs(global_z(4)-local_z(1)) > tol) &
        error stop 'endpoint-exclusive periodic grid embedding'

    print *, 'periodic embedding transition tests passed'
end program test_periodic_embedding
