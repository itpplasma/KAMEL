program test_periodic_embedding
    use QLBalance_kinds, only: dp
    use periodic_embedding_m, only: compact_transition
    implicit none

    real(dp), parameter :: lo = 2.0_dp, hi = 4.0_dp, width = 0.5_dp
    real(dp), parameter :: tol = 1.0e-14_dp
    real(dp) :: r, w_left, w_right, prev
    integer :: i

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

    print *, 'periodic embedding transition tests passed'
end program test_periodic_embedding
