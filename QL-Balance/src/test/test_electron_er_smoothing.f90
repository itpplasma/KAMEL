program test_electron_er_smoothing
    use QLBalance_kinds, only: dp
    use electron_er_smoothing_m, only: gauss_smooth_cells
    implicit none

    integer, parameter :: n = 401
    real(dp) :: x(n), y(n), y0(n)
    integer :: i

    do i = 1, n
        x(i) = real(i - 201, dp)
    end do

    ! Off: sigma <= 0 leaves the profile untouched.
    y = x**2
    y0 = y
    call gauss_smooth_cells(y, 0.0_dp)
    if (any(y /= y0)) error stop 'sigma = 0 must not change the profile'

    ! A linear profile is preserved in the interior (symmetric weights).
    y = 3.0_dp*x - 1.0_dp
    call gauss_smooth_cells(y, 5.0_dp)
    if (maxval(abs(y(30:n - 29) - (3.0_dp*x(30:n - 29) - 1.0_dp))) > 1.0e-10_dp) &
        error stop 'linear profile not preserved in the interior'

    ! A constant is preserved everywhere, including the truncated ends.
    y = 2.5_dp
    call gauss_smooth_cells(y, 20.0_dp)
    if (maxval(abs(y - 2.5_dp)) > 1.0e-12_dp) error stop 'constant not preserved'

    ! A single-cell spike of unit area keeps its area and gets the
    ! requested width: variance sigma^2 (the 4-sigma cut lowers it by ~1e-3).
    y = 0.0_dp
    y(201) = 1.0_dp
    call gauss_smooth_cells(y, 5.0_dp)
    if (abs(sum(y) - 1.0_dp) > 1.0e-12_dp) error stop 'spike area not conserved'
    if (abs(sum(y*x**2) - 25.0_dp) > 0.1_dp) error stop 'wrong smoothing width'

    write (*, *) 'test_electron_er_smoothing passed'
end program test_electron_er_smoothing
