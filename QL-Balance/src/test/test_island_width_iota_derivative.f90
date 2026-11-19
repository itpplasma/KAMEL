program test_island_width_iota_derivative
    !! The island width uses d(1/q)/dr on the grid that carries q (the balance boundary grid
    !! rb). With q = 1 + (r/a)^2 on a non-uniform grid, every point, including the shifted
    !! stencils at both ends, must reproduce the analytic -2 r / (a^2 (1 + (r/a)^2)^2).
    use magnetic_island_width_m, only: inverse_q_derivative
    implicit none

    integer, parameter :: dp = kind(1.0d0)
    integer, parameter :: n = 200
    real(dp), parameter :: a = 50.0_dp, tol = 1.0e-5_dp
    real(dp) :: r(n), q(n), diotadr(n), exact(n), rel_err
    integer :: i

    do i = 1, n
        r(i) = a*(real(i - 1, dp)/real(n - 1, dp))**1.5_dp
    end do
    q = 1.0_dp + (r/a)**2
    exact = -2.0_dp*r/(a**2*(1.0_dp + (r/a)**2)**2)

    call inverse_q_derivative(r, q, diotadr)

    rel_err = maxval(abs(diotadr - exact))/maxval(abs(exact))
    print *, 'max relative error of d(1/q)/dr:', rel_err
    if (rel_err > tol) error stop 'd(1/q)/dr does not match the analytic derivative'
    if (abs(diotadr(n) - exact(n)) > tol*abs(exact(n))) &
        error stop 'd(1/q)/dr at the outer boundary is not pointwise'
end program test_island_width_iota_derivative
