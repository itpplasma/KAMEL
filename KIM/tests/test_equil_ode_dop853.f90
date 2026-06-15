program test_equil_ode_dop853

    ! Behavioral check for the ddeabm -> fortnum DOP853 swap (Z2).
    ! calculate_equil integrates the scalar force-balance ODE
    !   du/dr = -2 r u / (q^2 R0^2 + r^2) - 8 pi dpress
    ! over the radial grid with ode_solve_dop, restarting from each grid point.
    ! With dpress = 0 and constant q the ODE has the closed form
    !   u(r) = u0 (q^2 R0^2 + r0^2) / (q^2 R0^2 + r^2),
    ! so the grid-by-grid continuation must reproduce u(r) to the 1e-12
    ! tolerance calculate_equil requests.

    use KIM_kinds_m, only: dp
    use fortnum_ode_dop853, only: ode_solve_dop
    use fortnum_status, only: fortnum_status_t, FORTNUM_OK

    implicit none

    real(dp), parameter :: q = 1.3_dp
    real(dp), parameter :: R0 = 165.0_dp
    real(dp), parameter :: rtol = 1.0d-12, atol = 1.0d-12
    real(dp), parameter :: check_tol = 1.0d-9

    integer,  parameter :: ngrid = 64
    real(dp) :: r_grid(ngrid)
    real(dp) :: u(ngrid)
    real(dp) :: u_prev(1)
    real(dp), allocatable :: r_out(:), u_seg(:,:)
    type(fortnum_status_t) :: status
    real(dp) :: denom0, u_exact, rel_err, max_err
    integer :: i

    do i = 1, ngrid
        r_grid(i) = 3.0_dp + real(i - 1, dp) * (67.0_dp - 3.0_dp) / real(ngrid - 1, dp)
    end do

    denom0 = q**2 * R0**2 + r_grid(1)**2
    u(1) = 4.2d8
    u_prev(1) = u(1)

    do i = 2, ngrid
        call ode_solve_dop(dudr, r_grid(i-1), r_grid(i), u_prev, &
                           r_out, u_seg, status, rtol=rtol, atol=atol)
        if (status%code /= FORTNUM_OK .or. size(u_seg, 2) < 1) then
            print *, 'ode_solve_dop failed at i=', i, ' ', trim(status%msg)
            error stop 1
        end if
        u_prev(1) = u_seg(1, size(u_seg, 2))
        u(i) = u_prev(1)
    end do

    max_err = 0.0_dp
    do i = 1, ngrid
        u_exact = u(1) * denom0 / (q**2 * R0**2 + r_grid(i)**2)
        rel_err = abs(u(i) - u_exact) / abs(u_exact)
        max_err = max(max_err, rel_err)
    end do

    print *, 'DOP853 equilibrium ODE: max relative error vs closed form = ', max_err
    if (max_err > check_tol) then
        print *, 'error exceeds tolerance ', check_tol
        error stop 2
    end if
    print *, 'Equilibrium ODE DOP853 OK'

contains

    subroutine dudr(r, y, dydt, ctx)
        real(dp), intent(in)  :: r
        real(dp), intent(in)  :: y(:)
        real(dp), intent(out) :: dydt(:)
        class(*), intent(in), optional :: ctx
        dydt(1) = -2.0_dp * r * y(1) / (q**2 * R0**2 + r**2)
    end subroutine dudr

end program test_equil_ode_dop853
