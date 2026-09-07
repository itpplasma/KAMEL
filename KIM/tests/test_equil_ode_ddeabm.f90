program test_equil_ode_ddeabm

    use KIM_kinds_m, only: dp
    use fortnum_ode_ddeabm, only: ddeabm_state_t, ddeabm_init, ddeabm_integrate_to
    use fortnum_status, only: fortnum_status_t, FORTNUM_OK

    implicit none

    real(dp), parameter :: q = 1.3_dp
    real(dp), parameter :: R0 = 165.0_dp
    real(dp), parameter :: rtol = 1.0d-12, atol = 1.0d-12
    real(dp), parameter :: check_tol = 1.0d-9

    integer, parameter :: ngrid = 64
    real(dp) :: r_grid(ngrid)
    real(dp) :: u(ngrid)
    real(dp), allocatable :: u_seg(:)
    type(ddeabm_state_t) :: ode_state
    type(fortnum_status_t) :: status
    real(dp) :: denom0, u_exact, rel_err, max_err
    integer :: i

    do i = 1, ngrid
        r_grid(i) = 3.0_dp + real(i - 1, dp) * (67.0_dp - 3.0_dp) / real(ngrid - 1, dp)
    end do

    denom0 = q**2 * R0**2 + r_grid(1)**2
    u(1) = 4.2d8
    call ddeabm_init(ode_state, 1, r_grid(1), [u(1)])

    do i = 2, ngrid
        call ddeabm_integrate_to(dudr, ode_state, r_grid(i), rtol, [atol], &
            u_seg, status, tstop=r_grid(ngrid))
        if (status%code /= FORTNUM_OK) then
            print *, 'ddeabm_integrate_to failed at i=', i, ' ', trim(status%msg)
            error stop 1
        end if
        if (.not. allocated(u_seg)) then
            print *, 'ddeabm_integrate_to did not return a solution at i=', i
            error stop 2
        end if
        u(i) = u_seg(1)
    end do

    max_err = 0.0_dp
    do i = 1, ngrid
        u_exact = u(1) * denom0 / (q**2 * R0**2 + r_grid(i)**2)
        rel_err = abs(u(i) - u_exact) / abs(u_exact)
        max_err = max(max_err, rel_err)
    end do

    print *, 'ddeabm equilibrium ODE: max relative error vs closed form = ', max_err
    if (max_err > check_tol) then
        print *, 'error exceeds tolerance ', check_tol
        error stop 3
    end if
    print *, 'Equilibrium ODE ddeabm OK'

contains

    subroutine dudr(r, y, dydt, ctx)
        real(dp), intent(in) :: r
        real(dp), intent(in) :: y(:)
        real(dp), intent(out) :: dydt(:)
        class(*), intent(in), optional :: ctx

        dydt(1) = -2.0_dp * r * y(1) / (q**2 * R0**2 + r**2)

    end subroutine dudr

end program test_equil_ode_ddeabm
