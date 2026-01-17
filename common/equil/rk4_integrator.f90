!> @file rk4_integrator.f90
!> @brief RK4 (Runge-Kutta 4th order) integrator module for field-line tracing
!>
!> Provides generic RK4 integration for ODEs, used primarily for field-line
!> tracing in equilibrium flux coordinate calculations.

module rk4_integrator_m

    implicit none
    private

    public :: rk4_step

    !> Module-level work arrays for RK4 integration
    !> These are allocated on first use and reused for efficiency
    integer, save :: nmax_alloc = 0
    double precision, dimension(:), allocatable, save :: dydx_work, yt_work, dyt_work, dym_work

contains

    !> Perform a single RK4 integration step
    !>
    !> @param[inout] y     State vector (dimension ndim)
    !> @param[in]    ndim  Dimension of state vector
    !> @param[inout] x     Independent variable (updated to x + h)
    !> @param[in]    h     Step size
    !> @param[in]    derivs Subroutine computing derivatives: derivs(ndim, x, y, dydx)
    subroutine rk4_step(y, ndim, x, h, derivs)
        integer, intent(in) :: ndim
        double precision, intent(inout) :: y(ndim)
        double precision, intent(inout) :: x
        double precision, intent(in) :: h
        external :: derivs

        double precision :: hh, h6, xh
        integer :: i

        ! Reallocate work arrays if needed
        if (ndim /= nmax_alloc) then
            if (allocated(dydx_work)) deallocate(dydx_work)
            if (allocated(yt_work)) deallocate(yt_work)
            if (allocated(dyt_work)) deallocate(dyt_work)
            if (allocated(dym_work)) deallocate(dym_work)
            nmax_alloc = ndim
            allocate(dydx_work(nmax_alloc))
            allocate(yt_work(nmax_alloc))
            allocate(dyt_work(nmax_alloc))
            allocate(dym_work(nmax_alloc))
        end if

        hh = h * 0.5d0
        h6 = h / 6.d0
        xh = x + hh

        ! First step: evaluate at initial point
        call derivs(ndim, x, y, dydx_work)

        ! Second step: evaluate at midpoint using first derivative
        do i = 1, ndim
            yt_work(i) = y(i) + hh * dydx_work(i)
        end do
        call derivs(ndim, xh, yt_work, dyt_work)

        ! Third step: evaluate at midpoint using second derivative
        do i = 1, ndim
            yt_work(i) = y(i) + hh * dyt_work(i)
        end do
        call derivs(ndim, xh, yt_work, dym_work)

        ! Fourth step: evaluate at endpoint
        do i = 1, ndim
            yt_work(i) = y(i) + h * dym_work(i)
            dym_work(i) = dyt_work(i) + dym_work(i)
        end do
        call derivs(ndim, x + h, yt_work, dyt_work)

        ! Combine to get final result
        do i = 1, ndim
            y(i) = y(i) + h6 * (dydx_work(i) + dyt_work(i) + 2.d0 * dym_work(i))
        end do

        x = x + h

    end subroutine rk4_step

end module rk4_integrator_m


!> Legacy interface for backward compatibility with existing code
!> This provides the RK4D subroutine interface used by fouriermodes
subroutine RK4D(y, n, x, h, derivs)
    use rk4_integrator_m, only: rk4_step

    implicit none
    integer, intent(in) :: n
    double precision, intent(inout) :: y(n)
    double precision, intent(inout) :: x
    double precision, intent(in) :: h
    external :: derivs

    call rk4_step(y, n, x, h, derivs)

end subroutine RK4D
