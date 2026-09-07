module flr2_tridiagonal_m
    ! Adapted from KiLCA-FLR2's progonka.f90. The module form adds explicit
    ! interfaces, validation/status returns, and defined endpoint values for g.
    use KIM_kinds_m, only: dp
    implicit none
    private

    integer, parameter, public :: FLR2_SOLVE_OK = 0
    integer, parameter, public :: FLR2_SOLVE_BAD_SIZE = 1
    integer, parameter, public :: FLR2_SOLVE_BAD_GRID = 2
    integer, parameter, public :: FLR2_SOLVE_SINGULAR = 3

    public :: solve_flr2_tridiagonal

contains

    subroutine solve_flr2_tridiagonal(isw_f, x, q, &
                                      a2_in, a2_out, a0, &
                                      b2_in, b2_out, b0, &
                                      c2_in, c2_out, c0, &
                                      d2_in, d2_out, d0, &
                                      f, g, stat)
        integer, intent(in) :: isw_f
        real(dp), intent(in) :: x(:)
        complex(dp), intent(in) :: q(:)
        complex(dp), intent(in) :: a2_in(:), a2_out(:), a0(:)
        complex(dp), intent(in) :: b2_in(:), b2_out(:), b0(:)
        complex(dp), intent(in) :: c2_in(:), c2_out(:), c0(:)
        complex(dp), intent(in) :: d2_in(:), d2_out(:), d0(:)
        complex(dp), intent(out) :: f(:), g(:)
        integer, intent(out) :: stat

        integer :: i, n
        real(dp) :: dxm, dxp, dxt
        real(dp), allocatable :: second_derivative_weights(:, :)
        complex(dp) :: denominator
        complex(dp), allocatable :: alpha(:), beta(:), equation(:, :), source(:)

        stat = FLR2_SOLVE_OK
        f = (0.0_dp, 0.0_dp)
        g = (0.0_dp, 0.0_dp)
        n = size(x)

        if (n < 3 .or. .not. arrays_have_size(n)) then
            stat = FLR2_SOLVE_BAD_SIZE
            return
        end if

        do i = 2, n
            if (x(i) <= x(i - 1)) then
                stat = FLR2_SOLVE_BAD_GRID
                return
            end if
        end do

        if (is_near_zero(a0(1)) .or. is_near_zero(a0(n))) then
            stat = FLR2_SOLVE_SINGULAR
            return
        end if

        allocate(second_derivative_weights(-1:1, n))
        allocate(equation(-1:1, n), alpha(n), beta(n), source(n))
        second_derivative_weights = 0.0_dp
        equation = (0.0_dp, 0.0_dp)
        alpha = (0.0_dp, 0.0_dp)
        beta = (0.0_dp, 0.0_dp)
        source = (0.0_dp, 0.0_dp)

        do i = 2, n - 1
            dxp = x(i + 1) - x(i)
            dxm = x(i) - x(i - 1)
            dxt = dxp + dxm
            second_derivative_weights(-1, i) = 2.0_dp/(dxm*dxt)
            second_derivative_weights(0, i) = -2.0_dp/(dxm*dxp)
            second_derivative_weights(1, i) = 2.0_dp/(dxp*dxt)

            source(i) = b2_out(i)*sum(q(i - 1:i + 1)*second_derivative_weights(:, i)) &
                        + sum(b2_in(i - 1:i + 1)*q(i - 1:i + 1)*second_derivative_weights(:, i)) &
                        + b0(i)*q(i)

            equation(:, i) = a2_out(i)*second_derivative_weights(:, i) &
                             + a2_in(i - 1:i + 1)*second_derivative_weights(:, i)
            equation(0, i) = equation(0, i) + a0(i)
        end do

        beta(2) = b0(1)*q(1)/a0(1)
        do i = 2, n - 1
            denominator = equation(0, i) + alpha(i)*equation(-1, i)
            if (is_near_zero(denominator)) then
                stat = FLR2_SOLVE_SINGULAR
                return
            end if
            alpha(i + 1) = -equation(1, i)/denominator
            beta(i + 1) = (source(i) - equation(-1, i)*beta(i))/denominator
        end do

        f(n) = b0(n)*q(n)/a0(n)
        do i = n - 1, 1, -1
            f(i) = alpha(i + 1)*f(i + 1) + beta(i + 1)
        end do
        f = f*real(isw_f, dp)

        do i = 1, n
            if (i == 1 .or. i == n) then
                g(i) = c0(i)*f(i) + d0(i)*q(i)
            else
                g(i) = c2_out(i)*sum(f(i - 1:i + 1)*second_derivative_weights(:, i)) &
                       + sum(c2_in(i - 1:i + 1)*f(i - 1:i + 1)*second_derivative_weights(:, i)) &
                       + c0(i)*f(i) &
                       + d2_out(i)*sum(q(i - 1:i + 1)*second_derivative_weights(:, i)) &
                       + sum(d2_in(i - 1:i + 1)*q(i - 1:i + 1)*second_derivative_weights(:, i)) &
                       + d0(i)*q(i)
            end if
        end do

    contains

        logical function arrays_have_size(expected_size)
            integer, intent(in) :: expected_size

            arrays_have_size = size(q) == expected_size &
                               .and. size(a2_in) == expected_size &
                               .and. size(a2_out) == expected_size &
                               .and. size(a0) == expected_size &
                               .and. size(b2_in) == expected_size &
                               .and. size(b2_out) == expected_size &
                               .and. size(b0) == expected_size &
                               .and. size(c2_in) == expected_size &
                               .and. size(c2_out) == expected_size &
                               .and. size(c0) == expected_size &
                               .and. size(d2_in) == expected_size &
                               .and. size(d2_out) == expected_size &
                               .and. size(d0) == expected_size &
                               .and. size(f) == expected_size &
                               .and. size(g) == expected_size
        end function arrays_have_size

        logical function is_near_zero(value)
            complex(dp), intent(in) :: value

            is_near_zero = abs(value) <= tiny(1.0_dp)
        end function is_near_zero

    end subroutine solve_flr2_tridiagonal

end module flr2_tridiagonal_m
