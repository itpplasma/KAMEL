program test_flr2_tridiagonal
    use KIM_kinds_m, only: dp
    use flr2_tridiagonal_m, only: solve_flr2_tridiagonal, FLR2_SOLVE_OK
    implicit none

    integer, parameter :: n = 5
    real(dp), parameter :: tol = 1.0e-12_dp
    real(dp) :: x(n)
    complex(dp) :: q(n), a2_in(n), a2_out(n), a0(n)
    complex(dp) :: b2_in(n), b2_out(n), b0(n)
    complex(dp) :: c2_in(n), c2_out(n), c0(n)
    complex(dp) :: d2_in(n), d2_out(n), d0(n)
    complex(dp) :: f(n), g(n)
    integer :: stat

    x = [0.0_dp, 0.5_dp, 1.5_dp, 3.0_dp, 5.0_dp]
    q = cmplx(x**2, 0.0_dp, dp)

    a2_in = (3.0_dp, 0.0_dp)
    a2_out = (2.0_dp, 0.0_dp)
    a0 = (1.0_dp, 0.0_dp)
    b2_in = (3.0_dp, 0.0_dp)
    b2_out = (2.0_dp, 0.0_dp)
    b0 = (1.0_dp, 0.0_dp)
    c2_in = (3.0_dp, 0.0_dp)
    c2_out = (2.0_dp, 0.0_dp)
    c0 = (4.0_dp, 0.0_dp)
    d2_in = (6.0_dp, 0.0_dp)
    d2_out = (5.0_dp, 0.0_dp)
    d0 = (7.0_dp, 0.0_dp)

    call solve_flr2_tridiagonal(1, x, q, a2_in, a2_out, a0, &
                                b2_in, b2_out, b0, c2_in, c2_out, c0, &
                                d2_in, d2_out, d0, f, g, stat)

    if (stat /= FLR2_SOLVE_OK) then
        write(*,*) 'FAIL: solver status = ', stat
        stop 1
    end if
    if (maxval(abs(f - q)) > tol) then
        write(*,*) 'FAIL: manufactured quadratic solution mismatch: ', f
        stop 1
    end if
    if (maxval(abs(g(2:n - 1) - (32.0_dp + 11.0_dp*q(2:n - 1)))) > tol) then
        write(*,*) 'FAIL: second-derivative response mismatch: ', g
        stop 1
    end if
    if (abs(g(1) - 11.0_dp*q(1)) > tol .or. abs(g(n) - 11.0_dp*q(n)) > tol) then
        write(*,*) 'FAIL: algebraic boundary response mismatch: ', g
        stop 1
    end if

    write(*,*) 'FLR2 tridiagonal test PASSED'
end program test_flr2_tridiagonal
