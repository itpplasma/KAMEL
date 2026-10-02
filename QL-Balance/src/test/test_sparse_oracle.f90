program test_sparse_oracle
    use, intrinsic :: iso_fortran_env, only: dp => real64
    use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
    use sparse_mod, only: sparse_solve, sparse_solve_method
    implicit none
    complex(dp) :: a(3, 3), b(3), x(3), expected(3), defect(3)
    real(dp) :: ar(3, 3), br(3), xr(3), expected_r(3), defect_r(3), err, residual

    ! Independently specified nonsymmetric systems, with a zero first pivot.
    ! Literal right-hand sides avoid sharing the conversion under test.
    a(1, :) = [cmplx(0, 0, dp), cmplx(2, 1, dp), cmplx(0, 0, dp)]
    a(2, :) = [cmplx(1, -1, dp), cmplx(4, 0, dp), cmplx(1, 0, dp)]
    a(3, :) = [cmplx(0, 0, dp), cmplx(2, 0, dp), cmplx(3, 1, dp)]
    b = [cmplx(-4.5_dp, -1.0_dp, dp), cmplx(-5.75_dp, 1.25_dp, dp), &
        cmplx(-2.5_dp, -1.0_dp, dp)]
    expected = [cmplx(1, 1, dp), cmplx(-2.0_dp, 0.5_dp, dp), &
        cmplx(0.25_dp, -0.75_dp, dp)]
    sparse_solve_method = 2
    x = b
    call sparse_solve(a, x, 0)
    if (.not. all(ieee_is_finite(real(x, dp)))) &
        error stop 'Nonfinite complex sparse solution real part'
    if (.not. all(ieee_is_finite(aimag(x)))) &
        error stop 'Nonfinite complex sparse solution imaginary part'
    err = maxval(abs(x - expected))/maxval(abs(expected))
    defect = matmul(a, x) - b
    residual = maxval(abs(defect))/ &
        (sqrt(sum(abs(a)**2))*sqrt(sum(abs(x)**2)) + sqrt(sum(abs(b)**2)))
    if (.not. ieee_is_finite(err)) error stop 'Nonfinite complex solution error'
    if (.not. ieee_is_finite(residual)) error stop 'Nonfinite complex residual'
    print *, 'Complex solution error, backward residual:', err, residual
    if (err > 1.0e-12_dp) error stop 'Complex sparse solution failed'
    if (residual > 1.0e-12_dp) error stop 'Complex sparse residual failed'

    ar = real(a, dp)
    br = [-4.0_dp, -6.75_dp, -3.25_dp]
    expected_r = [1.0_dp, -2.0_dp, 0.25_dp]
    xr = br
    call sparse_solve(ar, xr, 0)
    if (.not. all(ieee_is_finite(xr))) error stop 'Nonfinite real sparse solution'
    err = maxval(abs(xr - expected_r))/maxval(abs(expected_r))
    defect_r = matmul(ar, xr) - br
    residual = maxval(abs(defect_r))/ &
        (sqrt(sum(ar**2))*sqrt(sum(xr**2)) + sqrt(sum(br**2)))
    if (.not. ieee_is_finite(err)) error stop 'Nonfinite real solution error'
    if (.not. ieee_is_finite(residual)) error stop 'Nonfinite real residual'
    print *, 'Real solution error, backward residual:', err, residual
    if (err > 1.0e-12_dp) error stop 'Real sparse solution failed'
    if (residual > 1.0e-12_dp) error stop 'Real sparse residual failed'
end program test_sparse_oracle
