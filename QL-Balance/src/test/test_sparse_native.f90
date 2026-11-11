program test_sparse_native
    use, intrinsic :: iso_fortran_env, only: dp => real64
    use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
    use sparse_mod, only: sparse_solve, sparse_solve_method
    implicit none
    real(dp) :: a(3, 3), values(4), expected(3), rhs(3), x(3), block(3, 2)
    real(dp) :: expected_block(3, 2), original_block(3, 2)
    complex(dp) :: az(3, 3), values_z(4), expected_z(3), rhs_z(3), xz(3), block_z(3, 2)
    complex(dp) :: expected_block_z(3, 2), original_block_z(3, 2)
    integer :: row(4), pointer(4), columns(4), method
    ! Callers such as evolvestep pass work buffers longer than nz after
    ! compaction; only the first nz row/value entries define the matrix.
    integer :: row_long(6)
    real(dp) :: values_long(6)
    complex(dp) :: values_z_long(6)
    character(32) :: mode

    ! A zero leading diagonal demands row pivoting. Four nonzeros make
    ! size(irow)=size(pointer)=ncol+1: distinguish CSC from column IDs by content.
    a(1, :) = [0.0_dp, 2.0_dp, 0.0_dp]
    a(2, :) = [1.0_dp, 1.0_dp, 0.0_dp]
    a(3, :) = [0.0_dp, 0.0_dp, 3.0_dp]
    row = [2, 1, 2, 3]
    pointer = [1, 2, 4, 5]
    columns = [1, 2, 2, 3]
    values = [1.0_dp, 2.0_dp, 1.0_dp, 3.0_dp]
    expected = [1.0_dp, -2.0_dp, 0.25_dp]
    rhs = [-4.0_dp, -1.0_dp, 0.75_dp]
    az = cmplx(a, 0.0_dp, dp)
    az(2, 1) = cmplx(1.0_dp, -1.0_dp, dp)
    az(1, 2) = cmplx(2.0_dp, 1.0_dp, dp)
    az(3, 3) = cmplx(3.0_dp, 1.0_dp, dp)
    values_z = [az(2, 1), az(1, 2), az(2, 2), az(3, 3)]
    expected_z = [cmplx(1.0_dp, 1.0_dp, dp), cmplx(-2.0_dp, 0.5_dp, dp), &
                  cmplx(0.25_dp, -0.75_dp, dp)]
    ! Literal RHS independently checked against the declared dense matrix.
    rhs_z = [cmplx(-4.5_dp, -1.0_dp, dp), cmplx(0.0_dp, 0.5_dp, dp), &
             cmplx(1.5_dp, -2.0_dp, dp)]
    call get_command_argument(1, mode)
    select case (trim(mode))
    case ('singular')
        a(3, :) = 0.0_dp
        x = rhs
        call sparse_solve(a, x)
        stop 0
    case ('invalid_iopt')
        x = rhs
        call sparse_solve(a, x, 4)
        stop 0
    case ('method_change')
        sparse_solve_method = 2
        x = rhs
        call sparse_solve(a, x, 1)
        sparse_solve_method = 3
        call sparse_solve(a, x, 2)
        stop 0
    case ('invalid_method')
        sparse_solve_method = 1
        x = rhs
        call sparse_solve(a, x)
        stop 0
    case ('invalid_column')
        columns(4) = 4
        x = rhs
        call sparse_solve(3, 3, 4, row, columns, values, x)
        stop 0
    end select

    do method = 2, 3
        sparse_solve_method = method
        x = rhs
        call sparse_solve(a, x)
        call check_real(x, expected, rhs)
        xz = rhs_z
        call sparse_solve(az, xz)
        call check_complex(xz, expected_z, rhs_z)

        expected_block(:, 1) = expected
        expected_block(:, 2) = -3.0_dp*expected
        original_block(:, 1) = rhs
        original_block(:, 2) = -3.0_dp*rhs
        block = original_block
        call sparse_solve(a, block)
        call check_real(block(:, 1), expected_block(:, 1), original_block(:, 1))
        call check_real(block(:, 2), expected_block(:, 2), original_block(:, 2))
        expected_block_z(:, 1) = expected_z
        expected_block_z(:, 2) = cmplx(0.0_dp, 2.0_dp, dp)*expected_z
        original_block_z(:, 1) = rhs_z
        original_block_z(:, 2) = cmplx(0.0_dp, 2.0_dp, dp)*rhs_z
        block_z = original_block_z
        call sparse_solve(az, block_z)
        call check_complex(block_z(:, 1), expected_block_z(:, 1), original_block_z(:, 1))
        call check_complex(block_z(:, 2), expected_block_z(:, 2), original_block_z(:, 2))

        ! Factor-only leaves RHS intact; repeated solves use retained factors.
        x = rhs
        call sparse_solve(3, 3, 4, row, pointer, values, x, 1)
        if (maxval(abs(x - rhs)) > 0.0_dp) error stop 'factor-only changed RHS'
        call sparse_solve(3, 3, 4, row, pointer, values, x, 2)
        call check_real(x, expected, rhs)
        x = -2.0_dp*rhs
        call sparse_solve(3, 3, 4, row, pointer, values, x, 2)
        call check_real(x, -2.0_dp*expected, -2.0_dp*rhs)
        call sparse_solve(3, 3, 4, row, pointer, values, x, 3)
        ! iopt omitted must reset to zero after the earlier explicit iopt=3.
        x = rhs
        call sparse_solve(3, 3, 4, row, columns, values, x)
        call check_real(x, expected, rhs)
        row_long = [row, 3, 1]
        values_long = [values, 7.0_dp, -5.0_dp]
        x = rhs
        call sparse_solve(3, 3, 4, row_long, pointer, values_long, x)
        call check_real(x, expected, rhs)
        block = original_block
        call sparse_solve(3, 3, 4, row, pointer, values, block, 2)
        call check_real(block(:, 2), expected_block(:, 2), original_block(:, 2))
        call sparse_solve(3, 3, 4, row, pointer, values, block, 3)

        xz = rhs_z
        call sparse_solve(3, 3, 4, row, pointer, values_z, xz, 1)
        if (maxval(abs(xz - rhs_z)) > 0.0_dp) error stop 'complex factor-only changed RHS'
        call sparse_solve(3, 3, 4, row, pointer, values_z, xz, 2)
        call check_complex(xz, expected_z, rhs_z)
        call sparse_solve(3, 3, 4, row, pointer, values_z, xz, 3)
        xz = rhs_z
        call sparse_solve(3, 3, 4, row, columns, values_z, xz)
        call check_complex(xz, expected_z, rhs_z)
        values_z_long = [values_z, cmplx(7.0_dp, 1.0_dp, dp), cmplx(-5.0_dp, 0.0_dp, dp)]
        xz = rhs_z
        call sparse_solve(3, 3, 4, row_long, pointer, values_z_long, xz)
        call check_complex(xz, expected_z, rhs_z)
        block_z = original_block_z
        call sparse_solve(3, 3, 4, row, pointer, values_z, block_z, 2)
        call check_complex(block_z(:, 2), expected_block_z(:, 2), original_block_z(:, 2))
        call sparse_solve(3, 3, 4, row, pointer, values_z, block_z, 3)
    end do
contains
    subroutine check_real(solution, reference, b)
        real(dp), intent(in) :: solution(:), reference(:), b(:)
        real(dp) :: residual
        if (.not. all(ieee_is_finite(solution))) error stop 'nonfinite solution'
        if (maxval(abs(solution - reference)) > 1.0e-12_dp) error stop 'incorrect real solution'
        residual = maxval(abs(matmul(a, solution) - b))/ &
                   (maxval(sum(abs(a), dim=2))*maxval(abs(solution)) + maxval(abs(b)))
        if (residual > 1.0e-12_dp) error stop 'real original-matrix residual failed'
    end subroutine check_real
    subroutine check_complex(solution, reference, b)
        complex(dp), intent(in) :: solution(:), reference(:), b(:)
        real(dp) :: residual
        if (.not. all(ieee_is_finite(real(solution, dp)))) error stop 'nonfinite real part'
        if (.not. all(ieee_is_finite(aimag(solution)))) error stop 'nonfinite imaginary part'
        if (maxval(abs(solution - reference)) > 1.0e-12_dp) error stop 'incorrect complex solution'
        residual = maxval(abs(matmul(az, solution) - b))/ &
                   (maxval(sum(abs(az), dim=2))*maxval(abs(solution)) + maxval(abs(b)))
        if (residual > 1.0e-12_dp) error stop 'complex original-matrix residual failed'
    end subroutine check_complex
end program test_sparse_native
