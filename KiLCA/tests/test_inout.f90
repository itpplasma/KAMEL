program test_inout
    use iso_c_binding, only: c_int, c_double, c_null_char
    use kilca_inout_m, only: save_cmplx_matrix_
    implicit none
    real(c_double) :: x(2) = [1.0_c_double, 2.0_c_double]
    real(c_double) :: matrix(8) = [10d0, -10d0, 20d0, -20d0, 11d0, -11d0, 21d0, -21d0]
    real(c_double) :: row(3)
    integer :: rc, unit, col, k, ios
    character(len=32) :: filename

    rc = save_cmplx_matrix_(1_c_int, 2_c_int, 2_c_int, x, matrix, &
                            'matrix_output'//c_null_char)
    if (rc /= 0) error stop 'Matrix output failed'
    do col = 0, 1
        write (filename, '(a,i0,a)') 'matrix_output_', col, '.dat'
        open (newunit=unit, file=trim(filename), status='old')
        do k = 1, 2
            read (unit, *, iostat=ios) row
            if (ios /= 0) error stop 'Missing complex matrix columns'
            if (any(abs(row - [x(k), 10d0 * (col + 1) + k - 1, &
                               -10d0 * (col + 1) - k + 1]) > 1d-12)) &
                error stop 'Wrong complex matrix values or layout'
        end do
        close (unit, status='delete')
    end do
end program
