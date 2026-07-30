!> @file test_hdf5_tools.f90
!> @brief Unit tests for the KAMEL hdf5_tools compatibility shim
!>
!> Covers the two routines KAMEL adds on top of LIBNEO::hdf5_tools:
!>   - h5_add_float_1: writes a real(:) array as an on-disk float64 dataset
!>   - h5_append_double_1: appends a column into an unlimited 2D matrix

program test_hdf5_tools
    use hdf5, only: h5dopen_f, h5dclose_f, h5dget_type_f, h5tclose_f, &
                    h5tget_size_f, h5tget_class_f, H5T_FLOAT_F, SIZE_T, HSIZE_T
    use h5lt, only: h5ltread_dataset_double_f
    use KAMEL_hdf5_tools, only: HID_T, H5T_NATIVE_DOUBLE, h5_init, h5_deinit, &
                                h5_create, h5_open, h5_close, h5_get_double_1, &
                                h5_define_unlimited_matrix, h5_add_float_1, &
                                h5_append_double_1

    implicit none

    double precision, parameter :: tolerance = 1.0d-12
    character(len=*), parameter :: float_file = "test_hdf5_add_float.h5"
    character(len=*), parameter :: append_file = "test_hdf5_append.h5"

    logical :: all_passed

    all_passed = .true.

    print *, "========================================"
    print *, "Testing KAMEL hdf5_tools shim"
    print *, "========================================"

    call test_add_float_1(all_passed)
    call test_append_double_1(all_passed)

    print *, ""
    if (all_passed) then
        print *, "All hdf5_tools tests PASSED"
        stop 0
    else
        print *, "Some hdf5_tools tests FAILED"
        stop 1
    end if

contains

    !> h5_add_float_1 must round-trip real values and store them as float64.
    subroutine test_add_float_1(passed)
        logical, intent(inout) :: passed

        integer, parameter :: n = 3
        real, dimension(n) :: value
        double precision, dimension(n) :: got
        integer, dimension(1) :: lbounds, ubounds
        integer(HID_T) :: fid, dsetid, dtype
        integer(SIZE_T) :: type_size
        integer :: type_class, err, i

        print *, ""
        print *, "Test: h5_add_float_1"

        value = (/1.5, 2.25, -3.75/)
        lbounds = (/1/)
        ubounds = (/n/)

        call h5_init()
        call h5_create(float_file, fid)
        call h5_add_float_1(fid, "/floatvec", value, lbounds, ubounds)
        call h5_close(fid)

        call h5_open(float_file, fid)
        call h5_get_double_1(fid, "/floatvec", got)

        do i = 1, n
            if (abs(got(i) - dble(value(i))) > tolerance) then
                print *, "  FAIL: value mismatch at", i, got(i), dble(value(i))
                passed = .false.
            end if
        end do

        call h5dopen_f(fid, "/floatvec", dsetid, err)
        call h5dget_type_f(dsetid, dtype, err)
        call h5tget_size_f(dtype, type_size, err)
        call h5tget_class_f(dtype, type_class, err)
        call h5tclose_f(dtype, err)
        call h5dclose_f(dsetid, err)
        call h5_close(fid)
        call h5_deinit()

        if (type_class /= H5T_FLOAT_F) then
            print *, "  FAIL: dataset class is not floating point"
            passed = .false.
        end if
        if (type_size /= 8) then
            print *, "  FAIL: on-disk type size is", type_size, "expected 8 (float64)"
            passed = .false.
        end if

        if (passed) print *, "  PASS"
    end subroutine test_add_float_1

    !> h5_append_double_1 must place each vector into its target column.
    subroutine test_append_double_1(passed)
        logical, intent(inout) :: passed

        integer, parameter :: nrows = 4, ncols = 3
        double precision, dimension(nrows) :: col1, col2, col3
        double precision, dimension(nrows, ncols) :: got, expected
        integer(HSIZE_T), dimension(2) :: dims
        integer(HID_T) :: fid, dsetid
        integer :: err, i, j

        print *, ""
        print *, "Test: h5_append_double_1"

        col1 = (/1.d0, 2.d0, 3.d0, 4.d0/)
        col2 = (/10.d0, 20.d0, 30.d0, 40.d0/)
        col3 = (/100.d0, 200.d0, 300.d0, 400.d0/)
        expected(:, 1) = col1
        expected(:, 2) = col2
        expected(:, 3) = col3

        call h5_init()
        call h5_create(append_file, fid)
        call h5_define_unlimited_matrix(fid, "/mat", H5T_NATIVE_DOUBLE, &
                                        (/-1, ncols/), dsetid)
        call h5_append_double_1(dsetid, col1, 1)
        call h5_append_double_1(dsetid, col2, 2)
        call h5_append_double_1(dsetid, col3, 3)
        call h5_close(fid)

        call h5_open(append_file, fid)
        dims = (/int(nrows, kind=HSIZE_T), int(ncols, kind=HSIZE_T)/)
        call h5ltread_dataset_double_f(fid, "/mat", got, dims, err)
        call h5_close(fid)
        call h5_deinit()

        do j = 1, ncols
            do i = 1, nrows
                if (abs(got(i, j) - expected(i, j)) > tolerance) then
                    print *, "  FAIL: mismatch at", i, j, got(i, j), expected(i, j)
                    passed = .false.
                end if
            end do
        end do

        if (passed) print *, "  PASS"
    end subroutine test_append_double_1

end program test_hdf5_tools
