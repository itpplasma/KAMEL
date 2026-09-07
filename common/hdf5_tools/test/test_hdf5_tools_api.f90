program test_hdf5_tools_api
    use iso_fortran_env, only: real64
    use KAMEL_hdf5_tools, only: HID_T, h5_add, h5_close, h5_create, h5_deinit, &
        h5_get, h5_obj_exists, h5_open

    implicit none

    integer(HID_T) :: file_id
    logical :: matrix_exists, scalar_exists
    real(real64) :: matrix(2, 3), matrix_read(2, 3), scalar
    character(len=*), parameter :: path = "test_hdf5_tools_api.h5"

    matrix = reshape([1.0_real64, 2.0_real64, 3.0_real64, &
        4.0_real64, 5.0_real64, 6.0_real64], shape(matrix))

    call h5_create(path, file_id)
    call h5_add(file_id, "nested/group/matrix", matrix, [1, 1], [2, 3])
    call h5_add(file_id, "nested/group/scalar", 7.0_real64)
    call h5_close(file_id)
    call h5_deinit()

    call h5_open(path, file_id)
    call h5_obj_exists(file_id, "nested/group/matrix", matrix_exists)
    call h5_obj_exists(file_id, "nested/group/scalar", scalar_exists)
    if (.not. matrix_exists) error stop "2D HDF5 dataset was not written"
    if (.not. scalar_exists) error stop "nested scalar HDF5 dataset was not written"

    call h5_get(file_id, "nested/group/matrix", matrix_read)
    call h5_get(file_id, "nested/group/scalar", scalar)
    if (any(abs(matrix_read - matrix) > 1.0e-12_real64)) &
        error stop "2D HDF5 dataset values differ"
    if (abs(scalar - 7.0_real64) > 1.0e-12_real64) &
        error stop "nested scalar HDF5 value differs"

    call h5_close(file_id)
    call h5_deinit()
end program test_hdf5_tools_api
