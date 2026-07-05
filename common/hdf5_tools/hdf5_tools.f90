module KAMEL_hdf5_tools
  use hdf5, only: HID_T, HSIZE_T, H5S_SELECT_SET_F, H5T_NATIVE_DOUBLE, &
                  h5dset_extent_f, h5screate_simple_f, h5dget_space_f,  &
                  h5sselect_hyperslab_f, h5dwrite_f, h5sclose_f
  use hdf5_tools

  implicit none
  public

contains

  subroutine h5_add_float_1(h5id, dataset, value, lbounds, ubounds, comment, unit)
    integer(HID_T), intent(in)        :: h5id
    character(len=*), intent(in)      :: dataset
    real, dimension(:), intent(in)    :: value
    integer, dimension(:), intent(in) :: lbounds, ubounds
    character(len=*), optional        :: comment
    character(len=*), optional        :: unit

    ! Deliberate widening: the on-disk dataset is written as H5T_NATIVE_DOUBLE
    ! (float64) via the double-precision path, not float32. Widening is lossless.
    call h5_add_double_1(h5id, dataset, dble(value), lbounds, ubounds, comment, unit)
  end subroutine h5_add_float_1

  subroutine h5_append_double_1(dsetid, value, offset)
    integer(HID_T), intent(in)             :: dsetid
    double precision, dimension(:), intent(in) :: value
    integer, intent(in)                    :: offset

    integer(HSIZE_T), dimension(2)         :: dims
    integer(HSIZE_T), dimension(2)         :: extent
    integer(HID_T)                         :: memspace
    integer(HID_T)                         :: dspaceid
    integer                                :: rank
    integer                                :: nvalues
    integer(HSIZE_T), dimension(2)         :: offsetd

    rank = 2
    nvalues = size(value)
    extent = (/int(nvalues, kind=HSIZE_T), int(offset, kind=HSIZE_T)/)
    dims = (/int(nvalues, kind=HSIZE_T), 1_HSIZE_T/)
    offsetd = (/0_HSIZE_T, int(offset - 1, kind=HSIZE_T)/)

    call h5dset_extent_f(dsetid, extent, h5error)
    call h5_check()
    call h5screate_simple_f(rank, dims, memspace, h5error)
    call h5_check()
    call h5dget_space_f(dsetid, dspaceid, h5error)
    call h5_check()
    call h5sselect_hyperslab_f(dspaceid, H5S_SELECT_SET_F, offsetd, dims, h5error)
    call h5_check()
    call h5dwrite_f(dsetid, H5T_NATIVE_DOUBLE, value, dims, h5error, memspace, dspaceid)
    call h5_check()

    call h5sclose_f(memspace, h5error)
    call h5sclose_f(dspaceid, h5error)
  end subroutine h5_append_double_1
end module KAMEL_hdf5_tools
