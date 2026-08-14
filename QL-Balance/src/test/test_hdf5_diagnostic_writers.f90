program test_hdf5_diagnostic_writers
    use iso_fortran_env, only: real64
    use control_mod, only: data_verbosity, ihdf5IO
    use grid_mod, only: npoib
    use h5mod, only: h5_close, h5_create, h5_deinit, h5_id, h5_mode_groupname, &
        h5_obj_exists, h5_open, path2out
    use wave_code_data, only: Br, r

    implicit none

    integer, parameter :: nrad = 3
    integer, parameter :: nwrites = 70
    integer :: i
    logical :: dataset_exists
    complex(real64) :: brvac
    external :: write_Brvac

    ihdf5IO = 1
    data_verbosity = 2
    path2out = "test_hdf5_diagnostic_writers.h5"
    h5_mode_groupname = "test_mode"

    npoib = nrad
    allocate(r(nrad), Br(nrad))
    r = [1.0_real64, 2.0_real64, 3.0_real64]
    Br = [(cmplx(real(i, real64), -real(i, real64), real64), i = 1, nrad)]

    call h5_create(trim(path2out), h5_id)
    call h5_close(h5_id)
    call h5_deinit()

    do i = 1, nwrites
        brvac = cmplx(real(i, real64), -real(i, real64), real64)
        call write_Brvac(brvac)
    end do

    call h5_open(trim(path2out), h5_id)
    call h5_obj_exists(h5_id, "/test_mode/Brvac.dat", dataset_exists)
    if (.not. dataset_exists) error stop "Brvac diagnostic dataset was not written"
    call h5_close(h5_id)
    call h5_deinit()
end program test_hdf5_diagnostic_writers
