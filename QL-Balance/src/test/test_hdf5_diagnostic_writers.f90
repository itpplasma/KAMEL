program test_hdf5_diagnostic_writers
    use iso_fortran_env, only: real64
    use control_mod, only: data_verbosity, ihdf5IO
    use grid_mod, only: dery_equisource, dqle11, dqli11, npoi, npoib, rb, rc
    use h5mod, only: h5_create, h5_open, h5_close, h5_init, h5_deinit, h5_id, h5_mode_groupname, &
        h5_get, h5_obj_exists, path2out
    use QLbalance_diag, only: ind_dqle, ind_dqli, timscal_dql, timscal_dqli
    use wave_code_data, only: Br, r
    use writeData_m, only: writefort9999, writefort9999_stellarator

    implicit none

    integer, parameter :: nrad = 3
    integer, parameter :: nwrites = 70
    integer :: i
    logical :: dataset_exists
    complex(real64) :: brvac
    real(real64) :: dqle11_prev(nrad), dqli11_prev(nrad)
    real(real64) :: brvac_data(nrad, 2), expected_brvac(nrad, 2)
    real(real64) :: fort_data(nrad, 3), expected_fort(nrad, 3)
    real(real64) :: equisource_data(4, nrad), expected_equisource(4, nrad)
    external :: write_Brvac
    external :: write_balance_eqs_source_terms

    ihdf5IO = 1
    data_verbosity = 2
    path2out = "test_hdf5_diagnostic_writers.h5"
    h5_mode_groupname = "test_mode"

    npoib = nrad
    npoi = nrad
    allocate(r(nrad), Br(nrad), rb(nrad), rc(nrad), dqle11(nrad), dqli11(nrad))
    allocate(dery_equisource(4*nrad))
    r = [1.0_real64, 2.0_real64, 3.0_real64]
    Br = [(cmplx(real(i, real64), -real(i, real64), real64), i = 1, nrad)]
    rb = [10.0_real64, 20.0_real64, 30.0_real64]
    rc = rb
    dqle11 = [1.0_real64, 2.0_real64, 3.0_real64]
    dqli11 = [4.0_real64, 5.0_real64, 6.0_real64]
    dqle11_prev = [2.0_real64, 4.0_real64, 6.0_real64]
    dqli11_prev = [5.0_real64, 7.0_real64, 9.0_real64]
    dery_equisource = [(real(i, real64), i = 1, 4*nrad)]
    ind_dqle = [1]
    ind_dqli = [2]
    timscal_dql = 0.0_real64
    timscal_dqli = 0.0_real64

    expected_brvac(:, 1) = r
    expected_brvac(:, 2) = abs(Br)
    expected_fort(:, 1) = rb
    expected_fort(:, 2) = abs(dqle11_prev - dqle11)
    expected_fort(:, 3) = abs(dqli11_prev - dqli11)
    expected_equisource = reshape([(real(i, real64), i = 1, 4*nrad)], [4, nrad])

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
    call h5_get(h5_id, "/test_mode/Brvac.dat", brvac_data)
    call assert_close(brvac_data, expected_brvac, "Brvac.dat")
    call h5_close(h5_id)
    call h5_deinit()

    call writefort9999(dqle11_prev, dqli11_prev)
    call h5_open(trim(path2out), h5_id)
    call h5_get(h5_id, "/test_mode/fort.9999", fort_data)
    call assert_close(fort_data, expected_fort, "fort.9999")
    call h5_close(h5_id)
    call h5_deinit()

    call writefort9999_stellarator(dqle11_prev, dqli11_prev)
    call h5_open(trim(path2out), h5_id)
    call h5_get(h5_id, "/test_mode/fort.9999", fort_data)
    call assert_close(fort_data, expected_fort, "stellarator fort.9999")
    call h5_close(h5_id)
    call h5_deinit()

    call write_balance_eqs_source_terms
    call h5_open(trim(path2out), h5_id)
    call h5_get(h5_id, "/test_mode/equisource.dat", equisource_data)
    call assert_close(equisource_data, expected_equisource, "equisource.dat")
    call h5_close(h5_id)
    call h5_deinit()

contains

    subroutine assert_close(actual, expected, label)
        real(real64), intent(in) :: actual(:, :), expected(:, :)
        character(len=*), intent(in) :: label

        if (any(shape(actual) /= shape(expected))) then
            error stop "HDF5 diagnostic shape mismatch: "//trim(label)
        end if
        if (any(abs(actual - expected) > 1.0e-12_real64)) then
            error stop "HDF5 diagnostic values mismatch: "//trim(label)
        end if
    end subroutine assert_close

end program test_hdf5_diagnostic_writers
