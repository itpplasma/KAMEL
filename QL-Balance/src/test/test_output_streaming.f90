!> path2out lifecycle selected by prepare_output_file.
!>
!> Fresh run: an existing output file is replaced, and the QL-Balance pattern
!> of many open_rw / read back / delete / re-add / close cycles (growing
!> diagnostic history, float32 profile snapshots) gives the same content as
!> the in-memory writer. Restart (readfromtimestep > 0): the existing file is
!> extended.
program test_output_streaming
    use iso_fortran_env, only: real32, real64
    use control_mod, only: data_verbosity, ihdf5IO, readfromtimestep
    use h5mod, only: h5_create, h5_open, h5_open_rw, h5_close, h5_deinit, h5_id, &
        h5_mode_groupname, h5_get, h5_add, h5_add_float_1, h5_obj_exists, path2out, &
        prepare_output_file, h5_stream_write, h5_truncate_existing, &
        create_group_if_not_existent
    use grid_mod, only: npoib
    use wave_code_data, only: Br, r

    implicit none

    integer, parameter :: nrad = 3, nwrites = 40, nsnap = 8, nprof = 1000
    integer :: i, k
    logical :: exists
    complex(real64) :: brvac
    real(real64) :: brvac_data(nrad, 2), expected_brvac(nrad, 2)
    real(real64) :: profile(nprof)
    real(real32) :: snapshot(nprof)
    character(len=64) :: grp
    external :: write_Brvac

    ihdf5IO = 1
    data_verbosity = 2
    h5_mode_groupname = "test_mode"
    npoib = nrad
    allocate(r(nrad), Br(nrad))
    r = [1.0_real64, 2.0_real64, 3.0_real64]
    Br = [(cmplx(real(i, real64), -real(i, real64), real64), i = 1, nrad)]
    expected_brvac(:, 1) = r
    expected_brvac(:, 2) = abs(Br)
    profile = [(real(k, real64) / 7.0_real64, k = 1, nprof)]

    ! --- fresh run into a path that holds a stale file ---
    path2out = "test_output_streaming.h5"
    call write_stale_file()
    readfromtimestep = 0
    call prepare_output_file()
    if (.not. (h5_stream_write .and. h5_truncate_existing)) &
        error stop "fresh run did not select the streaming writer"

    do i = 1, nwrites
        brvac = cmplx(real(i, real64), -real(i, real64), real64)
        call write_Brvac(brvac)
        if (mod(i, nwrites / nsnap) == 0) then
            ! KinProfiles-like float32 snapshot in its own group
            write (grp, '("/test_mode/KinProfiles/",i0)') i
            call h5_open_rw(trim(path2out), h5_id)
            call create_group_if_not_existent(trim(grp))
            call h5_add_float_1(h5_id, trim(grp)//"/n", real(profile * i, real32), &
                lbound(profile), ubound(profile))
            call h5_close(h5_id)
            call h5_deinit()
        end if
    end do

    call h5_open(trim(path2out), h5_id)
    call h5_obj_exists(h5_id, "/stale", exists)
    if (exists) error stop "fresh run kept a dataset of the previous file"
    call h5_get(h5_id, "/test_mode/Brvac.dat", brvac_data)
    if (any(abs(brvac_data - expected_brvac) > 1.0e-12_real64)) &
        error stop "Brvac.dat differs after streaming writes"
    do i = nwrites / nsnap, nwrites, nwrites / nsnap
        write (grp, '("/test_mode/KinProfiles/",i0)') i
        call h5_get_float(trim(grp)//"/n", snapshot)
        if (any(snapshot /= real(profile * i, real32))) &
            error stop "float32 snapshot differs after streaming writes"
    end do
    call h5_close(h5_id)
    call h5_deinit()
    call check_h5dump(trim(path2out))

    ! --- restart: the existing file is extended ---
    path2out = "test_output_streaming_restart.h5"
    call write_stale_file()
    readfromtimestep = 5
    call prepare_output_file()
    if (h5_stream_write .or. h5_truncate_existing) &
        error stop "restart did not keep the in-memory writer"
    call write_Brvac(cmplx(1.0_real64, -1.0_real64, real64))
    call h5_open(trim(path2out), h5_id)
    call h5_obj_exists(h5_id, "/stale", exists)
    if (.not. exists) error stop "restart lost a dataset of the existing file"
    call h5_obj_exists(h5_id, "/test_mode/Brvac.dat", exists)
    if (.not. exists) error stop "restart did not write Brvac.dat"
    call h5_close(h5_id)
    call h5_deinit()
    call check_h5dump(trim(path2out))

contains

    !> Put a file from an earlier run at path2out. It is written under another
    !> name and copied: fortio keeps the writer state of every path this
    !> process wrote and would reuse it, which a real run never has for path2out.
    subroutine write_stale_file()
        integer :: status
        h5_stream_write = .false.
        h5_truncate_existing = .false.
        call h5_create(trim(path2out)//".stale", h5_id)
        call h5_add(h5_id, "/stale", 1)
        call h5_close(h5_id)
        call h5_deinit()
        call execute_command_line("cp "//trim(path2out)//".stale "//trim(path2out), &
            exitstat=status)
        if (status /= 0) error stop "could not place the stale output file"
    end subroutine write_stale_file

    subroutine h5_get_float(name, values)
        character(len=*), intent(in) :: name
        real(real32), intent(out) :: values(:)
        real(real64) :: buffer(size(values))
        call h5_get(h5_id, name, buffer)
        values = real(buffer, real32)
    end subroutine h5_get_float

    subroutine check_h5dump(path)
        character(len=*), intent(in) :: path
        integer :: status
        call execute_command_line("command -v h5dump > /dev/null", exitstat=status)
        if (status /= 0) return
        call execute_command_line("h5dump -H "//path//" > /dev/null", exitstat=status)
        if (status /= 0) error stop "h5dump rejects the output file"
    end subroutine check_h5dump

end program test_output_streaming
