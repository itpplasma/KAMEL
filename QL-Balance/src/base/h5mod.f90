
module h5mod

    use KAMEL_hdf5_tools
    use control_mod
    use logger_m, only: log_debug

    implicit none

    integer(HID_T) :: h5_id, group_id_1, group_id_2, group_id_3, dataset_id
    logical :: h5_exists_log
    integer :: mode_n, mode_m
    character(len=1024) :: path2inp ! path to hdf5 file with input data
    character(len=1024) :: path2time ! path to hdf5 file from which time evolved profiles are read
    character(len=1024) :: path2out ! path to hdf5 file where output is written
    character(len=1024) :: h5_mode_groupname
    character(len=1024) :: h5_currentgrp !> current hdf5 group string

    contains

    !> Select how fortio writes path2out. Call once after read_config, before
    !> the first write.
    !>
    !> QL-Balance opens and closes path2out many times per time step. With
    !> fortio's default in-memory writer every close rewrites the complete
    !> file image (float32 and int32 payloads one element per write call), so
    !> the cost per step grows with the number of stored steps and is dominated
    !> by file-system latency on network file systems. In streaming mode each
    !> dataset payload is written once when it is added and a close only
    !> rewrites the metadata.
    !>
    !> Streaming requires a fresh image: an existing path2out is replaced, not
    !> extended. A restart (readfromtimestep > 0) continues the histories
    !> already stored in path2out, so it keeps the in-memory writer and extends
    !> the file as before.
    subroutine prepare_output_file()

        use logger_m, only: log_info

        if (readfromtimestep > 0) then
            h5_stream_write = .false.
            h5_truncate_existing = .false.
            call log_info("    output: extending "//trim(path2out)// &
                " (restart; in-memory writer)")
            return
        end if
        h5_stream_write = .true.
        h5_truncate_existing = .true.
        call log_info("    output: new streaming file "//trim(path2out))

    end subroutine prepare_output_file

    subroutine write_reason_for_stop_to_h5(reason)

        implicit none

        character(*), intent(in) :: reason

        CALL h5_init()
        CALL h5_open_rw(path2out, h5_id)
        CALL h5_add_string(h5_id, trim(h5_mode_groupname)// &
        '/stopping_criterion', reason)
        CALL h5_close(h5_id)
        CALL h5_deinit()

    end subroutine

    subroutine create_group_if_not_existent(group_name)

        implicit none

        character(*), intent(in) :: group_name

        call log_debug("Creating group " // trim(group_name))

        call h5_obj_exists(h5_id, trim(group_name), h5_exists_log)
        if (.not. h5_exists_log) then
            call h5_define_group(h5_id, trim(group_name), group_id_1)
            call h5_close_group(group_id_1)
        end if

    end subroutine


end module h5mod
