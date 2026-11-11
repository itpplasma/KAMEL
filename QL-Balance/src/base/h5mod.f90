
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

    !> Record that |Br| at the resonance left its linear trend.
    subroutine write_br_discrepancy_reached_to_h5(time_ind, time)

        implicit none

        integer, intent(in) :: time_ind
        real(8), intent(in) :: time

        CALL h5_init()
        CALL h5_open_rw(path2out, h5_id)
        CALL h5_add_string(h5_id, trim(h5_mode_groupname)// &
            '/info', 'discrepancy to linearly predicted value of Br_abs_res > delta')
        CALL h5_add_double_1(h5_id, trim(h5_mode_groupname)// &
            '/discrep_time', (/time_ind*1.d0, time/), (/1/), (/2/))
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
