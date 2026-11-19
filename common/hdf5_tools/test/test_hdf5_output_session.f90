!> An output session keeps one file open across many open/close pairs, so a
!> time step that calls several writers pays for one open and one close.
!> Writers keep calling h5_open_rw, h5_close and h5_deinit as before.
program test_hdf5_output_session
    use iso_fortran_env, only: real64
    use KAMEL_hdf5_tools, only: HID_T, h5_add, h5_close, h5_create, h5_deinit, h5_get, &
        h5_obj_exists, h5_open, h5_open_rw, h5_begin_output_session, &
        h5_end_output_session, h5_output_session_active, h5_stream_write, &
        h5_truncate_existing

    implicit none

    character(len=*), parameter :: path = "test_hdf5_output_session.h5"
    character(len=*), parameter :: other_path = "test_hdf5_output_session_other.h5"
    character(len=*), parameter :: stop_path = "test_hdf5_output_session_stop.h5"
    character(len=256) :: mode

    h5_stream_write = .true.
    h5_truncate_existing = .true.

    call get_command_argument(1, mode)
    if (trim(mode) == "stop-inside-session") then
        call write_then_stop_inside_session
    end if

    call check_writers_share_one_handle
    call check_deinit_keeps_session_open
    call check_nested_sessions
    call check_other_files_unaffected
    call check_stop_inside_session_keeps_data
    print '(a)', "test_hdf5_output_session: all checks passed"

contains

    subroutine write_value(file_path, name, value)
        character(len=*), intent(in) :: file_path, name
        real(real64), intent(in) :: value
        integer(HID_T) :: h5id

        call h5_open_rw(file_path, h5id)
        call h5_add(h5id, name, value)
        call h5_close(h5id)
        call h5_deinit()
    end subroutine write_value

    subroutine require_value(file_path, name, expected, label)
        character(len=*), intent(in) :: file_path, name, label
        real(real64), intent(in) :: expected
        integer(HID_T) :: h5id
        logical :: exists
        real(real64) :: value

        call h5_open(file_path, h5id)
        call h5_obj_exists(h5id, name, exists)
        if (.not. exists) then
            print '(a,a,a,a)', "FAIL: ", label, ": missing ", name
            error stop 1
        end if
        call h5_get(h5id, name, value)
        call h5_close(h5id)
        call h5_deinit()
        if (value /= expected) then
            print '(a,a,a,es12.4,a,es12.4)', "FAIL: ", label, ": got ", value, &
                ", expected ", expected
            error stop 1
        end if
    end subroutine require_value

    subroutine fresh_file(file_path)
        character(len=*), intent(in) :: file_path
        integer(HID_T) :: h5id

        call h5_create(file_path, h5id)
        call h5_close(h5id)
        call h5_deinit()
    end subroutine fresh_file

    subroutine check_writers_share_one_handle
        integer(HID_T) :: first, second

        call fresh_file(path)
        call h5_begin_output_session(path)
        call h5_open_rw(path, first)
        call h5_close(first)
        call h5_open_rw(path, second)
        call h5_close(second)
        if (first /= second) then
            print '(a)', "FAIL: h5_open_rw inside a session did not return the session handle"
            error stop 1
        end if
        call write_value(path, "a", 1.0_real64)
        call write_value(path, "b", 2.0_real64)
        call h5_end_output_session()

        call require_value(path, "a", 1.0_real64, "first writer in session")
        call require_value(path, "b", 2.0_real64, "second writer in session")
    end subroutine check_writers_share_one_handle

    subroutine check_deinit_keeps_session_open
        call fresh_file(path)
        call h5_begin_output_session(path)
        call write_value(path, "before", 1.0_real64)
        call h5_deinit()
        if (.not. h5_output_session_active()) then
            print '(a)', "FAIL: h5_deinit ended the output session"
            error stop 1
        end if
        call write_value(path, "after", 2.0_real64)
        call h5_end_output_session()

        call require_value(path, "before", 1.0_real64, "write before h5_deinit")
        call require_value(path, "after", 2.0_real64, "write after h5_deinit")
    end subroutine check_deinit_keeps_session_open

    subroutine check_nested_sessions
        call fresh_file(path)
        call h5_begin_output_session(path)
        call h5_begin_output_session(path)
        call write_value(path, "inner", 1.0_real64)
        call h5_end_output_session()
        if (.not. h5_output_session_active()) then
            print '(a)', "FAIL: ending an inner session closed the outer one"
            error stop 1
        end if
        call write_value(path, "outer", 2.0_real64)
        call h5_end_output_session()
        if (h5_output_session_active()) then
            print '(a)', "FAIL: session still active after the outermost end"
            error stop 1
        end if

        call require_value(path, "inner", 1.0_real64, "nested session, inner write")
        call require_value(path, "outer", 2.0_real64, "nested session, outer write")
    end subroutine check_nested_sessions

    subroutine check_other_files_unaffected
        call fresh_file(path)
        call fresh_file(other_path)
        call h5_begin_output_session(path)
        call write_value(other_path, "x", 3.0_real64)
        call require_value(other_path, "x", 3.0_real64, "other file during session")
        call write_value(path, "y", 4.0_real64)
        call h5_end_output_session()

        call require_value(path, "y", 4.0_real64, "session file next to another file")
    end subroutine check_other_files_unaffected

    !> QL-Balance ends runs with stop from inside a step. The session must still
    !> write its file before the process exits.
    subroutine check_stop_inside_session_keeps_data
        character(len=4096) :: self
        integer :: exit_status

        call fresh_file(stop_path)
        call get_command_argument(0, self)
        call execute_command_line(trim(self)//" stop-inside-session", &
            exitstat=exit_status)
        if (exit_status /= 0) then
            print '(a,i0)', "FAIL: child process exited with ", exit_status
            error stop 1
        end if
        call require_value(stop_path, "last", 5.0_real64, "write before stop in session")
    end subroutine check_stop_inside_session_keeps_data

    subroutine write_then_stop_inside_session
        call h5_begin_output_session(stop_path)
        call write_value(stop_path, "last", 5.0_real64)
        stop
    end subroutine write_then_stop_inside_session

end program test_hdf5_output_session
