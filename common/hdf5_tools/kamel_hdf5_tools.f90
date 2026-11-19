!> KAMEL's HDF5 interface: fortio's hdf5_tools plus output sessions.
!>
!> With fortio, every close of a written file rewrites the metadata of all its
!> datasets. Code that opens and closes the same output file for each small
!> write (QL-Balance does this about 15 times per time step) then pays
!> O(datasets) per write. An output session keeps one file open: inside it,
!> h5_open_rw of that file returns the session handle, and h5_close of that
!> handle and h5_deinit do nothing. The file is closed, and its metadata
!> written, once when the outermost session ends.
!>
!> Writes inside a session reach the file at the end of the session. A session
!> still open when the program stops (stop or error stop) is closed by an exit
!> handler; only a killed process loses the writes of its open session.
module KAMEL_hdf5_tools
    use, intrinsic :: iso_c_binding, only: c_funptr, c_funloc, c_int
    use hdf5_tools, fortio_h5_open_rw => h5_open_rw, fortio_h5_close => h5_close, &
        fortio_h5_deinit => h5_deinit
    implicit none

    private :: fortio_h5_open_rw, fortio_h5_close, fortio_h5_deinit

    public :: h5_open_rw, h5_close, h5_deinit
    public :: h5_begin_output_session, h5_end_output_session, h5_output_session_active

    character(len=:), allocatable, save, private :: session_path
    integer(HID_T), save, private :: session_id = -1_HID_T
    integer, save, private :: session_depth = 0
    logical, save, private :: exit_handler_registered = .false.

    interface
        function c_atexit(handler) bind(C, name="atexit") result(status)
            import :: c_funptr, c_int
            type(c_funptr), value :: handler
            integer(c_int) :: status
        end function c_atexit
    end interface

contains

    !> Open filename for the session, or join the session already open for it.
    subroutine h5_begin_output_session(filename)
        character(len=*), intent(in) :: filename

        if (session_depth > 0) then
            if (trim(filename) /= session_path) &
                error stop "h5_begin_output_session: another output session is open"
            session_depth = session_depth + 1
            return
        end if
        call fortio_h5_open_rw(trim(filename), session_id)
        session_path = trim(filename)
        session_depth = 1
        if (.not. exit_handler_registered) then
            if (c_atexit(c_funloc(close_output_session_at_exit)) /= 0_c_int) &
                error stop "h5_begin_output_session: could not register the exit handler"
            exit_handler_registered = .true.
        end if
    end subroutine h5_begin_output_session

    !> Leave the session; the outermost end closes the file.
    subroutine h5_end_output_session()
        if (session_depth == 0) error stop "h5_end_output_session: no output session is open"
        session_depth = session_depth - 1
        if (session_depth > 0) return
        call fortio_h5_close(session_id)
        session_id = -1_HID_T
        ! Releases handles that writers left open, as their own h5_deinit would have.
        call fortio_h5_deinit()
    end subroutine h5_end_output_session

    logical function h5_output_session_active()
        h5_output_session_active = session_depth > 0
    end function h5_output_session_active

    subroutine close_output_session_at_exit() bind(C)
        if (session_depth == 0) return
        session_depth = 1
        call h5_end_output_session()
    end subroutine close_output_session_at_exit

    subroutine h5_open_rw(filename, h5id, opt_fileformat_version)
        character(len=*), intent(in) :: filename
        integer(HID_T), intent(out) :: h5id
        integer, intent(in), optional :: opt_fileformat_version

        if (session_depth > 0) then
            if (trim(filename) == session_path) then
                h5id = session_id
                return
            end if
        end if
        call fortio_h5_open_rw(filename, h5id, opt_fileformat_version)
    end subroutine h5_open_rw

    subroutine h5_close(h5id)
        integer(HID_T), intent(inout) :: h5id

        if (session_depth > 0 .and. h5id == session_id) then
            h5id = -1_HID_T
            return
        end if
        call fortio_h5_close(h5id)
    end subroutine h5_close

    subroutine h5_deinit()
        if (session_depth > 0) return
        call fortio_h5_deinit()
    end subroutine h5_deinit

end module KAMEL_hdf5_tools
