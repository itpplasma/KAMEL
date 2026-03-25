module logger_m
    implicit none
    private

    integer, parameter :: dp = kind(1.0d0)
    integer, parameter :: FMT_LEN = 256

    ! Log levels
    integer, public, parameter :: LVL_SILENT  = -1
    integer, public, parameter :: LVL_RESULT  =  0
    integer, public, parameter :: LVL_ERROR   =  1
    integer, public, parameter :: LVL_WARNING =  2
    integer, public, parameter :: LVL_INFO    =  3
    integer, public, parameter :: LVL_DEBUG   =  4
    integer, public, parameter :: LVL_TRACE   =  5

    integer :: current_level = LVL_INFO

    public :: set_log_level, get_log_level
    public :: log_error, log_warning, log_info, log_debug, log_trace, log_result

    interface fmt_val
        module procedure fmt_val_real, fmt_val_int, fmt_val_logical, fmt_val_str
    end interface fmt_val
    public :: fmt_val

contains

    subroutine set_log_level(level)
        integer, intent(in) :: level
        current_level = max(LVL_SILENT, min(LVL_TRACE, level))
    end subroutine set_log_level

    function get_log_level() result(level)
        integer :: level
        level = current_level
    end function get_log_level

    subroutine log_error(msg)
        character(*), intent(in) :: msg
        write(0, '(A)') '[ERROR] ' // trim(msg)
        error stop
    end subroutine log_error

    subroutine log_warning(msg)
        character(*), intent(in) :: msg
        if (current_level >= LVL_WARNING) write(0, '(A)') '[WARN ] ' // trim(msg)
    end subroutine log_warning

    subroutine log_info(msg)
        character(*), intent(in) :: msg
        if (current_level >= LVL_INFO) write(6, '(A)') '[INFO ] ' // trim(msg)
    end subroutine log_info

    subroutine log_debug(msg)
        character(*), intent(in) :: msg
        if (current_level >= LVL_DEBUG) write(6, '(A)') '[DEBUG] ' // trim(msg)
    end subroutine log_debug

    subroutine log_trace(msg)
        character(*), intent(in) :: msg
        if (current_level >= LVL_TRACE) write(6, '(A)') '[TRACE] ' // trim(msg)
    end subroutine log_trace

    subroutine log_result(msg)
        character(*), intent(in) :: msg
        if (current_level >= LVL_RESULT) write(6, '(A)') trim(msg)
    end subroutine log_result

    ! --- Format helpers ---

    function fmt_val_real(label, value, unit) result(s)
        character(*), intent(in) :: label
        real(dp), intent(in) :: value
        character(*), intent(in), optional :: unit
        character(len=FMT_LEN) :: s
        character(len=20) :: vbuf
        write(vbuf, '(ES15.8)') value
        if (present(unit)) then
            s = trim(label) // ' = ' // trim(adjustl(vbuf)) // ' ' // trim(unit)
        else
            s = trim(label) // ' = ' // trim(adjustl(vbuf))
        end if
    end function fmt_val_real

    function fmt_val_int(label, value, unit) result(s)
        character(*), intent(in) :: label
        integer, intent(in) :: value
        character(*), intent(in), optional :: unit
        character(len=FMT_LEN) :: s
        character(len=20) :: vbuf
        write(vbuf, '(I0)') value
        if (present(unit)) then
            s = trim(label) // ' = ' // trim(adjustl(vbuf)) // ' ' // trim(unit)
        else
            s = trim(label) // ' = ' // trim(adjustl(vbuf))
        end if
    end function fmt_val_int

    function fmt_val_logical(label, value, unit) result(s)
        character(*), intent(in) :: label
        logical, intent(in) :: value
        character(*), intent(in), optional :: unit
        character(len=FMT_LEN) :: s
        if (present(unit)) then
            s = trim(label) // ' = ' // merge('T', 'F', value) // ' ' // trim(unit)
        else
            s = trim(label) // ' = ' // merge('T', 'F', value)
        end if
    end function fmt_val_logical

    function fmt_val_str(label, value, unit) result(s)
        character(*), intent(in) :: label
        character(*), intent(in) :: value
        character(*), intent(in), optional :: unit
        character(len=FMT_LEN) :: s
        if (present(unit)) then
            s = trim(label) // ' = ' // trim(value) // ' ' // trim(unit)
        else
            s = trim(label) // ' = ' // trim(value)
        end if
    end function fmt_val_str

end module logger_m
