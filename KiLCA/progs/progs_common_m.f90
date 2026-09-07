!> Shared helpers for the KiLCA executable entry points (main_linear,
!> main_eig_param, main_post_proc), ported alongside them from the C++ drivers.
!> get_project_path reproduces the drivers' identical "argc==1 ? getcwd :
!> argv[1], then ensure a trailing slash" path setup. fmt_g reproduces C
!> printf "%.<prec>g" formatting (prec=6 for %lg, prec=15 for %.15lg), needed
!> to rebuild byte-for-byte the directory, glob and label strings the drivers
!> emitted with those conversions.
module kilca_progs_common_m
    use, intrinsic :: iso_c_binding, only: c_char, c_size_t, c_ptr, c_null_char, &
        c_int, c_associated
    use, intrinsic :: iso_fortran_env, only: dp => real64
    implicit none
    private

    public :: get_project_path, to_cstr, fmt_g, fmt_e
    public :: open_directory, read_directory, close_directory, clean_run, shell_quote

    interface
        function c_opendir(path) result(directory) bind(C, name="opendir")
            import :: c_char, c_ptr
            character(kind=c_char), intent(in) :: path(*)
            type(c_ptr) :: directory
        end function c_opendir

        function c_readdir_name(directory, name, capacity) result(status) &
                bind(C, name="kilca_readdir_name")
            import :: c_ptr, c_char, c_size_t, c_int
            type(c_ptr), value :: directory
            character(kind=c_char), intent(out) :: name(*)
            integer(c_size_t), value :: capacity
            integer(c_int) :: status
        end function c_readdir_name

        function c_closedir(directory) result(status) bind(C, name="closedir")
            import :: c_ptr, c_int
            type(c_ptr), value :: directory
            integer(c_int) :: status
        end function c_closedir

        function c_getcwd(buf, size) result(res) bind(C, name="getcwd")
            import :: c_char, c_size_t, c_ptr
            character(kind=c_char), intent(out) :: buf(*)
            integer(c_size_t), value :: size
            type(c_ptr) :: res
        end function c_getcwd
    end interface

contains

    function open_directory(path) result(directory)
        character(len=*), intent(in) :: path
        type(c_ptr) :: directory
        directory = c_opendir(to_cstr(path))
    end function open_directory

    !> Preserve readdir order, copying only d_name through the native C adapter.
    logical function read_directory(directory, name) result(found)
        type(c_ptr), intent(in) :: directory
        character(len=*), intent(out) :: name
        character(kind=c_char) :: buffer(len(name) + 1)
        integer :: i
        integer(c_int) :: status

        name = ''
        status = c_readdir_name(directory, buffer, int(size(buffer), c_size_t))
        if (status < 0) error stop 'Failed to read directory entry (or name too long)'
        found = status == 1
        if (.not. found) return
        do i = 1, len(name)
            if (buffer(i) == c_null_char) exit
            name(i:i) = buffer(i)
        end do
    end function read_directory

    subroutine close_directory(directory)
        type(c_ptr), intent(in) :: directory
        integer(c_int) :: status
        status = c_closedir(directory)
        if (status /= 0) error stop 'Failed to close directory'
    end subroutine close_directory

    !> Quote one shell argument, including embedded apostrophes.
    function shell_quote(value) result(quoted)
        character(len=*), intent(in) :: value
        character(len=:), allocatable :: quoted
        integer :: i
        quoted = "'"
        do i = 1, len_trim(value)
            if (value(i:i) == "'") then
                quoted = quoted//"'"//'"'//"'"//'"'//"'"
            else
                quoted = quoted//value(i:i)
            end if
        end do
        quoted = quoted//"'"
    end function shell_quote

    !> Retain the two bracketing eigenfunction folders and the legacy '...' entry.
    !> The production driver and regression test call this same cleanup routine.
    subroutine clean_run(fullpath, funct_first, funct_last)
        character(len=*), intent(in) :: fullpath
        complex(dp), intent(in) :: funct_first, funct_last
        character(len=:), allocatable :: base, dir_path, sub0, sub1
        character(len=1024) :: name
        type(c_ptr) :: directory

        if (len_trim(fullpath) == 0) return
        base = trim(fullpath)
        if (base(len(base):) /= '/') base = base//'/'
        dir_path = base//'linear-data/'
        sub0 = '['//fmt_g(real(funct_first, dp), 15)//','//fmt_g(aimag(funct_first), 15)//']'
        sub1 = '['//fmt_g(real(funct_last, dp), 15)//','//fmt_g(aimag(funct_last), 15)//']'

        directory = open_directory(dir_path)
        if (c_associated(directory)) then
            do while (read_directory(directory, name))
                if (len_trim(name) == 0) cycle
                if (name == '.' .or. name == '..' .or. name == '...') cycle
                if (index(name, sub0) /= 0 .or. index(name, sub1) /= 0) cycle
                call remove_tree(dir_path//trim(name))
            end do
            call close_directory(directory)
        else
            write (*, '(/,a,a,a)') 'clean_run: failed to open the directory ', dir_path, '.'
        end if
        call remove_tree(base//'dispersion-data')
        call remove_tree(base//'poincare-data')
    end subroutine clean_run

    subroutine remove_tree(path)
        character(len=*), intent(in) :: path
        integer :: command_status, exit_status
        if (len_trim(path) == 0) return
        call execute_command_line('rm -R -f -- '//shell_quote(path), wait=.true., &
            cmdstat=command_status, exitstat=exit_status)
        if (command_status /= 0) then
            write (*, '(a)') 'clean_run: failed to start cleanup command'
        else if (exit_status /= 0) then
            write (*, '(a,a)') 'clean_run: failed to remove ', path
        end if
    end subroutine remove_tree

    !> Mirrors the drivers' path bootstrap: no command-line argument means the
    !> current working directory, otherwise the first argument; a trailing '/'
    !> is appended when absent.
    function get_project_path() result(path)
        character(len=:), allocatable :: path
        character(kind=c_char) :: cbuf(1024)
        character(len=1024) :: arg
        type(c_ptr) :: res
        integer :: i, n, alen, stat

        if (command_argument_count() == 0) then
            res = c_getcwd(cbuf, 1024_c_size_t)
            n = 0
            do i = 1, 1024
                if (cbuf(i) == c_null_char) exit
                n = i
            end do
            allocate (character(len=n) :: path)
            do i = 1, n
                path(i:i) = cbuf(i)
            end do
        else
            call get_command_argument(1, arg, alen, stat)
            path = arg(1:alen)
        end if

        if (len(path) > 0) then
            if (path(len(path):len(path)) /= '/') path = path//'/'
        end if
    end function get_project_path

    function to_cstr(s) result(c)
        character(len=*), intent(in) :: s
        character(kind=c_char), allocatable :: c(:)
        integer :: i, n
        n = len_trim(s)
        allocate (c(n + 1))
        do i = 1, n
            c(i) = s(i:i)
        end do
        c(n + 1) = c_null_char
    end function to_cstr

    !> C99 %g of val at the given significant-digit precision: choose %e or %f
    !> by the decimal exponent, then strip trailing fraction zeros (and a bare
    !> trailing point). Matches glibc printf "%.<prec>g" byte for byte.
    function fmt_g(val, prec) result(s)
        real(dp), intent(in) :: val
        integer, intent(in) :: prec
        character(len=:), allocatable :: s
        integer :: p, x, ndec, ie
        character(len=64) :: buf, ebuf, fmtstr
        character(len=:), allocatable :: mant, expo

        p = prec
        if (p == 0) p = 1

        if (val == 0.0_dp) then
            if (sign(1.0_dp, val) < 0.0_dp) then
                s = '-0'
            else
                s = '0'
            end if
            return
        end if

        write (fmtstr, '(a,i0,a)') '(es45.', p - 1, 'e4)'
        write (ebuf, fmtstr) val
        ie = index(ebuf, 'E')
        read (ebuf(ie + 1:), *) x

        if (x >= -4 .and. x < p) then
            ndec = p - 1 - x
            if (ndec < 0) ndec = 0
            write (fmtstr, '(a,i0,a)') '(f0.', ndec, ')'
            write (buf, fmtstr) val
            s = trim(adjustl(buf))
            call strip_trailing_zeros(s)
            call ensure_leading_zero(s)
        else
            ndec = p - 1
            write (fmtstr, '(a,i0,a)') '(es45.', ndec, 'e2)'
            write (buf, fmtstr) val
            buf = adjustl(buf)
            ie = index(buf, 'E')
            mant = trim(buf(1:ie - 1))
            expo = trim(buf(ie + 1:))
            call strip_trailing_zeros(mant)
            call ensure_leading_zero(mant)
            s = mant//'e'//expo
        end if
    end function fmt_g

    !> C printf "%.<prec>e" of val: one digit before the point, prec after,
    !> lowercase 'e', signed exponent with at least two digits. Matches glibc.
    function fmt_e(val, prec) result(s)
        real(dp), intent(in) :: val
        integer, intent(in) :: prec
        character(len=:), allocatable :: s
        character(len=96) :: buf, fmtstr
        integer :: ie

        write (fmtstr, '(a,i0,a,i0,a)') '(es', prec + 12, '.', prec, 'e2)'
        write (buf, fmtstr) val
        s = trim(adjustl(buf))
        ie = index(s, 'E')
        if (ie > 0) s(ie:ie) = 'e'
    end function fmt_e

    subroutine strip_trailing_zeros(s)
        character(len=:), allocatable, intent(inout) :: s
        integer :: i
        if (index(s, '.') == 0) return
        i = len(s)
        do while (i > 1 .and. s(i:i) == '0')
            i = i - 1
        end do
        if (s(i:i) == '.') i = i - 1
        s = s(1:i)
    end subroutine strip_trailing_zeros

    subroutine ensure_leading_zero(s)
        character(len=:), allocatable, intent(inout) :: s
        if (len(s) == 0) return
        if (s(1:1) == '.') then
            s = '0'//s
        else if (len(s) >= 2) then
            if (s(1:1) == '-' .and. s(2:2) == '.') s = '-0'//s(2:)
        end if
    end subroutine ensure_leading_zero

end module kilca_progs_common_m
