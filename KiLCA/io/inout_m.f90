!> The file-I/O utilities ported from io/inout.h and io/inout.cpp: filename-based
!> savers/loaders (save_cmplx_matrix*, save_real_array, load_data_file,
!> count_lines_in_file and their with-comments variants) plus the FILE*-based
!> line readers (read_line_2get_*/read_line_2skip_it) called by the still-C++
!> post_processing and progs entry points. The line readers keep the C
!> "FILE *" handle opaque as type(c_ptr) and read through libc getline, so their
!> C++ callers pass the fopen'd handle through unchanged. The confirmed-dead
!> functions (save_real_matrix_to_one_file/save_complex_array/trim/load_profile/
!> load_and_alloc_profile/load_complex_profile) had zero callers across the
!> built targets and were dropped, completing the translation of inout.cpp.
module kilca_inout_m
    use, intrinsic :: iso_c_binding, only: c_int, c_double, c_char, c_null_char, &
        c_ptr, c_size_t, c_long, c_loc, c_f_pointer
    use constants, only: dp
    implicit none
    private

    public :: save_cmplx_matrix_, save_cmplx_matrix_to_one_file_
    public :: save_real_array_, load_data_file_, count_lines_in_file_
    public :: read_line_2get_double_, read_line_2get_complex_
    public :: read_line_2get_int_, read_line_2get_string_, read_line_2skip_it_
    public :: count_lines_in_file_with_comments_, load_data_file_with_comments_

    interface
        function c_getline(lineptr, n, stream) result(nread) &
            bind(C, name="getline")
            import :: c_ptr, c_long
            type(c_ptr), value :: lineptr
            type(c_ptr), value :: n
            type(c_ptr), value :: stream
            integer(c_long) :: nread
        end function c_getline

        function c_malloc(sz) result(p) bind(C, name="malloc")
            import :: c_ptr, c_size_t
            integer(c_size_t), value :: sz
            type(c_ptr) :: p
        end function c_malloc

        subroutine c_free(p) bind(C, name="free")
            import :: c_ptr
            type(c_ptr), value :: p
        end subroutine c_free
    end interface

contains

    !> Mirrors save_cmplx_matrix exactly: one file per column "i",
    !> path_name_i.dat, each row "x(k) re im re im ...".
    function save_cmplx_matrix_(Nrows, Ncols, Npoints, xgrid, arr, path_name) &
        result(ierr) bind(C, name="save_cmplx_matrix")
        integer(c_int), value :: Nrows, Ncols, Npoints
        real(c_double), intent(in) :: xgrid(0:Npoints - 1)
        real(c_double), intent(in) :: arr(0:2*Nrows*Ncols*Npoints - 1)
        character(kind=c_char), intent(in) :: path_name(*)
        integer(c_int) :: ierr

        character(len=1024) :: fname
        integer :: i, j, k, unit, ios

        do i = 0, Ncols - 1
            write (fname, '(a,a,i0,a)') c_string_to_fortran(path_name), '_', i, '.dat'
            open (newunit=unit, file=trim(fname), status='replace', action='write', iostat=ios)
            if (ios /= 0) then
                write (*, '(a,a)') 'Failed to open file ', trim(fname)
                cycle
            end if

            do k = 0, Npoints - 1
                write (unit, '(es24.16e3)', advance='no') xgrid(k)
                do j = 0, Nrows - 1
                    write (unit, '(a,es24.16e3,a,es24.16e3)', advance='no') char(9), &
                        arr(2*Nrows*Ncols*k + 2*Nrows*i + 2*j), char(9), &
                        arr(2*Nrows*Ncols*k + 2*Nrows*i + 2*j + 1)
                end do
                write (unit, *)
            end do
            close (unit)
        end do

        ierr = 0
    end function save_cmplx_matrix_

    !> Mirrors save_cmplx_matrix_to_one_file: a complex Nrows x Ncols matrix
    !> per point, columns-first (i outer, j inner), to one file.
    function save_cmplx_matrix_to_one_file_(Nrows, Ncols, Npoints, xgrid, arr, full_name) &
        result(ierr) bind(C, name="save_cmplx_matrix_to_one_file")
        integer(c_int), value :: Nrows, Ncols, Npoints
        real(c_double), intent(in) :: xgrid(0:Npoints - 1)
        real(c_double), intent(in) :: arr(0:2*Nrows*Ncols*Npoints - 1)
        character(kind=c_char), intent(in) :: full_name(*)
        integer(c_int) :: ierr

        integer :: i, j, k, unit, ios

        open (newunit=unit, file=trim(c_string_to_fortran(full_name)), status='replace', &
            action='write', iostat=ios)
        if (ios /= 0) then
            write (*, '(a,a)') 'Failed to open file ', trim(c_string_to_fortran(full_name))
            ierr = 1
            return
        end if

        do k = 0, Npoints - 1
            write (unit, '(es24.16e3)', advance='no') xgrid(k)
            do i = 0, Nrows - 1
                do j = 0, Ncols - 1
                    write (unit, '(a,es24.16e3,a,es24.16e3)', advance='no') char(9), &
                        arr(2*Nrows*Ncols*k + 2*Nrows*j + 2*i), char(9), &
                        arr(2*Nrows*Ncols*k + 2*Nrows*j + 2*i + 1)
                end do
            end do
            write (unit, *)
        end do
        close (unit)

        ierr = 0
    end function save_cmplx_matrix_to_one_file_

    !> Mirrors save_real_array: "x(k) arr(k)" per line.
    function save_real_array_(dim_, xgrid, arr, full_name) result(ierr) &
        bind(C, name="save_real_array")
        integer(c_int), value :: dim_
        real(c_double), intent(in) :: xgrid(0:dim_ - 1), arr(0:dim_ - 1)
        character(kind=c_char), intent(in) :: full_name(*)
        integer(c_int) :: ierr

        integer :: k, unit, ios

        open (newunit=unit, file=trim(c_string_to_fortran(full_name)), status='replace', &
            action='write', iostat=ios)
        if (ios /= 0) then
            write (*, '(a,a)') 'Failed to open file ', trim(c_string_to_fortran(full_name))
            ierr = 1
            return
        end if

        do k = 0, dim_ - 1
            write (unit, '(es24.16e3,a,es24.16e3)') xgrid(k), char(9), arr(k)
        end do
        close (unit)

        ierr = 0
    end function save_real_array_

    !> Mirrors count_lines_in_file: counts lines via a raw read loop
    !> (matching the oracle's getline-until-EOF convention, not a
    !> line-ending-aware text scan).
    function count_lines_in_file_(filename, flag_print) result(num_lines) &
        bind(C, name="count_lines_in_file")
        character(kind=c_char), intent(in) :: filename(*)
        integer(c_int), value :: flag_print
        integer(c_int) :: num_lines

        character(len=1024) :: fname
        integer :: unit, ios
        character(len=65536) :: line

        fname = c_string_to_fortran(filename)
        open (newunit=unit, file=trim(fname), status='old', action='read', iostat=ios)
        if (ios /= 0) then
            write (*, '(a,a)') 'count_lines_in_file: failed to open file ', trim(fname)
        end if

        num_lines = 0
        do
            read (unit, '(a)', iostat=ios) line
            if (ios /= 0) exit
            num_lines = num_lines + 1
        end do

        close (unit)

        if (flag_print /= 0) then
            write (*, '(a,a,a,i0,a)') 'file ', trim(fname), ' contains ', num_lines, ' non-empty lines'
        end if
    end function count_lines_in_file_

    !> Mirrors load_data_file exactly: dim must equal count_lines_in_file's
    !> result; each line is (x, y_1, ..., y_ncols), columns stored
    !> column-major in qgrid (qgrid(i + k*dim) for column k, row i).
    function load_data_file_(file_name, dim_, ncols, rgrid, qgrid) result(ierr) &
        bind(C, name="load_data_file")
        character(kind=c_char), intent(in) :: file_name(*)
        integer(c_int), value :: dim_, ncols
        real(c_double), intent(out) :: rgrid(0:dim_ - 1)
        real(c_double), intent(out) :: qgrid(0:dim_*ncols - 1)
        integer(c_int) :: ierr

        character(len=1024) :: fname
        integer :: unit, ios, i, k
        real(dp) :: vals(0:ncols)

        fname = c_string_to_fortran(file_name)

        if (dim_ /= count_lines_in_file_(file_name, 0_c_int)) then
            write (*, '(a,a)') &
                'error: load_data_file: read error or false input data dimension: ', trim(fname)
            stop 1
        end if

        open (newunit=unit, file=trim(fname), status='old', action='read', iostat=ios)
        if (ios /= 0) then
            write (*, '(a,a)') 'load_data_file: failed to open file ', trim(fname)
            stop 1
        end if

        ierr = 0
        do i = 0, dim_ - 1
            read (unit, *, iostat=ios) vals(0:ncols)
            if (ios /= 0) then
                write (*, '(a,a,a)') 'load_data_file: file ', trim(fname), ' reading error!'
                stop 1
            end if
            rgrid(i) = vals(0)
            do k = 0, ncols - 1
                qgrid(i + k*dim_) = vals(k + 1)
            end do
        end do

        close (unit)
    end function load_data_file_

    !> Mirrors read_line_2get_double: reads one line from the C FILE* handle,
    !> takes the text before '#' and parses it as a double (0 on no conversion,
    !> matching strtod).
    subroutine read_line_2get_double_(in, value) bind(C, name="read_line_2get_double")
        type(c_ptr), value :: in
        real(c_double), intent(out) :: value

        character(len=4096) :: sub
        logical :: ok
        integer :: ios

        call read_line_before_hash(in, sub, ok)
        if (.not. ok) write (*, '(a)') 'read_line_2get_double: file reading error'
        read (sub, *, iostat=ios) value
        if (ios /= 0) value = 0.0_dp
    end subroutine read_line_2get_double_

    !> Mirrors read_line_2get_complex: the text before '#' is "(re,im)", which is
    !> exactly Fortran list-directed complex input.
    subroutine read_line_2get_complex_(in, value) bind(C, name="read_line_2get_complex")
        type(c_ptr), value :: in
        complex(c_double), intent(out) :: value

        character(len=4096) :: sub
        logical :: ok
        integer :: ios

        call read_line_before_hash(in, sub, ok)
        if (.not. ok) write (*, '(a)') 'read_line_2get_complex: file reading error'
        read (sub, *, iostat=ios) value
        if (ios /= 0) value = (0.0_dp, 0.0_dp)
    end subroutine read_line_2get_complex_

    !> Mirrors read_line_2get_int: text before '#' parsed as an integer.
    subroutine read_line_2get_int_(in, value) bind(C, name="read_line_2get_int")
        type(c_ptr), value :: in
        integer(c_int), intent(out) :: value

        character(len=4096) :: sub
        logical :: ok
        integer :: ios

        call read_line_before_hash(in, sub, ok)
        if (.not. ok) write (*, '(a)') 'read_line_2get_int: file reading error'
        read (sub, *, iostat=ios) value
        if (ios /= 0) value = 0_c_int
    end subroutine read_line_2get_int_

    !> Mirrors read_line_2get_string: the first whitespace/tab-delimited token of
    !> the text before '#' is copied (null-terminated) into the caller's buffer.
    !> "value" is a C "char **"; *value is the destination char buffer.
    subroutine read_line_2get_string_(in, value) bind(C, name="read_line_2get_string")
        type(c_ptr), value :: in
        type(c_ptr), value :: value

        character(len=4096) :: sub
        logical :: ok
        integer :: i, s, e, n
        type(c_ptr), pointer :: destpp
        character(kind=c_char), pointer :: dest(:)

        call read_line_before_hash(in, sub, ok)
        if (.not. ok) write (*, '(a)') 'read_line_2get_string: file reading error'

        s = 0
        do i = 1, len(sub)
            if (sub(i:i) /= ' ') then
                if (sub(i:i) /= achar(9)) then
                    s = i
                    exit
                end if
            end if
        end do

        n = 0
        if (s > 0) then
            e = s
            do i = s, len(sub)
                if (sub(i:i) == ' ') exit
                if (sub(i:i) == achar(9)) exit
                e = i
            end do
            n = e - s + 1
        end if

        call c_f_pointer(value, destpp)
        call c_f_pointer(destpp, dest, [n + 1])
        do i = 1, n
            dest(i) = sub(s + i - 1:s + i - 1)
        end do
        dest(n + 1) = c_null_char
    end subroutine read_line_2get_string_

    !> Mirrors read_line_2skip_it: read and discard one line. "value" is unused,
    !> matching the oracle.
    subroutine read_line_2skip_it_(in, value) bind(C, name="read_line_2skip_it")
        type(c_ptr), value :: in
        type(c_ptr), value :: value

        character(len=4096) :: line
        integer(c_long) :: nread

        call read_raw_line(in, line, nread)
    end subroutine read_line_2skip_it_

    !> Mirrors count_lines_in_file_with_comments: counts lines that contain none
    !> of '%', '#', '!'.
    function count_lines_in_file_with_comments_(filename, flag_print) &
        result(num_lines) bind(C, name="count_lines_in_file_with_comments")
        character(kind=c_char), intent(in) :: filename(*)
        integer(c_int), value :: flag_print
        integer(c_int) :: num_lines

        character(len=1024) :: fname
        integer :: unit, ios
        character(len=65536) :: line

        fname = c_string_to_fortran(filename)
        open (newunit=unit, file=trim(fname), status='old', action='read', iostat=ios)
        if (ios /= 0) then
            write (*, '(a,a)') 'count_lines_in_file: failed to open file ', trim(fname)
        end if

        num_lines = 0
        do
            read (unit, '(a)', iostat=ios) line
            if (ios /= 0) exit
            if (index(line, '%') == 0 .and. index(line, '#') == 0 &
                .and. index(line, '!') == 0) then
                num_lines = num_lines + 1
            end if
        end do

        close (unit)

        if (flag_print /= 0) then
            write (*, '(a,a,a,i0,a)') 'file ', trim(fname), ' contains ', num_lines, &
                ' non-empty lines'
        end if
    end function count_lines_in_file_with_comments_

    !> Mirrors load_data_file_with_comments exactly: dim is the comment-aware line
    !> count; comment lines (containing '%', '#', or '!') are skipped. The loop
    !> bound dimf is the raw line count (count_lines_in_file), preserving the
    !> oracle's mix of the two counters.
    function load_data_file_with_comments_(file_name, dim_, ncols, rgrid, qgrid) &
        result(ierr) bind(C, name="load_data_file_with_comments")
        character(kind=c_char), intent(in) :: file_name(*)
        integer(c_int), value :: dim_, ncols
        real(c_double), intent(out) :: rgrid(0:dim_ - 1)
        real(c_double), intent(out) :: qgrid(0:dim_*ncols - 1)
        integer(c_int) :: ierr

        character(len=1024) :: fname
        integer :: unit, ios, i, k, ind, dimf
        real(dp) :: vals(0:ncols)
        character(len=65536) :: line

        fname = c_string_to_fortran(file_name)

        if (dim_ /= count_lines_in_file_with_comments_(file_name, 0_c_int)) then
            write (*, '(a,a)') 'error: load_data_file: read error or false input '// &
                'data dimension: ', trim(fname)
            stop 1
        end if

        open (newunit=unit, file=trim(fname), status='old', action='read', iostat=ios)
        if (ios /= 0) then
            write (*, '(a,a)') 'load_data_file: failed to open file ', trim(fname)
            stop 1
        end if

        dimf = count_lines_in_file_(file_name, 0_c_int)
        ind = 0
        do i = 0, dimf - 1
            read (unit, '(a)', iostat=ios) line
            if (ios /= 0) then
                write (*, '(a,a,a)') 'load_data_file: file ', trim(fname), &
                    ' reading error!'
                stop 1
            end if

            if (index(line, '%') /= 0 .or. index(line, '#') /= 0 &
                .or. index(line, '!') /= 0) cycle

            read (line, *, iostat=ios) vals(0:ncols)
            if (ios /= 0) then
                write (*, '(a,a,a)') 'error: load_data_file: ', trim(fname), &
                    ': read error or false data dimension!'
                stop 1
            end if

            rgrid(ind) = vals(0)
            do k = 0, ncols - 1
                qgrid(ind + k*dim_) = vals(k + 1)
            end do
            ind = ind + 1
        end do

        if (ind /= dim_) then
            write (*, '(a,a,a)') 'error: load_data_file: ', trim(fname), &
                ': read error or false data dimension!'
            stop 1
        end if

        close (unit)

        ierr = 0
    end function load_data_file_with_comments_

    !> Reads one line from a C FILE* via libc getline into a Fortran string,
    !> stripping the trailing newline. nread is getline's return (-1 at EOF).
    subroutine read_raw_line(stream, line, nread)
        type(c_ptr), value :: stream
        character(len=*), intent(out) :: line
        integer(c_long), intent(out) :: nread

        type(c_ptr), target :: buf
        integer(c_size_t), target :: nbytes
        character(kind=c_char), pointer :: cbuf(:)
        integer :: i, m

        buf = c_malloc(1024_c_size_t)
        nbytes = 1024_c_size_t
        nread = c_getline(c_loc(buf), c_loc(nbytes), stream)

        line = ''
        if (nread > 0) then
            m = int(nread)
            call c_f_pointer(buf, cbuf, [m])
            do i = 1, m
                if (cbuf(i) == c_null_char) exit
                if (cbuf(i) == achar(10)) exit
                if (i > len(line)) exit
                line(i:i) = cbuf(i)
            end do
        end if

        call c_free(buf)
    end subroutine read_raw_line

    !> Reads one line and returns the text before the first '#'. ok is false on a
    !> short/failed read (getline < 2), matching the oracle's warning condition.
    subroutine read_line_before_hash(stream, sub, ok)
        type(c_ptr), value :: stream
        character(len=*), intent(out) :: sub
        logical, intent(out) :: ok

        character(len=4096) :: line
        integer(c_long) :: nread
        integer :: p

        call read_raw_line(stream, line, nread)
        ok = nread >= 2
        p = index(line, '#')
        if (p > 0) then
            sub = line(1:p - 1)
        else
            sub = line
        end if
    end subroutine read_line_before_hash

    function c_string_to_fortran(cstr) result(fstr)
        character(kind=c_char), intent(in) :: cstr(*)
        character(len=1024) :: fstr
        integer :: i
        fstr = ''
        i = 0
        do
            if (cstr(i + 1) == c_null_char .or. i >= 1024) exit
            fstr(i + 1:i + 1) = cstr(i + 1)
            i = i + 1
        end do
    end function c_string_to_fortran

end module kilca_inout_m
