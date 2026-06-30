!> The live subset of io/inout.{h,cpp}'s file-I/O utilities: filename-based
!> functions called from already-Fortran kilca_lib modules (sysmat_profs_m,
!> disp_profs_m, flre_quants_m, background_data_m). The FILE*-taking
!> functions (read_line_2get_*/read_line_2skip_it/trim) and the confirmed-
!> dead functions (save_real_matrix_to_one_file/save_complex_array/
!> load_profile/load_and_alloc_profile/load_complex_profile) are used only
!> by post_processing/progs (still C++, S8 scope) and stay in inout.cpp for
!> now - this is a deliberate partial-file translation, not an oversight;
!> inout.cpp/inout.h are retired fully once those S8 consumers are ported.
module kilca_inout_m
    use, intrinsic :: iso_c_binding, only: c_int, c_double, c_char, c_null_char
    use constants, only: dp
    implicit none
    private

    public :: save_cmplx_matrix_, save_cmplx_matrix_to_one_file_
    public :: save_real_array_, load_data_file_, count_lines_in_file_

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
