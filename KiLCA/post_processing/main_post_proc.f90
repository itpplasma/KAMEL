!> Post-processing program to evaluate form-factors in various radial variables,
!> ported from post_processing/main_post_proc.cpp. Reads post_proc.in and
!> modes.in, loads the per-mode EB.dat (and, for the current-density output, the
!> zone_0_current_dens files) written by a prior KiLCA run, interpolates the
!> requested form-factor onto a surface-label grid, writes the human-readable
!> "<conf>.<quantity>.<label>" table and hands the packed complex form-factors to
!> the F77 writer save_form_facs_fortran for the field-line-tracing (.uff/.fff)
!> file. The library file I/O (read_line_2get_*, count_lines_in_file,
!> load_data_file*) and the polynomial interpolation (find_index_for_interp,
!> eval_neville_polynom) are the already-Fortran KiLCA library routines.
program main_post_proc
    use, intrinsic :: iso_c_binding, only: c_int, c_double, c_char, c_ptr, &
        c_null_char, c_null_ptr, c_loc, c_associated, c_f_pointer
    use, intrinsic :: iso_fortran_env, only: dp => real64
    use kilca_inout_m, only: read_line_2skip_it_, read_line_2get_string_, &
        read_line_2get_int_, read_line_2get_double_, read_line_2get_complex_, &
        count_lines_in_file_, load_data_file_, &
        count_lines_in_file_with_comments_, load_data_file_with_comments_
    use adaptive_grid_pol_m, only: find_index_for_interp, eval_neville_polynom
    use kilca_progs_common_m, only: fmt_g, to_cstr
    implicit none

    type :: darr_t
        real(dp), allocatable :: v(:)
    end type darr_t

    interface
        function c_fopen(path, mode) result(fp) bind(C, name="fopen")
            import :: c_char, c_ptr
            character(kind=c_char), intent(in) :: path(*), mode(*)
            type(c_ptr) :: fp
        end function c_fopen

        function c_fclose(fp) result(res) bind(C, name="fclose")
            import :: c_ptr, c_int
            type(c_ptr), value :: fp
            integer(c_int) :: res
        end function c_fclose
    end interface

    external :: save_form_facs_fortran

    character(len=32) :: name_out(0:2), name_label(0:6)
    character(len=1024) :: conf, path2projects(0:1), filename, fname_fort
    character(len=1024) :: eb_file
    type(c_ptr) :: in
    integer(c_int) :: dim, dma, imin, imax, m_min, m_max, n_min, n_max
    integer(c_int) :: pol_deg, output_flag, label_flag, format_flag
    real(c_double) :: label_min, label_max, rmini
    complex(c_double) :: flab
    real(dp), allocatable :: label(:), rout(:)
    integer, allocatable :: modes(:), modes_mn(:), modes_index(:)
    complex(dp), allocatable :: form_facs(:)
    complex(dp), allocatable :: ffpsi(:, :)
    type(darr_t) :: r_EB(0:1), prof_EB(0:1), r_jp(0:1), prof_jp(0:1)
    integer(c_int) :: dim_EB(0:1), dim_jp(0:1)
    integer, parameter :: ncols_EB = 12, ncols_jp = 2
    integer :: i, k, l, p, ind, m, n, unit, ios, nread_modes, modes_dim
    character(len=256) :: mline

    name_out(0) = 'Br_ff'
    name_out(1) = 'Jp_over_Br_over_r'
    name_out(2) = 'Tphi_ff'
    name_label(0) = 'cyl_radius'
    name_label(1) = 'sqrt_psi_pol'
    name_label(2) = 'psi_pol'
    name_label(3) = 'norm_psi_pol'
    name_label(4) = 'sqrt_psi_tor'
    name_label(5) = 'psi_tor'
    name_label(6) = 'norm_psi_tor'

    in = c_fopen(to_cstr('post_proc.in'), to_cstr('r'))
    if (.not. c_associated(in)) then
        write (*, '(/,a,a,a)') 'error: post_proc: failed to open file ', 'post_proc.in', achar(7)
        stop
    end if

    call read_line_2skip_it_(in, c_null_ptr)
    call read_string(in, conf)
    call read_string(in, path2projects(0))
    call read_string(in, path2projects(1))
    call read_line_2skip_it_(in, c_null_ptr)

    call read_line_2skip_it_(in, c_null_ptr)
    call read_line_2get_int_(in, dim)
    call read_line_2get_double_(in, label_min)
    call read_line_2get_double_(in, label_max)
    call read_line_2get_double_(in, rmini)
    call read_line_2skip_it_(in, c_null_ptr)

    call read_line_2skip_it_(in, c_null_ptr)
    call read_line_2get_complex_(in, flab)
    call read_line_2get_int_(in, dma)
    call read_line_2get_int_(in, imin)
    call read_line_2get_int_(in, imax)
    call read_line_2get_int_(in, m_min)
    call read_line_2get_int_(in, m_max)
    call read_line_2get_int_(in, n_min)
    call read_line_2get_int_(in, n_max)
    call read_line_2skip_it_(in, c_null_ptr)

    call read_line_2skip_it_(in, c_null_ptr)
    call read_line_2get_int_(in, pol_deg)
    call read_line_2skip_it_(in, c_null_ptr)

    call read_line_2skip_it_(in, c_null_ptr)
    call read_line_2get_int_(in, output_flag)
    call read_line_2skip_it_(in, c_null_ptr)

    call read_line_2skip_it_(in, c_null_ptr)
    call read_line_2get_int_(in, label_flag)
    call read_line_2skip_it_(in, c_null_ptr)

    call read_line_2skip_it_(in, c_null_ptr)
    call read_line_2get_int_(in, format_flag)
    call read_line_2skip_it_(in, c_null_ptr)

    ios = c_fclose(in)

    allocate (label(0:dim - 1), rout(0:dim - 1))
    do k = 0, dim - 1
        label(k) = label_min + k*(label_max - label_min)/(dim - 1)
    end do

    call transform_surface_label_to_cyl_radius(dim, rout, pol_deg, label_flag, label)

    open (newunit=unit, file='modes.in', status='old', action='read', iostat=ios)
    if (ios /= 0) then
        write (*, '(/,a,a,a)') 'error: failed to open file ', 'modes.in', achar(7)
        stop
    end if

    allocate (modes(0:2*dma - 1))
    nread_modes = 0
    do i = 0, dma - 1
        read (unit, '(a)', iostat=ios) mline
        if (ios /= 0) exit
        call parse_mode(mline, modes(2*i), modes(2*i + 1))
        nread_modes = i + 1
    end do
    close (unit)

    if (nread_modes /= dma) then
        write (*, '(/,a)') 'error: read error or false modes array dimension!'
        stop
    end if

    dim_EB = 0
    dim_jp = 0

    allocate (form_facs(0:(imax - imin + 1)*dim - 1))

    do i = imin, imax
        write (*, '(/,a,i0,a,i0,a)', advance='no') &
            '(m=', modes(2*i), ', n=', modes(2*i + 1), ') mode...'

        do p = 0, 1
            eb_file = trim(path2projects(p))//'linear-data/m_'// &
                itoa(modes(2*i))//'_n_'//itoa(modes(2*i + 1))//'_flab_['// &
                fmt_g(real(flab, dp), 6)//','//fmt_g(aimag(flab), 6)//']/EB.dat'
            dim_EB(p) = count_lines_in_file_(to_cstr(eb_file), 0_c_int)
            allocate (r_EB(p)%v(0:dim_EB(p) - 1))
            allocate (prof_EB(p)%v(0:ncols_EB*dim_EB(p) - 1))
            ios = load_data_file_(to_cstr(eb_file), dim_EB(p), ncols_EB, &
                                  r_EB(p)%v, prof_EB(p)%v)
        end do

        if (output_flag == 1) then
            p = 0
            eb_file = trim(path2projects(0))//'linear-data/m_'// &
                itoa(modes(2*i))//'_n_'//itoa(modes(2*i + 1))//'_flab_['// &
                fmt_g(real(flab, dp), 6)//','//fmt_g(aimag(flab), 6)// &
                ']/zone_0_current_dens_p_0_t_lab.dat'
            dim_jp(0) = count_lines_in_file_(to_cstr(eb_file), 0_c_int)
            allocate (r_jp(0)%v(0:dim_jp(0) - 1))
            allocate (prof_jp(0)%v(0:ncols_jp*dim_jp(0) - 1))
            ios = load_data_file_(to_cstr(eb_file), dim_jp(0), ncols_jp, &
                                  r_jp(0)%v, prof_jp(0)%v)
        end if

        do k = 0, dim - 1
            select case (output_flag)
            case (0)
                call evaluate_Br_formfactors(dim_EB(0), r_EB(0)%v, prof_EB(0)%v, &
                    dim_EB(1), r_EB(1)%v, prof_EB(1)%v, pol_deg, rout(k), &
                    form_facs(k + (i - imin)*dim))
            case (1)
                call evaluate_jp_over_Br(dim_EB(0), r_EB(0)%v, prof_EB(0)%v, &
                    dim_jp(0), r_jp(0)%v, prof_jp(0)%v, pol_deg, rout(k), &
                    form_facs(k + (i - imin)*dim))
            case (2)
                write (*, '(/,a)') 'error: not implemented!'
            case default
                write (*, '(/,a,i0)') 'error: unknown output flag = ', output_flag
                stop 1
            end select
        end do

        do l = 0, dim - 1
            if (rout(l) > rmini) exit
        end do

        do k = 0, l - 1
            form_facs(k + (i - imin)*dim) = form_facs(l + (i - imin)*dim) + &
                (rout(k) - rout(l))/(rout(l + 1) - rout(l))* &
                (form_facs(l + 1 + (i - imin)*dim) - form_facs(l + (i - imin)*dim))
        end do

        do p = 0, 1
            if (allocated(r_EB(p)%v)) deallocate (r_EB(p)%v)
            if (allocated(prof_EB(p)%v)) deallocate (prof_EB(p)%v)
            if (allocated(r_jp(p)%v)) deallocate (r_jp(p)%v)
            if (allocated(prof_jp(p)%v)) deallocate (prof_jp(p)%v)
        end do

        write (*, '(a)', advance='no') ' Ok.'
    end do

    write (*, '(/,a)', advance='no') 'Saving...'

    call transform_cyl_radius_to_surface_label(dim, rout, pol_deg, label_flag, label)

    filename = trim(conf)//'.'//trim(name_out(output_flag))//'.'//trim(name_label(label_flag))

    open (newunit=unit, file=trim(filename), status='replace', action='write', iostat=ios)
    if (ios /= 0) write (*, '(/,a,a,a)') 'Failed to open file ', trim(filename), achar(7)

    write (unit, '(a)', advance='no') '%modes:'
    do i = imin, imax
        write (unit, '(a,a,i0,a,i0,a)', advance='no') &
            char(9), '(', modes(2*i), ', ', modes(2*i + 1), ')'
    end do
    write (unit, '(a)') ''

    write (unit, '(a)', advance='no') '%r'
    do i = imin, imax
        write (unit, '(a,a)', advance='no') char(9), 'real(T_{mn}) imag(T_{mn})'
    end do
    write (unit, '(a)') ''

    do k = 0, dim - 1
        write (unit, '(es24.16e3)', advance='no') label(k)
        do i = imin, imax
            write (unit, '(a,es24.16e3,a,es24.16e3)', advance='no') char(9), &
                real(form_facs(k + (i - imin)*dim), dp), char(9), &
                aimag(form_facs(k + (i - imin)*dim))
        end do
        write (unit, '(a)') ''
    end do
    close (unit)

    modes_dim = imax - imin + 1
    allocate (modes_mn(0:2*modes_dim - 1))
    ind = 0
    do i = imin, imax
        modes_mn(2*ind) = modes(2*i)
        modes_mn(2*ind + 1) = modes(2*i + 1)
        ind = ind + 1
    end do

    allocate (modes_index(0:(m_max - m_min + 1)*(n_max - n_min + 1) - 1))
    do m = m_min, m_max
        do n = n_min, n_max
            modes_index(m - m_min + (m_max - m_min + 1)*(n - n_min)) = -1
        end do
    end do
    do i = 0, modes_dim - 1
        m = modes_mn(2*i)
        n = modes_mn(2*i + 1)
        modes_index(m - m_min + (m_max - m_min + 1)*(n - n_min)) = i + 1
    end do

    if (format_flag == 0) then
        fname_fort = trim(filename)//'.uff'
    else
        fname_fort = trim(filename)//'.fff'
    end if

    allocate (ffpsi(0:dim - 1, 0:modes_dim - 1))
    do i = imin, imax
        do k = 0, dim - 1
            ffpsi(k, i - imin) = form_facs(k + (i - imin)*dim)
        end do
    end do

    call save_form_facs_fortran(fname_fort, modes_dim, modes_mn, m_min, m_max, &
        n_min, n_max, modes_index, dim, label, ffpsi, format_flag)

    write (*, '(a)') ' Ok.'

contains

    !> Reads one whitespace/tab-delimited token before '#' from the C FILE*
    !> handle into a Fortran string, wrapping the library read_line_2get_string
    !> (whose C interface takes a char** destination).
    subroutine read_string(fh, dest)
        type(c_ptr), value :: fh
        character(len=*), intent(out) :: dest
        character(kind=c_char), target :: strbuf(1024)
        type(c_ptr), target :: strptr
        integer :: j

        strbuf = c_null_char
        strptr = c_loc(strbuf)
        call read_line_2get_string_(fh, c_loc(strptr))

        dest = ''
        do j = 1, 1024
            if (strbuf(j) == c_null_char) exit
            dest(j:j) = strbuf(j)
        end do
    end subroutine read_string

    !> Parses a "(m, n)" line into two integers, matching the oracle's
    !> fscanf "(%d, %d)".
    subroutine parse_mode(line, mval, nval)
        character(len=*), intent(in) :: line
        integer, intent(out) :: mval, nval
        character(len=len(line)) :: buf
        integer :: j
        buf = line
        do j = 1, len(buf)
            if (buf(j:j) == '(' .or. buf(j:j) == ')' .or. buf(j:j) == ',') buf(j:j) = ' '
        end do
        read (buf, *) mval, nval
    end subroutine parse_mode

    function itoa(v) result(s)
        integer, intent(in) :: v
        character(len=:), allocatable :: s
        character(len=32) :: b
        write (b, '(i0)') v
        s = trim(b)
    end function itoa

    subroutine transform_cyl_radius_to_surface_label(dim, r, pol_deg, label_flag, label)
        integer(c_int), intent(in) :: dim, pol_deg, label_flag
        real(c_double), intent(inout) :: r(0:dim - 1)
        real(c_double), intent(inout) :: label(0:dim - 1)
        real(dp), allocatable :: x_eq(:), label_eq(:)
        integer(c_int) :: dim_eq, ind
        integer :: k, u, io

        if (label_flag == 0) then
            do k = 0, dim - 1
                label(k) = r(k)
            end do
            return
        end if

        call load_equil_label(label_flag, x_eq, label_eq, dim_eq)

        ind = 0
        do k = 0, dim - 1
            call find_index_for_interp(pol_deg, r(k), dim_eq, x_eq, ind)
            call eval_neville_polynom(x_eq(ind), label_eq(ind), 1_c_int, pol_deg, r(k), label(k))
        end do

        open (newunit=u, file='label_grid.dat', status='replace', action='write', iostat=io)
        if (io /= 0) write (*, '(/,a,a,a)') 'Failed to open file ', 'label_grid.dat', achar(7)
        do k = 0, dim - 1
            write (u, '(es24.16e3,a,es24.16e3)') r(k), char(9), label(k)
        end do
        close (u)
    end subroutine transform_cyl_radius_to_surface_label

    subroutine transform_surface_label_to_cyl_radius(dim, r, pol_deg, label_flag, label)
        integer(c_int), intent(in) :: dim, pol_deg, label_flag
        real(c_double), intent(inout) :: r(0:dim - 1)
        real(c_double), intent(inout) :: label(0:dim - 1)
        real(dp), allocatable :: x_eq(:), label_eq(:)
        integer(c_int) :: dim_eq, ind
        integer :: k, u, io

        if (label_flag == 0) then
            do k = 0, dim - 1
                r(k) = label(k)
            end do
            return
        end if

        call load_equil_label(label_flag, x_eq, label_eq, dim_eq)

        ind = 0
        do k = 0, dim - 1
            call find_index_for_interp(pol_deg, label(k), dim_eq, label_eq, ind)
            call eval_neville_polynom(label_eq(ind), x_eq(ind), 1_c_int, pol_deg, label(k), r(k))
        end do

        open (newunit=u, file='radial_grid.dat', status='replace', action='write', iostat=io)
        if (io /= 0) write (*, '(/,a,a,a)') 'Failed to open file ', 'radial_grid.dat', achar(7)
        do k = 0, dim - 1
            write (u, '(es24.16e3,a,es24.16e3)') r(k), char(9), label(k)
        end do
        close (u)
    end subroutine transform_surface_label_to_cyl_radius

    !> Shared by both transforms: reads equil_r_q_psi.dat (5 columns,
    !> comment-aware) and builds the surface label array selected by label_flag.
    subroutine load_equil_label(label_flag, x_eq, label_eq, dim_eq)
        integer(c_int), intent(in) :: label_flag
        real(dp), allocatable, intent(out) :: x_eq(:), label_eq(:)
        integer(c_int), intent(out) :: dim_eq
        integer, parameter :: ncols = 5
        real(dp), allocatable :: prof_eq(:)
        real(dp) :: norm_coeff
        integer :: k, io

        dim_eq = count_lines_in_file_with_comments_(to_cstr('equil_r_q_psi.dat'), 1_c_int)
        allocate (x_eq(0:dim_eq - 1), prof_eq(0:ncols*dim_eq - 1), label_eq(0:dim_eq - 1))
        io = load_data_file_with_comments_(to_cstr('equil_r_q_psi.dat'), dim_eq, ncols, &
                                           x_eq, prof_eq)

        select case (label_flag)
        case (1)
            do k = 0, dim_eq - 1
                label_eq(k) = sqrt(abs(prof_eq(k + dim_eq)))
            end do
        case (2)
            do k = 0, dim_eq - 1
                label_eq(k) = prof_eq(k + dim_eq)
            end do
        case (3)
            norm_coeff = prof_eq(dim_eq - 1 + dim_eq)
            do k = 0, dim_eq - 1
                label_eq(k) = prof_eq(k + dim_eq)/norm_coeff
            end do
        case (4)
            do k = 0, dim_eq - 1
                label_eq(k) = sqrt(abs(prof_eq(k + 2*dim_eq)))
            end do
        case (5)
            do k = 0, dim_eq - 1
                label_eq(k) = prof_eq(k + 2*dim_eq)
            end do
        case (6)
            norm_coeff = prof_eq(dim_eq - 1 + 2*dim_eq)
            do k = 0, dim_eq - 1
                label_eq(k) = prof_eq(k + 2*dim_eq)/norm_coeff
            end do
        case default
            write (*, '(/,a,i0)') 'error: unknown grid flag = ', label_flag
            stop 1
        end select
    end subroutine load_equil_label

    !> C static-local quirk preserved: ind1/ind2 seed to dim/2 on the FIRST
    !> call only and then carry over between modes (the oracle's `static int
    !> ind = dim/2;`), so later modes start their index search from the
    !> previous mode's converged bracket.
    subroutine evaluate_Br_formfactors(dim1, r_EB1, prof_EB1, dim2, r_EB2, prof_EB2, &
        pol_deg, r, form_facs)
        integer(c_int), intent(in) :: dim1, dim2, pol_deg
        real(c_double), intent(in) :: r_EB1(0:*), prof_EB1(0:*), r_EB2(0:*), prof_EB2(0:*)
        real(c_double), intent(in) :: r
        complex(dp), intent(out) :: form_facs
        integer(c_int), save :: ind1, ind2
        logical, save :: ind1_init = .true., ind2_init = .true.
        integer, parameter :: re_Br_ind = 6, im_Br_ind = 7
        real(c_double) :: re_Q, im_Q
        complex(dp) :: Q(0:1)

        if (ind1_init) then
            ind1 = dim1/2
            ind1_init = .false.
        end if
        call find_index_for_interp(pol_deg, r, dim1, r_EB1, ind1)
        call eval_neville_polynom(r_EB1(ind1), prof_EB1(dim1*re_Br_ind + ind1), &
                                  1_c_int, pol_deg, r, re_Q)
        call eval_neville_polynom(r_EB1(ind1), prof_EB1(dim1*im_Br_ind + ind1), &
                                  1_c_int, pol_deg, r, im_Q)
        Q(0) = cmplx(re_Q, im_Q, dp)

        if (ind2_init) then
            ind2 = dim2/2
            ind2_init = .false.
        end if
        call find_index_for_interp(pol_deg, r, dim2, r_EB2, ind2)
        call eval_neville_polynom(r_EB2(ind2), prof_EB2(dim2*re_Br_ind + ind2), &
                                  1_c_int, pol_deg, r, re_Q)
        call eval_neville_polynom(r_EB2(ind2), prof_EB2(dim2*im_Br_ind + ind2), &
                                  1_c_int, pol_deg, r, im_Q)
        Q(1) = cmplx(re_Q, im_Q, dp)

        form_facs = Q(0)/Q(1)
    end subroutine evaluate_Br_formfactors

    subroutine evaluate_jp_over_Br(dim_EB, r_EB, prof_EB, dim_jp, r_jp, prof_jp, &
        pol_deg, r, form_facs)
        integer(c_int), intent(in) :: dim_EB, dim_jp, pol_deg
        real(c_double), intent(in) :: r_EB(0:*), prof_EB(0:*), r_jp(0:*), prof_jp(0:*)
        real(c_double), intent(in) :: r
        complex(dp), intent(out) :: form_facs
        integer(c_int), save :: ind1, ind2
        logical, save :: ind1_init = .true., ind2_init = .true.
        integer, parameter :: re_Br_ind = 6, im_Br_ind = 7
        integer, parameter :: re_jp_ind = 0, im_jp_ind = 1
        real(c_double) :: re_Q, im_Q
        complex(dp) :: Q(0:1)

        if (ind1_init) then
            ind1 = dim_EB/2
            ind1_init = .false.
        end if
        call find_index_for_interp(pol_deg, r, dim_EB, r_EB, ind1)
        call eval_neville_polynom(r_EB(ind1), prof_EB(dim_EB*re_Br_ind + ind1), &
                                  1_c_int, pol_deg, r, re_Q)
        call eval_neville_polynom(r_EB(ind1), prof_EB(dim_EB*im_Br_ind + ind1), &
                                  1_c_int, pol_deg, r, im_Q)
        Q(0) = cmplx(re_Q, im_Q, dp)

        if (ind2_init) then
            ind2 = dim_jp/2
            ind2_init = .false.
        end if
        call find_index_for_interp(pol_deg, r, dim_jp, r_jp, ind2)
        call eval_neville_polynom(r_jp(ind2), prof_jp(dim_jp*re_jp_ind + ind2), &
                                  1_c_int, pol_deg, r, re_Q)
        call eval_neville_polynom(r_jp(ind2), prof_jp(dim_jp*im_jp_ind + ind2), &
                                  1_c_int, pol_deg, r, im_Q)
        Q(1) = cmplx(re_Q, im_Q, dp)

        form_facs = Q(1)/Q(0)/r
    end subroutine evaluate_jp_over_Br

end program main_post_proc
