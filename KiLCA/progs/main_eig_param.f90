!> Additional driver program for parameter studies of eigenmodes, ported from
!> progs/main_eig_param.cpp. Walks a scan parameter from p1 to p2: at each step
!> it stages a per-value run directory (copying the *.in files, patching the
!> V_scale line of background.in and the starting-point line of eigmode.in),
!> runs the (now-Fortran) KiLCA core-data pipeline, then validates the root
!> written to roots.dat by ratio-of-determinant and initial-guess-accuracy
!> checks. A converged step is kept (its run directory pruned to the two
!> bracketing eigenfunction folders) and the root saved; a failed step halves
!> the parameter increment and retries with a linearly extrapolated guess.
!>
!> The scan builds directory names, glob patterns and file lines with C printf
!> %lg / %.15lg / %.20le, reproduced byte for byte via fmt_g / fmt_e. OS actions
!> the oracle performed with system() are issued through execute_command_line
!> (same /bin/sh -c semantics); the run-directory scan in clean_run keeps the
!> opendir/readdir walk but replaces fnmatch with the equivalent literal
!> substring / exact-name tests its known patterns reduce to.
program main_eig_param
    use, intrinsic :: iso_c_binding, only: c_int, c_intptr_t, c_double, c_char, &
        c_ptr, c_null_char, c_null_ptr, c_loc, c_associated, c_f_pointer
    use, intrinsic :: iso_fortran_env, only: dp => real64
    use kilca_inout_m, only: read_line_2skip_it_, read_line_2get_string_, &
        read_line_2get_double_, read_line_2get_complex_
    use kilca_shared_m, only: signum
    use kilca_core_data_m, only: core_data_create_, core_data_destroy_, &
        core_data_calc_and_set_mode_independent_, &
        core_data_calc_and_set_mode_dependent_antenna_, &
        core_data_calc_and_set_mode_dependent_eigmode_
    use kilca_progs_common_m, only: get_project_path, to_cstr, fmt_g, fmt_e
    implicit none

    interface
        subroutine set_core_data_in_core_module(cd) &
            bind(C, name="set_core_data_in_core_module_")
            import :: c_intptr_t
            integer(c_intptr_t), intent(in) :: cd
        end subroutine set_core_data_in_core_module

        integer(c_int) function get_antenna_flag_eigmode() &
            bind(C, name="get_antenna_flag_eigmode_")
            import :: c_int
        end function get_antenna_flag_eigmode

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

        function c_opendir(name) result(dp) bind(C, name="opendir")
            import :: c_char, c_ptr
            character(kind=c_char), intent(in) :: name(*)
            type(c_ptr) :: dp
        end function c_opendir

        function c_readdir(dp) result(ep) bind(C, name="readdir")
            import :: c_ptr
            type(c_ptr), value :: dp
            type(c_ptr) :: ep
        end function c_readdir

        function c_closedir(dp) result(r) bind(C, name="closedir")
            import :: c_ptr, c_int
            type(c_ptr), value :: dp
            integer(c_int) :: r
        end function c_closedir
    end interface

    integer, parameter :: max_iter_num = 10000

    character(len=:), allocatable :: path
    character(len=1024) :: config, param_name, func_name, filename, fullpath
    character(len=1024) :: file_set
    type(c_ptr) :: in
    real(c_double) :: p1, p2, dp_est, dp_min, dp_max, eps_det, eps_root
    complex(c_double) :: f_start
    real(dp) :: p_curr, p_prev, dpv
    complex(dp) :: f_prev, f_curr, f_est
    real(dp), allocatable :: param(:)
    complex(dp), allocatable :: funct(:)
    integer :: step, ind, status, ios
    character(len=1024) :: exe_name

    call get_command_argument(0, exe_name)
    open (newunit=ios, file='kilca_version', status='replace', action='write')
    write (ios, '(a,a)') 'KiLCA version used for the last run: ', trim(exe_name)
    close (ios)

    path = get_project_path()

    file_set = trim(path)//'/search.in'

    in = c_fopen(to_cstr(file_set), to_cstr('r'))
    if (.not. c_associated(in)) then
        write (*, '(/,a,a,a)') &
            'error: read_settings: failed to open file ', trim(file_set), achar(7)
        stop
    end if

    call read_line_2skip_it_(in, c_null_ptr)
    call read_string(in, config)
    call read_string(in, param_name)
    call read_string(in, func_name)
    call read_line_2get_double_(in, p1)
    call read_line_2get_double_(in, p2)
    call read_line_2get_double_(in, dp_min)
    call read_line_2get_double_(in, dp_est)
    call read_line_2get_double_(in, dp_max)
    call read_line_2get_double_(in, eps_det)
    call read_line_2get_double_(in, eps_root)
    call read_line_2get_complex_(in, f_start)
    call read_line_2skip_it_(in, c_null_ptr)
    ios = c_fclose(in)

    filename = trim(path)//trim(func_name)//'_'//trim(param_name)//'_'//trim(config)//'.dat'

    allocate (param(0:max_iter_num - 1), funct(0:max_iter_num - 1))

    p_prev = p1
    p_curr = p1
    dpv = dp_est

    f_prev = f_start
    f_curr = f_start
    f_est = f_start

    step = 0
    ind = -1

    do while (abs(p_curr) < abs(p2))
        step = step + 1

        if (ind == max_iter_num) then
            write (*, '(/,a)') 'error: maximum number of steps is reached!'
            stop 1
        end if

        fullpath = trim(path)//trim(param_name)//'_'//fmt_g(p_curr, 6)//'/'

        call make_run(fullpath, param_name, p_curr, f_start)

        status = get_root_value(fullpath, eps_det, eps_root, f_est, f_curr)

        if (status > 0) then
            call remove_run(fullpath)

            if (step == 1) then
                write (*, '(/,a,a,a)') &
                    'error: root approval failed for the first value: ', fmt_e(p_curr, 6), &
                    ', try another starting guess.'
                if (status /= 1) then
                    write (*, '(/,a,a,a)') 'root is found: ', cfmt(f_curr), &
                        ', use it as initial guess.'
                end if
                stop 1
            end if

            dpv = dpv/2.0_dp
            dpv = real(signum(dpv), dp)*min(abs(dpv), abs(dp_max))

            if (abs(dpv) < abs(dp_min)) then
                write (*, '(/,a,a)') 'warning: too small dp = ', fmt_g(dpv, 6)
            end if

            p_curr = p_prev + dpv

            if (ind > 1) then
                f_est = funct(ind) &
                        + ((funct(ind) - funct(ind - 1))/(param(ind) - param(ind - 1))) &
                        *(p_curr - param(ind))
            else
                f_est = f_prev
            end if

            write (*, '(/,a,a,a,a)') 'f_est = ', fmt_g(real(f_est, dp), 6), ' ', &
                fmt_g(aimag(f_est), 6)//'i'
            f_start = f_est
        else
            ind = ind + 1

            write (*, '(/,a)') 'message: succeeded to find a root.'
            write (*, '(/,a,i0,a,a,a,a,a,a)') 'ind = ', ind, ': at ', trim(param_name), &
                ' = ', fmt_g(p_curr, 6), ' found a root = ', cfmt(f_curr)

            call clean_run(fullpath, f_start, f_curr)

            param(ind) = p_curr
            funct(ind) = f_curr

            call save_result(filename, p_curr, f_curr)

            p_prev = p_curr
            f_prev = f_curr

            dpv = dpv*2.0_dp
            dpv = real(signum(dpv), dp)*min(abs(dpv), abs(dp_max))

            p_curr = p_prev + dpv

            if (ind > 1) then
                f_est = funct(ind) &
                        + ((funct(ind) - funct(ind - 1))/(param(ind) - param(ind - 1))) &
                        *(p_curr - param(ind))
            else
                f_est = f_prev
            end if

            write (*, '(/,a,a)') 'f_est = ', cfmt(f_est)
            f_start = f_est
        end if
    end do

contains

    !> Reads one token before '#' from the C FILE* handle, wrapping the library
    !> read_line_2get_string (whose C interface writes into a char**).
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

    !> C "%lg%+lgi": real part, then imaginary part with a forced sign and an
    !> 'i' suffix, as printed for the complex frequency diagnostics.
    function cfmt(z) result(s)
        complex(dp), intent(in) :: z
        character(len=:), allocatable :: s, im
        im = fmt_g(aimag(z), 6)
        if (im(1:1) /= '-') im = '+'//im
        s = fmt_g(real(z, dp), 6)//im//'i'
    end function cfmt

    subroutine make_run(fullpath, param_name, p_curr, f_start)
        character(len=*), intent(in) :: fullpath, param_name
        real(dp), intent(in) :: p_curr
        complex(dp), intent(in) :: f_start

        write (*, '(/,/,a,a,a,a,a,a)', advance='no') 'The next run is: ', trim(param_name), &
            ' = ', fmt_g(p_curr, 6), char(9)//'start = ', cfmt(f_start)

        call run_system('mkdir -p '//trim(fullpath))
        call run_system('cp -L -p *.in '//trim(fullpath))

        call edit_input_files_Vscale(fullpath, p_curr, f_start)

        call run_kilca(fullpath)
    end subroutine make_run

    !> Mirrors the oracle's system() error path: system() == -1 means the shell
    !> could not be started, printed as "error: system()!" and a hard exit.
    subroutine run_system(cmd)
        character(len=*), intent(in) :: cmd
        integer :: cstat
        character(len=256) :: cmsg

        call execute_command_line(cmd, wait=.true., cmdstat=cstat, cmdmsg=cmsg)
        if (cstat /= 0) then
            write (*, '(/,a)') 'error: system()!'
            stop 1
        end if
    end subroutine run_system

    subroutine run_kilca(runpath)
        character(len=*), intent(in) :: runpath
        integer(c_intptr_t) :: cd
        character(kind=c_char), allocatable :: cpath(:)

        cpath = to_cstr(runpath)
        cd = core_data_create_(cpath)
        call set_core_data_in_core_module(cd)

        call core_data_calc_and_set_mode_independent_(cd)

        if (get_antenna_flag_eigmode() == 0) then
            call core_data_calc_and_set_mode_dependent_antenna_(cd)
        else
            call core_data_calc_and_set_mode_dependent_eigmode_(cd)
        end if

        call core_data_destroy_(cd)
    end subroutine run_kilca

    subroutine edit_input_files_Vscale(fullpath, p_curr, f_prev)
        character(len=*), intent(in) :: fullpath
        real(dp), intent(in) :: p_curr
        complex(dp), intent(in) :: f_prev
        call edit_background_file_Vscale(fullpath, p_curr, f_prev)
        call edit_eigmode_file_Vscale(fullpath, p_curr, f_prev)
    end subroutine edit_input_files_Vscale

    !> Rewrites background.in in place, replacing line 12 with the V_scale line
    !> carrying the current parameter value.
    subroutine edit_background_file_Vscale(fullpath, p_curr, f_prev)
        character(len=*), intent(in) :: fullpath
        real(dp), intent(in) :: p_curr
        complex(dp), intent(in) :: f_prev
        character(len=1024) :: fin, fout
        character(len=4096) :: line
        integer :: uin, uout, io, count

        fin = trim(fullpath)//'background.in'
        fout = trim(fullpath)//'background.out'

        open (newunit=uin, file=trim(fin), status='old', action='read', iostat=io)
        if (io /= 0) write (*, '(/,a,a,a)') 'Failed to open file ', trim(fin), achar(7)
        open (newunit=uout, file=trim(fout), status='replace', action='write', iostat=io)
        if (io /= 0) write (*, '(/,a,a,a)') 'Failed to open file ', trim(fout), achar(7)

        count = 0
        do
            read (uin, '(a)', iostat=io) line
            if (io /= 0) exit
            count = count + 1
            if (count == 12) then
                write (uout, '(a,a)') fmt_g(p_curr, 6), &
                    char(9)//char(9)//'#V_scale: scale factor for the Vz velocity '// &
                    'profile: Vz = V_scale*Vz - V_gal_sys'
            else
                write (uout, '(a)') trim(line)
            end if
        end do

        close (uin)
        close (uout)

        call run_system('mv -f '//trim(fout)//' '//trim(fin))
    end subroutine edit_background_file_Vscale

    !> Rewrites eigmode.in in place, replacing line 34 with the starting-point
    !> line carrying the previous root as a complex seed.
    subroutine edit_eigmode_file_Vscale(fullpath, p_curr, f_prev)
        character(len=*), intent(in) :: fullpath
        real(dp), intent(in) :: p_curr
        complex(dp), intent(in) :: f_prev
        character(len=1024) :: fin, fout
        character(len=4096) :: line
        integer :: uin, uout, io, count

        fin = trim(fullpath)//'eigmode.in'
        fout = trim(fullpath)//'eigmode.out'

        open (newunit=uin, file=trim(fin), status='old', action='read', iostat=io)
        if (io /= 0) write (*, '(/,a,a,a)') 'Failed to open file ', trim(fin), achar(7)
        open (newunit=uout, file=trim(fout), status='replace', action='write', iostat=io)
        if (io /= 0) write (*, '(/,a,a,a)') 'Failed to open file ', trim(fout), achar(7)

        count = 0
        do
            read (uin, '(a)', iostat=io) line
            if (io /= 0) exit
            count = count + 1
            if (count == 34) then
                write (uout, '(a,a,a,a)') '(', fmt_g(real(f_prev, dp), 15), ', ', &
                    fmt_g(aimag(f_prev), 15)//')'//char(9)//char(9)//'#starting point'
            else
                write (uout, '(a)') trim(line)
            end if
        end do

        close (uin)
        close (uout)

        call run_system('mv -f '//trim(fout)//' '//trim(fin))
    end subroutine edit_eigmode_file_Vscale

    !> Reads roots.dat and decides whether the last converged frequency is an
    !> acceptable root: the data line preceding the "status" line gives the root
    !> and its determinant, line 2 gives the reference determinant. Returns 0 on
    !> acceptance, 1 if the determinant ratio is too large, 2 if the guess is too
    !> far. f_curr is set to the last parsed frequency.
    function get_root_value(fullpath, eps_det, eps_root, f_est, f_curr) result(status)
        character(len=*), intent(in) :: fullpath
        real(dp), intent(in) :: eps_det, eps_root
        complex(dp), intent(in) :: f_est
        complex(dp), intent(out) :: f_curr
        integer :: status
        character(len=1024) :: fname
        character(len=4096) :: str1, str2, buff
        integer :: unit, io, lineind, iter
        real(dp) :: freq_re, freq_im, det_re, det_im
        real(dp) :: freq_re_g, freq_im_g, det_re_g, det_im_g
        real(dp) :: ratio, err_re, err_im, err

        fname = trim(fullpath)//'roots.dat'

        freq_re = 0.0_dp
        freq_im = 0.0_dp
        freq_re_g = 0.0_dp
        freq_im_g = 0.0_dp
        det_re_g = 0.0_dp
        det_im_g = 0.0_dp

        status = 1

        open (newunit=unit, file=trim(fname), status='old', action='read', iostat=io)
        if (io /= 0) then
            write (*, '(/,a,a,a)') 'Failed to open file ', trim(fname), achar(7)
            f_curr = cmplx(freq_re, freq_im, dp)
            return
        end if

        str1 = ''
        str2 = ''
        lineind = 0
        do
            read (unit, '(a)', iostat=io) buff
            if (io /= 0) exit
            lineind = lineind + 1

            str1 = str2
            str2 = buff

            if (lineind == 2) then
                read (str2, *, iostat=io) iter, freq_re_g, freq_im_g, det_re_g, det_im_g
            end if

            if (index(buff, 'status') == 0) cycle

            read (str1, *, iostat=io) iter, freq_re, freq_im, det_re, det_im

            if (index(buff, 'success') /= 0) then
                write (*, '(/,a)') 'solver returned success.'
            else
                write (*, '(/,a)') 'warning: solver failed to locate a root.'
            end if

            ratio = abs(cmplx(det_re, det_im, dp))/abs(cmplx(det_re_g, det_im_g, dp))
            write (*, '(/,a,a)') 'check: ratio = ', fmt_g(ratio, 6)

            err_re = 0.0_dp
            err_im = abs(freq_im - aimag(f_est))
            if (abs(freq_im) > 1.0_dp) err_im = err_im/abs(freq_im)
            err = max(err_re, err_im)
            write (*, '(/,a,a)') 'check: accuracy of the initial guess = ', fmt_g(err, 6)

            if (ratio < eps_det .and. err < eps_root) then
                write (*, '(/,a)') 'The ratio of determinants and accuracy of the guess: PASSED.'
                status = 0
            else
                write (*, '(/,a)') 'The ratio of determinants and accuracy of the guess: FAILED.'
                if (ratio > eps_det) then
                    status = 1
                else if (ratio < eps_det) then
                    if (err > eps_root) status = 2
                end if
            end if

            if (status == 0) then
                write (*, '(/,a,a,a,a)') 'root: freq = ', cfmt_ei(freq_re, freq_im), &
                    char(9)//'det = ', cfmt_ei(det_re, det_im)
            else
                write (*, '(/,a)') 'failed to find a proper root.'
            end if
        end do

        close (unit)

        f_curr = cmplx(freq_re, freq_im, dp)
    end function get_root_value

    !> C "%le%+lei": real then imaginary (forced sign) with an 'i' suffix.
    function cfmt_ei(re, im) result(s)
        real(dp), intent(in) :: re, im
        character(len=:), allocatable :: s, ims
        ims = fmt_e(im, 6)
        if (ims(1:1) /= '-') ims = '+'//ims
        s = fmt_e(re, 6)//ims//'i'
    end function cfmt_ei

    subroutine remove_run(fullpath)
        character(len=*), intent(in) :: fullpath
        call run_system('rm -R -f '//trim(fullpath))
    end subroutine remove_run

    !> Prunes fullpath/linear-data to the two eigenfunction folders bracketing
    !> the accepted root, then drops the dispersion- and poincare-data trees.
    !> The oracle's five fnmatch patterns reduce to: keep any name containing
    !> "[re,im]" for the first or last frequency (%.15lg formatted), and keep
    !> ".", "..", "..." exactly; everything else is removed.
    subroutine clean_run(fullpath, funct_first, funct_last)
        character(len=*), intent(in) :: fullpath
        complex(dp), intent(in) :: funct_first, funct_last
        character(len=1024) :: dir_path
        character(len=:), allocatable :: sub0, sub1, dname
        type(c_ptr) :: dp_h, ep
        integer(c_intptr_t) :: addr
        character(kind=c_char), pointer :: cname(:)
        type(c_ptr) :: namep
        integer :: j, r
        logical :: keep

        dir_path = trim(fullpath)//'linear-data/'

        sub0 = '['//fmt_g(real(funct_first, dp), 15)//','//fmt_g(aimag(funct_first), 15)//']'
        sub1 = '['//fmt_g(real(funct_last, dp), 15)//','//fmt_g(aimag(funct_last), 15)//']'

        dp_h = c_opendir(to_cstr(dir_path))
        if (c_associated(dp_h)) then
            do
                ep = c_readdir(dp_h)
                if (.not. c_associated(ep)) exit

                addr = transfer(ep, addr) + 19_c_intptr_t
                namep = transfer(addr, namep)
                call c_f_pointer(namep, cname, [256])
                dname = ''
                do j = 1, 256
                    if (cname(j) == c_null_char) exit
                    dname = dname//cname(j)
                end do

                keep = .false.
                if (dname == '.' .or. dname == '..') keep = .true.
                if (dname == '...') keep = .true.
                if (index(dname, sub0) /= 0) keep = .true.
                if (index(dname, sub1) /= 0) keep = .true.
                if (keep) cycle

                call run_system_soft('rm -R -f '//trim(fullpath)//'linear-data/'//dname)
            end do
            r = c_closedir(dp_h)
        else
            write (*, '(/,a,a,a)') 'clean_run: faled to open the directory ', trim(dir_path), '.'
        end if

        call run_system_soft('rm -R -f '//trim(fullpath)//'dispersion-data')
        call run_system_soft('rm -R -f '//trim(fullpath)//'poincare-data')
    end subroutine clean_run

    !> system() variant matching clean_run's non-fatal error handling: a failed
    !> shell start only warns ("error: system()!"), the scan continues.
    subroutine run_system_soft(cmd)
        character(len=*), intent(in) :: cmd
        integer :: cstat
        character(len=256) :: cmsg
        call execute_command_line(cmd, wait=.true., cmdstat=cstat, cmdmsg=cmsg)
        if (cstat /= 0) write (*, '(/,a)') 'error: system()!'
    end subroutine run_system_soft

    subroutine save_result(filename, p_curr, f_curr)
        character(len=*), intent(in) :: filename
        real(dp), intent(in) :: p_curr
        complex(dp), intent(in) :: f_curr
        integer :: unit, io

        open (newunit=unit, file=trim(filename), status='unknown', position='append', &
            action='write', iostat=io)
        if (io /= 0) then
            write (*, '(/,a,a,a)') 'Failed to open file ', trim(filename), achar(7)
            return
        end if

        write (unit, '(a,a,a,a,a)') fmt_e(p_curr, 20), char(9), fmt_e(real(f_curr, dp), 20), &
            char(9), fmt_e(aimag(f_curr), 20)
        close (unit)
    end subroutine save_result

end program main_eig_param
