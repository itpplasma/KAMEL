program test_directory
    use, intrinsic :: iso_c_binding, only: c_ptr, c_associated
    use, intrinsic :: iso_fortran_env, only: dp => real64
    use kilca_progs_common_m, only: open_directory, read_directory, close_directory, &
                                    clean_run, shell_quote
    implicit none
    character(len=*), parameter :: root = "test directory's files"
    character(len=*), parameter :: runpath = root//'/run/'
    character(len=1024) :: name
    type(c_ptr) :: directory
    integer :: unit, count_zones

    ! A sibling sentinel catches accidental deletion through '..' or an empty name.
    call command('mkdir -p '//shell_quote(runpath//'linear-data/first[1,2]')//' '// &
                 shell_quote(runpath//'linear-data/last[3,-4]')//' '// &
                 shell_quote(runpath//"linear-data/rejected's folder")//' '// &
                 shell_quote(runpath//'linear-data/...')//' '// &
                shell_quote(runpath//'dispersion-data')//' '//shell_quote(runpath//'poincare-data'))
    call touch(root//'/sentinel')
    call touch(runpath//'zone_1.in')
    call touch(runpath//'other.in')
    call touch(runpath//'linear-data/first[1,2]/sentinel')
    call touch(runpath//'linear-data/last[3,-4]/sentinel')
    call touch(runpath//'linear-data/.../sentinel')

    directory = open_directory(runpath)
    if (.not. c_associated(directory)) error stop 'Cannot open test directory'
    count_zones = 0
    do while (read_directory(directory, name))
        if (trim(name) == 'zone_1.in') count_zones = count_zones + 1
    end do
    call close_directory(directory)
    if (count_zones /= 1) error stop 'Zone discovery did not return the native filename'

    call clean_run(runpath, cmplx(1.0_dp, 2.0_dp, dp), cmplx(3.0_dp, -4.0_dp, dp))
    call assert_exists(root//'/sentinel', .true.)
    call assert_exists(runpath//'zone_1.in', .true.)
    call assert_exists(runpath//'linear-data/first[1,2]/sentinel', .true.)
    call assert_exists(runpath//'linear-data/last[3,-4]/sentinel', .true.)
    call assert_exists(runpath//'linear-data/.../sentinel', .true.)
    call assert_exists(runpath//"linear-data/rejected's folder", .false.)
    call assert_exists(runpath//'dispersion-data', .false.)
    call assert_exists(runpath//'poincare-data', .false.)
    call command('rm -rf -- '//shell_quote(root))
    print *, 'PASS: native directory names and production eigenmode cleanup'
contains
    subroutine touch(path)
        character(len=*), intent(in) :: path
        open (newunit=unit, file=path, status='replace')
        close (unit)
    end subroutine touch

    subroutine assert_exists(path, expected)
        character(len=*), intent(in) :: path
        logical, intent(in) :: expected
        logical :: exists
        inquire (file=path, exist=exists)
        if (exists .neqv. expected) then
            print *, 'Unexpected existence: ', path, ' expected ', expected
            error stop 1
        end if
    end subroutine assert_exists

    subroutine command(cmd)
        character(len=*), intent(in) :: cmd
        integer :: status
        call execute_command_line(cmd, exitstat=status)
        if (status /= 0) error stop 'Test setup/teardown failed'
    end subroutine command
end program test_directory
