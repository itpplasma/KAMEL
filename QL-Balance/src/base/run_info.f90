!> Provenance of a run in its output file path2out:
!>   /run_info/start_time  ISO 8601 local time with UTC offset at which the run started
!>   /run_info/git_commit  git commit of the KAMEL build, "-dirty" appended if tracked
!>                         files differed from it at build time ("unknown" without git)
!>   /run_info/balancenml  namelist BALANCENML with every value as used by the run,
!>                         defaults included, as recorded by read_config
!> git_commit is compiled in from ql_balance_git_version.inc, which
!> cmake/WriteGitVersion.cmake regenerates on every build.
module run_info_m

    implicit none

    private
    public :: set_run_config, write_run_info

    include 'ql_balance_git_version.inc'

    character(len=:), allocatable :: config_nml

contains

    !> Store the records of a namelist write as text, one line per record.
    !> Drops the blanks gfortran pads values with, outside and at the end of strings.
    subroutine set_run_config(records)

        character(len=*), intent(in) :: records(:)
        character(len=:), allocatable :: line
        character :: c
        logical :: in_string
        integer :: i, j

        config_nml = ""
        do i = 1, size(records)
            line = ""
            in_string = .false.
            do j = 1, len_trim(records(i))
                c = records(i) (j:j)
                if (c == '"') then
                    if (in_string) line = trim(line)
                    in_string = .not. in_string
                else if (c == " " .and. .not. in_string) then
                    cycle
                end if
                line = line//c
            end do
            if (len(line) > 0) config_nml = config_nml//line//new_line("a")
        end do

    end subroutine set_run_config

    function start_time() result(stamp)

        character(len=25) :: stamp
        integer :: v(8)
        character :: sign

        call date_and_time(values=v)
        sign = merge("+", "-", v(4) >= 0)
        write (stamp, '(i4.4,"-",i2.2,"-",i2.2,"T",i2.2,":",i2.2,":",i2.2,a,i2.2,":",i2.2)') &
            v(1), v(2), v(3), v(5), v(6), v(7), sign, abs(v(4)) / 60, mod(abs(v(4)), 60)

    end function start_time

    !> Write /run_info to path2out. Called once after read_config, before the run
    !> writes its results.
    subroutine write_run_info

        use h5mod, only: h5_id, path2out, h5_init, h5_open_rw, h5_close, h5_deinit, h5_add, &
                         h5overwrite, create_group_if_not_existent
        use logger_m, only: log_info

        logical :: old_h5overwrite

        ! A rerun into an existing output file replaces the previous values.
        old_h5overwrite = h5overwrite
        h5overwrite = .true.
        call h5_init()
        call h5_open_rw(trim(path2out), h5_id)
        call create_group_if_not_existent("/run_info")
        call h5_add(h5_id, "/run_info/start_time", start_time())
        call h5_add(h5_id, "/run_info/git_commit", git_commit)
        call h5_add(h5_id, "/run_info/balancenml", config_nml)
        call h5_close(h5_id)
        call h5_deinit()
        h5overwrite = old_h5overwrite

        call log_info("    KAMEL commit "//git_commit)

    end subroutine write_run_info

end module run_info_m
