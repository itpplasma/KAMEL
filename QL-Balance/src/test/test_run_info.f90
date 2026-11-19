!> namelist_text turns the records of a namelist write into the text stored as
!> /run_info/balancenml. It must keep blanks inside strings, drop the padding,
!> and stay valid namelist input when gfortran wraps an array over several records.
program test_run_info
    use run_info_m, only: namelist_text
    implicit none

    character(len=1024) :: path
    integer :: modes(20)
    logical :: flag
    namelist /TESTNML/ path, modes, flag

    character(len=2048) :: records(64)
    character(len=2048), allocatable :: lines(:)
    character(len=:), allocatable :: text
    integer :: i, ios

    path = "/run dir/with spaces/out.h5"
    modes = [(7*i - 30, i = 1, size(modes))]
    flag = .true.

    records = ""
    write (records, nml=TESTNML, delim="quote", iostat=ios)
    call check(ios == 0, "namelist write fits the records")
    text = namelist_text(records)
    print '(a)', text

    call check(index(text, 'PATH="/run dir/with spaces/out.h5",') > 0, &
               "quoted path keeps its blanks and loses the padding")
    call check(index(text, "FLAG=T,") > 0, "logical value")
    call check(index(text, "MODES=-23,-16,") > 0, "array without padding")
    call split_lines(text, lines)
    ! &TESTNML, three keys and / on five lines unless the array wraps.
    call check(size(lines) > 5, "array wraps over several lines")
    call check(all([(index(trim(lines(i)), "  ") == 0, i = 1, size(lines))]), &
               "no padding left")

    ! The text must read back to the same values.
    path = ""
    modes = 0
    flag = .false.
    read (lines, nml=TESTNML, iostat=ios)
    call check(ios == 0, "text reads back as namelist")
    call check(path == "/run dir/with spaces/out.h5", "path read back")
    call check(all(modes == [(7*i - 30, i = 1, size(modes))]), "array read back")
    call check(flag, "logical read back")

contains

    subroutine split_lines(text, lines)
        character(len=*), intent(in) :: text
        character(len=2048), allocatable, intent(out) :: lines(:)
        integer :: start, newline, n

        allocate (lines(count([(text(i:i) == new_line("a"), i = 1, len(text))])))
        start = 1
        do n = 1, size(lines)
            newline = start - 1 + index(text(start:), new_line("a"))
            lines(n) = text(start:newline - 1)
            start = newline + 1
        end do
    end subroutine split_lines

    subroutine check(condition, label)
        logical, intent(in) :: condition
        character(len=*), intent(in) :: label

        if (.not. condition) then
            print '(a,a)', "FAIL: ", label
            error stop 1
        end if
    end subroutine check

end program test_run_info
