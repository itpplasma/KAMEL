! Minimal XERROR compatibility for QUADPACK when SLATEC error lib is absent
subroutine xerror(msg, nmsg, nerr, level)
    implicit none
    character(*), intent(in) :: msg
    integer, intent(in) :: nmsg, nerr, level
    integer, parameter :: MAX_CODE = 100
    integer, save :: counts(0:MAX_CODE) = 0
    integer :: idx
    character(:), allocatable :: m

    idx = nerr
    if (idx < 0) idx = 0
    if (idx > MAX_CODE) idx = MAX_CODE
    m = trim(msg(1:min(len(msg), nmsg)))

    ! Throttle: print first 3 occurrences per error code, then one suppression notice
    if (counts(idx) < 3) then
        !$omp critical(xerror_print)
        write(*,*) 'XERROR (nerr=', nerr, ', level=', level, '): ', m
        !$omp end critical(xerror_print)
    else if (counts(idx) == 3) then
        !$omp critical(xerror_print)
        write(*,*) 'XERROR (nerr=', nerr, '): further messages suppressed'
        !$omp end critical(xerror_print)
    end if
    counts(idx) = counts(idx) + 1
end subroutine xerror
