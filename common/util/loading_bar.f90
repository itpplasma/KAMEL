module loading_bar_m
    use iso_fortran_env, only: dp => real64

    implicit none

    private
    public :: updateLoadingBar, updateLoadingBarWithETA, drawLoadingBar

contains

    subroutine updateLoadingBar(current_step, total_steps)

        implicit none

        integer, intent(in) :: current_step, total_steps
        real(dp) :: percentage

        ! Calculate the percentage completion
        percentage = real(current_step, dp) / real(total_steps, dp) * 100.0_dp

        ! Clear the previous loading bar
        write (*, '(A)', advance='no') ACHAR(13)
        ! Display the loading bar
        write (*, '(A, F6.2, A)', advance='no') "Progress: ", percentage, "% ["
        call drawLoadingBar(percentage)
        write (*, '(A)', advance='no') "]"
    end subroutine updateLoadingBar

    subroutine updateLoadingBarWithETA(current_step, total_steps, start_count, count_rate)

        implicit none

        integer, intent(in) :: current_step, total_steps
        integer(kind=8), intent(in) :: start_count, count_rate
        real(dp) :: percentage, elapsed_time, estimated_total_time, eta
        integer(kind=8) :: current_count
        integer :: eta_hours, eta_minutes, eta_seconds

        ! Get current wall clock time
        call system_clock(current_count)

        ! Calculate the percentage completion
        percentage = real(current_step, dp) / real(total_steps, dp) * 100.0_dp

        ! Calculate ETA using wall time
        elapsed_time = real(current_count - start_count, dp) / real(count_rate, dp)
        if (current_step > 0) then
            estimated_total_time = elapsed_time * real(total_steps, dp) / real(current_step, dp)
            eta = estimated_total_time - elapsed_time

            ! Convert ETA to hours, minutes, seconds
            eta_hours = int(eta / 3600.0_dp)
            eta_minutes = int(mod(eta, 3600.0_dp) / 60.0_dp)
            eta_seconds = int(mod(eta, 60.0_dp))
        else
            eta_hours = 0
            eta_minutes = 0
            eta_seconds = 0
        end if

        ! Clear the previous loading bar
        write (*, '(A)', advance='no') ACHAR(13)

        ! Display the loading bar with ETA
        if (current_step > 0) then
            write (*, '(A, F6.2, A)', advance='no') "Progress: ", percentage, "% ["
            call drawLoadingBar(percentage)
            write (*, '(A, I0.2, A, I0.2, A, I0.2, A)', advance='no') "] ETA: ", &
                eta_hours, ":", eta_minutes, ":", eta_seconds, " "
        else
            write (*, '(A, F6.2, A)', advance='no') "Progress: ", percentage, "% ["
            call drawLoadingBar(percentage)
            write (*, '(A)', advance='no') "] ETA: --:--:-- "
        end if
    end subroutine updateLoadingBarWithETA

    subroutine drawLoadingBar(percentage)

        implicit none

        real(dp), intent(in) :: percentage
        integer :: num_blocks, i

        ! Calculate the number of blocks to display
        num_blocks = nint(percentage / 2.0_dp)

        ! Display blocks for the loading bar
        do i = 1, num_blocks
            write (*, '(A)', advance='no') "#"
        end do

        ! Display empty spaces for the loading bar
        do i = num_blocks + 1, 50
            write (*, '(A)', advance='no') " "
        end do
    end subroutine drawLoadingBar

end module loading_bar_m
