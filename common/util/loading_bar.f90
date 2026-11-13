module loading_bar_m
    use iso_c_binding, only: c_int
    use iso_fortran_env, only: dp => real64, output_unit

    implicit none

    private
    public :: updateLoadingBar, updateLoadingBarWithETA, drawLoadingBar

    integer, parameter :: default_bar_width = 50

    interface
        function get_terminal_width_c() bind(C, name="kamel_get_terminal_width") result(width)
            import c_int
            integer(c_int) :: width
        end function get_terminal_width_c
    end interface

contains

    subroutine updateLoadingBar(current_step, total_steps, label)

        implicit none

        integer, intent(in) :: current_step, total_steps
        character(len=*), intent(in), optional :: label
        real(dp) :: percentage
        character(len=32) :: display_label

        ! Set label (default: "Progress")
        if (present(label)) then
            display_label = label
        else
            display_label = "Progress"
        end if

        ! Calculate the percentage completion
        if (total_steps <= 0) then
            percentage = 0.0_dp
        else
            percentage = real(current_step, dp) / real(total_steps, dp) * 100.0_dp
        end if

        call renderLoadingBar(percentage, display_label, "]")
    end subroutine updateLoadingBar

    subroutine updateLoadingBarWithETA(current_step, total_steps, start_count, count_rate, label)

        implicit none

        integer, intent(in) :: current_step, total_steps
        integer(kind=8), intent(in) :: start_count, count_rate
        character(len=*), intent(in), optional :: label
        real(dp) :: percentage, elapsed_time, estimated_total_time, eta
        integer(kind=8) :: current_count
        integer :: eta_hours, eta_minutes, eta_seconds
        character(len=32) :: display_label
        character(len=64) :: eta_suffix

        ! Set label (default: "Progress")
        if (present(label)) then
            display_label = label
        else
            display_label = "Progress"
        end if

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

        if (current_step > 0) then
            write (eta_suffix, '(A, I0.2, A, I0.2, A, I0.2)') "] ETA: ", &
                eta_hours, ":", eta_minutes, ":", eta_seconds
        else
            eta_suffix = "] ETA: --:--:--"
        end if

        call renderLoadingBar(percentage, display_label, eta_suffix)
    end subroutine updateLoadingBarWithETA

    subroutine drawLoadingBar(percentage, bar_width)

        implicit none

        real(dp), intent(in) :: percentage
        integer, intent(in), optional :: bar_width
        integer :: width, num_blocks, i

        width = default_bar_width
        if (present(bar_width)) then
            width = max(0, bar_width)
        end if

        ! Calculate the number of blocks to display
        num_blocks = nint(percentage / 100.0_dp * real(width, dp))
        num_blocks = max(0, min(width, num_blocks))

        ! Display blocks for the loading bar
        do i = 1, num_blocks
            write (output_unit, '(A)', advance='no') "#"
        end do

        ! Display empty spaces for the loading bar
        do i = num_blocks + 1, width
            write (output_unit, '(A)', advance='no') " "
        end do
    end subroutine drawLoadingBar

    subroutine renderLoadingBar(percentage, display_label, suffix)

        implicit none

        real(dp), intent(in) :: percentage
        character(len=*), intent(in) :: display_label, suffix
        integer :: terminal_width, bar_width

        terminal_width = int(get_terminal_width_c())
        bar_width = calculateBarWidth(terminal_width, display_label, suffix)

        ! Clear the previous line when the output is a terminal.  The fallback
        ! keeps redirected output free of terminal escape sequences.
        write (output_unit, '(A)', advance='no') ACHAR(13)
        if (terminal_width > 0) then
            write (output_unit, '(A)', advance='no') ACHAR(27)//"[K"
        end if

        write (output_unit, '(A, A, F6.2, A)', advance='no') &
            trim(display_label), ": ", percentage, "% ["
        call drawLoadingBar(percentage, bar_width)
        write (output_unit, '(A)', advance='no') trim(suffix)
        flush (output_unit)
    end subroutine renderLoadingBar

    integer function calculateBarWidth(terminal_width, display_label, suffix)

        implicit none

        integer, intent(in) :: terminal_width
        character(len=*), intent(in) :: display_label, suffix
        integer :: fixed_width

        if (terminal_width <= 0) then
            calculateBarWidth = default_bar_width
            return
        end if

        fixed_width = len_trim(display_label) + len(": ") + 6 + len("% [") + len_trim(suffix)
        calculateBarWidth = max(0, terminal_width - fixed_width)
    end function calculateBarWidth

end module loading_bar_m
