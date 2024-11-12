module loading_bar

    implicit none

    contains

        subroutine updateLoadingBar(current_step, total_steps)

            implicit none 
            integer, intent(in) :: current_step, total_steps
            real(kind=8) :: percentage

            ! Calculate the percentage completion
            percentage = real(current_step) / real(total_steps) * 100.0

            ! Clear the previous loading bar
            write(*, '(A)', advance='no') ACHAR(13)!char(27)//"[2K"
            ! Display the loading bar
            write(*, '(A, F6.2, A)', advance='no') "Progress: ", percentage, "% ["
            call drawLoadingBar(percentage)
            write(*, '(A)', advance='no') "]"
        end subroutine updateLoadingBar

        subroutine drawLoadingBar(percentage)

            implicit none
            real(kind=8), intent(in) :: percentage
            integer :: num_blocks, i

            ! Calculate the number of blocks to display
            num_blocks = nint(percentage / 2.0)

            ! Display blocks for the loading bar
            do i = 1, num_blocks
                write(*, '(A)', advance='no') "#"
            end do

            ! Display empty spaces for the loading bar
            do i = num_blocks + 1, 50
                write(*, '(A)', advance='no') " "
            end do
        end subroutine drawLoadingBar

end module loading_bar