module plotting

    use KIM_kinds, only: dp

    implicit none

    contains

    subroutine write_array_to_file(filename, A, nx, ny)

        implicit none

        character(*), intent(in) :: filename
        real(dp), intent(in) :: A(nx, ny)
        integer, intent(in) :: nx, ny

        integer :: i, j
        open(unit=10, file=filename, status='replace')
        do j = 1, ny
            do i = 1, nx
                write(10, '(F12.6,1x,F12.6,1x,F12.6)') real(i), real(j), A(i,j)
            end do
            write(10, *) ! blank line to separate rows
        end do
        close(10)
    end subroutine write_array_to_file


    subroutine write_matrix(filename, A, nx, ny)
        implicit none
        character(*), intent(in) :: filename
        real(dp), intent(in) :: A(nx, ny)
        integer, intent(in) :: nx, ny

        integer :: i, j
        open(unit=10, file=filename, status='replace')

        do j = 1, ny
            write(10, *) (A(i,j), i = 1, nx)
        end do

        close(10)
    end subroutine write_matrix


    subroutine write_profile(x, y, n, filename)

        use KIM_kinds, only: dp

        implicit none

        integer, intent(in) :: n
        real(dp), intent(in) :: x(n), y(n)
        character(len=*), intent(in) :: filename

        ! Local variables
        integer :: i
        integer :: unit

        ! Choose a unit number (any unused unit number)
        unit = 10

        ! Open the file for writing
        open(unit=unit, file=filename, status='replace', action='write')

        ! Write data as two columns
        do i = 1, n
            write(unit, *) x(i), y(i)
        end do

        close(unit)

    end subroutine write_profile


    subroutine write_complex_profile(x, y, n, filename)

        use KIM_kinds, only: dp

        implicit none

        integer, intent(in) :: n
        real(dp), intent(in) :: x(n)
        complex(dp), intent(in) :: y(n)
        character(len=*), intent(in) :: filename

        ! Local variables
        integer :: i
        integer :: unit

        ! Choose a unit number (any unused unit number)
        unit = 10

        ! Open the file for writing
        open(unit=unit, file=filename, status='replace', action='write')

        ! Write data as two columns
        do i = 1, n
            write(unit, *) x(i), real(y(i)), dimag(y(i))
        end do

        close(unit)

    end subroutine write_complex_profile

    subroutine plot_2D(datafile)

        implicit none

        character(*), intent(in) :: datafile
        character(len=200) :: cmd

        ! Create command string
        write(cmd, '(A)') 'gnuplot -persist -e "splot '''//trim(datafile)//''' with lines"'

        call execute_command_line(trim(cmd))

    end subroutine plot_2D

    subroutine plot_1D(datafile)

        implicit none

        character(*), intent(in) :: datafile
        character(len=256) :: cmd

        write(cmd, '(A)') 'gnuplot -persist -e "plot '''//trim(datafile)//''' using 1:2 with lines"'

        call execute_command_line(trim(cmd))

    end subroutine

    subroutine plot_profile(x,y)

        use KIM_kinds, only: dp

        implicit none

        real(dp), intent(in) :: x(:), y(:)

        call write_profile(x, y, size(x), 'profile.dat')
        call plot_1D('profile.dat')
        call remove_file('profile.dat')

    end subroutine

    subroutine plot_complex_1D(datafile)

        implicit none

        character(*), intent(in) :: datafile
        character(len=256) :: cmd

        write(cmd, '(A)') 'gnuplot -persist -e "plot '''//trim(datafile)//''' using 1:2 with lines; plot '''//trim(datafile)//''' using 1:3 with lines"'

        call execute_command_line(trim(cmd))

    end subroutine

    subroutine plot_matrix(datafile)
        implicit none
        character(*), intent(in) :: datafile
        character(len=300) :: cmd

        write(cmd, '(A)') 'gnuplot -persist -e "plot '''//trim(datafile)//''' matrix with image"'
        call execute_command_line(trim(cmd))
    end subroutine plot_matrix

    subroutine plot_1D_labeled(datafile, xlabel, ylabel, title)

        implicit none

        character(*), intent(in) :: datafile
        character(*), intent(in) :: xlabel, ylabel, title
        character(len=1024) :: cmd
        character(len=1024) :: plot_cmd

        plot_cmd = 'set xlabel ''' // trim(xlabel) // '''; ' // &
                'set ylabel ''' // trim(ylabel) // '''; ' // &
                'set title '''  // trim(title)  // '''; ' // &
                'plot ''' // trim(datafile) // ''' using 1:2 with lines;'

        write(cmd, '(A)') 'gnuplot -persist -e "' // trim(plot_cmd) // '"'

        call execute_command_line(trim(cmd))

    end subroutine

    subroutine remove_file(filename)
        implicit none

        character(*), intent(in) :: filename
        character(len=300) :: cmd

        ! Build the command to remove the file
        write(cmd, '(A)') 'rm -f "' // trim(filename) // '"'

        ! Execute the command
        call execute_command_line(trim(cmd))

    end subroutine



end module