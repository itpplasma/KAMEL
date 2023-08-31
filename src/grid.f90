subroutine gengrid

    use plas_parameter
    use back_quants
    implicit none


end subroutine

subroutine generate_k_space_grid(write_out)
    
    use grid
    use setup
    use config, only: output_path, fstatus

    implicit none

    integer :: i
    logical, intent(in) :: write_out

    if (fstatus == 1) write(*,*) 'Status: Generating k-space grid, write out=', write_out

    allocate(kr(k_space_dim), krp(k_space_dim))

    do i=1, k_space_dim
        kr(i) = i
        krp(i) = i
    end do

    ! center around zero
    kr = kr - k_space_dim / 2d0
    krp = krp - k_space_dim / 2d0

    if (write_out) call write_k_space

    contains

        subroutine write_k_space

            implicit none

            open(unit = 78, file = trim(output_path)//'backs/kr.dat')
            open(unit = 79, file = trim(output_path)//'backs/krp.dat')
            do i = 1, k_space_dim
                write(78,*) kr(i)
                write(79,*) krp(i)
            end do
            close(78)
            close(79)

        end subroutine

end subroutine