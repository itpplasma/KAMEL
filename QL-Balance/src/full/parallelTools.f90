module parallelTools

    use mpi

    implicit none

    integer :: irank
    integer :: ierror
    integer :: np_num

    contains

    subroutine initMPI

        call MPI_INIT(ierror)
        call MPI_COMM_SIZE(MPI_COMM_WORLD, np_num, ierror)
        call MPI_COMM_RANK(MPI_COMM_WORLD, irank, ierror)


        if (irank .eq. 0) then
            write(*,*) ' '
            write(*,*) '******** MPI Init ************'
            write(*,*) 'number of processes:', np_num
            write(*,*) '              irank:', irank
            write(*,*) '******************************'
            write(*,*) ' '
        end if

    end subroutine

end module parallelTools