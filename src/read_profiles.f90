! Read the plasma profiles into global variables 
! defined in plas_parameter mod
subroutine read_profiles

    use plas_parameter
    use config
    implicit none

    if (hdf5_input) then
        ! read plasma profiles from hdf5 file
        call read_from_hdf5
    else
        ! read plasma profiles from text files
        call read_from_text
    endif

end subroutine


subroutine read_from_text

    use plas_parameter
    use config
    implicit none
    integer :: i

    if (fstatus == 1) write(*,*) 'Status: Reading profiles from text files'

    call find_file_length(trim(profile_location)//'n.dat', iprof_length)
    allocate(r_prof(iprof_length), n_prof(iprof_length), Te_prof(iprof_length), Ti_prof(iprof_length), &
    Er_prof(iprof_length), q_prof(iprof_length))

    open(11, file=trim(profile_location)//'n.dat')
    do i=1,iprof_length
        read(11, *) r_prof(i), n_prof(i)
    end do
    close(11)

    open(11, file=trim(profile_location)//'Te.dat')
    do i=1,iprof_length
        read(11, *) r_prof(i), Te_prof(i)
    end do
    close(11)

    open(11, file=trim(profile_location)//'Ti.dat')
    do i=1,iprof_length
        read(11, *) r_prof(i), Ti_prof(i)
    end do
    close(11)

    open(11, file=trim(profile_location)//'Er.dat')
    do i=1,iprof_length
        read(11, *) r_prof(i), Er_prof(i)
    end do
    close(11)

    open(11, file=trim(profile_location)//'q.dat')
    do i=1,iprof_length
        read(11, *) r_prof(i), q_prof(i)
    end do
    close(11)

    if (fstatus == 1) write(*,*) 'Status: Finished reading profiles from text files'

end subroutine


subroutine read_from_hdf5

    use plas_parameter
    use config
    implicit none

    if (fstatus == 1) write(*,*) 'Status: Reading profiles from hdf5 file'


end subroutine

! find the length of a profile file
subroutine find_file_length(filename, l)

    implicit none
    character(1024), intent(in) :: filename
    integer, intent(out) :: l
    integer :: ios = 0
    l = 0

    open(11, file=trim(filename))
    do while(ios == 0)
        read(11, *, iostat=ios)
        if (ios == 0) then
            l = l + 1
        end if
    end do
    close(11)

end subroutine