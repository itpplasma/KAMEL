! Read the plasma profiles into global variables 
! defined in plas_parameter mod
subroutine read_profiles(reduce)

    use plas_parameter
    use config
    use grid

    implicit none
    logical, intent(in) :: reduce ! reduce r dimension

    if (hdf5_input) then
        ! read plasma profiles from hdf5 file
        call read_from_hdf5
    else
        ! read plasma profiles from text files
        call read_from_text
    endif

    contains
    
    subroutine read_from_text

        implicit none
        integer :: i, sigma
        integer :: ierr
        integer :: ios
        character(256) :: fileloc, cwd

        if (fstatus == 1) write(*,*) 'Status: Reading profiles from text files, reduce=', reduce

        !call find_file_length(trim(profile_location)//'n.dat', iprof_length)

        ! find profile length
        iprof_length = 0

        open(99, file=trim(profile_location)//'n.dat')
        ios = 0
        do while(ios == 0)
            read(99, *, iostat=ios)
            if (ios == 0) then
                iprof_length = iprof_length + 1
            end if
        end do
        close(99)


        if (.not. allocated(r_prof)) allocate(r_prof(iprof_length), stat=ierr)
        if (ierr /= 0) print *, "array: Allocation request denied"
        
        !write(*,*) r_prof
        allocate(n_prof(iprof_length), Te_prof(iprof_length), Ti_prof(ispecies, iprof_length), &
        Er_prof(iprof_length), q_prof(iprof_length), ni_prof(ispecies, iprof_length))

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
            read(11, *) r_prof(i), Ti_prof(1, i)
        end do
        close(11)

        open(11, file=trim(profile_location)//'Er.dat')
        do i=1,iprof_length
            read(11, *) r_prof(i), Er_prof(i)
        end do
        close(11)

        open(11, file=trim(profile_location)//'q.dat')
        do i=1, iprof_length
            read(11, *) r_prof(i), q_prof(i)
        end do
        close(11)

        do i=2, ispecies
            Ti_prof(i,:) = Ti_prof(1,:)
        end do

        do i = 1, iprof_length
            do sigma = 1, ispecies
                ! ion density to fulfill quasineutrality
                ni_prof(sigma, i) = n_prof(i) * Zi(sigma) / sum(Zi)
            end do
        end do

        if (reduce) call reduce_dim

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
        character(256), intent(in) :: filename
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

    subroutine reduce_dim

        implicit none
        double precision :: step_h
        double precision, allocatable :: new_n_prof(:), new_Te_prof(:), new_Ti_prof(:,:), &
                                         new_ni_prof(:,:), new_Er_prof(:), new_q_prof(:), new_r_prof(:)
        integer :: i, sigma
        integer :: nlagr = 4
        integer :: nder = 0
        integer :: ibeg, iend, ir
        double precision, dimension(:,:), allocatable :: coef

        if (.not. allocated(coef)) allocate(coef(0:nder, nlagr))

        if (fstatus==1) write(*,*) 'Status: Reducing input profile r dimension'

        allocate(new_n_prof(reduced_r_dim), new_Te_prof(reduced_r_dim), &
                new_Ti_prof(ispecies, reduced_r_dim), new_ni_prof(ispecies, reduced_r_dim), &
                new_Er_prof(reduced_r_dim), new_q_prof(reduced_r_dim), new_r_prof(reduced_r_dim))


        step_h = (r_prof(iprof_length) - r_prof(1)) / reduced_r_dim ! new step size
        !write(*,*) 'step_h = ', step_h
        !write(*,*) 'r_prof = ', r_prof
        new_r_prof(1) = r_prof(1)

        do i = 2, reduced_r_dim
            new_r_prof(i) = new_r_prof(i-1) + step_h
        end do

        do i = 1, reduced_r_dim
            call binsrc(r_prof, 1, iprof_length, new_r_prof(i), ir) 
            ibeg = max(1, ir - nlagr/2)
            iend = ibeg + nlagr - 1
            if (iend .gt. iprof_length) then
                iend = iprof_length
                ibeg = iend -nlagr + 1
            end if

            call plag_coeff(nlagr, nder, new_r_prof(i), r_prof(ibeg:iend), coef)

            new_n_prof(i) = sum(coef(0,:) * n_prof(ibeg:iend))
            new_Te_prof(i) = sum(coef(0,:) * Te_prof(ibeg:iend))
            new_Er_prof(i) = sum(coef(0,:) * Er_prof(ibeg:iend))
            new_q_prof(i) = sum(coef(0,:) * q_prof(ibeg:iend))

            do sigma=1, ispecies 
                new_Ti_prof(sigma, i) = sum(coef(0,:) * Ti_prof(sigma, ibeg:iend))
                new_ni_prof(sigma, i) = sum(coef(0,:) * ni_prof(sigma, ibeg:iend))
            end do
        end do
            
        iprof_length = reduced_r_dim

        deallocate(r_prof, n_prof, Te_prof, Ti_prof, ni_prof, Er_prof, q_prof)
        allocate(r_prof(reduced_r_dim), n_prof(reduced_r_dim), Te_prof(reduced_r_dim), &
                 Ti_prof(ispecies, reduced_r_dim), ni_prof(ispecies, reduced_r_dim), &
                 Er_prof(reduced_r_dim), q_prof(reduced_r_dim))

        r_prof = new_r_prof
        n_prof = new_n_prof
        Te_prof = new_Te_prof
        Ti_prof = new_Ti_prof
        ni_prof = new_ni_prof
        Er_prof = new_Er_prof
        q_prof = new_q_prof
 
        deallocate(new_r_prof, new_n_prof, new_Te_prof, new_Ti_prof, new_ni_prof, new_Er_prof, new_q_prof)

        call write_profiles

    end subroutine

    subroutine write_profiles

        implicit none
        integer :: i

        if (fstatus == 1) write(*,*) 'Status: writing profiles to output_path'

        inquire(file=trim(output_path)//'profiles', exist=ex)
        if (.not. ex) then
            call system('mkdir -p '//trim(output_path)//'profiles')
        end if

        open(11, file=trim(output_path)//'profiles/n.dat')
        do i=1,iprof_length
            write(11, *) r_prof(i), n_prof(i)
        end do
        close(11)

        open(11, file=trim(output_path)//'profiles/Te.dat')
        do i=1,iprof_length
            write(11, *) r_prof(i), Te_prof(i)
        end do
        close(11)

        open(11, file=trim(output_path)//'profiles/Ti.dat')
        do i=1,iprof_length
            write(11, *) r_prof(i), Ti_prof(1, i)
        end do
        close(11)

        open(11, file=trim(output_path)//'profiles/Er.dat')
        do i=1,iprof_length
            write(11, *) r_prof(i), Er_prof(i)
        end do
        close(11)

        open(11, file=trim(output_path)//'profiles/q.dat')
        do i=1,iprof_length
            write(11, *) r_prof(i), q_prof(i)
        end do
        close(11)

        !do i=2, ispecies
        !    Ti_prof(i,:) = Ti_prof(1,:)
        !end do

    end subroutine

end subroutine


