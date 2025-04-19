module plasma_parameter

    use config, only: fstatus
    use KIM_kinds, only: dp

    implicit none

    integer :: set_profiles_constant = 1

    real(dp) :: r_plas
    integer :: iprof_length
    real(dp), allocatable :: r_prof(:)
    real(dp), allocatable :: n_prof(:)
    real(dp), allocatable :: Te_prof(:) 
    real(dp), allocatable :: Er_prof(:)
    real(dp), allocatable :: q_prof(:)

    real(dp), allocatable :: ni_prof(:, :)
    real(dp), allocatable :: Ti_prof(:, :)
    integer, dimension(:), allocatable :: Zi ! ion charge number
    integer, dimension(:), allocatable :: Ai ! ion mass number
    
    real(dp), allocatable :: dndr_prof(:)
    real(dp), allocatable :: dTedr_prof(:)
    real(dp), allocatable :: dTidr_prof(:,:)
    real(dp), allocatable :: dqdr_prof(:)
    real(dp), allocatable :: dnidr_prof(:, :)

    real(dp) :: rho_L

    contains 
        ! Read the plasma profiles into global variables 
        ! defined in plasma_parameter mod
        subroutine read_profiles(reduce)

            use config, only: hdf5_input            
            use grid, only: r_space_dim

            logical, intent(in) :: reduce ! reduce r dimension

            if (hdf5_input) then
                ! read plasma profiles from hdf5 file
                call read_from_hdf5
            else
                ! read plasma profiles from text files
                call read_from_text
            endif
            
            if (reduce) call reduce_dim

            r_space_dim = size(r_prof)

        end subroutine
    
        subroutine read_from_text

            use config, only: number_of_ion_species, profile_location
            use KIM_kinds, only: dp

            implicit none

            integer :: i, sigma
            integer :: ierr
            integer :: ios
            character(256) :: fileloc, cwd
            real(dp) :: r_temp

            if (fstatus == 1) write(*,*) 'Status: Reading profiles from text files'

            !call find_file_length(trim(profile_location)//'n.dat', iprof_length)

            ! find profile length
            iprof_length = 0

            open(99, file=trim(profile_location)//'n.dat')
            ios = 0
            do while(ios == 0)
                read(99, *, iostat=ios) r_temp
                if (r_temp < r_plas) then
                    iprof_length = iprof_length + 1
                else 
                    ios = 1
                end if
            end do
            close(99)

            if (.not. allocated(r_prof)) allocate(r_prof(iprof_length), stat=ierr)
            if (ierr /= 0) print *, "array: Allocation request denied"
        
            !write(*,*) r_prof
            allocate(n_prof(iprof_length), Te_prof(iprof_length), Ti_prof(number_of_ion_species, iprof_length), &
            Er_prof(iprof_length), q_prof(iprof_length), ni_prof(number_of_ion_species, iprof_length))

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

            do i=2, number_of_ion_species
                Ti_prof(i,:) = Ti_prof(1,:)
            end do

            do i = 1, iprof_length
                do sigma = 1, number_of_ion_species
                    ! ion density to fulfill quasineutrality
                    ni_prof(sigma, i) = n_prof(i) * Zi(sigma) / sum(Zi)
                end do
            end do


            if (set_profiles_constant == 1) then
                write(*,*) 'Info: Setting profiles to constant values'
                n_prof(:) = n_prof(1)
                Te_prof(:) = Te_prof(1)
                Er_prof(:) = Er_prof(1)
                do sigma = 1, number_of_ion_species
                    ni_prof(sigma, :) = ni_prof(sigma, 1)
                    Ti_prof(sigma, :) = Ti_prof(sigma, 1)
                end do
            end if

            if (fstatus == 1) write(*,*) 'Status: Finished reading profiles from text files'

        end subroutine


        subroutine read_from_hdf5

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

            use config, only: output_path, number_of_ion_species
            use grid, only: reduced_rg_dim
            use KIM_kinds, only: dp

            implicit none

            real(dp) :: step_h
            real(dp), allocatable :: new_n_prof(:), new_Te_prof(:), new_Ti_prof(:,:), &
                                            new_ni_prof(:,:), new_Er_prof(:), new_q_prof(:), new_r_prof(:)
            integer :: i, sigma
            integer :: nlagr = 4
            integer :: nder = 0
            integer :: ibeg, iend, ir
            real(dp), dimension(:,:), allocatable :: coef

            if (.not. allocated(coef)) allocate(coef(0:nder, nlagr))

            if (fstatus==1) write(*,*) 'Status: Reducing input profile r dimension'

            allocate(new_n_prof(reduced_rg_dim), new_Te_prof(reduced_rg_dim), &
                    new_Ti_prof(number_of_ion_species, reduced_rg_dim), new_ni_prof(number_of_ion_species, reduced_rg_dim), &
                    new_Er_prof(reduced_rg_dim), new_q_prof(reduced_rg_dim), new_r_prof(reduced_rg_dim))


            step_h = (r_prof(iprof_length) - r_prof(1)) / reduced_rg_dim ! new step size
            new_r_prof(1) = r_prof(1)

            do i = 2, reduced_rg_dim
                new_r_prof(i) = new_r_prof(i-1) + step_h
            end do

            do i = 1, reduced_rg_dim
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

                do sigma=1, number_of_ion_species 
                    new_Ti_prof(sigma, i) = sum(coef(0,:) * Ti_prof(sigma, ibeg:iend))
                    new_ni_prof(sigma, i) = sum(coef(0,:) * ni_prof(sigma, ibeg:iend))
                end do
            end do
            
            iprof_length = reduced_rg_dim

            deallocate(r_prof, n_prof, Te_prof, Ti_prof, ni_prof, Er_prof, q_prof)
            allocate(r_prof(reduced_rg_dim), n_prof(reduced_rg_dim), Te_prof(reduced_rg_dim), &
                    Ti_prof(number_of_ion_species, reduced_rg_dim), ni_prof(number_of_ion_species, reduced_rg_dim), &
                    Er_prof(reduced_rg_dim), q_prof(reduced_rg_dim))

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

            use config, only: output_path

            implicit none

            integer :: i
            logical :: ex

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

        end subroutine


end module plasma_parameter