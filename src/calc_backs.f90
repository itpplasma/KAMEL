! calculate background quantities, e.g. thermodynamic forces, thermal
! velocity, etc.
! input: write_out ... bool, if true, writes out background quantities
subroutine calc_backs(write_out)

    use back_quants
    use grid
    use constants
    use plas_parameter
    use config
    use setup

    implicit none

    logical, intent(in) :: write_out
    integer :: i  
    integer :: sigma
    double precision, dimension(:), allocatable :: Lee ! Coulomb log
    double precision, dimension(:,:), allocatable :: Lei ! Coulomb log

    if (fstatus == 1) write(*,*) 'Status: Calculating background quantities'

    allocate(Lee(iprof_length), Lei(ispecies, iprof_length))

    allocate(A1e(iprof_length), A2e(iprof_length), vTe(iprof_length), &
            omce(iprof_length), nue(iprof_length))
    allocate(A1i(ispecies, iprof_length), A2i(ispecies, iprof_length), &
            vTi(ispecies, iprof_length), omci(ispecies, iprof_length), &
            nui(ispecies, iprof_length))

    do i=1, iprof_length
        Lee(i) = 23.5d0 - log(sqrt(n_prof(i)) / Te_prof(i)**1.25) - &
        sqrt(1d-5 + (log(Te_prof(i)) -2.0)**2.0 / 16.0)

        vTe(i) = sqrt(Te_prof(i) * kB / e_mass)
        omce(i) = e_charge * btor / (e_mass * sol)
        nue(i) = 5.8e-6 * n_prof(i) * Lee(i) / Te_prof(i)**(3.0/2.0)
    end do

    do sigma=1, ispecies
        do i=1, iprof_length
            Lei(sigma, i) = 24.0d0 - log(sqrt(n_prof(i)) / Ti_prof(i))

            vTi(sigma, i) = sqrt(Ti_prof(i) * kB / (p_mass * Ai(sigma)))
            omci(sigma, i) = (e_charge * Zi(sigma)) * btor / (p_mass * Ai(sigma) * sol)

            nue(i) = nue(i) + 7.7d-6 * n_prof(i) * Lei(sigma, i) * Zi(sigma)**2 / Te_prof(i)**(3.0/2.0)

        end do
    end do

    if (write_out) then
        call write_backs
    endif

    contains
    ! write background quantities
    subroutine write_backs
        
        implicit none
        character(1024) :: filename

        if (fstatus == 1) write(*,*) 'Status: Writing background quants'

        if (hdf5_output) then
            ! write to hdf5
        else
            ! write to text files
            open(unit = 78, file = trim(output_path)//'backs/'//'vTe.dat')
            open(unit = 79, file = trim(output_path)//'backs/'//'nue.dat')
            do i=1, iprof_length
                write(78, *) r_prof(i), vTe(i)
                write(79, *) r_prof(i), nue(i)
            end do
            close(unit = 78)
            close(unit = 79)

            if (ispecies == 1) then
                open(unit = 78, file = trim(output_path)//'backs/'//'vTi.dat')
                do i=1, iprof_length
                    write(78, *) r_prof(i), vTi(1, i)
                end do
                close(unit = 78)
            else
                do sigma = 1, ispecies
                    write(filename, "(A4, I1, A4)") 'vT_', sigma, '.dat'
                    open(unit = 78, file = trim(output_path)//'backs/'//filename)
                    do i=1, iprof_length
                        write(78, *) r_prof(i), vTi(sigma, i)
                    end do
                    close(unit = 78)               
                end do
            end if

        endif
    end subroutine

end subroutine