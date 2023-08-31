! calculate background quantities, e.g. thermodynamic forces, thermal
! velocity, etc.
! input: write_out ... bool, if true, writes out background quantities
!
! Comments: Formulas from the Coulomb logarithms and the collision frequencies
! are from the 2019 NRL Formulary
subroutine calc_backs(write_out)

    use back_quants
    use grid
    use constants
    use plas_parameter
    use config
    use setup
    use equil, only: hz, hth

    implicit none

    logical, intent(in) :: write_out
    integer :: i  
    integer :: sigma, sigma_col
    double precision, dimension(:), allocatable :: Lee ! Coulomb log
    double precision, dimension(:,:), allocatable :: Lei ! Coulomb log
    double precision, dimension(:,:,:), allocatable :: Lii ! Coulomb log

    if (fstatus == 1) write(*,*) 'Status: Calculating background quantities'

    allocate(Lee(iprof_length), Lei(ispecies, iprof_length), &
            Lii(ispecies, ispecies, iprof_length))

    allocate(A1e(iprof_length), A2e(iprof_length), vTe(iprof_length), &
            omce(iprof_length), nue(iprof_length), lambda_De(iprof_length))
    allocate(A1i(ispecies, iprof_length), A2i(ispecies, iprof_length), &
            vTi(ispecies, iprof_length), omci(ispecies, iprof_length), &
            nui(ispecies, iprof_length), lambda_Di(ispecies, iprof_length))
    allocate(ks(iprof_length))

    call calc_plas_parameter_derivs

    do i=1, iprof_length
        ! Coulomb logarithm
        Lee(i) = 23.5d0 - log(sqrt(n_prof(i)) / Te_prof(i)**1.25) - &
        sqrt(1d-5 + (log(Te_prof(i)) -2.0)**2.0 / 16.0)
        ! Thermal velocity
        vTe(i) = sqrt(Te_prof(i) * kB / e_mass)
        ! cyclotron frequency
        omce(i) = e_charge * btor / (e_mass * sol)
        ! Collision frequency
        nue(i) = 5.8e-6 * n_prof(i) * Lee(i) / Te_prof(i)**(3.0/2.0)
        ! Debye length
        lambda_De(i) = sqrt(Te_prof(i) *kB/ (4*pi*n_prof(i) * e_charge**2))
        ! First thermodynamic force
        A1e(i) = dndr_prof(i) / n_prof(i) + e_charge/Te_prof(i) * Er_prof(i) - 3/(2*Te_prof(i)) * dTedr_prof(i)
        ! Second thermodynamic force
        A2e(i) = dTedr_prof(i) / Te_prof(i)

        ! "senkrecht" wavenumber
        ks(i) = (m_mode * hz(i) - n_mode * hth(i) / R0) / r_prof(i)
    end do

    do sigma=1, ispecies
        do i=1, iprof_length
            ! Coulomb logarithm electrons ions (= ions electrons)
            Lei(sigma, i) = 24.0d0 - log(sqrt(n_prof(i)) / Ti_prof(sigma, i))
            ! thermal velocity
            vTi(sigma, i) = sqrt(Ti_prof(sigma, i) * kB / (p_mass * Ai(sigma)))
            ! Cyclotron frequency
            omci(sigma, i) = (e_charge * Zi(sigma)) * btor / (p_mass * Ai(sigma) * sol)
            ! Collision frequency of electrons with ions
            nue(i) = nue(i) + 7.7d-6 * n_prof(i) * Lei(sigma, i) * Zi(sigma)**2 / Te_prof(i)**(3.0/2.0)
            ! Collision frequency ions with electrons
            nui(sigma, i) = 1.8d-7 * Ai(sigma)**(-1.0/2.0) * Ti_prof(sigma, i)**(-3.0/2.0) * n_prof(i) * &
                            Zi(sigma)**2 * Lei(sigma,i)

            do sigma_col=sigma, ispecies
                ! Coulomb logarithm ions - ions'
                Lii(sigma, sigma_col, i) = 23.0 - log(Zi(sigma) * Zi(sigma_col) * (Ai(sigma)+Ai(sigma_col)) /&
                                         (Ti_prof(sigma, i)*Ai(sigma) + Ti_prof(sigma_col, i) * Ai(sigma_col)) * &
                                          (ni_prof(sigma, i) * Zi(sigma)**2 / Ti_prof(sigma,i) + &
                                           ni_prof(sigma_col, i) * Zi(sigma_col)**2) / Ti_prof(sigma_col, i))
                ! Collision frequency ions - ions'
                nui(sigma, i) = nui(sigma, i) + 1.8d-7 * ni_prof(sigma_col, i) * Zi(sigma)**2 * &
                                Zi(sigma_col)**2 * Lii(sigma, sigma_col, i) * Ai(sigma)**(-1.0/2.0)&
                                * Ti_prof(sigma, i)**(-3.0/2.0)
            end do

            ! Debye length ions
            lambda_Di(sigma, i) = sqrt(Ti_prof(sigma, i) * kB/ (4*pi*ni_prof(sigma,i) * (e_charge*Zi(sigma))**2))
            ! First thermodynamic force
            A1i(sigma, i) = dnidr_prof(sigma, i) / ni_prof(sigma, i) - (e_charge*Zi(sigma))/Ti_prof(sigma, i) * Er_prof(i)&
                        - 3/(2*Ti_prof(sigma, i)) * dTidr_prof(sigma, i)
            ! Second thermodynamic force
            A2i(sigma, i) = dTidr_prof(sigma, i) / Ti_prof(sigma, i)
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
            open(unit = 80, file = trim(output_path)//'backs/'//'lambda_De.dat')
            open(unit = 81, file = trim(output_path)//'backs/'//'dndr.dat')
            open(unit = 82, file = trim(output_path)//'backs/'//'dTedr.dat')
            open(unit = 85, file = trim(output_path)//'backs/'//'dqdr.dat')
            open(unit = 86, file = trim(output_path)//'backs/'//'A1e.dat')
            open(unit = 87, file = trim(output_path)//'backs/'//'A2e.dat')
            open(unit = 88, file = trim(output_path)//'backs/'//'ks.dat')
            do i=1, iprof_length
                write(78, *) r_prof(i), vTe(i)
                write(79, *) r_prof(i), nue(i)
                write(80, *) r_prof(i), lambda_De(i)
                write(81, *) r_prof(i), dndr_prof(i)
                write(82, *) r_prof(i), dTedr_prof(i)
                write(85, *) r_prof(i), dqdr_prof(i)
                write(86, *) r_prof(i), A1e(i)
                write(87, *) r_prof(i), A2e(i)
                write(88, *) r_prof(i), ks(i)
            end do
            close(unit = 78)
            close(unit = 79)
            close(unit = 80)
            close(unit = 81)
            close(unit = 82)
            close(unit = 83)
            close(unit = 84)
            close(unit = 85)
            close(unit = 86)
            close(unit = 87)
            close(unit = 88)

            if (ispecies == 1) then
                open(unit = 78, file = trim(output_path)//'backs/'//'vTi.dat')
                open(unit = 79, file = trim(output_path)//'backs/'//'nui.dat')
                open(unit = 80, file = trim(output_path)//'backs/'//'lambda_Di.dat')
                open(unit = 81, file = trim(output_path)//'backs/'//'dTidr.dat')
                open(unit = 82, file = trim(output_path)//'backs/'//'dnidr.dat')
                open(unit = 83, file = trim(output_path)//'backs/'//'A1i.dat')
                open(unit = 84, file = trim(output_path)//'backs/'//'A2i.dat')
                do i=1, iprof_length
                    write(78, *) r_prof(i), vTi(1, i)
                    write(79, *) r_prof(i), nui(1, i)
                    write(80, *) r_prof(i), lambda_Di(1, i)
                    write(81, *) r_prof(i), dTidr_prof(1, i)
                    write(82, *) r_prof(i), dnidr_prof(1, i)
                    write(83, *) r_prof(i), A1i(1, i)
                    write(84, *) r_prof(i), A2i(1, i)
                end do
                close(unit = 78)
                close(unit = 79)
                close(unit = 80)
                close(unit = 81)
                close(unit = 82)
                close(unit = 83)
                close(unit = 84)
            else
                do sigma = 1, ispecies
                    ! thermal velocity
                    write(filename, "(A3, I1, A4)") 'vT_', sigma, '.dat'
                    open(unit = 78, file = trim(output_path)//'backs/'//filename)
                    do i=1, iprof_length
                        write(78, *) r_prof(i), vTi(sigma, i)
                    end do
                    close(unit = 78)               
                    ! Collision frequency
                    write(filename, "(A4, I1, A4)") 'nui_', sigma, '.dat'
                    open(unit = 78, file = trim(output_path)//'backs/'//filename)
                    do i=1, iprof_length
                        write(78, *) r_prof(i), nui(sigma, i)
                    end do
                    close(unit = 78)
                    ! Debye length
                    write(filename, "(A10, I1, A4)") 'lambda_Di_', sigma, '.dat'
                    open(unit = 78, file = trim(output_path)//'backs/'//filename)
                    do i=1, iprof_length
                        write(78, *) r_prof(i), lambda_Di(sigma, i)
                    end do
                    close(unit = 78)
                    ! density gradient
                    write(filename, "(A6, I1, A4)") 'dnidr_', sigma, '.dat'
                    open(unit = 78, file = trim(output_path)//'backs/'//filename)
                    do i=1, iprof_length
                        write(78, *) r_prof(i), dnidr_prof(sigma, i)
                    end do
                    close(unit = 78)
                    ! temperature gradient
                    write(filename, "(A6, I1, A4)") 'dTidr_', sigma, '.dat'
                    open(unit = 78, file = trim(output_path)//'backs/'//filename)
                    do i=1, iprof_length
                        write(78, *) r_prof(i), dTidr_prof(sigma, i)
                    end do
                    close(unit = 78)
                    ! first thermodynamic force
                    write(filename, "(A4, I1, A4)") 'A1i_', sigma, '.dat'
                    open(unit = 78, file = trim(output_path)//'backs/'//filename)
                    do i=1, iprof_length
                        write(78, *) r_prof(i), A1i(sigma, i)
                    end do
                    close(unit = 78)
                    ! second thermodynamic force
                    write(filename, "(A4, I1, A4)") 'A2i_', sigma, '.dat'
                    open(unit = 78, file = trim(output_path)//'backs/'//filename)
                    do i=1, iprof_length
                        write(78, *) r_prof(i), A2i(sigma, i)
                    end do
                    close(unit = 78)
                end do
            end if

        endif
    end subroutine

end subroutine


subroutine calc_plas_parameter_derivs

    use plas_parameter
    use grid
    use config

    implicit none

    integer :: i
    integer :: iend
    integer :: sigma

    if (.not. allocated(dndr_prof)) then
        write(*,*) 'dndr not allocated'
        allocate(dndr_prof(iprof_length), dTedr_prof(iprof_length), dTidr_prof(ispecies,iprof_length), &
        dqdr_prof(iprof_length), dnidr_prof(ispecies, iprof_length))
    end if

    iend = iprof_length - 1

    do i=1, iend
        dndr_prof(i) = (n_prof(i+1) - n_prof(i))/(r_prof(i+1)- r_prof(i))

        dTedr_prof(i) = (Te_prof(i+1) - Te_prof(i))/(r_prof(i+1)- r_prof(i))
        dqdr_prof(i) = (q_prof(i+1) - q_prof(i))/(r_prof(i+1)- r_prof(i))
        do sigma=1, ispecies
            dnidr_prof(sigma,i) = (ni_prof(sigma,i+1) - ni_prof(sigma,i))/(r_prof(i+1)- r_prof(i))
            dTidr_prof(sigma, i) = (Ti_prof(sigma, i+1) - Ti_prof(sigma, i))/(r_prof(i+1)- r_prof(i))
        end do
    end do


end subroutine
