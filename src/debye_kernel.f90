module debye_kernel
    ! Module used to check Debye shielding case
    ! This case should give a shielded Coulomb potential as a result of
    ! Poisson's equation

    implicit none

    double precision :: lambda_D = 0.0d0

    contains

    subroutine calculate_debye_length

        use config, only: ispecies
        use constants, only: pi, ev, e_charge
        use plasma_parameter, only: Te_prof, Ti_prof, n_prof, ni_prof
        implicit none
        integer :: i

        lambda_D = ((4.0d0 * pi * n_prof(1) * e_charge**2) / (Te_prof(1) * ev))
        do i=1, ispecies
            lambda_D = lambda_D + ((4.0d0 * pi * ni_prof(i,1) * e_charge**2) / (Ti_prof(i,1) * ev))
        end do

        lambda_D = 1/sqrt(lambda_D)

        write(*,*) "Debye length: ", lambda_D
        write(*,*) "n_prof(1): ", n_prof(1)
        write(*,*) "Te_prof(1): ", Te_prof(1)
        
    end subroutine

    double complex function func_debye_kernel(kr, krp)

        use constants, only: pi
        implicit none
        double precision, intent(in) :: kr, krp

        if (lambda_D == 0.0d0) then
            call calculate_debye_length
        end if

        if (kr == krp)then
            func_debye_kernel = -1.d0 / (4.0d0 * pi * lambda_D**2)
        else
            func_debye_kernel = cmplx(0.0d0,0.0d0)
        end if
    end function


end module