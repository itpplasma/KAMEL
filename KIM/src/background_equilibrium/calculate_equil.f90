module equilibrium_m

    use KIM_kinds_m, only: dp

    implicit none
    
    real(dp), allocatable :: B0z(:), B0th(:), B0(:), J0z(:), J0th(:)
    real(dp), allocatable :: hz(:), hth(:)
    
    integer :: ineq = 1 ! numbers of equations to be solved
    integer :: idid     ! indicator reporting what the code did
    real(dp) :: rtol = 1.0d-8 ! relative error tolerance
    real(dp) :: atol = 1.0d-8 ! absolute error tolerance
    integer, dimension(4) :: info      ! info vector to control solver
    integer, parameter :: lrw = 151
    integer, parameter :: liw = 51
    real(dp), dimension(lrw) :: rwork
    integer, dimension(liw) :: iwork
    real(dp) :: rpar
    integer :: ipar

    integer :: i, sigma
    real(dp) :: r1, radius0
    real(dp) :: u0
    real(dp), allocatable :: u(:)
    real(dp), allocatable :: dpress_prof(:)
    real(dp), allocatable :: press_prof(:)

    integer :: nlagr = 4
    integer :: nder = 0
    integer :: ibeg, iend, ir
    real(dp), dimension(:,:), allocatable :: coef


    contains

        subroutine calculate_equil(write_out)
        ! Calculate the magnetic field and current equilibrium from the input profiles
        ! To solve the force balance equation for B0z, which is a first order ODE, we use
        ! the ddeabm subroutine from the SLATEC library. This method uses the Adams-Bashforth
        ! method.

            use species_m, only: plasma, calc_plasma_parameter_derivs
            use constants_m, only: ev, pi, sol
            use setup_m, only: btor, R0, m_mode, n_mode
            use config_m

            implicit none

            logical, intent(in) :: write_out

            if(.not. allocated(coef)) allocate(coef(0:nder, nlagr))

            allocate(u(plasma%grid_size), B0z(plasma%grid_size), dpress_prof(plasma%grid_size), &
                    B0th(plasma%grid_size), B0(plasma%grid_size), hz(plasma%grid_size), hth(plasma%grid_size))

            if (fstatus == 1) write(*,*) 'Status: Calculating equilibrium, write_out=', write_out
            if (.not. allocated(plasma%spec(0)%dndr)) then
                call calc_plasma_parameter_derivs
            end if
    
            dpress_prof = 0.0d0
            do i=1, plasma%grid_size
                do sigma = 0, number_of_ion_species
                    dpress_prof(i) = dpress_prof(i) + plasma%spec(sigma)%dndr(i) * plasma%spec(sigma)%T(i) &
                        + plasma%spec(sigma)%n(i) * plasma%spec(sigma)%dTdr(i)
                end do
            end do
            dpress_prof = dpress_prof * ev

            ! configuration of the solver:
            ! info(1) = 0: initialization, i.e. tell code it is a new problem
            ! info(2) = 0: input scalars for rtol and atol (instead of vectors)
            ! info(3) = 0: solution is only given at TOUT (or r1 in this implementation)
            info = 0 
            ! info(4) = 1: the integration can NOT be carried out without any restrictions
            ! on the independent variable T
            info(4) = 1

            rwork = 0.0d0
            rwork(1) = plasma%r_grid(plasma%grid_size) ! rwork(1) has to be set to the r stopping point

            radius0 = plasma%r_grid(1)
            u0 = btor**2.0d0 * (1.0d0 + radius0**2.0d0 / (R0**2.0d0 * plasma%q(1)**2.0d0)) ! initial value
            u(1) = u0

            do i=2, plasma%grid_size
                r1 = plasma%r_grid(i)
                call ddeabm(dudr, ineq, radius0, u0, r1, info, rtol, atol, idid, rwork, lrw, &
                            iwork, liw, rpar, ipar)

                if (idid .lt. 1) write(*,*) 'Warning: calculate_equil: r=', r1, ' idid=', idid

                u(i) = u0
                info(1) = 1
            end do

            
            allocate(plasma%B0(plasma%grid_size))
            allocate(plasma%ks(plasma%grid_size))
            allocate(plasma%kp(plasma%grid_size))
            allocate(plasma%om_E(plasma%grid_size))

            do i=1, plasma%grid_size
                B0z(i) = sign(1d0, btor) * sqrt(u(i) /(1d0 + plasma%r_grid(i)**2d0/(R0**2d0 * plasma%q(i)**2)))
                B0th(i) = B0z(i) * plasma%r_grid(i) /(plasma%q(i) * R0)
                B0(i) = sqrt(B0th(i)**2d0 + B0z(i)**2d0)
                plasma%B0(i) = B0(i)

                hz(i) = B0z(i) / B0(i)
                hth(i) = B0th(i) / B0(i)

                ! "senkrecht" wavenumber
                plasma%ks(i) = (m_mode * hz(i) - n_mode * hth(i) / R0) / plasma%r_grid(i)
                ! parallel wavenumber
                plasma%kp(i) = m_mode/(plasma%r_grid(i)) * hth(i) + n_mode / R0 * hz(i)
                ! ExB rotation frequency
                plasma%om_E(i) = - sol * plasma%ks(i) * plasma%Er(i) / B0(i)

            end do


            if (write_out) call write_equil

            contains

                subroutine dudr(r, u, du)

                    implicit none

                    real(dp), intent(in) :: r
                    real(dp), intent(in) :: u
                    real(dp), intent(out) :: du

                    real(dp) :: q, dpress, g

                    ! interpolate q and pressure profiles at radial variable
                    call binsrc(plasma%r_grid, 1, plasma%grid_size, r, ir)
                    ibeg = max(1, ir - nlagr/2)
                    iend = ibeg + nlagr - 1
                    if (iend .gt. plasma%grid_size) then
                        iend = plasma%grid_size
                        ibeg = iend - nlagr + 1
                    end if

                    call plag_coeff(nlagr, nder, r, plasma%r_grid(ibeg:iend), coef)

                    q = sum(coef(0,:) * plasma%q(ibeg:iend))
                    dpress = sum(coef(0, :) * dpress_prof(ibeg:iend))
                    
                    du = -2.0d0 * r * u / (q**2.0d0 * R0**2.0d0 + r**2.0d0) - 8.0d0 * pi * dpress

                end subroutine
    
                subroutine write_equil

                    implicit none
                    character(1024) :: filename
                    logical :: ex

                    if(fstatus == 1) write(*,*) 'Status: Writing equilibrium'

                    if (hdf5_output) then

                    else

                        inquire(file=trim(output_path)//'backs', exist=ex)
                        if (.not. ex) then
                            call system('mkdir -p '//trim(output_path)//'backs')
                        end if
                        open(unit = 78, file = trim(output_path)//'backs/'//'B0z.dat')
                        open(unit = 80, file = trim(output_path)//'backs/'//'B0th.dat')
                        open(unit = 81, file = trim(output_path)//'backs/'//'B0.dat')
                        open(unit = 82, file = trim(output_path)//'backs/'//'hz.dat')
                        open(unit = 83, file = trim(output_path)//'backs/'//'hth.dat')
                        open(unit = 79, file = trim(output_path)//'backs/'//'dpress.dat')
                        open(unit = 87, file = trim(output_path)//'backs/'//'p_tot.dat')

                        do i = 1, plasma%grid_size
                            write(78, *) plasma%r_grid(i), B0z(i)
                            write(80, *) plasma%r_grid(i), B0th(i)
                            write(81, *) plasma%r_grid(i), B0(i)
                            write(82, *) plasma%r_grid(i), hz(i)
                            write(83, *) plasma%r_grid(i), hth(i)
                            write(79, *) plasma%r_grid(i), dpress_prof(i)
                        end do
                        close(78)
                        close(79)
                        close(80)
                        close(82)
                        close(81)
                        close(83)
                        close(87)

                    end if


                end subroutine

        end subroutine

end module equilibrium_m
