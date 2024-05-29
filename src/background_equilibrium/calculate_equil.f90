module equilibrium

    implicit none
    
    double precision, allocatable :: B0z(:), B0th(:), B0(:), J0z(:), J0th(:)
    double precision, allocatable :: hz(:), hth(:)
    
    integer :: ineq = 1 ! numbers of equations to be solved
    integer :: idid     ! indicator reporting what the code did
    double precision :: rtol = 1.0d-8 ! relative error tolerance
    double precision :: atol = 1.0d-8 ! absolute error tolerance
    integer, dimension(4) :: info      ! info vector to control solver
    integer, parameter :: lrw = 151
    integer, parameter :: liw = 51
    double precision, dimension(lrw) :: rwork
    integer, dimension(liw) :: iwork
    double precision :: rpar
    integer :: ipar

    integer :: i, sigma
    double precision :: r1, radius0
    double precision :: u0
    double precision, allocatable :: u(:)
    double precision, allocatable :: dpress_prof(:)
    double precision, allocatable :: press_prof(:)

    integer :: nlagr = 4
    integer :: nder = 0
    integer :: ibeg, iend, ir
    double precision, dimension(:,:), allocatable :: coef


    contains

        subroutine calculate_equil(write_out)
        ! Calculate the magnetic field and current equilibrium from the input profiles
        ! To solve the force balance equation for B0z, which is a first order ODE, we use
        ! the ddeabm subroutine from the SLATEC library. This method uses the Adams-Bashforth
        ! method.

            use plasma_parameter
            use constants, only: ev, pi
            use setup, only: btor, R0
            use config

            implicit none

            logical, intent(in) :: write_out

            if(.not. allocated(coef)) allocate(coef(0:nder, nlagr))

            allocate(u(iprof_length), B0z(iprof_length), dpress_prof(iprof_length), &
                    B0th(iprof_length), B0(iprof_length), hz(iprof_length), hth(iprof_length))

            if (fstatus == 1) write(*,*) 'Status: Calculating equilibrium, write_out=', write_out
            if (.not. allocated(dndr_prof)) then
                call calc_plasma_parameter_derivs
            end if

    
            do i=1, iprof_length
                dpress_prof(i) = ev *(dndr_prof(i) * Te_prof(i) + n_prof(i) * dTedr_prof(i))
            end do

            do sigma=1, ispecies
                press_prof = ev * n_prof * (Te_prof + Ti_prof(sigma, :))
                do i=1, iprof_length
                    dpress_prof(i) = dpress_prof(i) + ev * (dnidr_prof(sigma,i) * &
                        Ti_prof(sigma, i) + ni_prof(sigma,i) * dTidr_prof(sigma,i))
                end do
            end do

            ! configuration of the solver:
            ! info(1) = 0: initialization, i.e. tell code it is a new problem
            ! info(2) = 0: input scalars for rtol and atol (instead of vectors)
            ! info(3) = 0: solution is only given at TOUT (or r1 in this implementation)
            info = 0 
            ! info(4) = 1: the integration can NOT be carried out without any restrictions
            ! on the independent variable T
            info(4) = 1

            rwork = 0.0d0
            rwork(1) = r_prof(iprof_length) ! rwork(1) has to be set to the r stopping point

            radius0 = r_prof(1)
            u0 = btor**2.0d0 * (1.0d0 + radius0**2.0d0 / (R0**2.0d0 * q_prof(1)**2.0d0)) ! initial value
            u(1) = u0

            do i=2, iprof_length
                r1 = r_prof(i)
                call ddeabm(dudr, ineq, radius0, u0, r1, info, rtol, atol, idid, rwork, lrw, &
                            iwork, liw, rpar, ipar)

        if (idid .lt. 1) write(*,*) 'Warning: calculate_equil: r=', r1, ' idid=', idid

                u(i) = u0
                info(1) = 1
            end do

            do i=1, iprof_length
                B0z(i) = sign(1d0, btor) * sqrt(u(i) /(1d0 + r_prof(i)**2d0/(R0**2d0 * q_prof(i)**2)))
                B0th(i) = B0z(i) * r_prof(i) /(q_prof(i) * R0)
                B0(i) = sqrt(B0th(i)**2d0 + B0z(i)**2d0)

                hz(i) = B0z(i) / B0(i)
                hth(i) = B0th(i) / B0(i)
            end do


            if (write_out) call write_equil

            contains

                subroutine dudr(r, u, du)
                    implicit none
                    double precision, intent(in) :: r
                    double precision, intent(in) :: u
                    double precision, intent(out) :: du

                    double precision :: q, dpress, g
            
                    ! interpolate q and pressure profiles at radial variable
                    call binsrc(r_prof, 1, iprof_length, r, ir)
                    ibeg = max(1, ir - nlagr/2)
                    iend = ibeg + nlagr - 1
                    if (iend .gt. iprof_length) then
                        iend = iprof_length
                        ibeg = iend - nlagr + 1
                    end if

                    call plag_coeff(nlagr, nder, r, r_prof(ibeg:iend), coef)

                    q = sum(coef(0,:) * q_prof(ibeg:iend))
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

                        do i = 1, iprof_length
                            write(78, *) r_prof(i), B0z(i)
                            write(80, *) r_prof(i), B0th(i)
                            write(81, *) r_prof(i), B0(i)
                            write(82, *) r_prof(i), hz(i)
                            write(83, *) r_prof(i), hth(i)
                            write(79, *) r_prof(i), dpress_prof(i)
                            write(87, *) r_prof(i), press_prof(i)
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

end module equilibrium
