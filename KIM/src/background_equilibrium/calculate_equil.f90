module equilibrium_m

    use KIM_kinds_m, only: dp

    implicit none

    real(dp), allocatable :: B0z(:), B0th(:), B0(:), J0z(:), J0th(:)
    real(dp), allocatable :: hz(:), hth(:)
    real(dp), allocatable :: equil_grid(:)

    real(dp) :: rtol = 1.0d-12 ! relative error tolerance
    real(dp) :: atol = 1.0d-12 ! absolute error tolerance

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
        ! Calculate the magnetic field and current equilibrium from the input profiles.
        ! The force-balance equation for B0z is a scalar first-order ODE du/dr;
        ! it is integrated over the monotone radial grid with the fortnum
        ! variable-order Adams integrator (ddeabm), a clean-room equivalent of
        ! SLATEC ddeabm. A single ddeabm_init seeds the state at r_grid(1); each
        ! grid value u(i) is then produced by continuing the same re-entrant
        ! state to r_grid(i) (the SLATEC INFO(1)=1 restart), with the integration
        ! barred from stepping past r_grid(end) (the SLATEC RWORK(1)/INFO(4)=1
        ! tstop bound).

            use species_m, only: plasma, calc_plasma_parameter_derivs
            use constants_m, only: ev, pi, sol
            use setup_m, only: btor, R0, m_mode, n_mode
            use config_m, only: number_of_ion_species, output_path, hdf5_output
            use logger_m, only: log_info, log_warning
            use fortnum_ode_ddeabm, only: ddeabm_state_t, ddeabm_init, &
                ddeabm_integrate_to
            use fortnum_status, only: fortnum_status_t, FORTNUM_OK

            implicit none

            logical, intent(in) :: write_out

            real(dp), allocatable :: u_seg(:)
            real(dp) :: rstop
            type(ddeabm_state_t) :: ode_state
            type(fortnum_status_t) :: ode_status

            if(.not. allocated(coef)) allocate(coef(0:nder, nlagr))

            if (allocated(u)) deallocate(u)
            if (allocated(B0z)) deallocate(B0z)
            if (allocated(dpress_prof)) deallocate(dpress_prof)
            if (allocated(B0th)) deallocate(B0th)
            if (allocated(B0)) deallocate(B0)
            if (allocated(hz)) deallocate(hz)
            if (allocated(hth)) deallocate(hth)
            if (allocated(equil_grid)) deallocate(equil_grid)
            allocate(u(plasma%grid_size), &
                    B0z(plasma%grid_size), &
                    dpress_prof(plasma%grid_size), &
                    B0th(plasma%grid_size), &
                    B0(plasma%grid_size), &
                    hz(plasma%grid_size), &
                    hth(plasma%grid_size),&
                    equil_grid(plasma%grid_size))

            equil_grid = plasma%r_grid ! for future modifications

            call log_info('Calculating equilibrium')
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

            radius0 = plasma%r_grid(1)
            u0 = btor**2.0d0 * (1.0d0 + radius0**2.0d0 / (R0**2.0d0 * plasma%q(1)**2.0d0)) ! initial value
            u(1) = u0

            ! Seed the re-entrant integrator at the inner grid point, then carry
            ! the same state forward to every outer point. rstop bars the Adams
            ! step from overshooting the outer boundary.
            rstop = plasma%r_grid(plasma%grid_size)
            call ddeabm_init(ode_state, 1, radius0, [u0])

            do i=2, plasma%grid_size
                r1 = plasma%r_grid(i)
                call ddeabm_integrate_to(dudr, ode_state, r1, rtol, [atol], &
                                         u_seg, ode_status, tstop=rstop)

                if (ode_status%code /= FORTNUM_OK .or. .not. allocated(u_seg)) then
                    block
                        character(len=160) :: wbuf
                        write(wbuf, '(A,ES12.5,A,A)') &
                            'calculate_equil: r=', r1, ' ddeabm: ', trim(ode_status%msg)
                        call log_warning(trim(wbuf))
                    end block
                    u(i) = u(i-1)
                else
                    u(i) = u_seg(1)
                end if
            end do

            if (allocated(plasma%B0)) deallocate(plasma%B0)
            if (allocated(plasma%ks)) deallocate(plasma%ks)
            if (allocated(plasma%kp)) deallocate(plasma%kp)
            if (allocated(plasma%om_E)) deallocate(plasma%om_E)
            allocate(plasma%B0(plasma%grid_size))
            allocate(plasma%ks(plasma%grid_size))
            allocate(plasma%kp(plasma%grid_size))
            allocate(plasma%om_E(plasma%grid_size))

            do i = 1, plasma%grid_size
                ! covariant components of the magnetic field vector
                B0z(i) = sign(1d0, btor) * sqrt(u(i) /(1d0 + plasma%r_grid(i)**2d0/(R0**2d0 * plasma%q(i)**2)))
                B0th(i) = B0z(i) * plasma%r_grid(i) /(plasma%q(i) * R0)

                B0(i) = sqrt(B0th(i)**2d0 + B0z(i)**2d0)

                plasma%B0(i) = B0(i)

                ! covariant components of the magnetic field unit vector
                hz(i) = B0z(i) / B0(i)
                hth(i) = B0th(i) / B0(i)

                ! "senkrecht" wavenumber
                plasma%ks(i) = (m_mode * hz(i) - n_mode * hth(i) / R0) / plasma%r_grid(i)
                ! parallel wavenumber
                plasma%kp(i) = (m_mode/(plasma%r_grid(i)) * hth(i) + n_mode / R0 * hz(i))
                ! ExB rotation frequency
                plasma%om_E(i) = - sol * plasma%ks(i) * plasma%Er(i) / plasma%B0(i)

            end do

            if (write_out) call write_equil

            contains

                subroutine dudr(r, u, du, ctx)
                    ! fortnum ode_rhs_t: scalar force-balance RHS du/dr.
                    ! u and du are size-1; the q-profile and dpress data ride on
                    ! host association, so ctx is unused.

                    implicit none

                    real(dp), intent(in)  :: r
                    real(dp), intent(in)  :: u(:)
                    real(dp), intent(out) :: du(:)
                    class(*), intent(in), optional :: ctx

                    real(dp) :: q, dpress

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

                    du(1) = -2.0d0 * r * u(1) / (q**2.0d0 * R0**2.0d0 + r**2.0d0) - 8.0d0 * pi * dpress

                end subroutine

                subroutine write_equil

                    use IO_collection_m, only: write_profile

                    implicit none

                    call log_info('Writing equilibrium')

                    call write_profile(equil_grid, B0z, size(equil_grid), 'backs/B0z', &
                                        'z component of equilibrium magnetic field', 'G')
                    call write_profile(equil_grid, B0th, size(equil_grid), 'backs/B0th', &
                                        'Poloidal component of equilibrium magnetic field', 'G')
                    call write_profile(equil_grid, B0, size(equil_grid), 'backs/B0', &
                                        'Magnitude of equilibrium magnetic field', 'G')
                    call write_profile(equil_grid, hz, size(equil_grid), 'backs/hz', &
                                        'z direction of equilibrium magnetic field', '1')
                    call write_profile(equil_grid, hth, size(equil_grid), 'backs/hth', &
                                        'Poloidal direction of equilibrium magnetic field', '1')

                end subroutine


        end subroutine


        subroutine interpolate_equil(grid)

        use KIM_kinds_m, only: dp

        implicit none

        real(dp), intent(in) :: grid(:)
        integer :: i
        integer :: nlagr = 4
        integer :: nder = 0
        integer :: ibeg, iend, ir
        real(dp), dimension(:,:), allocatable :: coef
        real(dp), allocatable :: B0z_new(:), B0th_new(:), B0_new(:)
        real(dp), allocatable :: hz_new(:), hth_new(:)

        if (.not. allocated(coef)) allocate(coef(0:nder, nlagr))

        allocate(B0z_new(size(grid)))
        allocate(B0th_new(size(grid)))
        allocate(B0_new(size(grid)))
        allocate(hz_new(size(grid)))
        allocate(hth_new(size(grid)))

        do i = 1, size(grid)
            call binsrc(equil_grid, 1, size(equil_grid), grid(i), ir)
            ibeg = max(1, ir - nlagr/2)
            iend = ibeg + nlagr - 1
            if (iend .gt. size(equil_grid)) then
                iend = size(equil_grid)
                ibeg = iend - nlagr + 1
            end if

            call plag_coeff(nlagr, nder, grid(i), equil_grid(ibeg:iend), coef)

            B0z_new(i)   = sum(coef(0,:) * B0z(ibeg:iend))
            B0th_new(i)  = sum(coef(0,:) * B0th(ibeg:iend))
            B0_new(i)    = sum(coef(0,:) * B0(ibeg:iend))
            hz_new(i)    = sum(coef(0,:) * hz(ibeg:iend))
            hth_new(i)   = sum(coef(0,:) * hth(ibeg:iend))

        end do

        call move_alloc(B0z_new, B0z)
        call move_alloc(B0th_new, B0th)
        call move_alloc(B0_new, B0)
        call move_alloc(hz_new, hz)
        call move_alloc(hth_new, hth)

        equil_grid = grid

    end subroutine

end module equilibrium_m
