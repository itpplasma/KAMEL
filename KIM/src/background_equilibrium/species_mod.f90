module species

    use KIM_kinds, only: dp

    implicit none

    type :: plasma_t
        type(species_t), allocatable :: spec(:)
        integer :: n_species
        real(dp), allocatable :: om_E(:)
        real(dp), allocatable :: ks(:)
        real(dp), allocatable :: kp(:)
        real(dp), allocatable :: q(:) ! safety factor profile
        real(dp), allocatable :: dqdr(:)
        real(dp), allocatable :: Er(:)
        real(dp), allocatable :: r_grid(:)
    end type

    type :: species_t
        character(len=2) :: name
        integer :: Aspec ! mass number
        integer :: Zspec ! charge number
        real(dp) :: mass ! mass in g
        real(dp), allocatable :: n(:) ! density
        real(dp), allocatable :: dndr(:) ! density gradient
        real(dp), allocatable :: T(:) ! temperature
        real(dp), allocatable :: dTdr(:) ! temperature gradient
        real(dp), allocatable :: A1(:) ! thermodynamic force
        real(dp), allocatable :: A2(:) ! thermodynamic force
        real(dp), allocatable :: nu(:) ! collision frequency
        real(dp), allocatable :: vT(:) ! thermal velocity
        real(dp), allocatable :: omega_c(:) ! cyclotron frequency
        real(dp), allocatable :: lambda_D(:) ! Debye length
        real(dp), allocatable :: rho_L(:) ! Larmor radius
        real(dp), allocatable :: z0(:) ! Larmor radius
    end type

    type(plasma_t) :: plasma

    contains

    subroutine init_deuterium_plasma(plasma)

        implicit none

        type(plasma_t), intent(inout) :: plasma

        plasma%n_species = 2
        allocate(plasma%spec(0:plasma%n_species-1))
        call init_electron_species(plasma%spec(0))
        call init_deuterium_species(plasma%spec(1))

    end subroutine

    subroutine init_deuterium_species(deut)
    
        use constants, only: p_mass

        implicit none

        type(species_t), intent(inout) :: deut

        deut%name = 'i'
        deut%Aspec = 2
        deut%Zspec = 1
        deut%mass = p_mass * deut%Aspec

    end subroutine

    subroutine init_electron_species(elec)

        use constants, only: e_mass
        implicit none

        type(species_t), intent(inout) :: elec

        elec%name = 'e'
        elec%Aspec = 1
        elec%Zspec = -1
        elec%mass = e_mass

    end subroutine

    subroutine set_species_profiles(spec, n_prof, dndr_prof, T_prof, dTdr_prof, prof_length)

        implicit none

        type(species_t), intent(inout) :: spec
        integer, intent(in) :: prof_length
        real(dp), intent(in) :: n_prof(prof_length), dndr_prof(prof_length)
        real(dp), intent(in) :: T_prof(prof_length), dTdr_prof(prof_length)

        allocate(spec%n(prof_length), spec%dndr(prof_length))
        allocate(spec%T(prof_length), spec%dTdr(prof_length))

        spec%n = n_prof
        spec%T = T_prof
        spec%dndr = dndr_prof
        spec%dTdr = dTdr_prof

    end subroutine


    subroutine set_deuterium_plasma(plasma)

        use plasma_parameter, only: iprof_length, n_prof, Te_prof, Ti_prof, dTidr_prof, &
            dndr_prof, dTedr_prof, dnidr_prof, dqdr_prof, r_prof, q_prof, Er_prof

        implicit none

        type(plasma_t), intent(inout) :: plasma

        allocate(plasma%q(iprof_length))
        allocate(plasma%dqdr(iprof_length))
        allocate(plasma%Er(iprof_length))
        allocate(plasma%ks(iprof_length))
        allocate(plasma%kp(iprof_length))
        allocate(plasma%om_E(iprof_length))

        plasma%q = q_prof
        plasma%Er = Er_prof
        plasma%dqdr = dqdr_prof
        plasma%r_grid = r_prof

        call calc_plasma_parameter_derivs

        call set_species_profiles(plasma%spec(0), n_prof, dndr_prof, Te_prof, dTedr_prof, iprof_length)
        call set_species_profiles(plasma%spec(1), n_prof, dnidr_prof(1,:), Ti_prof, dTidr_prof(1,:), iprof_length)
        
        call calculate_plasma_backs(plasma)

    end subroutine

    subroutine calculate_plasma_backs(plasma)

        use constants, only: sol, e_charge, ev, pi, com_unit
        use equilibrium, only: hz, hth, B0
        use plasma_parameter, only: iprof_length, r_prof
        use setup, only: m_mode, n_mode, omega, R0
        use grid, only: rg_grid
        use plotting, only: plot_profile

        implicit none

        type(plasma_t), intent(inout) :: plasma
        integer :: i, sp

        do i = 1, iprof_length
            ! "senkrecht" wavenumber
            plasma%ks(i) = (m_mode * hz(i) - n_mode * hth(i) / R0) / r_prof(i)
            ! parallel wavenumber
            plasma%kp(i) = m_mode/(r_prof(i)) * hth(i) + n_mode / R0 * hz(i)
            ! ExB rotation frequency
            plasma%om_E(i) = - sol * plasma%ks(i) * plasma%Er(i) / B0(i)
        end do

        call plot_profile(r_prof, plasma%kp)

        do sp = 0, plasma%n_species-1

            allocate(plasma%spec(sp)%rho_L(iprof_length))
            allocate(plasma%spec(sp)%z0(iprof_length))
            allocate(plasma%spec(sp)%vT(iprof_length))
            allocate(plasma%spec(sp)%lambda_D(iprof_length))
            allocate(plasma%spec(sp)%A1(iprof_length))
            allocate(plasma%spec(sp)%A2(iprof_length))
            allocate(plasma%spec(sp)%nu(iprof_length))
            allocate(plasma%spec(sp)%omega_c(iprof_length))

            do i=1, iprof_length
                plasma%spec(sp)%vT(i) = sqrt(plasma%spec(sp)%T(i) * ev / (plasma%spec(sp)%mass))
                plasma%spec(sp)%omega_c(i) = plasma%spec(sp)%Zspec * e_charge * abs(B0(i)) &
                    / (plasma%spec(sp)%mass * sol)
                !nue(i) = 5.8e-6 * n_prof(i) * Lee(i) / Te_prof(i)**(3.0/2.0)
                plasma%spec(sp)%nu(i) = 0.0d0
                plasma%spec(sp)%lambda_D(i) = sqrt(plasma%spec(sp)%T(i) *ev / (4.0d0*pi* plasma%spec(sp)%n(i) &
                    * (plasma%spec(sp)%Zspec * e_charge)**2.0d0))
                plasma%spec(sp)%A1(i) = plasma%spec(sp)%dndr(i) / plasma%spec(sp)%n(i) - plasma%spec(sp)%Zspec *e_charge&
                    /(plasma%spec(sp)%T(i) * ev) * plasma%Er(i) - 3/(2*plasma%spec(sp)%T(i)) * plasma%spec(sp)%dTdr(i)
                plasma%spec(sp)%A2(i) = plasma%spec(sp)%dTdr(i) / plasma%spec(sp)%T(i)
                plasma%spec(sp)%rho_L(i) = plasma%spec(sp)%vT(i) / abs(plasma%spec(sp)%omega_c(i))
                plasma%spec(sp)%z0(i) = - (plasma%om_E(i) - omega - com_unit * plasma%spec(sp)%nu(i)) &
                    / (abs(plasma%kp(i)) * sqrt(2d0) * plasma%spec(sp)%vT(i) )
            end do

        end do

        

    end subroutine

    subroutine interpolate_plasma_backs(plasma_in, grid)

        use KIM_kinds, only: dp

        implicit none

        type(plasma_t), intent(inout) :: plasma_in
        type(plasma_t), allocatable :: plasma_temp
        real(dp), intent(in) :: grid(:)
        integer :: i, sp
        integer :: nlagr = 4
        integer :: nder = 0
        integer :: ibeg, iend, ir
        real(dp), dimension(:,:), allocatable :: coef

        if (.not. allocated(coef)) allocate(coef(0:nder, nlagr))

        plasma_temp = plasma_in

        do sp = 0, plasma_temp%n_species-1
            deallocate(plasma_temp%spec(sp)%n)
            allocate(plasma_temp%spec(sp)%n(size(grid)))
            deallocate(plasma_temp%spec(sp)%dndr)
            allocate(plasma_temp%spec(sp)%dndr(size(grid)))
            deallocate(plasma_temp%spec(sp)%T)
            allocate(plasma_temp%spec(sp)%T(size(grid)))
            deallocate(plasma_temp%spec(sp)%dTdr)
            allocate(plasma_temp%spec(sp)%dTdr(size(grid)))
            deallocate(plasma_temp%spec(sp)%A1)
            allocate(plasma_temp%spec(sp)%A1(size(grid)))
            deallocate(plasma_temp%spec(sp)%A2)
            allocate(plasma_temp%spec(sp)%A2(size(grid)))
            deallocate(plasma_temp%spec(sp)%nu)
            allocate(plasma_temp%spec(sp)%nu(size(grid)))
            deallocate(plasma_temp%spec(sp)%vT)
            allocate(plasma_temp%spec(sp)%vT(size(grid)))
            deallocate(plasma_temp%spec(sp)%omega_c)
            allocate(plasma_temp%spec(sp)%omega_c(size(grid)))
            deallocate(plasma_temp%spec(sp)%lambda_D)
            allocate(plasma_temp%spec(sp)%lambda_D(size(grid)))
            deallocate(plasma_temp%spec(sp)%rho_L)
            allocate(plasma_temp%spec(sp)%rho_L(size(grid)))
            deallocate(plasma_temp%spec(sp)%z0)
            allocate(plasma_temp%spec(sp)%z0(size(grid)))

            deallocate(plasma_temp%kp)
            allocate(plasma_temp%kp(size(grid)))
            deallocate(plasma_temp%ks)
            allocate(plasma_temp%ks(size(grid)))
            deallocate(plasma_temp%om_E)
            allocate(plasma_temp%om_E(size(grid)))

        end do

        do sp = 0, plasma_temp%n_species-1
            do i = 1, size(grid)
                call binsrc(plasma_temp%r_grid, 1, size(plasma_temp%r_grid), grid, ir) 
                ibeg = max(1, ir - nlagr/2)
                iend = ibeg + nlagr - 1
                if (iend .gt. size(plasma_temp%r_grid)) then
                    iend = size(plasma_temp%r_grid)
                    ibeg = iend -nlagr + 1
                end if

                call plag_coeff(nlagr, nder, grid(i), plasma_temp%r_grid(ibeg:iend), coef)

                plasma_temp%spec(sp)%n(i) = sum(coef(0,:) * plasma_in%spec(sp)%n(ibeg:iend))
                plasma_temp%spec(sp)%dndr(i) = sum(coef(0,:) * plasma_in%spec(sp)%dndr(ibeg:iend))
                plasma_temp%spec(sp)%T(i) = sum(coef(0,:) * plasma_in%spec(sp)%T(ibeg:iend))
                plasma_temp%spec(sp)%dTdr(i) = sum(coef(0,:) * plasma_in%spec(sp)%dTdr(ibeg:iend))
                plasma_temp%spec(sp)%A1(i) = sum(coef(0,:) * plasma_in%spec(sp)%A1(ibeg:iend))
                plasma_temp%spec(sp)%A2(i) = sum(coef(0,:) * plasma_in%spec(sp)%A2(ibeg:iend))
                plasma_temp%spec(sp)%nu(i) = sum(coef(0,:) * plasma_in%spec(sp)%nu(ibeg:iend))
                plasma_temp%spec(sp)%vT(i) = sum(coef(0,:) * plasma_in%spec(sp)%vT(ibeg:iend))
                plasma_temp%spec(sp)%omega_c(i) = sum(coef(0,:) * plasma_in%spec(sp)%omega_c(ibeg:iend))
                plasma_temp%spec(sp)%lambda_D(i) = sum(coef(0,:) * plasma_in%spec(sp)%lambda_D(ibeg:iend))
                plasma_temp%spec(sp)%rho_L(i) = sum(coef(0,:) * plasma_in%spec(sp)%rho_L(ibeg:iend))
                plasma_temp%spec(sp)%z0(i) = sum(coef(0,:) * plasma_in%spec(sp)%z0(ibeg:iend))

                plasma_temp%ks(i) = sum(coef(0,:) * plasma_in%ks(ibeg:iend))
                plasma_temp%kp(i) = sum(coef(0,:) * plasma_in%kp(ibeg:iend))
                plasma_temp%om_E(i) = sum(coef(0,:) * plasma_in%om_E(ibeg:iend))
                plasma_temp%q(i) = sum(coef(0,:) * plasma_in%q(ibeg:iend))
                plasma_temp%dqdr(i) = sum(coef(0,:) * plasma_in%dqdr(ibeg:iend))
                plasma_temp%Er(i) = sum(coef(0,:) * plasma_in%Er(ibeg:iend))

            end do
        end do

        plasma_temp%r_grid = grid
        plasma_in = plasma_temp

        deallocate(plasma_temp)

    end subroutine


    subroutine plot_species(spec)

        use plotting, only: plot_1D_labeled, write_profile, remove_file
        use grid, only: rg_grid

        implicit none

        type(species_t), intent(in) :: spec

        call write_profile(rg_grid%xb, spec%rho_L, rg_grid%npts_b, 'rho_L.dat')
        call plot_1D_labeled('rho_L.dat', ' r [cm] ', ' rho_L [cm] ', '')
        call remove_file('rho_L.dat')

    end subroutine

end module