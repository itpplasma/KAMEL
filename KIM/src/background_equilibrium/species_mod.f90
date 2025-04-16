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
                plasma%spec(sp)%vT(i) = 4.19d7 * sqrt(plasma%spec(sp)%T(i))
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

end module