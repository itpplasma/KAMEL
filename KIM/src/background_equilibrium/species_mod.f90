module species

    use KIM_kinds, only: dp

    implicit none

    type :: plasma_t
        type(species_t), allocatable :: spec(:)
        integer :: n_species
        real(dp), allocatable :: om_E(:) ! ExB rotation frequency
        real(dp), allocatable :: ks(:) ! "senkrecht" wavenumber
        real(dp), allocatable :: kp(:) ! parallel wavenumber
        real(dp), allocatable :: q(:) ! safety factor profile
        real(dp), allocatable :: dqdr(:) ! safety factor gradient
        real(dp), allocatable :: Er(:) ! radial electric field
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
        real(dp), allocatable :: x1(:) ! normalized distance to res. surface
        real(dp), allocatable :: x2(:) ! normalized inverse collisionality
        complex(dp), allocatable :: z0(:) ! = x2/(sqrt(2) x1), argument of plasma dispersion function
        complex(dp), dimension(:,:), allocatable :: symbI ! susceptibility function used for Fokker-Planck collision model
        ! susceptibility function profiles
        complex(dp), allocatable :: I00(:)
        complex(dp), allocatable :: I20(:)
        complex(dp), allocatable :: I01(:)
        complex(dp), allocatable :: I21(:)
    end type

    type(plasma_t) :: plasma

    integer, parameter :: nmmax = 3

    contains

    subroutine read_species_from_nml(plasma_in)

        use KIM_kinds, only: dp
        use config, only: number_of_ion_species, nml_config_path
        use constants, only: p_mass

        implicit none

        type(plasma_t), intent(inout) :: plasma_in
        integer :: ai(16), zi(16)
        integer :: i


        namelist /KIM_species/ ai, zi

        open(unit=77, file=nml_config_path)
        read(unit=77, nml=KIM_species)
        close(unit=77)

        do i = 1, number_of_ion_species
            plasma_in%spec(i)%Aspec = ai(i)
            plasma_in%spec(i)%Zspec = zi(i)
            plasma_in%spec(i)%name = 'i'
            plasma_in%spec(i)%mass = p_mass * plasma_in%spec(i)%Aspec
        end do

    end subroutine

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
        use grid, only: rg_grid

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

        call interpolate_plasma_backs(plasma, rg_grid%xb)

        call write_species_backs(plasma%spec(0), plasma%r_grid)
        call write_species_backs(plasma%spec(1), plasma%r_grid)
        call write_plasma_backs(plasma, plasma%r_grid)

    end subroutine

    subroutine calculate_plasma_backs(plasma)

        use constants, only: sol, e_charge, ev, pi, com_unit
        use equilibrium, only: hz, hth, B0
        use plasma_parameter, only: iprof_length, r_prof
        use setup, only: m_mode, n_mode, omega, R0, collisions_off
        use config, only: number_of_ion_species

        implicit none

        type(plasma_t), intent(inout) :: plasma
        integer :: i, sp, sp_col
        real(dp) :: Lee(iprof_length), Lei(number_of_ion_species, iprof_length), &
            Lii(number_of_ion_species, number_of_ion_species, iprof_length),&
            nue(iprof_length), nui(number_of_ion_species, iprof_length)

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
            allocate(plasma%spec(sp)%x1(iprof_length))
            allocate(plasma%spec(sp)%x2(iprof_length))
            allocate(plasma%spec(sp)%vT(iprof_length))
            allocate(plasma%spec(sp)%lambda_D(iprof_length))
            allocate(plasma%spec(sp)%A1(iprof_length))
            allocate(plasma%spec(sp)%A2(iprof_length))
            allocate(plasma%spec(sp)%I00(iprof_length))
            allocate(plasma%spec(sp)%I20(iprof_length))
            allocate(plasma%spec(sp)%I01(iprof_length))
            allocate(plasma%spec(sp)%I21(iprof_length))
            allocate(plasma%spec(sp)%nu(iprof_length))
            allocate(plasma%spec(sp)%omega_c(iprof_length))

            do i=1, iprof_length
                plasma%spec(sp)%vT(i) = sqrt(plasma%spec(sp)%T(i) * ev / (plasma%spec(sp)%mass))
                plasma%spec(sp)%omega_c(i) = plasma%spec(sp)%Zspec * e_charge * abs(B0(i)) &
                    / (plasma%spec(sp)%mass * sol)
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

        do i=1, iprof_length
            ! Coulomb logarithm
            Lee(i) = 23.5d0 - log(sqrt(plasma%spec(0)%n(i)) / plasma%spec(0)%T(i)**1.25) - &
                sqrt(1d-5 + (log(plasma%spec(0)%T(i)) -2.0)**2.0 / 16.0)
            nue(i) = 5.8e-6 * plasma%spec(0)%n(i) * Lee(i) / plasma%spec(0)%T(i)**(3.0/2.0)
        end do

        do sp=1, number_of_ion_species
            do i=1, iprof_length
                ! Coulomb logarithm electrons ions (= ions electrons)
                Lei(sp, i) = 24.0d0 - log(sqrt(plasma%spec(0)%n(i)) / plasma%spec(sp)%T(i))

                nue(i) = nue(i) + 7.7d-6 * plasma%spec(sp)%n(i) * Lei(sp, i) * plasma%spec(sp)%Zspec**2 / plasma%spec(0)%T(i)**(3.0/2.0)

                nui(sp, i) = 1.8d-7 * plasma%spec(sp)%Aspec**(-1.0/2.0) * plasma%spec(sp)%T(i)**(-3.0/2.0) * plasma%spec(0)%n(i) * &
                                plasma%spec(sp)%Zspec**2 * Lei(sp,i)

                do sp_col=sp, number_of_ion_species
                    ! Coulomb logarithm ions - ions'
                    Lii(sp, sp_col, i) = 23.0d0 - log(plasma%spec(sp)%Zspec * plasma%spec(sp_col)%Zspec &
                                        * (plasma%spec(sp)%Aspec + plasma%spec(sp_col)%Aspec)&
                                        / (plasma%spec(sp)%T(i) * plasma%spec(sp_col)%Aspec + plasma%spec(sp_col)%T(i) * plasma%spec(sp)%Aspec)&
                                        * (plasma%spec(sp)%n(i) * plasma%spec(sp)%Zspec**2 / plasma%spec(sp)%T(i) &
                                        + plasma%spec(sp_col)%n(i) * plasma%spec(sp_col)%Zspec**2) / plasma%spec(sp_col)%T(i))
                    ! Collision frequency ions - ions'
                    nui(sp, i) = nui(sp, i) + 1.8d-7 * plasma%spec(sp_col)%n(i) * plasma%spec(sp)%Zspec**2 &
                                    * plasma%spec(sp_col)%Zspec**2 * Lii(sp, sp_col, i) * plasma%spec(sp)%Aspec**(-1.0/2.0)&
                                    * plasma%spec(sp)%T(i)**(-3.0/2.0)
                end do
            end do
        end do

        do i = 1, iprof_length
            plasma%spec(0)%nu(i) = nue(i)
        end do

        do sp = 1, plasma%n_species-1
            do i = 1, iprof_length
                plasma%spec(sp)%nu(i) = nui(sp, i)
            end do
        end do

        do sp =0, plasma%n_species-1
            do i = 1,iprof_length
                plasma%spec(sp)%x1(i) = plasma%kp(i) * plasma%spec(sp)%vT(i) / plasma%spec(sp)%nu(i)
                plasma%spec(sp)%x2(i) = - (plasma%om_E(i) - omega) / plasma%spec(sp)%nu(i)
                if (collisions_off .eqv. .true.)then
                    nue(i) = 0.0d0
                    nui(sp+1, i) = 0.0d0
                end if
            end do
        end do

            
        do sp=0, plasma%n_species-1
            call calculate_susc_funcs_profiles(plasma%spec(sp))
        end do

    end subroutine

    subroutine write_plasma_backs(plasma, r_grid)

        use IO_collection, only: write_profile
        use config, only: output_path

        implicit none

        type(plasma_t), intent(in) :: plasma
        real(dp), intent(in) :: r_grid(:)
        logical :: ex

        inquire(file=trim(output_path)//'profiles', exist=ex)
        if (.not. ex) then
            call system('mkdir -p '//trim(output_path)//'backs/')
        end if

        call write_profile(r_grid, plasma%ks, size(r_grid), trim(output_path)//'backs/'//'ks.dat')
        call write_profile(r_grid, plasma%kp, size(r_grid), trim(output_path)//'backs/'//'kp.dat')
        call write_profile(r_grid, plasma%om_E, size(r_grid), trim(output_path)//'backs/'//'om_E.dat')
        call write_profile(r_grid, plasma%q, size(r_grid), trim(output_path)//'backs/'//'q.dat')
        call write_profile(r_grid, plasma%dqdr, size(r_grid), trim(output_path)//'backs/'//'dqdr.dat')
        call write_profile(r_grid, plasma%Er, size(r_grid), trim(output_path)//'backs/'//'Er.dat')

    end subroutine

    subroutine write_species_backs(spec, r_grid)

        use KIM_kinds, only: dp
        use IO_collection, only: plot_1D_labeled, write_profile, remove_file, write_complex_profile
        use config, only: output_path

        implicit none

        type(species_t), intent(in) :: spec
        real(dp), intent(in) :: r_grid(:)
        logical :: ex

        inquire(file=trim(output_path)//'profiles', exist=ex)
        if (.not. ex) then
            call system('mkdir -p '//trim(output_path)//'backs/'//trim(spec%name))
        end if

        call write_profile(r_grid, spec%lambda_D, size(r_grid), trim(output_path)//'backs/'//trim(spec%name)//'/lambda_D.dat')
        call write_profile(r_grid, spec%rho_L, size(r_grid), trim(output_path)//'backs/'//trim(spec%name)//'/rho_L.dat')
        call write_profile(r_grid, spec%vT, size(r_grid), trim(output_path)//'backs/'//trim(spec%name)//'/vT.dat')
        call write_profile(r_grid, spec%omega_c, size(r_grid), trim(output_path)//'backs/'//trim(spec%name)//'/omega_c.dat')
        call write_profile(r_grid, spec%nu, size(r_grid), trim(output_path)//'backs/'//trim(spec%name)//'/nu.dat')
        call write_profile(r_grid, spec%A1, size(r_grid), trim(output_path)//'backs/'//trim(spec%name)//'/A1.dat')
        call write_profile(r_grid, spec%A2, size(r_grid), trim(output_path)//'backs/'//trim(spec%name)//'/A2.dat')
        call write_complex_profile(r_grid, spec%z0, size(r_grid), trim(output_path)//'backs/'//trim(spec%name)//'/z0.dat')

        call write_complex_profile(r_grid, spec%I00, size(r_grid), trim(output_path)//'backs/'//trim(spec%name)//'/I00.dat')
        call write_complex_profile(r_grid, spec%I20, size(r_grid), trim(output_path)//'backs/'//trim(spec%name)//'/I20.dat')
        call write_complex_profile(r_grid, spec%I01, size(r_grid), trim(output_path)//'backs/'//trim(spec%name)//'/I01.dat')
        call write_complex_profile(r_grid, spec%I21, size(r_grid), trim(output_path)//'backs/'//trim(spec%name)//'/I21.dat')
        call write_profile(r_grid, spec%x1, size(r_grid), trim(output_path)//'backs/'//trim(spec%name)//'/x1.dat')
        call write_profile(r_grid, spec%x2, size(r_grid), trim(output_path)//'backs/'//trim(spec%name)//'/x2.dat')

        call write_profile(r_grid, spec%dTdr, size(r_grid), trim(output_path)//'backs/'//trim(spec%name)//'/dTdr.dat')
        call write_profile(r_grid, spec%dndr, size(r_grid), trim(output_path)//'backs/'//trim(spec%name)//'/dndr.dat')
        call write_profile(r_grid, spec%T, size(r_grid), trim(output_path)//'backs/'//trim(spec%name)//'/T.dat')
        call write_profile(r_grid, spec%n, size(r_grid), trim(output_path)//'backs/'//trim(spec%name)//'/n.dat')

    end subroutine

    subroutine interpolate_plasma_backs(plasma_in, grid)

        use KIM_kinds, only: dp
        use IO_collection, only: plot_profile

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
            call allocate_species_fields(plasma_temp%spec(sp), size(grid))
            call allocate_plasma_fields(plasma_temp, size(grid))
            plasma_temp%r_grid = grid
        end do

        do sp = 0, plasma_temp%n_species-1
            do i = 1, size(grid)
                call binsrc(plasma_in%r_grid, 1, size(plasma_in%r_grid), grid(i), ir) 
                ibeg = max(1, ir - nlagr/2)
                iend = ibeg + nlagr - 1
                if (iend .gt. size(plasma_in%r_grid)) then
                    iend = size(plasma_in%r_grid)
                    ibeg = iend -nlagr + 1
                end if

                call plag_coeff(nlagr, nder, grid(i), plasma_in%r_grid(ibeg:iend), coef)

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

                plasma_temp%spec(sp)%x1(i) = sum(coef(0,:) * plasma_in%spec(sp)%x1(ibeg:iend))
                plasma_temp%spec(sp)%x2(i) = sum(coef(0,:) * plasma_in%spec(sp)%x2(ibeg:iend))

                plasma_temp%spec(sp)%I00(i) = sum(coef(0,:) * plasma_in%spec(sp)%I00(ibeg:iend))
                plasma_temp%spec(sp)%I01(i) = sum(coef(0,:) * plasma_in%spec(sp)%I01(ibeg:iend))
                plasma_temp%spec(sp)%I20(i) = sum(coef(0,:) * plasma_in%spec(sp)%I20(ibeg:iend))
                plasma_temp%spec(sp)%I21(i) = sum(coef(0,:) * plasma_in%spec(sp)%I21(ibeg:iend))

                plasma_temp%ks(i) = sum(coef(0,:) * plasma_in%ks(ibeg:iend))
                plasma_temp%kp(i) = sum(coef(0,:) * plasma_in%kp(ibeg:iend))
                plasma_temp%om_E(i) = sum(coef(0,:) * plasma_in%om_E(ibeg:iend))
                plasma_temp%q(i) = sum(coef(0,:) * plasma_in%q(ibeg:iend))
                plasma_temp%dqdr(i) = sum(coef(0,:) * plasma_in%dqdr(ibeg:iend))
                plasma_temp%Er(i) = sum(coef(0,:) * plasma_in%Er(ibeg:iend))

            end do
            
        end do
        
        plasma_in = plasma_temp

        deallocate(plasma_temp)

    end subroutine

    subroutine allocate_species_fields(spec, grid_size)

        implicit none

        integer, intent(in) :: grid_size

        type(species_t), intent(inout) :: spec

        call reallocate(spec%n, grid_size)
        call reallocate(spec%dndr, grid_size)
        call reallocate(spec%T, grid_size)
        call reallocate(spec%dTdr, grid_size)
        call reallocate(spec%A1, grid_size)
        call reallocate(spec%A2, grid_size)
        call reallocate(spec%nu, grid_size)
        call reallocate(spec%vT, grid_size)
        call reallocate(spec%omega_c, grid_size)
        call reallocate(spec%lambda_D, grid_size)
        call reallocate(spec%rho_L, grid_size)
        call reallocate_complex(spec%z0, grid_size)
        call reallocate(spec%x1, grid_size)
        call reallocate(spec%x2, grid_size)

        call reallocate_complex(spec%I00, grid_size)
        call reallocate_complex(spec%I01, grid_size)
        call reallocate_complex(spec%I20, grid_size)
        call reallocate_complex(spec%I21, grid_size)

    end subroutine

    subroutine allocate_plasma_fields(plasma_in, grid_size)

        implicit none

        integer, intent(in) :: grid_size
        type(plasma_t), intent(inout) :: plasma_in

        call reallocate(plasma_in%ks, grid_size)
        call reallocate(plasma_in%kp, grid_size)
        call reallocate(plasma_in%om_E, grid_size)
        call reallocate(plasma_in%Er, grid_size)
        call reallocate(plasma_in%q, grid_size)
        call reallocate(plasma_in%dqdr, grid_size)
        call reallocate(plasma_in%r_grid, grid_size)

    end subroutine

    subroutine reallocate(array, n)

        use KIM_kinds, only: dp

        implicit none

        real(dp), allocatable, intent(inout) :: array(:)
        integer, intent(in) :: n

        if (allocated(array)) deallocate(array)
        allocate(array(n))

    end subroutine

    subroutine reallocate_complex(array, n)

        use KIM_kinds, only: dp

        implicit none

        complex(dp), allocatable, intent(inout) :: array(:)
        integer, intent(in) :: n

        if (allocated(array)) deallocate(array)
        allocate(array(n))

    end subroutine

    subroutine plot_species(spec)

        use IO_collection, only: plot_1D_labeled, write_profile, remove_file
        use grid, only: rg_grid

        implicit none

        type(species_t), intent(in) :: spec

        call write_profile(rg_grid%xb, spec%rho_L, rg_grid%npts_b, 'rho_L.dat')
        call plot_1D_labeled('rho_L.dat', ' r [cm] ', ' rho_L [cm] ', '')
        call remove_file('rho_L.dat')

    end subroutine


    subroutine calculate_susc_funcs_profiles(spec)

        use plasma_parameter, only: iprof_length, r_prof
        use resonances_mod, only: width_res, r_res

        implicit none

        type(species_t), intent(inout) :: spec
        integer :: j

        if (.not. allocated(spec%symbI)) allocate(spec%symbI(0:nmmax, 0:nmmax))
        spec%symbI = 0.0d0
        do j = 1, iprof_length
            if (.false.) then !(r_prof(j) .lt. r_res - 5.d0*width_res .or. r_prof(j) .gt. r_res + 5.d0 * width_res) then
                spec%symbI = 0.0d0
            else
                call getIfunc(spec%x1(j), spec%x2(j), spec%symbI)
            end if
            spec%I00(j) = spec%symbI(0, 0)
            spec%I20(j) = spec%symbI(2, 0)
            spec%I01(j) = spec%symbI(0, 1)
            spec%I21(j) = spec%symbI(2, 1)
        end do

    end subroutine

end module