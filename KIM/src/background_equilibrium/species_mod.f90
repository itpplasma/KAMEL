module species_m

    use KIM_kinds_m, only: dp

    implicit none

    type :: plasma_t
        type(species_t), allocatable :: spec(:)
        integer :: n_species
        integer :: grid_size
        real(dp), allocatable :: om_E(:) ! ExB rotation frequency
        real(dp), allocatable :: ks(:) ! "senkrecht" wavenumber
        real(dp), allocatable :: kp(:) ! parallel wavenumber
        real(dp), allocatable :: B0(:) ! equilibrium magnetic field magnitude
        real(dp), allocatable :: q(:) ! safety factor profile
        real(dp), allocatable :: dqdr(:) ! safety factor gradient
        real(dp), allocatable :: Er(:) ! radial electric field
        real(dp), allocatable :: r_grid(:)
        ! Cell-centered quantities on rg_grid%xc (size rg_grid%npts_c)
        real(dp), allocatable :: ks_cc(:) ! wavenumber at cell centers
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
        !
        ! susceptibility function profiles
        complex(dp), allocatable :: I00(:)
        complex(dp), allocatable :: I11(:)
        complex(dp), allocatable :: I13(:)
        complex(dp), allocatable :: I20(:)
        complex(dp), allocatable :: I02(:)
        complex(dp), allocatable :: I01(:)
        complex(dp), allocatable :: I21(:)
        complex(dp), allocatable :: I22(:)

        ! Cell-centered quantities on rg_grid%xc (size rg_grid%npts_c)
        real(dp), allocatable :: n_cc(:)
        real(dp), allocatable :: dndr_cc(:)
        real(dp), allocatable :: T_cc(:)
        real(dp), allocatable :: dTdr_cc(:)
        real(dp), allocatable :: A1_cc(:)
        real(dp), allocatable :: A2_cc(:)
        real(dp), allocatable :: nu_cc(:)
        real(dp), allocatable :: vT_cc(:)
        real(dp), allocatable :: omega_c_cc(:)
        real(dp), allocatable :: lambda_D_cc(:)
        real(dp), allocatable :: rho_L_cc(:)
        complex(dp), allocatable :: I00_cc(:)
        complex(dp), allocatable :: I01_cc(:)
        complex(dp), allocatable :: I20_cc(:)
        complex(dp), allocatable :: I21_cc(:)
        complex(dp), allocatable :: I22_cc(:)
        complex(dp), allocatable :: I02_cc(:)
    end type

    type(plasma_t) :: plasma

    integer, parameter :: nmmax = 3 ! max order of susceptibility functions to be calculated

    contains

    subroutine init_plasma(plasma_in)

        use config_m, only: read_species_from_namelist, plasma_type

        implicit none

        type(plasma_t), intent(inout) :: plasma_in

        if (read_species_from_namelist .eqv. .true.) then
            call read_species_from_nml(plasma_in)
            call init_electron_species(plasma_in%spec(0))
        else
            if (plasma_type == 'H') then
                call init_hydrogen_plasma(plasma_in)
            else if (plasma_type == 'D') then
                call init_deuterium_plasma(plasma_in)
            else
                print *, "Error: Unknown plasma type. Please choose 'H' or 'D'."
                stop
            end if
            !call init_deuterium_plasma(plasma_in)
        end if

        !call check_quasineutrality(plasma_in)

    end subroutine

    subroutine allocate_plasma

        use config_m, only: number_of_ion_species

        implicit none

        plasma%n_species = number_of_ion_species+1
        allocate(plasma%spec(0:plasma%n_species-1))

    end subroutine

    subroutine read_species_from_nml(plasma_in)

        use KIM_kinds_m, only: dp
        use config_m, only: number_of_ion_species, nml_config_path
        use constants_m, only: p_mass

        implicit none

        type(plasma_t), intent(inout) :: plasma_in
        integer :: ai(16), zi(16)
        integer :: i


        namelist /KIM_species/ ai, zi

        open(unit=77, file=nml_config_path)
        read(unit=77, nml=KIM_species)
        close(unit=77)

        print *,  "Reading species from namelist:"
        print *, "Number of ion species: ", number_of_ion_species
        print *, "Mass numbers: ", ai(1:number_of_ion_species)
        print *, "Charge numbers: ", zi(1:number_of_ion_species)

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

        call init_electron_species(plasma%spec(0))
        call init_deuterium_species(plasma%spec(1))

    end subroutine

    subroutine init_deuterium_species(deut)
    
        use constants_m, only: p_mass

        implicit none

        type(species_t), intent(inout) :: deut

        deut%name = 'i'
        deut%Aspec = 2
        deut%Zspec = 1
        deut%mass = p_mass * deut%Aspec

    end subroutine

    subroutine init_hydrogen_plasma(plasma)

        implicit none

        type(plasma_t), intent(inout) :: plasma

        call init_electron_species(plasma%spec(0))
        call init_hydrogen_species(plasma%spec(1))

    end subroutine

    subroutine init_hydrogen_species(deut)
    
        use constants_m, only: p_mass

        implicit none

        type(species_t), intent(inout) :: deut

        deut%name = 'i'
        deut%Aspec = 1
        deut%Zspec = 1
        deut%mass = p_mass * deut%Aspec

    end subroutine

    subroutine init_electron_species(elec)

        use constants_m, only: e_mass
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


    subroutine set_plasma_quantities(plasma)

        use grid_m, only: rg_grid

        implicit none

        type(plasma_t), intent(inout) :: plasma

        call calc_plasma_parameter_derivs
        
        call calculate_plasma_backs(plasma)

        call interpolate_plasma_backs(plasma, rg_grid%xb)
        call compute_rg_cell_centers(plasma)

        call write_species_backs(plasma%spec(0), plasma%r_grid)
        call write_species_backs(plasma%spec(1), plasma%r_grid)
        call write_plasma_backs(plasma, plasma%r_grid)

    end subroutine

    subroutine compute_rg_cell_centers(plasma_in)

        use grid_m, only: rg_grid

        implicit none

        type(plasma_t), intent(inout) :: plasma_in
        integer :: sp, j

        if (.not. allocated(plasma_in%ks_cc)) allocate(plasma_in%ks_cc(rg_grid%npts_c))

        do sp = 0, plasma_in%n_species-1
            if (.not. allocated(plasma_in%spec(sp)%n_cc)) then
                allocate(plasma_in%spec(sp)%n_cc(rg_grid%npts_c))
                allocate(plasma_in%spec(sp)%dndr_cc(rg_grid%npts_c))
                allocate(plasma_in%spec(sp)%T_cc(rg_grid%npts_c))
                allocate(plasma_in%spec(sp)%dTdr_cc(rg_grid%npts_c))
                allocate(plasma_in%spec(sp)%A1_cc(rg_grid%npts_c))
                allocate(plasma_in%spec(sp)%A2_cc(rg_grid%npts_c))
                allocate(plasma_in%spec(sp)%nu_cc(rg_grid%npts_c))
                allocate(plasma_in%spec(sp)%vT_cc(rg_grid%npts_c))
                allocate(plasma_in%spec(sp)%omega_c_cc(rg_grid%npts_c))
                allocate(plasma_in%spec(sp)%lambda_D_cc(rg_grid%npts_c))
                allocate(plasma_in%spec(sp)%rho_L_cc(rg_grid%npts_c))
                allocate(plasma_in%spec(sp)%I00_cc(rg_grid%npts_c))
                allocate(plasma_in%spec(sp)%I01_cc(rg_grid%npts_c))
                allocate(plasma_in%spec(sp)%I20_cc(rg_grid%npts_c))
                allocate(plasma_in%spec(sp)%I21_cc(rg_grid%npts_c))
                allocate(plasma_in%spec(sp)%I22_cc(rg_grid%npts_c))
                allocate(plasma_in%spec(sp)%I02_cc(rg_grid%npts_c))
            end if
        end do

        do j = 1, rg_grid%npts_c
            plasma_in%ks_cc(j) = 0.5d0 * (plasma_in%ks(j) + plasma_in%ks(j+1))
            do sp = 0, plasma_in%n_species-1
                plasma_in%spec(sp)%n_cc(j)        = 0.5d0 * (plasma_in%spec(sp)%n(j)        + plasma_in%spec(sp)%n(j+1))
                plasma_in%spec(sp)%dndr_cc(j)     = 0.5d0 * (plasma_in%spec(sp)%dndr(j)     + plasma_in%spec(sp)%dndr(j+1))
                plasma_in%spec(sp)%T_cc(j)        = 0.5d0 * (plasma_in%spec(sp)%T(j)        + plasma_in%spec(sp)%T(j+1))
                plasma_in%spec(sp)%dTdr_cc(j)     = 0.5d0 * (plasma_in%spec(sp)%dTdr(j)     + plasma_in%spec(sp)%dTdr(j+1))
                plasma_in%spec(sp)%A1_cc(j)       = 0.5d0 * (plasma_in%spec(sp)%A1(j)       + plasma_in%spec(sp)%A1(j+1))
                plasma_in%spec(sp)%A2_cc(j)       = 0.5d0 * (plasma_in%spec(sp)%A2(j)       + plasma_in%spec(sp)%A2(j+1))
                plasma_in%spec(sp)%nu_cc(j)       = 0.5d0 * (plasma_in%spec(sp)%nu(j)       + plasma_in%spec(sp)%nu(j+1))
                plasma_in%spec(sp)%vT_cc(j)       = 0.5d0 * (plasma_in%spec(sp)%vT(j)       + plasma_in%spec(sp)%vT(j+1))
                plasma_in%spec(sp)%omega_c_cc(j)  = 0.5d0 * (plasma_in%spec(sp)%omega_c(j)  + plasma_in%spec(sp)%omega_c(j+1))
                plasma_in%spec(sp)%lambda_D_cc(j) = 0.5d0 * (plasma_in%spec(sp)%lambda_D(j) + plasma_in%spec(sp)%lambda_D(j+1))
                plasma_in%spec(sp)%rho_L_cc(j)    = 0.5d0 * (plasma_in%spec(sp)%rho_L(j)    + plasma_in%spec(sp)%rho_L(j+1))
                plasma_in%spec(sp)%I00_cc(j)      = 0.5d0 * (plasma_in%spec(sp)%I00(j)      + plasma_in%spec(sp)%I00(j+1))
                plasma_in%spec(sp)%I01_cc(j)      = 0.5d0 * (plasma_in%spec(sp)%I01(j)      + plasma_in%spec(sp)%I01(j+1))
                plasma_in%spec(sp)%I20_cc(j)      = 0.5d0 * (plasma_in%spec(sp)%I20(j)      + plasma_in%spec(sp)%I20(j+1))
                plasma_in%spec(sp)%I21_cc(j)      = 0.5d0 * (plasma_in%spec(sp)%I21(j)      + plasma_in%spec(sp)%I21(j+1))
                plasma_in%spec(sp)%I22_cc(j)      = 0.5d0 * (plasma_in%spec(sp)%I22(j)      + plasma_in%spec(sp)%I22(j+1))
                plasma_in%spec(sp)%I02_cc(j)      = 0.5d0 * (plasma_in%spec(sp)%I02(j)      + plasma_in%spec(sp)%I02(j+1))
            end do
        end do

    end subroutine compute_rg_cell_centers

    subroutine calculate_plasma_backs(plasma_in)

        use constants_m, only: sol, e_charge, ev, pi, com_unit
        use setup_m, only: omega, collisions_off
        use config_m, only: number_of_ion_species

        implicit none

        type(plasma_t), intent(inout) :: plasma_in
        integer :: i, sp, sp_col
        real(dp) :: Lee(plasma%grid_size), Lei(number_of_ion_species, plasma%grid_size), &
            Lii(number_of_ion_species, number_of_ion_species, plasma%grid_size),&
            nue(plasma%grid_size), nui(number_of_ion_species, plasma%grid_size)


        do sp = 0, plasma%n_species-1

            allocate(plasma_in%spec(sp)%rho_L(plasma_in%grid_size))
            allocate(plasma_in%spec(sp)%z0(plasma_in%grid_size))
            allocate(plasma_in%spec(sp)%x1(plasma_in%grid_size))
            allocate(plasma_in%spec(sp)%x2(plasma_in%grid_size))
            allocate(plasma_in%spec(sp)%vT(plasma_in%grid_size))
            allocate(plasma_in%spec(sp)%lambda_D(plasma_in%grid_size))
            allocate(plasma_in%spec(sp)%A1(plasma_in%grid_size))
            allocate(plasma_in%spec(sp)%A2(plasma_in%grid_size))
            allocate(plasma_in%spec(sp)%I00(plasma_in%grid_size))
            allocate(plasma_in%spec(sp)%I20(plasma_in%grid_size))
            allocate(plasma_in%spec(sp)%I01(plasma_in%grid_size))
            allocate(plasma_in%spec(sp)%I02(plasma_in%grid_size))
            allocate(plasma_in%spec(sp)%I21(plasma_in%grid_size))
            allocate(plasma_in%spec(sp)%I22(plasma_in%grid_size))
            allocate(plasma_in%spec(sp)%I11(plasma_in%grid_size))
            allocate(plasma_in%spec(sp)%I13(plasma_in%grid_size))
            allocate(plasma_in%spec(sp)%nu(plasma_in%grid_size))
            allocate(plasma_in%spec(sp)%omega_c(plasma_in%grid_size))

            do i=1, plasma_in%grid_size
                plasma_in%spec(sp)%vT(i) = sqrt(plasma_in%spec(sp)%T(i) * ev / (plasma_in%spec(sp)%mass))
                plasma_in%spec(sp)%omega_c(i) = plasma_in%spec(sp)%Zspec * e_charge * abs(plasma_in%B0(i)) &
                    / (plasma%spec(sp)%mass * sol)

                plasma_in%spec(sp)%rho_L(i) = abs(plasma_in%spec(sp)%vT(i) / (plasma_in%spec(sp)%omega_c(i)))

                plasma_in%spec(sp)%lambda_D(i) = sqrt(plasma_in%spec(sp)%T(i) *ev / (4.0d0*pi* plasma_in%spec(sp)%n(i) &
                    * (plasma_in%spec(sp)%Zspec * e_charge)**2.0d0))

                plasma_in%spec(sp)%A1(i) = plasma_in%spec(sp)%dndr(i) / plasma_in%spec(sp)%n(i) - plasma_in%spec(sp)%Zspec *e_charge&
                    /(plasma_in%spec(sp)%T(i) * ev) * plasma_in%Er(i) - 3.0d0/(2.0d0 * plasma_in%spec(sp)%T(i)) * plasma_in%spec(sp)%dTdr(i)
                plasma_in%spec(sp)%A2(i) = plasma_in%spec(sp)%dTdr(i) / plasma_in%spec(sp)%T(i)

                plasma_in%spec(sp)%z0(i) = - (plasma_in%om_E(i) - omega - com_unit * plasma_in%spec(sp)%nu(i)) &
                    / (abs(plasma_in%kp(i)) * sqrt(2d0) * plasma_in%spec(sp)%vT(i) )
            end do
        end do

        do i=1, plasma_in%grid_size
            ! Coulomb logarithm
            Lee(i) = 23.5d0 - log(sqrt(plasma_in%spec(0)%n(i)) / plasma_in%spec(0)%T(i)**1.25d0) - &
                sqrt(1d-5 + (log(plasma_in%spec(0)%T(i)) -2.0d0)**2.0d0 / 16.0d0)
            nue(i) = 5.8e-6 * plasma_in%spec(0)%n(i) * Lee(i) / plasma_in%spec(0)%T(i)**(1.5d0)
            Lei(:, i) = 24.0d0 - log(sqrt(plasma_in%spec(0)%n(i)) / plasma_in%spec(0)%T(i))
        end do

        do sp=1, number_of_ion_species
            do i=1, plasma_in%grid_size
                ! Coulomb logarithm electrons ions (= ions electrons)
                

                nue(i) = nue(i) + 7.7d-6 * plasma_in%spec(sp)%n(i) * Lei(sp, i) * plasma_in%spec(sp)%Zspec**2 / plasma_in%spec(0)%T(i)**(1.5d0)

                nui(sp, i) = 1.8d-7 * plasma_in%spec(sp)%Aspec**(-1.0/2.0) * plasma_in%spec(sp)%T(i)**(-3.0/2.0) * plasma_in%spec(0)%n(i) * &
                                plasma_in%spec(sp)%Zspec**2 * Lei(sp,i)

                do sp_col=sp, number_of_ion_species
                    ! Coulomb logarithm ions - ions'
                    Lii(sp, sp_col, i) = 23.0d0 - log(plasma_in%spec(sp)%Zspec * plasma_in%spec(sp_col)%Zspec &
                                        * (plasma_in%spec(sp)%Aspec + plasma_in%spec(sp_col)%Aspec)&
                                        / (plasma_in%spec(sp)%T(i) * plasma_in%spec(sp_col)%Aspec + plasma_in%spec(sp_col)%T(i) * plasma_in%spec(sp)%Aspec)&
                                        * (plasma_in%spec(sp)%n(i) * plasma_in%spec(sp)%Zspec**2 / plasma_in%spec(sp)%T(i) &
                                        + plasma_in%spec(sp_col)%n(i) * plasma_in%spec(sp_col)%Zspec**2) / plasma_in%spec(sp_col)%T(i))
                    ! Collision frequency ions - ions'
                    nui(sp, i) = nui(sp, i) + 1.8d-7 * plasma_in%spec(sp_col)%n(i) * plasma_in%spec(sp)%Zspec**2 &
                                    * plasma_in%spec(sp_col)%Zspec**2 * Lii(sp, sp_col, i) * plasma_in%spec(sp)%Aspec**(-1.0/2.0)&
                                    * plasma_in%spec(sp)%T(i)**(-3.0/2.0)
                end do
            end do
        end do

        do i = 1, plasma_in%grid_size
            plasma_in%spec(0)%nu(i) = nue(i)
        end do

        do sp = 1, plasma_in%n_species-1
            do i = 1, plasma_in%grid_size
                plasma_in%spec(sp)%nu(i) = nui(sp, i)
            end do
        end do

        do sp =0, plasma_in%n_species-1
            do i = 1,plasma_in%grid_size
                plasma_in%spec(sp)%x1(i) = plasma_in%kp(i) * plasma_in%spec(sp)%vT(i) / plasma_in%spec(sp)%nu(i)
                plasma_in%spec(sp)%x2(i) = - (plasma_in%om_E(i) - omega) / plasma_in%spec(sp)%nu(i)
                if (collisions_off .eqv. .true.)then
                    plasma_in%spec(sp)%nu(i) = 0.0d0
                end if
            end do
        end do

        do sp=0, plasma_in%n_species-1
            call calculate_susc_funcs_profiles(plasma_in%spec(sp))
        end do

    end subroutine

    subroutine write_plasma_backs(plasma, r_grid)

        use IO_collection_m, only: write_profile
        use config_m, only: output_path

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

        use KIM_kinds_m, only: dp
        use IO_collection_m, only: plot_1D_labeled, write_profile, remove_file, write_complex_profile
        use config_m, only: output_path

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

        use KIM_kinds_m, only: dp
        use IO_collection_m, only: plot_profile

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
                plasma_temp%spec(sp)%I22(i) = sum(coef(0,:) * plasma_in%spec(sp)%I22(ibeg:iend))
                plasma_temp%spec(sp)%I11(i) = sum(coef(0,:) * plasma_in%spec(sp)%I11(ibeg:iend))
                plasma_temp%spec(sp)%I13(i) = sum(coef(0,:) * plasma_in%spec(sp)%I13(ibeg:iend))
                plasma_temp%spec(sp)%I02(i) = sum(coef(0,:) * plasma_in%spec(sp)%I02(ibeg:iend))

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

        use KIM_kinds_m, only: dp

        implicit none

        real(dp), allocatable, intent(inout) :: array(:)
        integer, intent(in) :: n

        if (allocated(array)) deallocate(array)
        allocate(array(n))

    end subroutine

    subroutine reallocate_complex(array, n)

        use KIM_kinds_m, only: dp

        implicit none

        complex(dp), allocatable, intent(inout) :: array(:)
        integer, intent(in) :: n

        if (allocated(array)) deallocate(array)
        allocate(array(n))

    end subroutine

    subroutine plot_species(spec)

        use IO_collection_m, only: plot_1D_labeled, write_profile, remove_file
        use grid_m, only: rg_grid

        implicit none

        type(species_t), intent(in) :: spec

        call write_profile(rg_grid%xb, spec%rho_L, rg_grid%npts_b, 'rho_L.dat')
        call plot_1D_labeled('rho_L.dat', ' r [cm] ', ' rho_L [cm] ', '')
        call remove_file('rho_L.dat')

    end subroutine


    subroutine calculate_susc_funcs_profiles(spec)

        use resonances_mod, only: r_res
        use grid_m, only: width_res

        implicit none

        type(species_t), intent(inout) :: spec
        integer :: j

        if (.not. allocated(spec%symbI)) allocate(spec%symbI(0:nmmax, 0:nmmax))
        spec%symbI = 0.0d0
        do j = 1, plasma%grid_size
            if (.false.) then 
            !if (plasma%r_grid(j) .lt. r_res - 7.d0*width_res .or. plasma%r_grid(j) .gt. r_res + 7.d0 * width_res) then
                spec%symbI = 0.0d0
            else
                call getIfunc(spec%x1(j), spec%x2(j), spec%symbI)
            end if
            spec%I00(j) = spec%symbI(0, 0)
            spec%I20(j) = spec%symbI(2, 0)
            spec%I02(j) = spec%symbI(0, 2)
            spec%I01(j) = spec%symbI(0, 1)
            spec%I21(j) = spec%symbI(2, 1)
            spec%I22(j) = spec%symbI(2, 2)
            spec%I11(j) = spec%symbI(1, 1)
            spec%I13(j) = spec%symbI(1, 3)
        end do

    end subroutine

    subroutine check_quasineutrality(plasma_in)

        use KIM_kinds_m, only: dp
        use config_m, only: number_of_ion_species

        implicit none

        type(plasma_t), intent(in) :: plasma_in
        real(dp) :: n_zero
        logical :: check_succeeded = .true.
        integer :: i, sp

        print *, "Checking quasineutrality..." 
        do i = 1, size(plasma_in%spec(0)%n)
            n_zero = plasma_in%spec(0)%n(i)
            do sp=1, number_of_ion_species
                n_zero = n_zero - plasma_in%spec(sp)%Zspec * plasma_in%spec(sp)%n(i)
            end do
            if (.not. (abs(n_zero) <1.0)) then
                check_succeeded = .false.
                print *, "Warning: quasineutality check failed for r = ", plasma_in%r_grid(i), &
                    " with n_zero = ", n_zero
            end if
        end do

        if (check_succeeded) then
            print *, "Quasineutrality check succeeded."
        else
            print *, "Quasineutrality check failed. Please check your profiles."
        end if

    end subroutine

    subroutine calc_plasma_parameter_derivs

        use grid_m
        use config_m, only: number_of_ion_species

        implicit none

        integer :: i, sigma
        real(dp) :: h1, h2, c1, c2, c3

        ! Allocate derivative arrays if not already allocated
        if (.not. allocated(plasma%spec(0)%dndr)) then
            do i = 0, number_of_ion_species
                allocate(plasma%spec(i)%dndr(plasma%grid_size), &
                        plasma%spec(i)%dTdr(plasma%grid_size))
            end do
            allocate(plasma%dqdr(plasma%grid_size))
        end if

        ! Handle the edge case of too few grid points
        if (plasma%grid_size < 3) then
            ! For very small grids, use simple first-order differences
            if (plasma%grid_size == 1) then
                do sigma = 0, number_of_ion_species
                    plasma%spec(sigma)%dndr(1) = 0.0d0
                    plasma%spec(sigma)%dTdr(1) = 0.0d0
                end do
                plasma%dqdr(1) = 0.0d0
            else if (plasma%grid_size == 2) then
                h1 = plasma%r_grid(2) - plasma%r_grid(1)
                do sigma = 0, number_of_ion_species
                    plasma%spec(sigma)%dndr(1) = (plasma%spec(sigma)%n(2) - plasma%spec(sigma)%n(1)) / h1
                    plasma%spec(sigma)%dTdr(1) = (plasma%spec(sigma)%T(2) - plasma%spec(sigma)%T(1)) / h1
                    plasma%spec(sigma)%dndr(2) = plasma%spec(sigma)%dndr(1)
                    plasma%spec(sigma)%dTdr(2) = plasma%spec(sigma)%dTdr(1)
                end do
                plasma%dqdr(1) = (plasma%q(2) - plasma%q(1)) / h1
                plasma%dqdr(2) = plasma%dqdr(1)
            end if
            return
        end if

        ! LEFT BOUNDARY: Second-order forward differences
        ! f'(x0) = (-3*f(x0) + 4*f(x1) - f(x2)) / (2*h) for uniform grid
        ! For non-uniform grid, we need to account for varying spacing
        i = 1
        h1 = plasma%r_grid(2) - plasma%r_grid(1)
        h2 = plasma%r_grid(3) - plasma%r_grid(2)
        
        ! Coefficients for non-uniform grid
        c1 = -(2.0d0*h1 + h2) / (h1*(h1 + h2))
        c2 = (h1 + h2) / (h1*h2)
        c3 = -h1 / (h2*(h1 + h2))
        
        do sigma = 0, number_of_ion_species
            plasma%spec(sigma)%dndr(i) = c1*plasma%spec(sigma)%n(1) + c2*plasma%spec(sigma)%n(2) + c3*plasma%spec(sigma)%n(3)
            plasma%spec(sigma)%dTdr(i) = c1*plasma%spec(sigma)%T(1) + c2*plasma%spec(sigma)%T(2) + c3*plasma%spec(sigma)%T(3)
        end do
        plasma%dqdr(i) = c1*plasma%q(1) + c2*plasma%q(2) + c3*plasma%q(3)

        ! INTERIOR POINTS: Central differences
        ! f'(xi) = (f(xi+1) - f(xi-1)) / (2*h) for uniform grid
        do i = 2, plasma%grid_size - 1
            h1 = plasma%r_grid(i) - plasma%r_grid(i-1)
            h2 = plasma%r_grid(i+1) - plasma%r_grid(i)
            
            ! Coefficients for non-uniform grid central differences
            c1 = -h2 / (h1*(h1 + h2))
            c2 = (h2 - h1) / (h1*h2)
            c3 = h1 / (h2*(h1 + h2))
            
            do sigma = 0, number_of_ion_species
                plasma%spec(sigma)%dndr(i) = c1*plasma%spec(sigma)%n(i-1) + c2*plasma%spec(sigma)%n(i) + c3*plasma%spec(sigma)%n(i+1)
                plasma%spec(sigma)%dTdr(i) = c1*plasma%spec(sigma)%T(i-1) + c2*plasma%spec(sigma)%T(i) + c3*plasma%spec(sigma)%T(i+1)
            end do
            plasma%dqdr(i) = c1*plasma%q(i-1) + c2*plasma%q(i) + c3*plasma%q(i+1)
        end do

        ! RIGHT BOUNDARY: Second-order backward differences
        ! f'(xn) = (f(xn-2) - 4*f(xn-1) + 3*f(xn)) / (2*h) for uniform grid
        i = plasma%grid_size
        h1 = plasma%r_grid(i-1) - plasma%r_grid(i-2)
        h2 = plasma%r_grid(i) - plasma%r_grid(i-1)
        
        ! Coefficients for non-uniform grid
        c1 = h2 / (h1*(h1 + h2))
        c2 = -(h1 + h2) / (h1*h2)
        c3 = (2.0d0*h2 + h1) / (h2*(h1 + h2))
        
        do sigma = 0, number_of_ion_species
            plasma%spec(sigma)%dndr(i) = c1*plasma%spec(sigma)%n(i-2) + c2*plasma%spec(sigma)%n(i-1) + c3*plasma%spec(sigma)%n(i)
            plasma%spec(sigma)%dTdr(i) = c1*plasma%spec(sigma)%T(i-2) + c2*plasma%spec(sigma)%T(i-1) + c3*plasma%spec(sigma)%T(i)
        end do
        plasma%dqdr(i) = c1*plasma%q(i-2) + c2*plasma%q(i-1) + c3*plasma%q(i)

    end subroutine

    subroutine read_profiles()

        use config_m, only: hdf5_input            
        use grid_m, only: r_space_dim

        if (hdf5_input) then
            ! read plasma profiles from hdf5 file
            call read_from_hdf5
        else
            ! read plasma profiles from text files
            call read_from_text
        endif
            
        r_space_dim = plasma%grid_size

    end subroutine
    
    subroutine read_from_text

        use config_m, only: number_of_ion_species, profile_location, fstatus
        use KIM_kinds_m, only: dp
        use setup_m, only: set_profiles_constant
        use grid_m, only: r_plas

        implicit none

        integer :: i, sigma
        integer :: ierr
        integer :: ios
        integer :: total_Z
        real(dp) :: r_temp

        if (fstatus == 1) write(*,*) 'Status: Reading profiles from text files'

        ! find profile length
        plasma%grid_size = 0

        open(99, file=trim(profile_location)//'n.dat')
        ios = 0
        do while(ios == 0)
            read(99, *, iostat=ios) r_temp
            if (r_temp < r_plas) then
                plasma%grid_size = plasma%grid_size + 1
            else 
                ios = 1
            end if
        end do
        close(99)

        ierr = 0
        if (.not. allocated(plasma%r_grid)) allocate(plasma%r_grid(plasma%grid_size), stat=ierr)
        if (ierr /= 0) print *, "array: Allocation request denied"
        
        do sigma = 0, number_of_ion_species
            allocate(plasma%spec(sigma)%n(plasma%grid_size),&
                plasma%spec(sigma)%T(plasma%grid_size))
        end do
        
        allocate(plasma%q(plasma%grid_size), &
                plasma%Er(plasma%grid_size))

        open(11, file=trim(profile_location)//'n.dat')
        do i=1, plasma%grid_size
            read(11, *) plasma%r_grid(i), plasma%spec(0)%n(i)
        end do
        close(11)

        open(11, file=trim(profile_location)//'Te.dat')
        do i=1, plasma%grid_size
            read(11, *) r_temp, plasma%spec(0)%T(i)
        end do
        close(11)

        open(11, file=trim(profile_location)//'Ti.dat')
        do i=1, plasma%grid_size
            read(11, *) r_temp, plasma%spec(1)%T(i)
            do sigma = 2, number_of_ion_species
                plasma%spec(sigma)%T(i) = plasma%spec(1)%T(i)
            end do
        end do
        close(11)

        open(11, file=trim(profile_location)//'Er.dat')
        do i=1, plasma%grid_size
            read(11, *) r_temp, plasma%Er(i)
        end do
        close(11)

        open(11, file=trim(profile_location)//'q.dat')
        do i=1, plasma%grid_size
            read(11, *) r_temp, plasma%q(i)
        end do
        close(11)

        total_Z = 0
        do sigma = 1, number_of_ion_species
            total_Z = total_Z + plasma%spec(sigma)%Zspec
        end do

        do i = 1, plasma%grid_size
            do sigma = 1, number_of_ion_species
                ! ion density to fulfill quasineutrality
                plasma%spec(sigma)%n(i) = plasma%spec(0)%n(i) * plasma%spec(sigma)%Zspec / total_Z
            end do
        end do

        if (set_profiles_constant == 1) then
            write(*,*) 'Info: Setting profiles to constant values'
            plasma%Er(:) = plasma%Er(1)
            do sigma = 0, number_of_ion_species
                plasma%spec(sigma)%n(:) = plasma%spec(sigma)%n(1)
                plasma%spec(sigma)%T(:) = plasma%spec(sigma)%T(1)
            end do
        end if

        if (set_profiles_constant == 2) then
            write(*,*) 'Info: Setting profiles to constant values'
            do sigma = 0, number_of_ion_species
                plasma%spec(sigma)%T(:) = plasma%spec(sigma)%T(1)
            end do
        end if

        if (fstatus == 1) write(*,*) 'Status: Finished reading profiles from text files'

    end subroutine


    subroutine read_from_hdf5

        use config_m, only: fstatus

        implicit none

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

    subroutine write_profiles

        use config_m, only: output_path, fstatus, number_of_ion_species
        use IO_collection_m, only: write_profile

        implicit none

        integer :: sigma
        logical :: ex

        if (fstatus == 1) write(*,*) 'Status: writing profiles to output_path'

        inquire(file=trim(output_path)//'profiles', exist=ex)
        if (.not. ex) then
            call system('mkdir -p '//trim(output_path)//'profiles')
        end if

        do sigma = 0, number_of_ion_species
            call write_profile(plasma%r_grid, plasma%spec(sigma)%n, plasma%grid_size, trim(output_path)//'profiles/n.dat')
            call write_profile(plasma%r_grid, plasma%spec(sigma)%T, plasma%grid_size, trim(output_path)//'profiles/T.dat')
        end do
        call write_profile(plasma%r_grid, plasma%Er, plasma%grid_size, trim(output_path)//'profiles/Er.dat')
        call write_profile(plasma%r_grid, plasma%q, plasma%grid_size, trim(output_path)//'profiles/q.dat')

    end subroutine


end module
