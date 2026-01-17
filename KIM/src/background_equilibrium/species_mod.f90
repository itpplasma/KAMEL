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
        real(dp), allocatable :: Er_cc(:) ! radial electric field at cell centers
        real(dp), allocatable :: om_E_cc(:) ! ExB rotation frequency at cell centers
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
        real(dp), allocatable :: x2(:, :) ! normalized inverse collisionality, second index for cyclotron harmonic
        complex(dp), allocatable :: z0(:) ! = x2/(sqrt(2) x1), argument of plasma dispersion function
        complex(dp), dimension(:,:), allocatable :: symbI ! susceptibility function used for Fokker-Planck collision model
        !
        ! susceptibility function profiles
        complex(dp), allocatable :: I00(:, :)
        complex(dp), allocatable :: I11(:, :)
        complex(dp), allocatable :: I13(:, :)
        complex(dp), allocatable :: I20(:, :)
        complex(dp), allocatable :: I02(:, :)
        complex(dp), allocatable :: I01(:, :)
        complex(dp), allocatable :: I21(:, :)
        complex(dp), allocatable :: I22(:, :)

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
        real(dp), allocatable :: x1_cc(:)
        real(dp), allocatable :: x2_cc(:, :)
        complex(dp), allocatable :: I00_cc(:, :)
        complex(dp), allocatable :: I01_cc(:, :)
        complex(dp), allocatable :: I10_cc(:, :)
        complex(dp), allocatable :: I20_cc(:, :)
        complex(dp), allocatable :: I12_cc(:, :)
        complex(dp), allocatable :: I21_cc(:, :)
        complex(dp), allocatable :: I22_cc(:, :)
        complex(dp), allocatable :: I02_cc(:, :)
        complex(dp), allocatable :: I13_cc(:, :)
        complex(dp), allocatable :: I11_cc(:, :)
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

    subroutine init_hydrogen_species(hydro)
    
        use constants_m, only: p_mass

        implicit none

        type(species_t), intent(inout) :: hydro

        hydro%name = 'i'
        hydro%Aspec = 1
        hydro%Zspec = 1
        hydro%mass = p_mass * hydro%Aspec

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

        integer :: sp

        call calc_plasma_parameter_derivs

        call calculate_plasma_backs(plasma)

        call interpolate_plasma_backs(plasma, rg_grid%xb)
        call compute_rg_cell_centers(plasma)

        call calculate_thermodynamic_forces_and_susc(plasma)

        do sp = 0, plasma%n_species-1
            call write_species_backs(plasma%spec(sp), plasma%r_grid)
            call write_species_cc_quantities(plasma%spec(sp), rg_grid%xc)
        end do

        call write_plasma_backs(plasma, plasma%r_grid)

    end subroutine

    subroutine compute_rg_cell_centers(plasma_in)
        ! Compute cell-centered values of all plasma background profiles on rg_grid%xc
        ! Simple averaging of boundary values for smooth quantities

        use grid_m, only: rg_grid

        implicit none

        type(plasma_t), intent(inout) :: plasma_in
        integer :: sp, j

        if (.not. allocated(plasma_in%ks_cc)) allocate(plasma_in%ks_cc(rg_grid%npts_c))
        if (.not. allocated(plasma_in%Er_cc)) allocate(plasma_in%Er_cc(rg_grid%npts_c))
        if (.not. allocated(plasma_in%om_E_cc)) allocate(plasma_in%om_E_cc(rg_grid%npts_c))

        do sp = 0, plasma_in%n_species-1
            if (.not. allocated(plasma_in%spec(sp)%n_cc)) then
                allocate(plasma_in%spec(sp)%n_cc(rg_grid%npts_c))
                allocate(plasma_in%spec(sp)%dndr_cc(rg_grid%npts_c))
                allocate(plasma_in%spec(sp)%T_cc(rg_grid%npts_c))
                allocate(plasma_in%spec(sp)%dTdr_cc(rg_grid%npts_c))
                allocate(plasma_in%spec(sp)%nu_cc(rg_grid%npts_c))
                allocate(plasma_in%spec(sp)%vT_cc(rg_grid%npts_c))
                allocate(plasma_in%spec(sp)%omega_c_cc(rg_grid%npts_c))
                allocate(plasma_in%spec(sp)%lambda_D_cc(rg_grid%npts_c))
                allocate(plasma_in%spec(sp)%rho_L_cc(rg_grid%npts_c))
            end if
        end do

        do j = 1, rg_grid%npts_c
            plasma_in%ks_cc(j) = 0.5d0 * (plasma_in%ks(j) + plasma_in%ks(j+1))
            plasma_in%Er_cc(j) = 0.5d0 * (plasma_in%Er(j) + plasma_in%Er(j+1))
            plasma_in%om_E_cc(j) = 0.5d0 * (plasma_in%om_E(j) + plasma_in%om_E(j+1))
            do sp = 0, plasma_in%n_species-1
                plasma_in%spec(sp)%n_cc(j)        = 0.5d0 * (plasma_in%spec(sp)%n(j)        + plasma_in%spec(sp)%n(j+1))
                plasma_in%spec(sp)%dndr_cc(j)     = 0.5d0 * (plasma_in%spec(sp)%dndr(j)     + plasma_in%spec(sp)%dndr(j+1))
                plasma_in%spec(sp)%T_cc(j)        = 0.5d0 * (plasma_in%spec(sp)%T(j)        + plasma_in%spec(sp)%T(j+1))
                plasma_in%spec(sp)%dTdr_cc(j)     = 0.5d0 * (plasma_in%spec(sp)%dTdr(j)     + plasma_in%spec(sp)%dTdr(j+1))
                plasma_in%spec(sp)%nu_cc(j)       = 0.5d0 * (plasma_in%spec(sp)%nu(j)       + plasma_in%spec(sp)%nu(j+1))
                plasma_in%spec(sp)%vT_cc(j)       = 0.5d0 * (plasma_in%spec(sp)%vT(j)       + plasma_in%spec(sp)%vT(j+1))
                plasma_in%spec(sp)%omega_c_cc(j)  = 0.5d0 * (plasma_in%spec(sp)%omega_c(j)  + plasma_in%spec(sp)%omega_c(j+1))
                plasma_in%spec(sp)%lambda_D_cc(j) = 0.5d0 * (plasma_in%spec(sp)%lambda_D(j) + plasma_in%spec(sp)%lambda_D(j+1))
                plasma_in%spec(sp)%rho_L_cc(j)    = 0.5d0 * (plasma_in%spec(sp)%rho_L(j)    + plasma_in%spec(sp)%rho_L(j+1))
            end do
        end do

    end subroutine compute_rg_cell_centers

    subroutine calculate_thermodynamic_forces_and_susc(plasma_in)
        ! Calculate x1, x2, A1, A2, and susceptibility functions on rg_grid
        ! This is called AFTER interpolation to rg_grid to avoid grid-dependent aliasing
        ! Calculates BOTH boundary values (for FLR2 asymptotics) AND cell-center values (for kernels)

        use constants_m, only: e_charge, ev
        use setup_m, only: omega, mphi_max
        use grid_m, only: rg_grid
        use KIM_kinds_m, only: dp
        use config_m, only: ion_flr_scale_factor

        implicit none

        type(plasma_t), intent(inout) :: plasma_in
        integer :: sp, j, mphi

        ! Allocate arrays
        do sp = 0, plasma_in%n_species-1
            if (.not. allocated(plasma_in%spec(sp)%symbI)) allocate(plasma_in%spec(sp)%symbI(0:nmmax, 0:nmmax))

            ! Boundary arrays (for FLR2 asymptotics and Krook model)
            if (.not. allocated(plasma_in%spec(sp)%A1)) then
                allocate(plasma_in%spec(sp)%A1(rg_grid%npts_b))
                allocate(plasma_in%spec(sp)%A2(rg_grid%npts_b))
                allocate(plasma_in%spec(sp)%x1(rg_grid%npts_b))
                allocate(plasma_in%spec(sp)%x2(rg_grid%npts_b, -mphi_max:mphi_max))
                allocate(plasma_in%spec(sp)%I00(rg_grid%npts_b, -mphi_max:mphi_max))
                allocate(plasma_in%spec(sp)%I01(rg_grid%npts_b, -mphi_max:mphi_max))
                allocate(plasma_in%spec(sp)%I20(rg_grid%npts_b, -mphi_max:mphi_max))
                allocate(plasma_in%spec(sp)%I21(rg_grid%npts_b, -mphi_max:mphi_max))
                allocate(plasma_in%spec(sp)%I22(rg_grid%npts_b, -mphi_max:mphi_max))
                allocate(plasma_in%spec(sp)%I02(rg_grid%npts_b, -mphi_max:mphi_max))
                allocate(plasma_in%spec(sp)%I11(rg_grid%npts_b, -mphi_max:mphi_max))
                allocate(plasma_in%spec(sp)%I13(rg_grid%npts_b, -mphi_max:mphi_max))
            end if

            ! Cell-center arrays (for FP kernels - computed on cell centers to avoid aliasing)
            if (.not. allocated(plasma_in%spec(sp)%A1_cc)) then
                allocate(plasma_in%spec(sp)%A1_cc(rg_grid%npts_c))
                allocate(plasma_in%spec(sp)%A2_cc(rg_grid%npts_c))
                allocate(plasma_in%spec(sp)%x1_cc(rg_grid%npts_c))
                allocate(plasma_in%spec(sp)%x2_cc(rg_grid%npts_c, -mphi_max:mphi_max))
                allocate(plasma_in%spec(sp)%I00_cc(rg_grid%npts_c, -mphi_max:mphi_max))
                allocate(plasma_in%spec(sp)%I01_cc(rg_grid%npts_c, -mphi_max:mphi_max))
                allocate(plasma_in%spec(sp)%I10_cc(rg_grid%npts_c, -mphi_max:mphi_max))
                allocate(plasma_in%spec(sp)%I20_cc(rg_grid%npts_c, -mphi_max:mphi_max))
                allocate(plasma_in%spec(sp)%I21_cc(rg_grid%npts_c, -mphi_max:mphi_max))
                allocate(plasma_in%spec(sp)%I12_cc(rg_grid%npts_c, -mphi_max:mphi_max))
                allocate(plasma_in%spec(sp)%I22_cc(rg_grid%npts_c, -mphi_max:mphi_max))
                allocate(plasma_in%spec(sp)%I02_cc(rg_grid%npts_c, -mphi_max:mphi_max))
                allocate(plasma_in%spec(sp)%I13_cc(rg_grid%npts_c, -mphi_max:mphi_max))
                allocate(plasma_in%spec(sp)%I11_cc(rg_grid%npts_c, -mphi_max:mphi_max))
            end if
        end do

        ! Calculate on boundary points (npts_b)
        do sp = 0, plasma_in%n_species-1
            plasma_in%spec(sp)%symbI = 0.0d0

            do j = 1, rg_grid%npts_b

                plasma_in%spec(sp)%A1(j) = plasma_in%spec(sp)%dndr(j) / plasma_in%spec(sp)%n(j) &
                    - plasma_in%spec(sp)%Zspec * e_charge / (plasma_in%spec(sp)%T(j) * ev) * plasma_in%Er(j) &
                    - 3.0d0 / (2.0d0 * plasma_in%spec(sp)%T(j)) * plasma_in%spec(sp)%dTdr(j)
                plasma_in%spec(sp)%A2(j) = plasma_in%spec(sp)%dTdr(j) / plasma_in%spec(sp)%T(j)


                plasma_in%spec(sp)%x1(j) = plasma_in%kp(j) * plasma_in%spec(sp)%vT(j) / plasma_in%spec(sp)%nu(j)
                do mphi = -mphi_max, mphi_max
                    plasma_in%spec(sp)%x2(j, mphi) = - (plasma_in%om_E(j) & !* ion_flr_scale_factor & !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!1
                                                    + mphi * plasma%spec(sp)%omega_c(j)  - omega) &
                                                    / plasma_in%spec(sp)%nu(j)

                    call getIfunc(plasma_in%spec(sp)%x1(j), plasma_in%spec(sp)%x2(j, mphi), plasma_in%spec(sp)%symbI)
                    plasma_in%spec(sp)%I00(j, mphi) = plasma_in%spec(sp)%symbI(0, 0)
                    plasma_in%spec(sp)%I20(j, mphi) = plasma_in%spec(sp)%symbI(2, 0)
                    plasma_in%spec(sp)%I02(j, mphi) = plasma_in%spec(sp)%symbI(0, 2)
                    plasma_in%spec(sp)%I01(j, mphi) = plasma_in%spec(sp)%symbI(0, 1)
                    plasma_in%spec(sp)%I21(j, mphi) = plasma_in%spec(sp)%symbI(2, 1)
                    plasma_in%spec(sp)%I22(j, mphi) = plasma_in%spec(sp)%symbI(2, 2)
                    plasma_in%spec(sp)%I11(j, mphi) = plasma_in%spec(sp)%symbI(1, 1)
                    plasma_in%spec(sp)%I13(j, mphi) = plasma_in%spec(sp)%symbI(1, 3)
                end do
            end do
        end do

        ! Calculate on cell centers (npts_c) 
        do sp = 0, plasma_in%n_species-1
            plasma_in%spec(sp)%symbI = 0.0d0

            do j = 1, rg_grid%npts_c
                plasma_in%spec(sp)%A1_cc(j) = plasma_in%spec(sp)%dndr_cc(j) / plasma_in%spec(sp)%n_cc(j) &
                    - plasma_in%spec(sp)%Zspec * e_charge / (plasma_in%spec(sp)%T_cc(j) * ev) * plasma_in%Er_cc(j) &
                    - 3.0d0 / (2.0d0 * plasma_in%spec(sp)%T_cc(j)) * plasma_in%spec(sp)%dTdr_cc(j)
                plasma_in%spec(sp)%A2_cc(j) = plasma_in%spec(sp)%dTdr_cc(j) / plasma_in%spec(sp)%T_cc(j)

                plasma_in%spec(sp)%x1_cc(j) = 0.5d0 * (plasma_in%kp(j) + plasma_in%kp(j+1)) &
                    * plasma_in%spec(sp)%vT_cc(j) / plasma_in%spec(sp)%nu_cc(j)

                do mphi = -mphi_max, mphi_max
                    plasma_in%spec(sp)%x2_cc(j, mphi) = - (plasma_in%om_E_cc(j)  & !* ion_flr_scale_factor & !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
                                                        + mphi * plasma_in%spec(sp)%omega_c(j) - omega) &
                                                        / plasma_in%spec(sp)%nu_cc(j)

                    call getIfunc(plasma_in%spec(sp)%x1_cc(j), plasma_in%spec(sp)%x2_cc(j, mphi), plasma_in%spec(sp)%symbI)
                    plasma_in%spec(sp)%I00_cc(j, mphi) = plasma_in%spec(sp)%symbI(0, 0)
                    plasma_in%spec(sp)%I20_cc(j, mphi) = plasma_in%spec(sp)%symbI(2, 0)
                    plasma_in%spec(sp)%I10_cc(j, mphi) = plasma_in%spec(sp)%symbI(1, 0)
                    plasma_in%spec(sp)%I12_cc(j, mphi) = plasma_in%spec(sp)%symbI(1, 2)
                    plasma_in%spec(sp)%I02_cc(j, mphi) = plasma_in%spec(sp)%symbI(0, 2)
                    plasma_in%spec(sp)%I01_cc(j, mphi) = plasma_in%spec(sp)%symbI(0, 1)
                    plasma_in%spec(sp)%I21_cc(j, mphi) = plasma_in%spec(sp)%symbI(2, 1)
                    plasma_in%spec(sp)%I22_cc(j, mphi) = plasma_in%spec(sp)%symbI(2, 2)
                    plasma_in%spec(sp)%I11_cc(j, mphi) = plasma_in%spec(sp)%symbI(1, 1)
                    plasma_in%spec(sp)%I13_cc(j, mphi) = plasma_in%spec(sp)%symbI(1, 3)
                end do
            end do
        end do

    end subroutine calculate_thermodynamic_forces_and_susc

    subroutine calculate_plasma_backs(plasma_in)
        ! Calculate basic parameters of the plasma: thermal velocity, collision frequency,
        ! Larmor radius, Debye length, z0. Does NOT calculate x1, x2, A1, A2, or susceptibility
        ! functions - those are calculated AFTER interpolation to avoid grid-dependent aliasing.

        use constants_m, only: sol, e_charge, ev, pi, com_unit
        use setup_m, only: omega, collisions_off
        use config_m, only: number_of_ion_species, rescale_density, number_density_rescale, ion_flr_scale_factor

        implicit none

        type(plasma_t), intent(inout) :: plasma_in
        integer :: i, sp, sp_col
        real(dp) :: Lee(plasma%grid_size)
        real(dp) :: Lei(number_of_ion_species, plasma%grid_size)
        real(dp) :: Lii(number_of_ion_species, number_of_ion_species, plasma%grid_size)
        real(dp) :: nue(plasma%grid_size)
        real(dp) :: nui(number_of_ion_species, plasma%grid_size)


        do sp = 0, plasma%n_species-1

            allocate(plasma_in%spec(sp)%rho_L(plasma_in%grid_size))
            allocate(plasma_in%spec(sp)%z0(plasma_in%grid_size))
            allocate(plasma_in%spec(sp)%vT(plasma_in%grid_size))
            allocate(plasma_in%spec(sp)%lambda_D(plasma_in%grid_size))
            allocate(plasma_in%spec(sp)%nu(plasma_in%grid_size))
            allocate(plasma_in%spec(sp)%omega_c(plasma_in%grid_size))

            do i=1, plasma_in%grid_size
                plasma_in%spec(sp)%vT(i) = sqrt(plasma_in%spec(sp)%T(i) * ev / (plasma_in%spec(sp)%mass))
                plasma_in%spec(sp)%omega_c(i) = plasma_in%spec(sp)%Zspec * e_charge * abs(plasma_in%B0(i)) &
                    / (plasma%spec(sp)%mass * sol)

                if (sp > 0)then
                    plasma_in%spec(sp)%rho_L(i) = abs(plasma_in%spec(sp)%vT(i) / (plasma_in%spec(sp)%omega_c(i))) * ion_flr_scale_factor
                else
                    plasma_in%spec(sp)%rho_L(i) = abs(plasma_in%spec(sp)%vT(i) / (plasma_in%spec(sp)%omega_c(i)))
                end if

                plasma_in%spec(sp)%lambda_D(i) = sqrt((plasma_in%spec(sp)%T(i) * ev) / (4.0d0*pi* plasma_in%spec(sp)%n(i) &
                    * (plasma_in%spec(sp)%Zspec * e_charge)**2.0d0))
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
                plasma_in%spec(sp)%z0(i) = - (plasma_in%om_E(i) - omega - com_unit * plasma_in%spec(sp)%nu(i)) &
                    / (abs(plasma_in%kp(i)) * sqrt(2d0) * plasma_in%spec(sp)%vT(i) )
                if (collisions_off .eqv. .true.)then
                    plasma_in%spec(sp)%nu(i) = 0.0d0
                end if
            end do
        end do

        if (rescale_density .eqv. .true.) then
            print *, " "
            print *, " ! ! ! "
            print *, " ! ! ! Rescaling density (after collision frequency calculation) ! ! !"
            print *, " ! ! ! "
            print *, " "
            do sp = 0, plasma_in%n_species-1
                do i = 1, plasma_in%grid_size
                    ! density rescaling only affects lambda_D (in A_1 the rescaling cancels out)
                    plasma_in%spec(sp)%n(i) = plasma_in%spec(sp)%n(i) * number_density_rescale
                    plasma_in%spec(sp)%dndr(i) = plasma_in%spec(sp)%dndr(i) * number_density_rescale
                    plasma_in%spec(sp)%lambda_D(i) = sqrt(plasma_in%spec(sp)%T(i) * ev / (4.0d0 * pi * plasma_in%spec(sp)%n(i) &
                        * (plasma_in%spec(sp)%Zspec * e_charge)**2.0d0))
                end do
            end do
        end if

    end subroutine

    subroutine write_plasma_backs(plasma, r_grid)

        use IO_collection_m, only: write_profile
        use config_m, only: output_path, hdf5_output

        implicit none

        type(plasma_t), intent(in) :: plasma
        real(dp), intent(in) :: r_grid(:)
        logical :: ex

        if (.not. hdf5_output) then
            inquire(file=trim(output_path)//'profiles', exist=ex)
            if (.not. ex) then
                call system('mkdir -p '//trim(output_path)//'backs/')
            end if
        end if

        call write_profile(r_grid, plasma%ks, size(r_grid), 'backs/ks', &
            'Perpendicular ("senkrecht") wave number', '1/cm')
        call write_profile(r_grid, plasma%kp, size(r_grid), 'backs/kp', &
            'Parallel wave number', '1/cm')
        call write_profile(r_grid, plasma%om_E, size(r_grid), 'backs/om_E', &
            'E x B drift frequency', 'rad/s')
        call write_profile(r_grid, plasma%q, size(r_grid), 'backs/q', &
            'Safety factor', '1')
        call write_profile(r_grid, plasma%dqdr, size(r_grid), 'backs/dqdr', &
            'Radial derivative of safety factor', '1/cm')
        call write_profile(r_grid, plasma%Er, size(r_grid), 'backs/E0r', &
            'Equilibrium radial electric field', 'statV/cm')

    end subroutine

    subroutine write_species_backs(spec, r_grid)

        use KIM_kinds_m, only: dp
        use IO_collection_m, only: plot_1D_labeled, write_profile, remove_file, write_complex_profile, h5id
        use config_m, only: output_path, hdf5_output
        use KAMEL_hdf5_tools, only: h5_define_group, h5_obj_exists, HID_T

        implicit none

        type(species_t), intent(in) :: spec
        real(dp), intent(in) :: r_grid(:)
        logical :: ex
        integer(HID_T) :: h5grpid

        if (hdf5_output) then
            call h5_obj_exists(h5id, 'backs/'//trim(spec%name), ex)
            if (.not. ex) then
                call h5_define_group(h5id, 'backs/'//trim(spec%name), h5grpid)
            end if

            call write_profile(r_grid, r_grid, size(r_grid), 'backs/'//trim(spec%name)//'/r', &
                'Effective radius (rg space)', 'cm')
        else
            inquire(file=trim(output_path)//'backs'//trim(spec%name), exist=ex)
            if (.not. ex) then
                call system('mkdir -p '//trim(output_path)//'backs/'//trim(spec%name))
            end if
        end if

        call write_profile(r_grid, spec%lambda_D, size(r_grid), 'backs/'//trim(spec%name)//'/lambda_D', &
            'Debye length', 'cm')
        call write_profile(r_grid, spec%rho_L, size(r_grid), 'backs/'//trim(spec%name)//'/rho_L', &
            'Larmor radius', 'cm')
        call write_profile(r_grid, spec%vT, size(r_grid), 'backs/'//trim(spec%name)//'/vT', &
            'Thermal velocity (w/o factor 2)', 'cm/s')
        call write_profile(r_grid, spec%omega_c, size(r_grid), 'backs/'//trim(spec%name)//'/omega_c', &
            'Cyclotron frequency', 'rad/s')
        call write_profile(r_grid, spec%nu, size(r_grid), 'backs/'//trim(spec%name)//'/nu', &
            '(Perpendicular) collision frequency', '1/s')
        call write_complex_profile(r_grid, spec%z0, size(r_grid), 'backs/'//trim(spec%name)//'/z0', &
            'Normalized complex frequency z0 (= x2/(sqrt(2)x1))', '1')

        call write_profile(r_grid, spec%dTdr, size(r_grid), 'backs/'//trim(spec%name)//'/dTdr', &
            'Temperature gradient dT/dr', 'eV/cm')
        call write_profile(r_grid, spec%dndr, size(r_grid), 'backs/'//trim(spec%name)//'/dndr', &
            'Particle density gradient dn/dr', '1/cm^4')
        call write_profile(r_grid, spec%T, size(r_grid), 'backs/'//trim(spec%name)//'/T', &
            'Temperature T', 'eV')
        call write_profile(r_grid, spec%n, size(r_grid), 'backs/'//trim(spec%name)//'/n', &
            'Particle density n', '1/cm^3')

    end subroutine

    subroutine write_species_cc_quantities(spec, r_grid_cc)
        ! Write cell-center quantities (A1, A2, x1, x2, I00, I20, etc.)
        ! These are computed after interpolation to avoid grid-dependent aliasing

        use KIM_kinds_m, only: dp
        use IO_collection_m, only: write_profile, write_complex_profile, h5id, itoa
        use config_m, only: output_path, hdf5_output
        use KAMEL_hdf5_tools, only: h5_define_group, h5_obj_exists, HID_T
        use setup_m, only: mphi_max

        implicit none

        type(species_t), intent(in) :: spec
        real(dp), intent(in) :: r_grid_cc(:)
        logical :: ex
        integer(HID_T) :: h5grpid
        integer :: mphi

        if (hdf5_output) then
            call h5_obj_exists(h5id, 'backs/'//trim(spec%name), ex)
            if (.not. ex) then
                call h5_define_group(h5id, 'fields/'//trim(spec%name), h5grpid)
            end if

            call write_profile(r_grid_cc, r_grid_cc, size(r_grid_cc), 'backs/'//trim(spec%name)//'/r_c', &
                'Effective radius of cell centers (rg space)', 'cm')
        else
            inquire(file=trim(output_path)//'profiles', exist=ex)
            if (.not. ex) then
                call system('mkdir -p '//trim(output_path)//'backs/'//trim(spec%name))
            end if
        end if

        if (allocated(spec%A1_cc)) then
            call write_profile(r_grid_cc, spec%A1_cc, size(r_grid_cc), &
                'backs/'//trim(spec%name)//'/A1_cc', &
                'Thermodynamic force A1 at cell centers', '1/cm')
            call write_profile(r_grid_cc, spec%A2_cc, size(r_grid_cc), &
                'backs/'//trim(spec%name)//'/A2_cc', &
                'Thermodynamic force A2 at cell centers', '1/cm')
            call write_profile(r_grid_cc, spec%x1_cc, size(r_grid_cc), &
                'backs/'//trim(spec%name)//'/x1_cc', &
                'Normalized distance to resonance x1 at cell centers', '1')

            do mphi = -mphi_max, mphi_max
                call write_profile(r_grid_cc, spec%x2_cc(:, mphi), size(r_grid_cc), &
                    'backs/'//trim(spec%name)//'/x2_cc_mphi_'//trim(adjustl(itoa(mphi))), &
                    'Normalized collision frequency x2 at cell centers, mphi='//trim(adjustl(itoa(mphi))), '1')

                call write_complex_profile(r_grid_cc, spec%I00_cc(:, mphi), size(r_grid_cc), &
                    'backs/'//trim(spec%name)//'/I00_cc_mphi_'//trim(adjustl(itoa(mphi))), &
                    'Susceptibility function I00 at cell centers, mphi='//trim(adjustl(itoa(mphi))), '1')
                call write_complex_profile(r_grid_cc, spec%I20_cc(:, mphi), size(r_grid_cc), &
                    'backs/'//trim(spec%name)//'/I20_cc_mphi_'//trim(adjustl(itoa(mphi))), &
                    'Susceptibility function I20 at cell centers, mphi='//trim(adjustl(itoa(mphi))), '1')
                call write_complex_profile(r_grid_cc, spec%I01_cc(:, mphi), size(r_grid_cc), &
                    'backs/'//trim(spec%name)//'/I01_cc_mphi_'//trim(adjustl(itoa(mphi))), &
                    'Susceptibility function I01 at cell centers, mphi='//trim(adjustl(itoa(mphi))), '1')
                call write_complex_profile(r_grid_cc, spec%I21_cc(:, mphi), size(r_grid_cc), &
                    'backs/'//trim(spec%name)//'/I21_cc_mphi_'//trim(adjustl(itoa(mphi))), &
                    'Susceptibility function I21 at cell centers, mphi='//trim(adjustl(itoa(mphi))), '1')
                call write_complex_profile(r_grid_cc, spec%I22_cc(:, mphi), size(r_grid_cc), &
                    'backs/'//trim(spec%name)//'/I22_cc_mphi_'//trim(adjustl(itoa(mphi))), &
                    'Susceptibility function I22 at cell centers, mphi='//trim(adjustl(itoa(mphi))), '1')
                call write_complex_profile(r_grid_cc, spec%I02_cc(:, mphi), size(r_grid_cc), &
                    'backs/'//trim(spec%name)//'/I02_cc_mphi_'//trim(adjustl(itoa(mphi))), &
                    'Susceptibility function I02 at cell centers, mphi='//trim(adjustl(itoa(mphi))), '1')
                call write_complex_profile(r_grid_cc, spec%I13_cc(:, mphi), size(r_grid_cc), &
                    'backs/'//trim(spec%name)//'/I13_cc_mphi_'//trim(adjustl(itoa(mphi))), &
                    'Susceptibility function I13 at cell centers, mphi='//trim(adjustl(itoa(mphi))), '1')
                call write_complex_profile(r_grid_cc, spec%I11_cc(:, mphi), size(r_grid_cc), &
                    'backs/'//trim(spec%name)//'/I11_cc_mphi_'//trim(adjustl(itoa(mphi))), &
                    'Susceptibility function I11 at cell centers, mphi='//trim(adjustl(itoa(mphi))), '1')
                call write_complex_profile(r_grid_cc, spec%I10_cc(:, mphi), size(r_grid_cc), &
                    'backs/'//trim(spec%name)//'/I10_cc_mphi_'//trim(adjustl(itoa(mphi))), &
                    'Susceptibility function I10 at cell centers, mphi='//trim(adjustl(itoa(mphi))), '1')
                call write_complex_profile(r_grid_cc, spec%I12_cc(:, mphi), size(r_grid_cc), &
                    'backs/'//trim(spec%name)//'/I12_cc_mphi_'//trim(adjustl(itoa(mphi))), &
                    'Susceptibility function I12 at cell centers, mphi='//trim(adjustl(itoa(mphi))), '1')
            end do

        end if

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
        end do
        call allocate_plasma_fields(plasma_temp, size(grid))
        plasma_temp%grid_size = size(grid)
        plasma_temp%r_grid = grid

        do i = 1, size(grid)
            call binsrc(plasma_in%r_grid, 1, size(plasma_in%r_grid), grid(i), ir)
            ibeg = max(1, ir - nlagr/2)
            iend = ibeg + nlagr - 1
            if (iend .gt. size(plasma_in%r_grid)) then
                iend = size(plasma_in%r_grid)
                ibeg = iend - nlagr + 1
            end if

            call plag_coeff(nlagr, nder, grid(i), plasma_in%r_grid(ibeg:iend), coef)

            do sp = 0, plasma_temp%n_species-1
                plasma_temp%spec(sp)%n(i)        = sum(coef(0,:) * plasma_in%spec(sp)%n(ibeg:iend))
                plasma_temp%spec(sp)%dndr(i)     = sum(coef(0,:) * plasma_in%spec(sp)%dndr(ibeg:iend))
                plasma_temp%spec(sp)%T(i)        = sum(coef(0,:) * plasma_in%spec(sp)%T(ibeg:iend))
                plasma_temp%spec(sp)%dTdr(i)     = sum(coef(0,:) * plasma_in%spec(sp)%dTdr(ibeg:iend))
                plasma_temp%spec(sp)%nu(i)       = sum(coef(0,:) * plasma_in%spec(sp)%nu(ibeg:iend))
                plasma_temp%spec(sp)%vT(i)       = sum(coef(0,:) * plasma_in%spec(sp)%vT(ibeg:iend))
                plasma_temp%spec(sp)%omega_c(i)  = sum(coef(0,:) * plasma_in%spec(sp)%omega_c(ibeg:iend))
                plasma_temp%spec(sp)%lambda_D(i) = sum(coef(0,:) * plasma_in%spec(sp)%lambda_D(ibeg:iend))
                plasma_temp%spec(sp)%rho_L(i)    = sum(coef(0,:) * plasma_in%spec(sp)%rho_L(ibeg:iend))
                plasma_temp%spec(sp)%z0(i)       = sum(coef(0,:) * plasma_in%spec(sp)%z0(ibeg:iend))
            end do

            plasma_temp%B0(i)   = sum(coef(0,:) * plasma_in%B0(ibeg:iend))
            plasma_temp%ks(i)   = sum(coef(0,:) * plasma_in%ks(ibeg:iend))
            plasma_temp%kp(i)   = sum(coef(0,:) * plasma_in%kp(ibeg:iend))
            plasma_temp%om_E(i) = sum(coef(0,:) * plasma_in%om_E(ibeg:iend))
            plasma_temp%q(i)    = sum(coef(0,:) * plasma_in%q(ibeg:iend))
            plasma_temp%dqdr(i) = sum(coef(0,:) * plasma_in%dqdr(ibeg:iend))
            plasma_temp%Er(i)   = sum(coef(0,:) * plasma_in%Er(ibeg:iend))

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
        call reallocate(spec%nu, grid_size)
        call reallocate(spec%vT, grid_size)
        call reallocate(spec%omega_c, grid_size)
        call reallocate(spec%lambda_D, grid_size)
        call reallocate(spec%rho_L, grid_size)
        call reallocate_complex(spec%z0, grid_size)

    end subroutine

    subroutine allocate_plasma_fields(plasma_in, grid_size)

        implicit none

        integer, intent(in) :: grid_size
        type(plasma_t), intent(inout) :: plasma_in

        call reallocate(plasma_in%ks, grid_size)
        call reallocate(plasma_in%B0, grid_size)
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

        call write_profile(rg_grid%xb, spec%rho_L, rg_grid%npts_b, 'rho_L')
        call plot_1D_labeled('rho_L', ' r [cm] ', ' rho_L [cm] ', '')
        call remove_file('rho_L.dat')

    end subroutine


    subroutine calculate_susc_funcs_profiles(spec, mphi)

        use resonances_mod, only: r_res
        use grid_m, only: width_res

        implicit none

        type(species_t), intent(inout) :: spec
        integer, intent(in) :: mphi
        integer :: j

        if (.not. allocated(spec%symbI)) allocate(spec%symbI(0:nmmax, 0:nmmax))
        spec%symbI = 0.0d0
        do j = 1, plasma%grid_size
            call getIfunc(spec%x1(j), spec%x2(j, mphi), spec%symbI)
            spec%I00(j, mphi) = spec%symbI(0, 0)
            spec%I20(j, mphi) = spec%symbI(2, 0)
            spec%I02(j, mphi) = spec%symbI(0, 2)
            spec%I01(j, mphi) = spec%symbI(0, 1)
            spec%I21(j, mphi) = spec%symbI(2, 1)
            spec%I22(j, mphi) = spec%symbI(2, 2)
            spec%I11(j, mphi) = spec%symbI(1, 1)
            spec%I13(j, mphi) = spec%symbI(1, 3)
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

        ! Read Er.dat with interpolation if grid doesn't match
        call read_and_interpolate_profile(trim(profile_location)//'Er.dat', &
                                          plasma%r_grid, plasma%Er, plasma%grid_size)

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

        ! Validate units - density should be in CGS (1/cm^3), typically 10^12 to 10^15
        call validate_profile_units(plasma)

        if (fstatus == 1) write(*,*) 'Status: Finished reading profiles from text files'

    end subroutine

    subroutine validate_profile_units(plasma)
        !> Check that density and temperature are in expected CGS units
        !> Density: 1/cm^3, typically 10^12 to 10^15 for fusion plasmas
        !> Temperature: eV, typically 10 to 20000 eV
        !> Also checks q vs m_mode sign consistency
        use setup_m, only: m_mode
        implicit none

        type(plasma_t), intent(in) :: plasma
        real(dp) :: n_max, T_max, q_mean
        integer :: sigma

        ! Check electron density
        n_max = maxval(plasma%spec(0)%n)

        if (n_max > 1.0d17) then
            write(*,*) ''
            write(*,*) '╔══════════════════════════════════════════════════════════════════╗'
            write(*,*) '║                    ERROR: DENSITY UNITS                          ║'
            write(*,*) '╠══════════════════════════════════════════════════════════════════╣'
            write(*,*) '║  Density appears to be in SI units (1/m^3) instead of CGS!       ║'
            write(*,*) '╚══════════════════════════════════════════════════════════════════╝'
            write(*,*) ''
            write(*,*) '  Maximum density found: ', n_max, ' 1/cm^3'
            write(*,*) '  Expected range (CGS):  1e12 to 1e15 1/cm^3'
            write(*,*) ''
            write(*,*) '  Your density values suggest SI units (1/m^3).'
            write(*,*) '  Please convert to CGS by dividing by 1e6.'
            write(*,*) ''
            write(*,*) '  Example: n_SI = 4.67e19 1/m^3  -->  n_CGS = 4.67e13 1/cm^3'
            write(*,*) ''
            stop 1
        else if (n_max < 1.0d10) then
            write(*,*) ''
            write(*,*) 'WARNING: Density values appear very low'
            write(*,*) '  Maximum density: ', n_max, ' 1/cm^3'
            write(*,*) '  Expected range:  1e12 to 1e15 1/cm^3'
            write(*,*) ''
        end if

        ! Check temperatures (should be in eV)
        ! Note: spec array is 0:n_species-1 but read_profiles uses number_of_ion_species from config
        do sigma = 0, min(plasma%n_species, size(plasma%spec)-1)
            if (.not. allocated(plasma%spec(sigma)%T)) cycle
            T_max = maxval(plasma%spec(sigma)%T)
            if (T_max > 1.0d6) then
                write(*,*) ''
                write(*,*) 'WARNING: Temperature appears very high'
                write(*,*) '  Species ', sigma, ' max T = ', T_max, ' eV'
                write(*,*) '  Expected range: 10 to 20000 eV'
                write(*,*) '  Check if temperature is in Kelvin instead of eV'
                write(*,*) ''
            end if
        end do

        ! Check q vs m_mode sign consistency
        ! For positive q, m should be negative (and vice versa) for proper helicity
        if (allocated(plasma%q) .and. size(plasma%q) > 0) then
            q_mean = sum(plasma%q) / size(plasma%q)
            if (q_mean > 0.0_dp .and. m_mode > 0) then
                write(*,*) ''
                write(*,*) '╔══════════════════════════════════════════════════════════════════╗'
                write(*,*) '║               WARNING: q AND m_mode SIGN MISMATCH                ║'
                write(*,*) '╠══════════════════════════════════════════════════════════════════╣'
                write(*,*) '║  Safety factor q > 0 typically requires m < 0 for resonance      ║'
                write(*,*) '╚══════════════════════════════════════════════════════════════════╝'
                write(*,*) ''
                write(*,*) '  Mean safety factor q = ', q_mean
                write(*,*) '  Poloidal mode number m = ', m_mode
                write(*,*) ''
                write(*,*) '  No resonant surfaces will be found with this configuration.'
                write(*,*) '  For resonant behavior, use m = ', -abs(m_mode)
                write(*,*) ''
            else if (q_mean < 0.0_dp .and. m_mode < 0) then
                write(*,*) ''
                write(*,*) '╔══════════════════════════════════════════════════════════════════╗'
                write(*,*) '║               WARNING: q AND m_mode SIGN MISMATCH                ║'
                write(*,*) '╠══════════════════════════════════════════════════════════════════╣'
                write(*,*) '║  Safety factor q < 0 typically requires m > 0 for resonance      ║'
                write(*,*) '╚══════════════════════════════════════════════════════════════════╝'
                write(*,*) ''
                write(*,*) '  Mean safety factor q = ', q_mean
                write(*,*) '  Poloidal mode number m = ', m_mode
                write(*,*) ''
                write(*,*) '  No resonant surfaces will be found with this configuration.'
                write(*,*) '  For resonant behavior, use m = ', abs(m_mode)
                write(*,*) ''
            end if
        end if

    end subroutine validate_profile_units

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

    subroutine read_and_interpolate_profile(filename, target_grid, profile_out, n_target)
        !> Read a profile file and interpolate to target grid if grids don't match
        !> Uses linear interpolation when source grid differs from target grid
        implicit none

        character(*), intent(in) :: filename
        real(dp), intent(in) :: target_grid(:)
        real(dp), intent(out) :: profile_out(:)
        integer, intent(in) :: n_target

        real(dp), allocatable :: r_src(:), val_src(:)
        real(dp) :: r_temp, val_temp, r_first_src, r_first_tgt
        integer :: n_src, i, j, ios, iunit
        real(dp) :: t
        logical :: grids_match
        real(dp), parameter :: GRID_TOL = 0.01_dp  ! 1% tolerance for grid matching

        ! First pass: count lines in source file
        n_src = 0
        open(newunit=iunit, file=trim(filename), status='old', action='read')
        do
            read(iunit, *, iostat=ios) r_temp, val_temp
            if (ios /= 0) exit
            n_src = n_src + 1
        end do
        close(iunit)

        if (n_src == 0) then
            write(*,*) 'ERROR: Empty profile file: ', trim(filename)
            stop 1
        end if

        ! Allocate and read source data
        allocate(r_src(n_src), val_src(n_src))
        open(newunit=iunit, file=trim(filename), status='old', action='read')
        do i = 1, n_src
            read(iunit, *) r_src(i), val_src(i)
        end do
        close(iunit)

        ! Check if grids match (same number of points and similar first value)
        grids_match = .false.
        if (n_src == n_target) then
            r_first_src = r_src(1)
            r_first_tgt = target_grid(1)
            if (abs(r_first_src - r_first_tgt) < GRID_TOL * abs(r_first_tgt + 1.0d-10)) then
                grids_match = .true.
            end if
        end if

        if (grids_match) then
            ! Grids match - direct copy
            profile_out(1:n_target) = val_src(1:n_target)
        else
            ! Grids don't match - interpolate
            write(*,*) 'Note: Interpolating ', trim(filename), ' to match grid'
            write(*,*) '  Source: ', n_src, ' points, r = [', r_src(1), ', ', r_src(n_src), ']'
            write(*,*) '  Target: ', n_target, ' points, r = [', target_grid(1), ', ', target_grid(n_target), ']'

            do i = 1, n_target
                r_temp = target_grid(i)

                ! Find bracketing indices in source grid
                if (r_temp <= r_src(1)) then
                    ! Extrapolate below (use first value)
                    profile_out(i) = val_src(1)
                else if (r_temp >= r_src(n_src)) then
                    ! Extrapolate above (use last value)
                    profile_out(i) = val_src(n_src)
                else
                    ! Find j such that r_src(j) <= r_temp < r_src(j+1)
                    j = 1
                    do while (j < n_src .and. r_src(j+1) < r_temp)
                        j = j + 1
                    end do
                    ! Linear interpolation
                    t = (r_temp - r_src(j)) / (r_src(j+1) - r_src(j))
                    profile_out(i) = val_src(j) + t * (val_src(j+1) - val_src(j))
                end if
            end do
        end if

        deallocate(r_src, val_src)

    end subroutine read_and_interpolate_profile

end module
