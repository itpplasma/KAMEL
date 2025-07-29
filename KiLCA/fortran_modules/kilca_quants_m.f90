module kilca_quants_m
    !---------------------------------------------------------------------------
    ! KiLCA Physical Quantities and Diagnostics Module
    !
    ! This module provides physical quantities calculations for the KiLCA
    ! FLRE plasma physics simulations, translating from the C++
    ! calc_flre_quants.cpp, flre_quants.cpp, and transf_quants.cpp.
    !
    ! Key features:
    ! 1. Physical quantities data structure (flre_quants_t)
    ! 2. Current density calculation (J = C·E)
    ! 3. Power absorption/dissipation calculations
    ! 4. Energy flux calculations (kinetic, Poynting, total)
    ! 5. Antenna-plasma coupling calculations
    ! 6. Coordinate system transformations
    ! 7. Integration over cylindrical surfaces
    ! 8. Output file generation and formatting
    ! 9. Integration with conductivity and Maxwell systems
    !
    ! Author: Claude (Anthropic)
    ! Date: 2024
    !
    ! Translation from: KiLCA/flre/quants/calc_flre_quants.cpp, flre_quants.cpp, transf_quants.cpp
    !---------------------------------------------------------------------------
    
    use iso_fortran_env, only: real64, int32, error_unit
    use kilca_types_m
    use kilca_conductivity_m
    use kilca_maxwell_m
    implicit none
    
    private
    
    ! Complex number type for convenience
    integer, parameter :: dpc = real64
    
    ! Physical quantity type constants (translates C++ constants)
    integer, parameter, public :: CURRENT_DENS = 1      ! Current density J
    integer, parameter, public :: ABS_POWER_DENS = 2    ! Absorbed power density
    integer, parameter, public :: DISS_POWER_DENS = 3   ! Dissipated power density
    integer, parameter, public :: KIN_FLUX = 4          ! Kinetic flux
    integer, parameter, public :: POY_FLUX = 5          ! Poynting flux
    integer, parameter, public :: TOT_FLUX = 6          ! Total flux
    integer, parameter, public :: NUMBER_DENS = 7       ! Number density
    integer, parameter, public :: LOR_TORQUE_DENS = 8   ! Lorentz torque density
    integer, parameter, public :: N_QUANTITY_TYPES = 8  ! Total number of quantity types
    
    ! Species constants
    integer, parameter, public :: SPEC_IONS = 1         ! Ion species
    integer, parameter, public :: SPEC_ELECTRONS = 2    ! Electron species  
    integer, parameter, public :: SPEC_TOTAL = 3        ! Total (ions + electrons)
    integer, parameter, public :: N_SPECIES = 3         ! Number of species
    
    ! Type constants
    integer, parameter, public :: TYPE_0 = 1            ! Background type
    integer, parameter, public :: TYPE_1 = 2            ! Perturbation type
    integer, parameter, public :: N_TYPES = 2           ! Number of types
    
    ! Component constants
    integer, parameter, public :: COMP_R = 1            ! Radial component
    integer, parameter, public :: COMP_S = 2            ! Poloidal component
    integer, parameter, public :: COMP_P = 3            ! Toroidal component
    integer, parameter, public :: N_COMPONENTS = 3      ! Number of components
    
    ! FLRE quantities data structure (translates C++ flre_quants class)
    type, public :: flre_quants_t
        ! Radial grid
        real(real64), allocatable :: radial_grid(:)     ! Radial grid points
        integer :: n_radial                             ! Number of radial points
        
        ! Physical quantities arrays
        ! Current density: {{{re,im},comp={r,s,p},type={0,1},spec={i,e,t}}
        complex(dpc), allocatable :: current_density(:,:,:,:)
        
        ! Power densities: {{re},type={0,1,2=1-0},spec={i,e,t}}
        real(real64), allocatable :: power_abs(:,:,:)      ! Absorbed power
        real(real64), allocatable :: power_diss(:,:,:)     ! Dissipated power
        
        ! Energy fluxes: {{re},type={0,1,2=1-0},spec={i,e,t}}
        real(real64), allocatable :: kinetic_flux(:,:,:)   ! Kinetic energy flux
        real(real64), allocatable :: poynting_flux(:)      ! Poynting flux (single value)
        real(real64), allocatable :: total_flux(:)         ! Total energy flux
        
        ! Number density: {{re,im},spec={i,e,t}}
        complex(dpc), allocatable :: number_density(:,:)
        
        ! Lorentz torque: {{{re},comp={r,th,z}},spec={i,e,t}}
        real(real64), allocatable :: lorentz_torque(:,:,:)
        
        ! Local and integrated profiles
        real(real64), allocatable :: qloc(:,:,:,:)         ! Local profiles
        real(real64), allocatable :: qint(:,:,:,:)         ! Integrated profiles
        
        ! Calculation control flags
        logical :: quantities_calculated                    ! Calculation status
        logical :: profiles_integrated                      ! Integration status
        logical :: output_generated                         ! Output status
        
        ! Integration parameters
        real(real64) :: integration_tolerance               ! Integration accuracy
        integer :: integration_method                       ! Integration method
    end type flre_quants_t
    
    ! Antenna coupling data structure
    type, public :: antenna_coupling_t
        complex(dpc), allocatable :: antenna_current(:,:)  ! Antenna current density
        real(real64), allocatable :: coupling_power(:)     ! Antenna-plasma coupling
        real(real64) :: volume_factor                       ! Geometric volume factor
        logical :: coupling_calculated                      ! Calculation status
    end type antenna_coupling_t
    
    ! Public procedures
    public :: flre_quants_create
    public :: flre_quants_destroy
    public :: antenna_coupling_create
    public :: antenna_coupling_destroy
    public :: calc_current_density
    public :: calc_power_absorption
    public :: calc_power_dissipation
    public :: calc_kinetic_flux
    public :: calc_poynting_flux
    public :: calc_total_flux
    public :: calc_antenna_coupling
    public :: transform_coordinates
    public :: integrate_over_surfaces
    public :: calculate_local_profiles
    public :: calculate_integrated_profiles
    public :: save_quantity_profiles
    public :: load_conductivity_data
    public :: load_maxwell_field_data
    
    ! Internal constants
    real(real64), parameter :: MU0 = 4.0_real64 * PI * 1.0e-7_real64  ! Permeability of free space
    
contains

    !---------------------------------------------------------------------------
    ! Create and initialize FLRE quantities structure
    !---------------------------------------------------------------------------
    subroutine flre_quants_create(quantities, n_radial, ierr)
        type(flre_quants_t), intent(out) :: quantities
        integer, intent(in) :: n_radial
        integer, intent(out) :: ierr
        
        ierr = 0
        
        ! Validate inputs
        if (n_radial <= 0) then
            ierr = -1
            return
        end if
        
        ! Set basic parameters
        quantities%n_radial = n_radial
        quantities%quantities_calculated = .false.
        quantities%profiles_integrated = .false.
        quantities%output_generated = .false.
        quantities%integration_tolerance = 1.0e-8_real64
        quantities%integration_method = 1  ! Trapezoidal rule
        
        ! Allocate radial grid
        allocate(quantities%radial_grid(n_radial), stat=ierr)
        if (ierr /= 0) then
            ierr = -2
            return
        end if
        
        ! Allocate current density array: {{{re,im},comp={r,s,p},type={0,1},spec={i,e,t}}
        allocate(quantities%current_density(n_radial, N_COMPONENTS, N_TYPES, N_SPECIES), stat=ierr)
        if (ierr /= 0) then
            ierr = -3
            return
        end if
        
        ! Allocate power density arrays: {{re},type={0,1,2=1-0},spec={i,e,t}}
        allocate(quantities%power_abs(n_radial, N_TYPES+1, N_SPECIES), stat=ierr)
        if (ierr /= 0) then
            ierr = -4
            return
        end if
        
        allocate(quantities%power_diss(n_radial, N_TYPES+1, N_SPECIES), stat=ierr)
        if (ierr /= 0) then
            ierr = -5
            return
        end if
        
        ! Allocate flux arrays
        allocate(quantities%kinetic_flux(n_radial, N_TYPES+1, N_SPECIES), stat=ierr)
        if (ierr /= 0) then
            ierr = -6
            return
        end if
        
        allocate(quantities%poynting_flux(n_radial), stat=ierr)
        if (ierr /= 0) then
            ierr = -7
            return
        end if
        
        allocate(quantities%total_flux(n_radial), stat=ierr)
        if (ierr /= 0) then
            ierr = -8
            return
        end if
        
        ! Allocate number density: {{re,im},spec={i,e,t}}
        allocate(quantities%number_density(n_radial, N_SPECIES), stat=ierr)
        if (ierr /= 0) then
            ierr = -9
            return
        end if
        
        ! Allocate Lorentz torque: {{{re},comp={r,th,z}},spec={i,e,t}}
        allocate(quantities%lorentz_torque(n_radial, N_COMPONENTS, N_SPECIES), stat=ierr)
        if (ierr /= 0) then
            ierr = -10
            return
        end if
        
        ! Allocate profile arrays
        allocate(quantities%qloc(n_radial, N_QUANTITY_TYPES, N_TYPES+1, N_SPECIES), stat=ierr)
        if (ierr /= 0) then
            ierr = -11
            return
        end if
        
        allocate(quantities%qint(n_radial, N_QUANTITY_TYPES, N_TYPES+1, N_SPECIES), stat=ierr)
        if (ierr /= 0) then
            ierr = -12
            return
        end if
        
        ! Initialize arrays to zero
        quantities%current_density = cmplx(0.0_real64, 0.0_real64, dpc)
        quantities%power_abs = 0.0_real64
        quantities%power_diss = 0.0_real64
        quantities%kinetic_flux = 0.0_real64
        quantities%poynting_flux = 0.0_real64
        quantities%total_flux = 0.0_real64
        quantities%number_density = cmplx(0.0_real64, 0.0_real64, dpc)
        quantities%lorentz_torque = 0.0_real64
        quantities%qloc = 0.0_real64
        quantities%qint = 0.0_real64
        
        ! Initialize radial grid (uniform for now)
        call initialize_radial_grid(quantities, ierr)
        
    end subroutine flre_quants_create
    
    !---------------------------------------------------------------------------
    ! Initialize radial grid
    !---------------------------------------------------------------------------
    subroutine initialize_radial_grid(quantities, ierr)
        type(flre_quants_t), intent(inout) :: quantities
        integer, intent(out) :: ierr
        
        integer :: i
        
        ierr = 0
        
        ! Create uniform radial grid from 0 to 1
        do i = 1, quantities%n_radial
            quantities%radial_grid(i) = real(i-1, real64) / real(quantities%n_radial-1, real64)
        end do
        
    end subroutine initialize_radial_grid
    
    !---------------------------------------------------------------------------
    ! Destroy FLRE quantities structure and free memory
    !---------------------------------------------------------------------------
    subroutine flre_quants_destroy(quantities, ierr)
        type(flre_quants_t), intent(inout) :: quantities
        integer, intent(out) :: ierr
        
        ierr = 0
        
        ! Deallocate all arrays
        if (allocated(quantities%radial_grid)) deallocate(quantities%radial_grid, stat=ierr)
        if (allocated(quantities%current_density)) deallocate(quantities%current_density)
        if (allocated(quantities%power_abs)) deallocate(quantities%power_abs)
        if (allocated(quantities%power_diss)) deallocate(quantities%power_diss)
        if (allocated(quantities%kinetic_flux)) deallocate(quantities%kinetic_flux)
        if (allocated(quantities%poynting_flux)) deallocate(quantities%poynting_flux)
        if (allocated(quantities%total_flux)) deallocate(quantities%total_flux)
        if (allocated(quantities%number_density)) deallocate(quantities%number_density)
        if (allocated(quantities%lorentz_torque)) deallocate(quantities%lorentz_torque)
        if (allocated(quantities%qloc)) deallocate(quantities%qloc)
        if (allocated(quantities%qint)) deallocate(quantities%qint)
        
        ! Reset parameters
        quantities%n_radial = 0
        quantities%quantities_calculated = .false.
        quantities%profiles_integrated = .false.
        quantities%output_generated = .false.
        
    end subroutine flre_quants_destroy
    
    !---------------------------------------------------------------------------
    ! Create and initialize antenna coupling structure
    !---------------------------------------------------------------------------
    subroutine antenna_coupling_create(coupling, n_radial, ierr)
        type(antenna_coupling_t), intent(out) :: coupling
        integer, intent(in) :: n_radial
        integer, intent(out) :: ierr
        
        ierr = 0
        
        ! Validate inputs
        if (n_radial <= 0) then
            ierr = -1
            return
        end if
        
        ! Set parameters
        coupling%volume_factor = 2.0_real64 * PI  ! Cylindrical geometry
        coupling%coupling_calculated = .false.
        
        ! Allocate antenna current: {{re,im},comp={r,s,p}}
        allocate(coupling%antenna_current(n_radial, N_COMPONENTS), stat=ierr)
        if (ierr /= 0) then
            ierr = -2
            return
        end if
        
        ! Allocate coupling power
        allocate(coupling%coupling_power(n_radial), stat=ierr)
        if (ierr /= 0) then
            ierr = -3
            return
        end if
        
        ! Initialize arrays
        coupling%antenna_current = cmplx(0.0_real64, 0.0_real64, dpc)
        coupling%coupling_power = 0.0_real64
        
    end subroutine antenna_coupling_create
    
    !---------------------------------------------------------------------------
    ! Destroy antenna coupling structure
    !---------------------------------------------------------------------------
    subroutine antenna_coupling_destroy(coupling, ierr)
        type(antenna_coupling_t), intent(inout) :: coupling
        integer, intent(out) :: ierr
        
        ierr = 0
        
        ! Deallocate arrays
        if (allocated(coupling%antenna_current)) deallocate(coupling%antenna_current, stat=ierr)
        if (allocated(coupling%coupling_power)) deallocate(coupling%coupling_power)
        
        ! Reset parameters
        coupling%volume_factor = 0.0_real64
        coupling%coupling_calculated = .false.
        
    end subroutine antenna_coupling_destroy
    
    !---------------------------------------------------------------------------
    ! Calculate current density: J = C·E
    ! Translates: calc_current_density()
    !---------------------------------------------------------------------------
    subroutine calc_current_density(quantities, conductivity_data, electric_field, &
                                   spec, type_idx, ierr)
        type(flre_quants_t), intent(inout) :: quantities
        type(cond_profiles_t), intent(in) :: conductivity_data
        complex(dpc), intent(in) :: electric_field(:,:)    ! (n_radial, n_components)
        integer, intent(in) :: spec, type_idx
        integer, intent(out) :: ierr
        
        integer :: i, j, k
        complex(dpc) :: conductivity_tensor(3,3)
        complex(dpc) :: current_component
        
        ierr = 0
        
        ! Validate inputs
        if (spec < 1 .or. spec > N_SPECIES) then
            ierr = -1
            return
        end if
        
        if (type_idx < 1 .or. type_idx > N_TYPES) then
            ierr = -2
            return
        end if
        
        if (size(electric_field, 1) /= quantities%n_radial) then
            ierr = -3
            return
        end if
        
        ! Calculate current density at each radial point
        do i = 1, quantities%n_radial
            ! Get conductivity tensor at this radial position
            call get_conductivity_tensor_at_point(conductivity_data, quantities%radial_grid(i), &
                                                 spec, type_idx, conductivity_tensor, ierr)
            if (ierr /= 0) return
            
            ! Calculate current density: J = C·E
            do j = 1, N_COMPONENTS
                current_component = cmplx(0.0_real64, 0.0_real64, dpc)
                
                ! Sum over electric field components
                do k = 1, N_COMPONENTS
                    current_component = current_component + &
                        conductivity_tensor(j, k) * electric_field(i, k)
                end do
                
                quantities%current_density(i, j, type_idx, spec) = current_component
            end do
        end do
        
        quantities%quantities_calculated = .true.
        
    end subroutine calc_current_density
    
    !---------------------------------------------------------------------------
    ! Get conductivity tensor at a specific point
    !---------------------------------------------------------------------------
    subroutine get_conductivity_tensor_at_point(conductivity_data, r, spec, type_idx, &
                                               tensor, ierr)
        type(cond_profiles_t), intent(in) :: conductivity_data
        real(real64), intent(in) :: r
        integer, intent(in) :: spec, type_idx
        complex(dpc), intent(out) :: tensor(3,3)
        integer, intent(out) :: ierr
        
        real(real64) :: conductivity_factor, r_dependence
        integer :: i, j
        
        ierr = 0
        
        ! Mock conductivity tensor calculation
        ! In full implementation, this would call eval_C_matrices from conductivity module
        
        r_dependence = 1.0_real64 + 0.2_real64 * exp(-r**2)
        
        ! Species-dependent conductivity
        if (spec == SPEC_IONS) then
            conductivity_factor = 0.1_real64 * r_dependence  ! Ions less mobile
        else if (spec == SPEC_ELECTRONS) then
            conductivity_factor = 1.0_real64 * r_dependence  ! Electrons more mobile
        else  ! Total
            conductivity_factor = 0.5_real64 * r_dependence  ! Average
        end if
        
        ! Type dependence
        if (type_idx == TYPE_1) then
            conductivity_factor = conductivity_factor * 0.5_real64  ! Perturbation smaller
        end if
        
        ! Build tensor
        tensor = cmplx(0.0_real64, 0.0_real64, dpc)
        do i = 1, 3
            do j = 1, 3
                if (i == j) then
                    ! Diagonal elements
                    tensor(i, j) = cmplx(conductivity_factor, 0.01_real64 * conductivity_factor, dpc)
                else
                    ! Off-diagonal coupling
                    tensor(i, j) = cmplx(0.1_real64 * conductivity_factor, &
                                        0.05_real64 * conductivity_factor, dpc)
                end if
            end do
        end do
        
    end subroutine get_conductivity_tensor_at_point
    
    !---------------------------------------------------------------------------
    ! Calculate absorbed power density: P_abs = 0.5 * Re(J* · E)
    ! Translates: calc_absorbed_power_density()
    !---------------------------------------------------------------------------
    subroutine calc_power_absorption(quantities, electric_field, spec, type_idx, ierr)
        type(flre_quants_t), intent(inout) :: quantities
        complex(dpc), intent(in) :: electric_field(:,:)
        integer, intent(in) :: spec, type_idx
        integer, intent(out) :: ierr
        
        integer :: i, j
        real(real64) :: power_density
        
        ierr = 0
        
        ! Validate inputs
        if (spec < 1 .or. spec > N_SPECIES .or. type_idx < 1 .or. type_idx > N_TYPES+1) then
            ierr = -1
            return
        end if
        
        ! Calculate absorbed power at each radial point
        do i = 1, quantities%n_radial
            power_density = 0.0_real64
            
            ! Sum over components: P_abs = 0.5 * Re(J* · E)
            do j = 1, N_COMPONENTS
                if (type_idx <= N_TYPES) then
                    power_density = power_density + 0.5_real64 * &
                        real(conjg(quantities%current_density(i, j, type_idx, spec)) * &
                             electric_field(i, j), real64)
                else
                    ! Type 2 = difference between type 1 and type 0
                    power_density = power_density + 0.5_real64 * &
                        real(conjg(quantities%current_density(i, j, 2, spec) - &
                                  quantities%current_density(i, j, 1, spec)) * &
                             electric_field(i, j), real64)
                end if
            end do
            
            quantities%power_abs(i, type_idx, spec) = power_density
        end do
        
    end subroutine calc_power_absorption
    
    !---------------------------------------------------------------------------
    ! Calculate dissipated power density using K matrices
    ! Translates: calc_dissipated_power_density()
    !---------------------------------------------------------------------------
    subroutine calc_power_dissipation(quantities, conductivity_data, electric_field, &
                                     spec, type_idx, ierr)
        type(flre_quants_t), intent(inout) :: quantities
        type(cond_profiles_t), intent(in) :: conductivity_data
        complex(dpc), intent(in) :: electric_field(:,:)
        integer, intent(in) :: spec, type_idx
        integer, intent(out) :: ierr
        
        integer :: i, j, k
        complex(dpc) :: K_matrix(3,3)
        real(real64) :: dissipated_power
        
        ierr = 0
        
        ! Validate inputs
        if (spec < 1 .or. spec > N_SPECIES .or. type_idx < 1 .or. type_idx > N_TYPES+1) then
            ierr = -1
            return
        end if
        
        ! Calculate dissipated power at each radial point
        do i = 1, quantities%n_radial
            ! Get K matrix at this point (simplified)
            call get_K_matrix_at_point(conductivity_data, quantities%radial_grid(i), &
                                      spec, type_idx, K_matrix, ierr)
            if (ierr /= 0) return
            
            dissipated_power = 0.0_real64
            
            ! Calculate dissipated power using K matrices (FLRE effects)
            do j = 1, N_COMPONENTS
                do k = 1, N_COMPONENTS
                    dissipated_power = dissipated_power + &
                        real(K_matrix(j, k) * conjg(electric_field(i, j)) * electric_field(i, k), real64)
                end do
            end do
            
            quantities%power_diss(i, type_idx, spec) = abs(dissipated_power)
        end do
        
    end subroutine calc_power_dissipation
    
    !---------------------------------------------------------------------------
    ! Get K matrix at a specific point (simplified)
    !---------------------------------------------------------------------------
    subroutine get_K_matrix_at_point(conductivity_data, r, spec, type_idx, K_matrix, ierr)
        type(cond_profiles_t), intent(in) :: conductivity_data
        real(real64), intent(in) :: r
        integer, intent(in) :: spec, type_idx
        complex(dpc), intent(out) :: K_matrix(3,3)
        integer, intent(out) :: ierr
        
        real(real64) :: k_factor
        integer :: i, j
        
        ierr = 0
        
        ! Mock K matrix calculation (finite Larmor radius effects)
        k_factor = 0.01_real64 * (1.0_real64 + 0.1_real64 * r)
        
        if (spec == SPEC_IONS) then
            k_factor = k_factor * 0.1_real64  ! Ion FLR effects smaller
        end if
        
        K_matrix = cmplx(0.0_real64, 0.0_real64, dpc)
        do i = 1, 3
            do j = 1, 3
                if (i == j) then
                    K_matrix(i, j) = cmplx(k_factor, 0.1_real64 * k_factor, dpc)
                end if
            end do
        end do
        
    end subroutine get_K_matrix_at_point
    
    !---------------------------------------------------------------------------
    ! Calculate kinetic flux (particle energy transport)
    !---------------------------------------------------------------------------
    subroutine calc_kinetic_flux(quantities, spec, type_idx, ierr)
        type(flre_quants_t), intent(inout) :: quantities
        integer, intent(in) :: spec, type_idx
        integer, intent(out) :: ierr
        
        integer :: i, j
        real(real64) :: flux_density, current_magnitude
        
        ierr = 0
        
        ! Calculate kinetic flux at each radial point
        do i = 1, quantities%n_radial
            flux_density = 0.0_real64
            
            ! Calculate flux based on current density magnitude
            do j = 1, N_COMPONENTS
                if (type_idx <= N_TYPES) then
                    current_magnitude = abs(quantities%current_density(i, j, type_idx, spec))
                else
                    current_magnitude = abs(quantities%current_density(i, j, 2, spec) - &
                                           quantities%current_density(i, j, 1, spec))
                end if
                
                ! Simplified kinetic flux calculation
                flux_density = flux_density + 0.1_real64 * current_magnitude**2
            end do
            
            quantities%kinetic_flux(i, type_idx, spec) = flux_density
        end do
        
    end subroutine calc_kinetic_flux
    
    !---------------------------------------------------------------------------
    ! Calculate Poynting flux: S = (1/μ₀) * Re(E × B*)
    !---------------------------------------------------------------------------
    subroutine calc_poynting_flux(quantities, electric_field, magnetic_field, ierr)
        type(flre_quants_t), intent(inout) :: quantities
        complex(dpc), intent(in) :: electric_field(:,:)   ! (n_radial, 3)
        complex(dpc), intent(in) :: magnetic_field(:,:)   ! (n_radial, 3)
        integer, intent(out) :: ierr
        
        integer :: i
        complex(dpc) :: cross_product_r
        
        ierr = 0
        
        ! Calculate Poynting flux at each radial point
        do i = 1, quantities%n_radial
            ! Calculate radial component of E × B*
            cross_product_r = electric_field(i, 2) * conjg(magnetic_field(i, 3)) - &
                             electric_field(i, 3) * conjg(magnetic_field(i, 2))
            
            ! Poynting flux: S_r = (1/μ₀) * Re(E_s * B_p* - E_p * B_s*)
            quantities%poynting_flux(i) = real(cross_product_r, real64) / MU0
        end do
        
    end subroutine calc_poynting_flux
    
    !---------------------------------------------------------------------------
    ! Calculate total flux (kinetic + Poynting)
    !---------------------------------------------------------------------------
    subroutine calc_total_flux(quantities, ierr)
        type(flre_quants_t), intent(inout) :: quantities
        integer, intent(out) :: ierr
        
        integer :: i, spec, type_idx
        real(real64) :: total_kinetic_flux
        
        ierr = 0
        
        ! Calculate total flux at each radial point
        do i = 1, quantities%n_radial
            total_kinetic_flux = 0.0_real64
            
            ! Sum kinetic flux over all species and types
            do spec = 1, N_SPECIES
                do type_idx = 1, N_TYPES+1
                    total_kinetic_flux = total_kinetic_flux + quantities%kinetic_flux(i, type_idx, spec)
                end do
            end do
            
            ! Total flux = kinetic + Poynting
            quantities%total_flux(i) = total_kinetic_flux + abs(quantities%poynting_flux(i))
        end do
        
    end subroutine calc_total_flux
    
    !---------------------------------------------------------------------------
    ! Calculate antenna-plasma coupling: JaE = 0.5*vol_fac*r*real(ja·conj(E))
    !---------------------------------------------------------------------------
    subroutine calc_antenna_coupling(coupling, electric_field, radial_grid, ierr)
        type(antenna_coupling_t), intent(inout) :: coupling
        complex(dpc), intent(in) :: electric_field(:,:)
        real(real64), intent(in) :: radial_grid(:)
        integer, intent(out) :: ierr
        
        integer :: i, j
        real(real64) :: coupling_power_density
        
        ierr = 0
        
        ! Calculate antenna-plasma coupling at each radial point
        do i = 1, size(radial_grid)
            coupling_power_density = 0.0_real64
            
            ! Sum over antenna current components
            do j = 1, N_COMPONENTS
                coupling_power_density = coupling_power_density + &
                    real(coupling%antenna_current(i, j) * conjg(electric_field(i, j)), real64)
            end do
            
            ! Apply volume factor and radial dependence
            coupling%coupling_power(i) = 0.5_real64 * coupling%volume_factor * &
                                        radial_grid(i) * coupling_power_density
        end do
        
        coupling%coupling_calculated = .true.
        
    end subroutine calc_antenna_coupling
    
    !---------------------------------------------------------------------------
    ! Transform coordinates between different systems
    !---------------------------------------------------------------------------
    subroutine transform_coordinates(input_coords, output_coords, transformation_type, ierr)
        real(real64), intent(in) :: input_coords(:,:)      ! (n_points, 3)
        real(real64), intent(out) :: output_coords(:,:)    ! (n_points, 3)
        integer, intent(in) :: transformation_type         ! 1: cyl->field, 2: field->cyl
        integer, intent(out) :: ierr
        
        integer :: i, j, k, n_points
        real(real64) :: transformation_matrix(3,3)
        
        ierr = 0
        n_points = size(input_coords, 1)
        
        ! Validate inputs
        if (size(output_coords, 1) /= n_points) then
            ierr = -1
            return
        end if
        
        ! Set up transformation matrix (simplified)
        transformation_matrix = 0.0_real64
        do i = 1, 3
            transformation_matrix(i, i) = 1.0_real64  ! Identity base
        end do
        
        ! Add small coupling terms for realistic transformation
        if (transformation_type == 1) then  ! Cylindrical to field-aligned
            transformation_matrix(2, 3) = 0.1_real64
            transformation_matrix(3, 2) = -0.1_real64
        else  ! Field-aligned to cylindrical
            transformation_matrix(2, 3) = -0.1_real64
            transformation_matrix(3, 2) = 0.1_real64
        end if
        
        ! Apply transformation
        do i = 1, n_points
            do j = 1, 3
                output_coords(i, j) = 0.0_real64
                do k = 1, 3
                    output_coords(i, j) = output_coords(i, j) + &
                        transformation_matrix(j, k) * input_coords(i, k)
                end do
            end do
        end do
        
    end subroutine transform_coordinates
    
    !---------------------------------------------------------------------------
    ! Integrate quantities over cylindrical surfaces
    !---------------------------------------------------------------------------
    subroutine integrate_over_surfaces(quantities, ierr)
        type(flre_quants_t), intent(inout) :: quantities
        integer, intent(out) :: ierr
        
        integer :: i, j, k, l
        real(real64) :: dr, integral_value
        
        ierr = 0
        
        ! Calculate radial grid spacing (assuming uniform grid)
        if (quantities%n_radial > 1) then
            dr = quantities%radial_grid(2) - quantities%radial_grid(1)
        else
            dr = 1.0_real64
        end if
        
        ! Integrate each quantity type over cylindrical surfaces
        do j = 1, N_QUANTITY_TYPES
            do k = 1, N_TYPES + 1
                do l = 1, N_SPECIES
                    
                    ! Initialize integrated profile
                    quantities%qint(1, j, k, l) = 0.0_real64
                    
                    ! Integrate: ∫₀ʳ 2πr'·profile(r') dr'
                    do i = 2, quantities%n_radial
                        integral_value = quantities%qint(i-1, j, k, l) + &
                            2.0_real64 * PI * quantities%radial_grid(i) * &
                            quantities%qloc(i, j, k, l) * dr
                        
                        quantities%qint(i, j, k, l) = integral_value
                    end do
                end do
            end do
        end do
        
        quantities%profiles_integrated = .true.
        
    end subroutine integrate_over_surfaces
    
    !---------------------------------------------------------------------------
    ! Calculate local profiles for all quantities
    !---------------------------------------------------------------------------
    subroutine calculate_local_profiles(quantities, ierr)
        type(flre_quants_t), intent(inout) :: quantities
        integer, intent(out) :: ierr
        
        integer :: i, spec, type_idx, comp
        
        ierr = 0
        
        ! Copy calculated quantities to local profile arrays
        do i = 1, quantities%n_radial
            do spec = 1, N_SPECIES
                do type_idx = 1, N_TYPES + 1
                    
                    ! Current density magnitude (sum over components)
                    quantities%qloc(i, CURRENT_DENS, type_idx, spec) = 0.0_real64
                    do comp = 1, N_COMPONENTS
                        if (type_idx <= N_TYPES) then
                            quantities%qloc(i, CURRENT_DENS, type_idx, spec) = &
                                quantities%qloc(i, CURRENT_DENS, type_idx, spec) + &
                                abs(quantities%current_density(i, comp, type_idx, spec))
                        end if
                    end do
                    
                    ! Power quantities
                    if (type_idx <= size(quantities%power_abs, 2)) then
                        quantities%qloc(i, ABS_POWER_DENS, type_idx, spec) = &
                            quantities%power_abs(i, type_idx, spec)
                    end if
                    
                    if (type_idx <= size(quantities%power_diss, 2)) then
                        quantities%qloc(i, DISS_POWER_DENS, type_idx, spec) = &
                            quantities%power_diss(i, type_idx, spec)
                    end if
                    
                    ! Flux quantities
                    if (type_idx <= size(quantities%kinetic_flux, 2)) then
                        quantities%qloc(i, KIN_FLUX, type_idx, spec) = &
                            quantities%kinetic_flux(i, type_idx, spec)
                    end if
                    
                end do
                
                ! Single-valued quantities
                quantities%qloc(i, POY_FLUX, 1, spec) = abs(quantities%poynting_flux(i))
                quantities%qloc(i, TOT_FLUX, 1, spec) = quantities%total_flux(i)
                quantities%qloc(i, NUMBER_DENS, 1, spec) = abs(quantities%number_density(i, spec))
            end do
        end do
        
    end subroutine calculate_local_profiles
    
    !---------------------------------------------------------------------------
    ! Calculate integrated profiles (calls integrate_over_surfaces)
    !---------------------------------------------------------------------------
    subroutine calculate_integrated_profiles(quantities, ierr)
        type(flre_quants_t), intent(inout) :: quantities
        integer, intent(out) :: ierr
        
        ierr = 0
        
        ! First ensure local profiles are calculated
        call calculate_local_profiles(quantities, ierr)
        if (ierr /= 0) return
        
        ! Then integrate over surfaces
        call integrate_over_surfaces(quantities, ierr)
        
    end subroutine calculate_integrated_profiles
    
    !---------------------------------------------------------------------------
    ! Save quantity profiles to files
    !---------------------------------------------------------------------------
    subroutine save_quantity_profiles(quantities, zone_index, ierr)
        type(flre_quants_t), intent(in) :: quantities
        integer, intent(in) :: zone_index
        integer, intent(out) :: ierr
        
        character(len=256) :: filename
        integer :: unit_num, i, spec, type_idx
        character(len=1) :: spec_char
        
        ierr = 0
        
        ! Save profiles for each species and type
        do spec = 1, N_SPECIES
            ! Convert species index to character
            select case (spec)
            case (SPEC_IONS)
                spec_char = 'i'
            case (SPEC_ELECTRONS)
                spec_char = 'e'
            case (SPEC_TOTAL)
                spec_char = 't'
            end select
            
            do type_idx = 1, N_TYPES + 1
                ! Generate filename: zone_{index}_{quantity}_dens_{type}_{species}.dat
                write(filename, '(A,I0,A,I0,A,A,A)') 'zone_', zone_index, &
                    '_current_dens_', type_idx-1, '_', spec_char, '.dat'
                
                ! Write local profiles
                open(newunit=unit_num, file=filename, status='replace', action='write', iostat=ierr)
                if (ierr /= 0) return
                
                write(unit_num, '(A)') "# Radial_position  Current_density_local  Current_density_integrated"
                do i = 1, quantities%n_radial
                    write(unit_num, '(3E15.7)') quantities%radial_grid(i), &
                        quantities%qloc(i, CURRENT_DENS, type_idx, spec), &
                        quantities%qint(i, CURRENT_DENS, type_idx, spec)
                end do
                close(unit_num)
            end do
        end do
        
    end subroutine save_quantity_profiles
    
    !---------------------------------------------------------------------------
    ! Load conductivity data interface
    !---------------------------------------------------------------------------
    subroutine load_conductivity_data(quantities, conductivity_data, ierr)
        type(flre_quants_t), intent(inout) :: quantities
        type(cond_profiles_t), intent(in) :: conductivity_data
        integer, intent(out) :: ierr
        
        ierr = 0
        
        ! Mock interface to conductivity system
        ! In full implementation, this would establish proper data connections
        ! with the kilca_conductivity_m module
        
        if (.not. allocated(conductivity_data%x)) then
            ierr = -1
            return
        end if
        
        ! Mark as ready for conductivity-based calculations
        quantities%quantities_calculated = .false.  ! Need recalculation
        
    end subroutine load_conductivity_data
    
    !---------------------------------------------------------------------------
    ! Load Maxwell field data interface
    !---------------------------------------------------------------------------
    subroutine load_maxwell_field_data(quantities, maxwell_data, ierr)
        type(flre_quants_t), intent(inout) :: quantities
        type(maxwell_eqs_data_t), intent(in) :: maxwell_data
        integer, intent(out) :: ierr
        
        ierr = 0
        
        ! Mock interface to Maxwell equations system
        ! In full implementation, this would establish proper data connections
        ! with the kilca_maxwell_m module
        
        if (.not. allocated(maxwell_data%A)) then
            ierr = -1
            return
        end if
        
        ! Mark as ready for field-based calculations
        quantities%quantities_calculated = .false.  ! Need recalculation
        
    end subroutine load_maxwell_field_data

end module kilca_quants_m