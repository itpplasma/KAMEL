module config_m

    use KIM_kinds_m, only: dp

    implicit none

    ! KIM_CONFIG namelist variables
    integer :: number_of_ion_species ! number of ion species
    logical :: read_species_from_namelist ! read species from namelist or use deuterium plasma
    character(100) :: type_of_run
    character(100) :: collision_model ! type of collision model
    character(20) :: ion_collision_model = 'FokkerPlanck' ! FokkerPlanck or collisionless Krook/Hamiltonian
    real(dp) :: collisionless_kpar_epsilon = -1.0_dp ! causal +i epsilon pole [1/cm]; required for collisionless ions
    real(dp) :: ion_fp_collision_scale = 1.0_dp ! multiplier applied only to computed FP ion collision frequencies
    integer :: artificial_debye_case
    logical :: turn_off_ions ! if true, only the first species (electrons) is considered in calculations
    logical :: turn_off_electrons
    character(100) :: plasma_type ! type of plasma ('H' for hydrogen, 'D' for deuterium)
    integer, parameter :: IFUNC_MODEL_INHERIT = -1
    integer, parameter :: IFUNC_MODEL_NUMBER = 0
    integer, parameter :: IFUNC_MODEL_ENERGY = 1
    integer, parameter :: IFUNC_MODEL_MOMENTUM = 2
    integer, parameter :: IFUNC_MODEL_ENERGY_MOMENTUM = 3
    ! KiLCA-compatible per-species conservation models: 0=N, 1=N+E,
    ! 2=N+P_parallel, 3=N+E+P_parallel. The -1 default inherits the legacy
    ! boole_energy_conservation setting.
    integer :: electron_ifunc_conservation_model = IFUNC_MODEL_INHERIT
    integer :: ion_ifunc_conservation_model = IFUNC_MODEL_INHERIT
    integer :: resolved_electron_ifunc_conservation_model = IFUNC_MODEL_ENERGY
    integer :: resolved_ion_ifunc_conservation_model = IFUNC_MODEL_ENERGY
    ! Deprecated compatibility setting. It is consulted only when a species'
    ! integer conservation model is left at IFUNC_MODEL_INHERIT.
    logical :: boole_energy_conservation = .true.

    ! WKB_DISPERSION namelist variables
    character(20) :: WKB_dispersion_mode = 'KIM'      ! 'KIM' (full Bessel) or 'FLRE' (finite Larmor radius)
    character(20) :: WKB_dispersion_solver = 'Muller' ! 'Muller'
    logical :: WKB_solve_for_kr_squared = .false.     ! If true, solve for kr^2; if false, solve for kr
    integer :: WKB_max_tracked_branches = 4           ! Max branches for per-branch root tracking
    real(dp) :: WKB_branch_search_halfwidth = 1.5_dp  ! Search window half-width for tracking
    real(dp) :: WKB_broad_search_halfwidth = 5.0_dp   ! Broad search window half-width
    integer :: WKB_broad_search_interval = 0          ! Run broad search every N points (0=only at start)
    real(dp) :: WKB_root_tolerance = 1.0d-6           ! |f(z)| tolerance for valid roots
    logical :: WKB_verbose = .false.                  ! Verbose dispersion solver output

    ! KIM_PERIODIC namelist variables (forced-periodicity electrostatic run-type).
    ! Window half-widths and Fourier cutoff are set in units of the resonant-
    ! surface Larmor radius rho_L(rm): dx_asis = scale*rho_L, dx_tr = scale*rho_L,
    ! k_max = scale/rho_L. Defaults apply when the &KIM_PERIODIC group is absent.
    real(dp) :: periodic_dr_asis_scale = 5.0_dp   ! as-is half-width / rho_L(rm)
    real(dp) :: periodic_dr_tr_scale   = 10.0_dp  ! transition width / rho_L(rm)
    real(dp) :: periodic_kmax_scale    = 5.0_dp   ! k_max * rho_L(rm)
    integer  :: periodic_n_rg          = 96       ! window grid boundary points
    ! Reproduce the global FEM kernel approximations for direct comparisons:
    ! drop k_s^2 from Bessel arguments and take electrons in the zero-FLR limit.
    ! The forced-periodic solver uses its full Fourier kernel by default.
    logical :: periodic_match_global_kernel_approximations = .false.

    ! KIM_FLR2 namelist variables. These switches only select terms in the
    ! imported FLR2 response; all profiles and susceptibilities come from KIM.
    logical :: flr2_electron_flr = .true.
    logical :: flr2_ion_flr = .true.
    logical :: flr2_electron_potential = .true.
    logical :: flr2_ion_potential = .true.
    logical :: flr2_electron_current = .true.
    logical :: flr2_ion_current = .true.
    logical :: flr2_include_potential_in_current = .true.

    ! KIM_IO namelist variables
    character(256) :: profile_location ! path to profile directory
    character(256) :: output_path         ! path to output directory
    character(256) :: dispersion_output_path  ! path to dispersion output subdirectory
    character(256) :: h5_out_file ! file name of the hdf5 output file
    logical :: hdf5_input, hdf5_output
    integer :: log_level = 3         ! maps to LVL_INFO
    integer :: data_verbosity = 1    ! standard output
    logical :: calculate_asymptotics ! enable/disable asymptotic calculations
    ! Write kim_diagnostics.dat (dqle22, Ipar, Ipar_e) regardless of
    ! hdf5_output; this flat file is the scan driver's interface.
    logical :: write_diagnostics_dat = .false.

    character(256) :: nml_config_path = "./KIM_config.nml" ! path to the namelist file

    logical :: profiles_in_memory = .false.  ! When true, skip file-based profile reading

    logical :: rescale_density
    real(dp) :: number_density_rescale
    real(dp) :: ion_flr_scale_factor
    real(dp) :: collision_frequency_scale = 1.0_dp

    ! KIM_PROFILES namelist variables
    character(20) :: coord_type = 'auto'           ! 'auto', 'sqrt_psiN', or 'r_eff'
    character(256) :: input_profile_dir = './'     ! Directory for raw input profiles
    character(256) :: equil_file = ''              ! Path to equil_r_q_psi.dat (empty = compute from geqdsk)
    character(256) :: geqdsk_file = ''             ! Path to GEQDSK g-file for equilibrium calculation

    ! Input profile filenames (for sqrt_psiN coordinate input)
    character(256) :: n_input_file = 'n_of_psiN.dat'
    character(256) :: Te_input_file = 'Te_of_psiN.dat'
    character(256) :: Ti_input_file = 'Ti_of_psiN.dat'
    character(256) :: Vz_input_file = 'Vz_of_psiN.dat'

    ! Output profile filenames (in profile_location)
    character(256) :: n_file = 'n.dat'
    character(256) :: Te_file = 'Te.dat'
    character(256) :: Ti_file = 'Ti.dat'
    character(256) :: Vz_file = 'Vz.dat'
    character(256) :: Er_file = 'Er.dat'
    character(256) :: q_file = 'q.dat'

contains

    pure logical function ifunc_model_is_valid(model)

        integer, intent(in) :: model

        ifunc_model_is_valid = model >= IFUNC_MODEL_INHERIT &
            .and. model <= IFUNC_MODEL_ENERGY_MOMENTUM

    end function ifunc_model_is_valid

    pure integer function resolve_ifunc_model(model, legacy_energy_conservation)

        integer, intent(in) :: model
        logical, intent(in) :: legacy_energy_conservation

        if (model == IFUNC_MODEL_INHERIT) then
            resolve_ifunc_model = merge(IFUNC_MODEL_ENERGY, IFUNC_MODEL_NUMBER, &
                legacy_energy_conservation)
        else
            resolve_ifunc_model = model
        end if

    end function resolve_ifunc_model

    pure function ifunc_model_name(model) result(name)

        integer, intent(in) :: model
        character(len=16) :: name

        select case (model)
        case (IFUNC_MODEL_NUMBER)
            name = '0 (N)'
        case (IFUNC_MODEL_ENERGY)
            name = '1 (N+E)'
        case (IFUNC_MODEL_MOMENTUM)
            name = '2 (N+P)'
        case (IFUNC_MODEL_ENERGY_MOMENTUM)
            name = '3 (N+E+P)'
        case (IFUNC_MODEL_INHERIT)
            name = '-1 (inherit)'
        case default
            name = 'invalid'
        end select

    end function ifunc_model_name

    integer function ifunc_model_for_species(species_index)

        integer, intent(in) :: species_index

        if (species_index == 0) then
            ifunc_model_for_species = resolved_electron_ifunc_conservation_model
        else
            ifunc_model_for_species = resolved_ion_ifunc_conservation_model
        end if

    end function ifunc_model_for_species

end module config_m
