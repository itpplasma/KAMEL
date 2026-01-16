module config_m

    use KIM_kinds_m, only: dp

    implicit none

    ! KIM_CONFIG namelist variables
    integer :: number_of_ion_species ! number of ion species
    logical :: read_species_from_namelist ! read species from namelist or use deuterium plasma
    character(100) :: type_of_run
    character(100) :: collision_model ! type of collision model
    integer :: artificial_debye_case
    logical :: turn_off_ions ! if true, only the first species (electrons) is considered in calculations
    logical :: turn_off_electrons
    character(100) :: plasma_type ! type of plasma ('H' for hydrogen, 'D' for deuterium)

    ! WKB_DISPERSION namelist variables
    character(20) :: WKB_dispersion_mode = 'KIM'      ! 'KIM' (full Bessel) or 'FLRE' (finite Larmor radius)
    character(20) :: WKB_dispersion_solver = 'Muller' ! 'Muller' or 'ZEAL'
    integer :: WKB_max_tracked_branches = 4           ! Max branches for ZEAL per-branch tracking
    real(dp) :: WKB_branch_search_halfwidth = 1.5_dp  ! Search window half-width for tracking
    real(dp) :: WKB_broad_search_halfwidth = 5.0_dp   ! Broad search window half-width
    integer :: WKB_broad_search_interval = 0          ! Run broad search every N points (0=only at start)
    real(dp) :: WKB_root_tolerance = 1.0d-6           ! |f(z)| tolerance for valid roots
    logical :: WKB_verbose = .true.                   ! ZEAL verbose output

    ! KIM_IO namelist variables
    character(256) :: profile_location ! path to profile directory
    character(256) :: output_path         ! path to output directory
    character(256) :: h5_out_file ! file name of the hdf5 output file
    logical :: hdf5_input, hdf5_output
    integer :: fdebug, fstatus, fdiagnostics
    logical :: calculate_asymptotics ! enable/disable asymptotic calculations

    character(256) :: nml_config_path = "./KIM_config.nml" ! path to the namelist file

    logical :: rescale_density
    real(dp) :: number_density_rescale
    real(dp) :: ion_flr_scale_factor

    ! KIM_PROFILES namelist variables
    character(20) :: coord_type = 'auto'           ! 'auto', 'sqrt_psiN', or 'r_eff'
    character(256) :: input_profile_dir = './'     ! Directory for raw input profiles
    character(256) :: equil_file = ''              ! Path to equilibrium file (empty = compute)

end module
