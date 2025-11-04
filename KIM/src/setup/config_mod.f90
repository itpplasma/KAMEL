module config_m

    use KIM_kinds_m, only: dp
    
    implicit none

    ! KIM_CONFIG namelist variables
    integer :: number_of_ion_species ! number of ion species
    logical :: read_species_from_namelist ! read species from namelist or use deuterium plasma
    character(100) :: type_of_run         
    character(100) :: collision_model ! type of collision model
    integer :: artificial_debye_case
    logical :: kernel_debye_case
    logical :: turn_off_ions ! if true, only the first species (electrons) is considered in calculations
    logical :: turn_off_electrons
    character(100) :: plasma_type ! type of plasma ('H' for hydrogen, 'D' for deuterium)

    ! KIM_IO namelist variables
    character(256) :: profile_location ! path to profile directory
    character(256) :: output_path         ! path to output directory
    logical :: hdf5_input, hdf5_output
    integer :: fdebug, fstatus, fdiagnostics
    logical :: calculate_asymptotics ! enable/disable asymptotic calculations

    character(256) :: nml_config_path = "./KIM_config.nml" ! path to the namelist file

    logical :: rescale_density
    real(dp) :: number_density_rescale

end module