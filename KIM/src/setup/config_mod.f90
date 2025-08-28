module config_m

    implicit none

    ! KIM_CONFIG namelist variables
    integer :: number_of_ion_species ! number of ion species
    logical :: read_species_from_namelist ! read species from namelist or use deuterium plasma
    character(100) :: type_of_run         
    character(100) :: collision_model ! type of collision model
    logical :: artificial_debye_case
    logical :: kernel_debye_case

    ! KIM_IO namelist variables
    character(256) :: profile_location ! path to profile directory
    character(256) :: output_path         ! path to output directory
    logical :: hdf5_input, hdf5_output
    integer :: fdebug, fstatus
    logical :: calculate_asymptotics ! enable/disable asymptotic calculations

    character(256) :: nml_config_path = "./KIM_config.nml" ! path to the namelist file

end module