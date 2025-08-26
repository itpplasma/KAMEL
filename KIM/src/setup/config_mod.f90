module config_m

    implicit none

    character(256) :: profile_location ! path to profile directory
    character(256) :: output_path         ! path to output directory
    character(100) :: type_of_run         
    character(100) :: collision_model ! type of collision model
    logical :: hdf5_input, hdf5_output
    integer :: fdebug, fstatus
    integer :: number_of_ion_species ! number of ion species
    logical :: read_species_from_namelist ! read species from namelist or use deuterium plasma
    logical :: artificial_debye_case
    logical :: kernel_debye_case
    character(256) :: nml_config_path = "./KIM_config.nml" ! path to the namelist file

end module