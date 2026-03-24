subroutine kim_init

    use species_m, only: read_profiles, allocate_plasma, init_plasma, plasma
    use IO_collection_m, only: initialize_hdf5_output, write_KIM_namelist_to_hdf5
    use config_m, only: hdf5_output, profiles_in_memory
    use profile_input_m, only: prepare_profiles

    implicit none

    call kim_read_config

    if (hdf5_output) then
        call initialize_hdf5_output()
        call write_KIM_namelist_to_hdf5()
    end if

    if (.not. profiles_in_memory) then
        call prepare_profiles()
    end if

    call allocate_plasma
    call init_plasma(plasma)

    if (.not. profiles_in_memory) then
        call read_profiles()
    end if

end subroutine
