subroutine kim_init

    use species_m, only: read_profiles, allocate_plasma, init_plasma, plasma
    use IO_collection_m, only: initialize_hdf5_output, write_KIM_namelist_to_hdf5
    use config_m, only: hdf5_output
    use profile_input_m, only: prepare_profiles

    implicit none

    call read_config

    if (hdf5_output) then
        call initialize_hdf5_output()
        call write_KIM_namelist_to_hdf5()
    end if

    call prepare_profiles()

    call allocate_plasma
    call init_plasma(plasma)
    call read_profiles()

end subroutine
