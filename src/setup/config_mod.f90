module config

    implicit none

    character(1024) :: profile_location ! path to profile directory
    character(1024) :: output_path         ! path to output directory
    logical :: hdf5_input, hdf5_output
    integer :: fdebug, fstatus
    integer :: ispecies ! number of ion species
    logical :: artificial_debye_case

end module