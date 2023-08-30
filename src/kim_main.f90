program kim_main

    use plas_parameter
    use config
    use constants
    use setup
    use use_libcerf

    implicit none

    double complex :: z, res

    namelist /KIM_CONFIG/ profile_location, hdf5_input, hdf5_output, &
                          fdebug, fstatus, ispecies, output_path

    namelist /KIM_SETUP/ btor, R0, m_mode, n_mode, Zi, Ai

    open(unit = 77, file = './nmls/KIM_config.nml')
    read(unit = 77, nml = KIM_CONFIG)
    allocate(Zi(ispecies), Ai(ispecies))
    read(unit = 77, nml = KIM_SETUP)
    close(unit = 77)

    write(*,*) ''
    write(*,*) ' - - - - - - - - - - - - - - - - - -'
    write(*,*) 'Configuration namelist'
    write(*,*) '  profile_location = ', trim(profile_location)
    write(*,*) '  output_path         = ', trim(output_path)
    write(*,*) '  hdf5_input       = ', hdf5_input
    write(*,*) '  hdf5_output      = ', hdf5_output
    write(*,*) '  fdebug           = ', fdebug
    write(*,*) '  fstatus          = ', fstatus
    write(*,*) '  ispecies         = ', ispecies
    write(*,*) ' - - - - - - - - - - - - - - - - - -'

    write(*,*) ''
    write(*,*) ' - - - - - - - - - - - - - - - - - -'
    write(*,*) 'Setup namelist'
    write(*,*) '  btor             = ', btor
    write(*,*) '  R0               = ', R0
    write(*,*) '  m_mode           = ', m_mode
    write(*,*) '  n_mode           = ', n_mode
    write(*,*) '  Zi               = ', Zi
    write(*,*) '  Ai               = ', Ai
    write(*,*) ' - - - - - - - - - - - - - - - - - -'

    call read_profiles
    call calc_backs(.true.)
    call kernel_rho_phi

    z = (1d0, 0d0)

    res = cerfc_F(z)
    write(*,*) res

end program