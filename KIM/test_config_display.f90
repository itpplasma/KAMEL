program test_config_display
    use config_display, only: display_kim_configuration
    use config
    use setup
    use grid
    use species
    
    implicit none
    
    ! Set some test values
    profile_location = './profiles/'
    output_path = './output/'
    hdf5_input = .true.
    hdf5_output = .false.
    fdebug = 1
    fstatus = 1
    number_of_ion_species = 2
    type_of_run = 'electrostatic'
    collision_model = 'Krook'
    collisions_off = .false.
    artificial_debye_case = .false.
    
    btor = -17977.413d0
    r0 = 165.0d0
    r_plas = 67.0d0
    m_mode = -6
    n_mode = 2
    omega = 0.0d0
    eps_reg = 0.01d0
    
    k_space_dim = 100
    l_space_dim = 1000
    grid_spacing = 3
    Larmor_skip_factor = 5.0d0
    gauss_int_nodes_Nx = 6
    gauss_int_nodes_Nxp = 7
    gauss_int_nodes_Ntheta = 7
    
    ! Display the configuration
    call display_kim_configuration()
    
end program test_config_display