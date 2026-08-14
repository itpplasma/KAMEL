program test_periodic_workflow_validation
    use QLBalance_kinds, only: dp
    use periodic_workflow_validation_m, only: validate_periodic_workflow
    use time_evolution, only: write_periodic_workflow_provenance
    use control_mod, only: wave_code, kim_run_type, type_of_run, kim_profiles_from_balance, &
        kim_n_modes, kim_m_list, kim_n_list, kim_electron_transport_model, &
        kim_ion_transport_model, kim_bparallel_source, kim_benchmark_mode, ihdf5IO
    use wave_code_data, only: I_par_toroidal
    use h5mod, only: path2out, h5_mode_groupname
    use KAMEL_hdf5_tools, only: HID_T, h5_init, h5_create, h5_open, h5_close, h5_deinit, &
        h5_obj_exists
    implicit none
    integer :: modes_m(2), modes_n(2)
    integer(HID_T) :: h5id
    integer :: file_unit
    logical :: exists

    modes_m = [6, -6]
    modes_n = [2, 2]
    call validate_periodic_workflow('KIM', 'electrostatic_periodic', 'TimeEvolution', .true., &
        2, modes_m, modes_n, 4.0_dp, 'conductivity', 'drift_kinetic', 'integral', &
        'periodic', 'none')
    call validate_periodic_workflow('KiLCA', 'electrostatic_periodic', 'SingleStep', .true., &
        1, modes_m, modes_n, 0.0_dp, 'conductivity', 'ignored', 'ignored', 'ignored', 'ignored')
    call validate_periodic_workflow('KIM', 'future_nonperiodic_solver', 'SingleStep', .false., &
        0, modes_m, modes_n, 0.0_dp, 'conductivity', 'ignored', 'ignored', 'ignored', 'ignored')

    path2out = 'test_periodic_workflow_provenance.h5'
    h5_mode_groupname = 'multi_mode'
    wave_code = 'KIM'
    kim_run_type = 'electrostatic_periodic'
    type_of_run = 'TimeEvolution'
    kim_profiles_from_balance = .true.
    kim_n_modes = 2
    kim_m_list(1:2) = modes_m
    kim_n_list(1:2) = modes_n
    kim_electron_transport_model = 'drift_kinetic'
    kim_ion_transport_model = 'integral'
    kim_bparallel_source = 'periodic'
    kim_benchmark_mode = 'none'
    I_par_toroidal = 4.0_dp
    ihdf5IO = 1
    call h5_init()
    call h5_create(trim(path2out), h5id)
    call h5_close(h5id)
    call h5_deinit()
    call write_periodic_workflow_provenance()
    call h5_init()
    call h5_open(trim(path2out), h5id)
    call h5_obj_exists(h5id, '/multi_mode/periodic_workflow/mode_m', exists)
    if (.not. exists) error stop 'periodic provenance writer did not create mode metadata'
    call h5_close(h5id)
    call h5_deinit()
    open(newunit=file_unit, file=trim(path2out), status='old')
    close(file_unit, status='delete')
    print *, 'periodic workflow validation tests passed'
end program test_periodic_workflow_validation
