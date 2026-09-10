program test_periodic_workflow_validation
    use QLBalance_kinds, only: dp
    use periodic_workflow_validation_m, only: validate_periodic_workflow, periodic_benchmark_enabled
    use time_evolution, only: write_periodic_workflow_provenance
    use control_mod, only: wave_code, kim_run_type, type_of_run, kim_profiles_from_balance, &
        kim_n_modes, kim_m_list, kim_n_list, kim_electron_transport_model, &
        kim_ion_transport_model, kim_bparallel_source, kim_benchmark_mode, ihdf5IO, &
        kim_config_sha256
    use wave_code_data, only: I_par_toroidal
    use h5mod, only: path2out, h5_mode_groupname
    use KAMEL_hdf5_tools, only: HID_T, h5_init, h5_create, h5_open, h5_close, h5_deinit, &
        h5_obj_exists
    implicit none
    integer :: modes_m(2), modes_n(2)
    integer(HID_T) :: h5id
    integer :: file_unit
    character(32) :: argument
    logical :: exists

    modes_m = [6, -6]
    modes_n = [2, 2]
    call get_command_argument(1, argument)
    select case (trim(argument))
    case ('duplicate')
        call validate_periodic_workflow('KIM', 'electrostatic_periodic', 'TimeEvolution', .true., &
            2, [6, 6], [2, 2], 4.0_dp, 'conductivity', 'drift_kinetic', &
            'finite_larmor_radius', 'disabled', 'none')
        error stop 'duplicate modes were accepted'
    case ('shape')
        call validate_periodic_workflow('KIM', 'electrostatic_periodic', 'TimeEvolution', .true., &
            2, modes_m, [2], 4.0_dp, 'conductivity', 'drift_kinetic', &
            'finite_larmor_radius', 'disabled', 'none')
        error stop 'mismatched mode arrays were accepted'
    case ('bparallel')
        call validate_periodic_workflow('KIM', 'electrostatic_periodic', 'TimeEvolution', .true., &
            2, modes_m, modes_n, 4.0_dp, 'conductivity', 'drift_kinetic', &
                                        'finite_larmor_radius', 'self_consistent', 'none')
        error stop 'unsupported Bparallel source was accepted'
    end select
    call validate_periodic_workflow('KIM', 'electrostatic_periodic', 'TimeEvolution', .true., &
        2, modes_m, modes_n, 4.0_dp, 'conductivity', 'drift_kinetic', &
        'finite_larmor_radius', &
        'disabled', 'none')
    call validate_periodic_workflow('KIM', 'electrostatic_periodic', 'SingleStep', .true., &
        2, modes_m, modes_n, 0.0_dp, 'conductivity', 'drift_kinetic', 'drift_kinetic', &
        'disabled', 'none')
    call validate_periodic_workflow('KIM', 'electrostatic_periodic', 'SingleStep', .true., &
                                    2, modes_m, modes_n, 0.0_dp, 'conductivity', 'drift_kinetic', &
                                    'finite_larmor_radius', 'prescribed_zero_mode', 'none')
    call validate_periodic_workflow('KiLCA', 'electrostatic_periodic', 'SingleStep', .true., &
        1, modes_m, modes_n, 0.0_dp, 'conductivity', 'ignored', 'ignored', 'ignored', 'ignored')
    call validate_periodic_workflow('KIM', 'future_nonperiodic_solver', 'SingleStep', .false., &
        0, modes_m, modes_n, 0.0_dp, 'conductivity', 'ignored', 'ignored', 'ignored', 'ignored')
    if (periodic_benchmark_enabled('none')) error stop 'none benchmark mode was enabled'
    if (.not. periodic_benchmark_enabled('drift_kinetic_limit')) &
        error stop 'drift-kinetic benchmark mode was disabled'

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
    kim_ion_transport_model = 'finite_larmor_radius'
    kim_bparallel_source = 'disabled'
    kim_benchmark_mode = 'none'
    kim_config_sha256 = 'test-config-digest'
    I_par_toroidal = 4.0_dp
    ihdf5IO = 1
    call h5_init()
    call h5_create(trim(path2out), h5id)
    call h5_close(h5id)
    call h5_deinit()
    call write_periodic_workflow_provenance()
    I_par_toroidal = 5.0_dp
    call write_periodic_workflow_provenance()
    call h5_init()
    call h5_open(trim(path2out), h5id)
    call h5_obj_exists(h5id, '/multi_mode/periodic_workflow/mode_m', exists)
    if (.not. exists) error stop 'periodic provenance writer did not create mode metadata'
    call h5_obj_exists(h5id, '/multi_mode/periodic_workflow/kim_config_sha256', exists)
    if (.not. exists) error stop 'periodic provenance writer omitted the KIM config digest'
    call h5_obj_exists(h5id, '/multi_mode/periodic_workflow/periodic_kmax_scale', exists)
    if (.not. exists) error stop 'periodic provenance writer omitted resolved KIM settings'
    call h5_close(h5id)
    call h5_deinit()
    open(newunit=file_unit, file=trim(path2out), status='old')
    close(file_unit, status='delete')
    print *, 'periodic workflow validation tests passed'
end program test_periodic_workflow_validation
