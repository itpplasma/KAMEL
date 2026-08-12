program test_periodic_workflow_validation
    use QLBalance_kinds, only: dp
    use periodic_workflow_validation_m, only: validate_periodic_workflow
    implicit none
    integer :: modes_m(2), modes_n(2)

    modes_m = [6, -6]
    modes_n = [2, 2]
    call validate_periodic_workflow('KIM', 'electrostatic_periodic', 'TimeEvolution', .true., &
        2, modes_m, modes_n, 4.0_dp, 'conductivity', 'drift_kinetic', 'integral', &
        'periodic', 'none')
    call validate_periodic_workflow('KiLCA', 'electrostatic_periodic', 'SingleStep', .true., &
        1, modes_m, modes_n, 0.0_dp, 'conductivity', 'ignored', 'ignored', 'ignored', 'ignored')
    print *, 'periodic workflow validation tests passed'
end program test_periodic_workflow_validation
