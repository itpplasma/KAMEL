program test_periodic_kim_selector
    use control_mod, only: kim_run_type
    implicit none

    if (trim(kim_run_type) /= 'electrostatic_periodic') &
        error stop 'periodic KIM coupling is not the production default'
    print *, 'periodic KIM selector test passed'
end program test_periodic_kim_selector
