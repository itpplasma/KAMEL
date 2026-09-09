program test_periodic_kim_selector
    use control_mod, only: kim_run_type, kim_ion_transport_model, &
        ion_transport_model_id, ION_TRANSPORT_FLR, &
        ION_TRANSPORT_DRIFT_KINETIC, ION_TRANSPORT_INVALID
    implicit none

    if (trim(kim_run_type) /= 'electromagnetic') &
        error stop 'legacy KIM field model changed without explicit selection'
    if (trim(kim_ion_transport_model) /= 'finite_larmor_radius') &
        error stop 'finite-Larmor-radius ion transport is not the production default'
    if (ion_transport_model_id('finite_larmor_radius') /= ION_TRANSPORT_FLR) &
        error stop 'finite-Larmor-radius ion transport selector is broken'
    if (ion_transport_model_id('integral') /= ION_TRANSPORT_INVALID) &
        error stop 'obsolete integral ion transport selector is still accepted'
    if (ion_transport_model_id('drift_kinetic') /= ION_TRANSPORT_DRIFT_KINETIC) &
        error stop 'drift-kinetic ion transport selector is broken'
    if (ion_transport_model_id('unknown') /= ION_TRANSPORT_INVALID) &
        error stop 'unknown ion transport model is not rejected'
    print *, 'periodic KIM selector test passed'
end program test_periodic_kim_selector
