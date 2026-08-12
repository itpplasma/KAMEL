program test_periodic_ion_tensor
    use KIM_kinds_m, only: dp
    use kim_solver_m, only: kim_results_t
    use rt_electrostatic_periodic_m, only: compute_periodic_ion_tensor
    implicit none

    type(kim_results_t) :: result
    complex(dp) :: fields(3)
    real(dp) :: tensor_bpar(2,2), tensor_no_bpar(2,2)
    real(dp), parameter :: ks = 0.2_dp, kr = 0.0_dp, kpar = 0.15_dp
    real(dp), parameter :: vti = 2.0e7_dp, nui = 4.0e5_dp
    real(dp), parameter :: omega_ci = 1.0e7_dp, omega_mode = 0.0_dp
    real(dp), parameter :: om_e = 2.0e5_dp, b0 = 1.0e4_dp

    allocate(result%D_ion(2,2,3))
    result%D_ion = 0.0_dp
    if (size(result%D_ion, 1) /= 2 .or. size(result%D_ion, 2) /= 2) &
        error stop 'periodic result ion tensor shape'

    fields = [(0.0_dp, 0.0_dp), (0.0_dp, 0.0_dp), (1.0_dp, 0.0_dp)]
    call compute_periodic_ion_tensor(fields, ks, kr, kpar, vti, nui, omega_ci, &
                                      omega_mode, om_e, b0, tensor_bpar)
    fields(3) = (0.0_dp, 0.0_dp)
    call compute_periodic_ion_tensor(fields, ks, kr, kpar, vti, nui, omega_ci, &
                                      omega_mode, om_e, b0, tensor_no_bpar)

    if (maxval(abs(tensor_bpar - tensor_no_bpar)) <= 1.0e-30_dp) &
        error stop 'Bparallel did not contribute to periodic ion tensor'
    if (min(tensor_bpar(1,1), tensor_bpar(2,2)) < -1.0e-12_dp) &
        error stop 'periodic ion tensor lost positive diagonal'
    print *, 'periodic ion tensor contract tests passed'
end program test_periodic_ion_tensor
