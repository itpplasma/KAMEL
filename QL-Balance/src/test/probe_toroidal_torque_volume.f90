program probe_toroidal_torque_volume
    use iso_fortran_env, only: dp => real64
    use baseparam_mod, only: pi, rtor
    use grid_mod, only: T_EM_phi_e, T_EM_phi_i, T_tot_phi_e, T_tot_phi_i
    use wave_code_data, only: r

    implicit none

    real(dp) :: full_cylinder

    rtor = 3.0_dp
    allocate (r(3), T_EM_phi_e(3), T_EM_phi_i(3))
    allocate (T_tot_phi_e(1), T_tot_phi_i(1))

    r = [0.0_dp, 0.5_dp, 1.0_dp]
    T_EM_phi_e = 1.0_dp
    T_EM_phi_i = -2.0_dp

    call calculate_total_toroidal_torque(1)

    full_cylinder = 4.0_dp*pi**2*rtor*0.5_dp
    write (*, '(a,es24.16)') 'ql_e=', T_tot_phi_e(1)
    write (*, '(a,es24.16)') 'ql_i=', T_tot_phi_i(1)
    write (*, '(a,es24.16)') 'full_e=', full_cylinder
    write (*, '(a,es24.16)') 'ql_to_full=', T_tot_phi_e(1)/full_cylinder
end program probe_toroidal_torque_volume
