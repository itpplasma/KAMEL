subroutine calculate_total_toroidal_torque

    use iso_fortran_env, only: dp => real64
    use baseparam_mod, only: pi
    use grid_mod, only: T_EM_phi_e, T_EM_phi_i, T_tot_phi_e, T_tot_phi_i
    use integration, only: simpson_nonequi
    use wave_code_data, only: r

    implicit none

    integer :: n

    n = size(r)
    if (mod(n - 1, 2) /= 0) then
        n = n - 1 ! Ensure n is even for Simpson's rule
    end if

    call simpson_nonequi(T_tot_phi_e, r(1:n), r(1:n) * T_EM_phi_e(1:n))
    call simpson_nonequi(T_tot_phi_i, r(1:n), r(1:n) * T_EM_phi_i(1:n))

    T_tot_phi_e = T_tot_phi_e * 2.0_dp * pi
    T_tot_phi_i = T_tot_phi_i * 2.0_dp * pi

end subroutine calculate_total_toroidal_torque
