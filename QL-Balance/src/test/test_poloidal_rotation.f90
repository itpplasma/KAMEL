program test_poloidal_rotation
    !! Neoclassical poloidal rotation V_pol = C dT_i/dr (poloidal_rotation_m) and the
    !! effective poloidal velocity used for E0r (rhs_balance_m), Markl 2023 eq. (26).
    use QLBalance_kinds, only: dp
    use poloidal_rotation_m, only: init_poloidal_rotation_coefficient, poloidal_velocity, &
                                   poloidal_rotation_ready, poloidal_rotation_coefficient, &
                                   reset_poloidal_rotation
    use control_mod, only: poloidal_rotation_model
    use wave_code_data, only: Vth
    use rhs_balance_m, only: effective_poloidal_velocity, store_poloidal_velocity
    implicit none

    real(dp), parameter :: tol = 1.0e-14_dp
    real(dp) :: vpol_in(4), dTi0(4), v(4)

    vpol_in = [2.0_dp, -4.0_dp, 0.0_dp, 5.0_dp]
    dTi0 = [1.0_dp, 2.0_dp, 0.0_dp, -5.0_dp]

    ! Module: C = Vth / dTi at the reference state; zero where dTi = 0.
    call reset_poloidal_rotation()
    if (poloidal_rotation_ready()) error stop 'coefficient must start unset'
    call init_poloidal_rotation_coefficient(vpol_in, dTi0)
    if (any(abs(poloidal_rotation_coefficient() - [2.0_dp, -2.0_dp, 0.0_dp, -1.0_dp]) > tol)) &
        error stop 'wrong coefficient'
    ! Reproduces the input at the reference state ...
    if (any(abs(poloidal_velocity(dTi0) - vpol_in) > tol)) error stop 'input not reproduced'
    ! ... and is proportional to dT_i/dr afterwards (k, geometry fixed).
    if (any(abs(poloidal_velocity(0.5_dp*dTi0) - 0.5_dp*vpol_in) > tol)) &
        error stop 'not proportional to dT_i/dr'

    ! Helper, 'fixed': the input Vth regardless of dT_i/dr.
    allocate (Vth(4))
    Vth = vpol_in
    call reset_poloidal_rotation()
    poloidal_rotation_model = 'fixed'
    v = effective_poloidal_velocity(3.0_dp*dTi0)
    if (any(abs(v - vpol_in) > tol)) error stop 'fixed model must keep Vth'
    if (poloidal_rotation_ready()) error stop 'fixed model must not set the coefficient'

    ! Helper, 'neoclassical': first call (initial state) sets C from Vth and reproduces it;
    ! later states scale with dT_i/dr; the stored Vth follows.
    poloidal_rotation_model = 'neoclassical'
    v = effective_poloidal_velocity(dTi0)
    if (any(abs(v - vpol_in) > tol)) error stop 'initial state must reproduce Vth'
    call store_poloidal_velocity(v)
    v = effective_poloidal_velocity(2.0_dp*dTi0)
    if (any(abs(v - 2.0_dp*vpol_in) > tol)) error stop 'V_pol must follow dT_i/dr'
    call store_poloidal_velocity(v)
    if (any(abs(Vth - 2.0_dp*vpol_in) > tol)) error stop 'stored Vth must follow'
    ! A restart from the stored state recovers the same coefficient.
    call reset_poloidal_rotation()
    v = effective_poloidal_velocity(2.0_dp*dTi0)
    if (any(abs(poloidal_rotation_coefficient() - [2.0_dp, -2.0_dp, 0.0_dp, -1.0_dp]) > tol)) &
        error stop 'restart must recover the coefficient'

    write (*, *) 'test_poloidal_rotation passed'
end program test_poloidal_rotation
