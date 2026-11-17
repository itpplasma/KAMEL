program test_flr2_response
    use, intrinsic :: ieee_arithmetic, only: ieee_is_finite, ieee_quiet_nan, ieee_value
    use KIM_kinds_m, only: dp
    use constants_m, only: e_charge
    use species_m, only: plasma_t
    use flr2_response_m, only: flr2_options_t, solve_flr2_response, &
                              FLR2_RESPONSE_BAD_INPUT, FLR2_RESPONSE_OK
    implicit none

    integer, parameter :: n = 5
    real(dp), parameter :: tol = 1.0e-11_dp
    real(dp), parameter :: major_radius = 100.0_dp
    type(plasma_t) :: background
    type(flr2_options_t) :: options
    complex(dp) :: bpsi_over_bphi(n), phi(n), current(n), jpar_e(n), jpar_i(n)
    complex(dp) :: current_with_ions(n)
    complex(dp) :: phi_reference(n)
    complex(dp) :: principal_coefficient(n)
    complex(dp) :: factor_of_phi, expected_phi, expected_current
    real(dp) :: f_electron, d0_electron
    real(dp) :: current_scale, closure_scale
    integer :: i, stat

    call make_constant_background(background)
    bpsi_over_bphi = cmplx(1.0_dp, 0.0_dp, dp)
    options%ion_potential = .false.
    options%ion_current = .false.

    call solve_flr2_response(background, 1, 1, major_radius, background%B0, &
                             bpsi_over_bphi, options, phi, current, stat, &
                             jpar_e, jpar_i, principal_coefficient)
    if (stat /= FLR2_RESPONSE_OK) then
        write(*,*) 'FAIL: response status = ', stat
        stop 1
    end if

    factor_of_phi = cmplx(0.0_dp, 1.0_dp, dp)*(1.0_dp + 2.0_dp) &
                    /(2.0_dp*(-10.0_dp))
    expected_phi = -bpsi_over_bphi(1)/factor_of_phi
    if (maxval(abs(phi - expected_phi)) > tol) then
        write(*,*) 'FAIL: aligned potential mismatch: ', phi
        stop 1
    end if
    closure_scale = max(maxval(abs(current)), maxval(abs(jpar_e)), &
                        maxval(abs(jpar_i)), 1.0e-14_dp)
    if (maxval(abs(jpar_i)) > tol .or. &
        maxval(abs(current-jpar_e-jpar_i)) > 1.0e-8_dp*closure_scale) then
        write(*,*) 'FAIL: species current attribution mismatch: ', jpar_e, jpar_i
        stop 1
    end if
    if (.not. all(ieee_is_finite(real(principal_coefficient))) .or. &
        .not. all(ieee_is_finite(aimag(principal_coefficient)))) then
        write(*,*) 'FAIL: principal coefficient profile was not finite'
        stop 1
    end if
    phi_reference = phi

    options%include_potential_in_current = .false.
    call solve_flr2_response(background, 1, 1, major_radius, background%B0, &
                             bpsi_over_bphi, options, phi, current, stat, &
                             jpar_e, jpar_i)
    if (stat /= FLR2_RESPONSE_OK) then
        write(*,*) 'FAIL: magnetic-only response status = ', stat
        stop 1
    end if

    f_electron = -e_charge
    d0_electron = -f_electron/(background%B0(1)*major_radius)
    expected_current = cmplx(d0_electron, 0.0_dp, dp)*bpsi_over_bphi(1)
    current_scale = abs(expected_current)
    if (maxval(abs(phi)) > tol) then
        write(*,*) 'FAIL: disabled potential should be zero: ', phi
        stop 1
    end if
    if (maxval(abs(current - expected_current)) > &
        1.0e-8_dp*max(current_scale, 1.0e-14_dp)) then
        write(*,*) 'FAIL: magnetic-only current mismatch: ', current
        stop 1
    end if
    closure_scale = max(maxval(abs(current)), maxval(abs(jpar_e)), &
                        maxval(abs(jpar_i)), 1.0e-14_dp)
    if (maxval(abs(current-jpar_e-jpar_i)) > 1.0e-8_dp*closure_scale) then
        write(*,*) 'FAIL: magnetic-only species current attribution mismatch: ', &
                   current, jpar_e, jpar_i
        stop 1
    end if
    if (maxval(abs(current)) <= max(1.0e-3_dp*current_scale, 1.0e-14_dp)) then
        write(*,*) 'FAIL: magnetic-only current was numerically zero: ', current
        stop 1
    end if

    options%include_potential_in_current = .true.
    options%electron_current = .false.
    call solve_flr2_response(background, 1, 1, major_radius, background%B0, &
                             bpsi_over_bphi, options, phi, current, stat, &
                             jpar_e, jpar_i)
    if (stat /= FLR2_RESPONSE_OK .or. maxval(abs(phi-phi_reference)) > tol) then
        write(*,*) 'FAIL: electron current switch changed solved Phi: ', stat, phi
        stop 1
    end if
    closure_scale = max(maxval(abs(current)), maxval(abs(jpar_e)), &
                        maxval(abs(jpar_i)), 1.0e-14_dp)
    if (maxval(abs(jpar_e)) > tol .or. &
        maxval(abs(current-jpar_e-jpar_i)) > 1.0e-8_dp*closure_scale) then
        write(*,*) 'FAIL: electron current switch attribution mismatch: ', jpar_e, jpar_i
        stop 1
    end if

    options%electron_current = .true.
    options%ion_current = .false.
    call solve_flr2_response(background, 1, 1, major_radius, background%B0, &
                             bpsi_over_bphi, options, phi, current, stat, &
                             jpar_e, jpar_i)
    if (stat /= FLR2_RESPONSE_OK .or. maxval(abs(phi-phi_reference)) > tol) then
        write(*,*) 'FAIL: ion current switch changed solved Phi: ', stat, phi
        stop 1
    end if
    closure_scale = max(maxval(abs(current)), maxval(abs(jpar_e)), &
                        maxval(abs(jpar_i)), 1.0e-14_dp)
    if (maxval(abs(jpar_i)) > tol .or. &
        maxval(abs(current-jpar_e-jpar_i)) > 1.0e-8_dp*closure_scale) then
        write(*,*) 'FAIL: ion current switch attribution mismatch: ', jpar_e, jpar_i
        stop 1
    end if
    if (maxval(abs(jpar_i)) > max(1.0e-3_dp*current_scale, 1.0e-14_dp)) then
        write(*,*) 'FAIL: ion-disabled attribution was not zero: ', jpar_i
        stop 1
    end if

    ! Use a nonconstant magnetic drive for the ion attribution check so the
    ! second-derivative current channel is observable.
    do i = 1, n
        bpsi_over_bphi(i) = cmplx(1.0_dp + 0.1_dp*real(i - 1, dp), 0.0_dp, dp)
    end do
    call solve_flr2_response(background, 1, 1, major_radius, background%B0, &
                             bpsi_over_bphi, options, phi, current, stat, &
                             jpar_e, jpar_i)
    if (stat /= FLR2_RESPONSE_OK) then
        write(*,*) 'FAIL: nonconstant-drive ion baseline status = ', stat
        stop 1
    end if
    phi_reference = phi

    options%ion_current = .true.
    options%include_potential_in_current = .true.
    call solve_flr2_response(background, 1, 1, major_radius, background%B0, &
                             bpsi_over_bphi, options, phi, current, stat, &
                             jpar_e, jpar_i)
    if (stat /= FLR2_RESPONSE_OK .or. maxval(abs(phi-phi_reference)) > tol) then
        write(*,*) 'FAIL: enabling ion current changed solved Phi: ', stat, phi
        stop 1
    end if
    closure_scale = max(maxval(abs(current)), maxval(abs(jpar_e)), &
                        maxval(abs(jpar_i)), 1.0e-14_dp)
    if (maxval(abs(current-jpar_e-jpar_i)) > 1.0e-8_dp*closure_scale) then
        write(*,*) 'FAIL: ion-enabled species current attribution mismatch: ', &
                   current, jpar_e, jpar_i
        stop 1
    end if
    if (maxval(abs(jpar_i)) <= max(1.0e-3_dp*current_scale, 1.0e-14_dp)) then
        write(*,*) 'FAIL: ion-enabled baseline was numerically zero: ', jpar_i
        stop 1
    end if
    current_with_ions = current

    call solve_flr2_response(background, 1, 1, major_radius, -background%B0, &
                             bpsi_over_bphi, options, phi, current, stat)
    if (stat /= FLR2_RESPONSE_OK &
        .or. maxval(abs(current + current_with_ions)) &
             > 1.0e-8_dp*max(maxval(abs(current_with_ions)), 1.0e-14_dp)) then
        write(*,*) 'FAIL: signed toroidal field did not reverse current: ', current
        stop 1
    end if

    background%spec(0)%n(3) = ieee_value(0.0_dp, ieee_quiet_nan)
    call solve_flr2_response(background, 1, 1, major_radius, background%B0, &
                             bpsi_over_bphi, options, phi, current, stat)
    if (stat /= FLR2_RESPONSE_BAD_INPUT) then
        write(*,*) 'FAIL: non-finite background was accepted, status = ', stat
        stop 1
    end if

    write(*,*) 'FLR2 response test PASSED'

contains

    subroutine make_constant_background(plasma)
        type(plasma_t), intent(out) :: plasma
        integer :: sp

        plasma%n_species = 2
        plasma%grid_size = n
        allocate(plasma%r_grid(n), plasma%q(n), plasma%B0(n), plasma%Er(n), plasma%om_E(n))
        allocate(plasma%spec(0:1))

        plasma%r_grid = [0.0_dp, 0.5_dp, 1.5_dp, 3.0_dp, 5.0_dp]
        plasma%q = 2.0_dp
        plasma%B0 = 1.0_dp
        plasma%Er = 10.0_dp
        plasma%om_E = 3.0_dp

        do sp = 0, 1
            allocate(plasma%spec(sp)%n(n), plasma%spec(sp)%vT(n))
            allocate(plasma%spec(sp)%nu(n), plasma%spec(sp)%rho_L(n))
            allocate(plasma%spec(sp)%A1(n), plasma%spec(sp)%A2(n))
            allocate(plasma%spec(sp)%I11(n, 0:0), plasma%spec(sp)%I13(n, 0:0))
            plasma%spec(sp)%n = 1.0_dp
            plasma%spec(sp)%vT = 1.0_dp
            plasma%spec(sp)%nu = 1.0_dp
            plasma%spec(sp)%rho_L = 0.1_dp
            plasma%spec(sp)%A1 = 1.0_dp
            plasma%spec(sp)%A2 = 0.0_dp
            plasma%spec(sp)%I11(:, 0) = cmplx(1.0_dp, 0.0_dp, dp)
            plasma%spec(sp)%I13(:, 0) = cmplx(0.0_dp, 0.0_dp, dp)
        end do
        plasma%spec(0)%Zspec = -1
        plasma%spec(1)%Zspec = 1
    end subroutine make_constant_background

end program test_flr2_response
