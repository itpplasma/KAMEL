program test_flr2_response
    use, intrinsic :: ieee_arithmetic, only: ieee_quiet_nan, ieee_value
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
    complex(dp) :: bpsi_over_bphi(n), phi(n), current(n)
    complex(dp) :: factor_of_phi, expected_phi, expected_current
    real(dp) :: f_electron, d0_electron
    integer :: stat

    call make_constant_background(background)
    bpsi_over_bphi = cmplx(1.0_dp, 0.0_dp, dp)
    options%ion_potential = .false.
    options%ion_current = .false.

    call solve_flr2_response(background, 1, 1, major_radius, background%B0, &
                             bpsi_over_bphi, options, phi, current, stat)
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
    if (maxval(abs(current)) > tol) then
        write(*,*) 'FAIL: aligned potential should cancel parallel current: ', current
        stop 1
    end if

    options%include_potential_in_current = .false.
    call solve_flr2_response(background, 1, 1, major_radius, background%B0, &
                             bpsi_over_bphi, options, phi, current, stat)
    if (stat /= FLR2_RESPONSE_OK) then
        write(*,*) 'FAIL: magnetic-only response status = ', stat
        stop 1
    end if

    f_electron = -e_charge
    d0_electron = -f_electron/(background%B0(1)*major_radius)
    expected_current = cmplx(d0_electron, 0.0_dp, dp)*bpsi_over_bphi(1)
    if (maxval(abs(phi)) > tol) then
        write(*,*) 'FAIL: disabled potential should be zero: ', phi
        stop 1
    end if
    if (maxval(abs(current - expected_current)) > tol*max(1.0_dp, abs(expected_current))) then
        write(*,*) 'FAIL: magnetic-only current mismatch: ', current
        stop 1
    end if

    call solve_flr2_response(background, 1, 1, major_radius, -background%B0, &
                             bpsi_over_bphi, options, phi, current, stat)
    if (stat /= FLR2_RESPONSE_OK &
        .or. maxval(abs(current + expected_current)) &
             > tol*max(1.0_dp, abs(expected_current))) then
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
