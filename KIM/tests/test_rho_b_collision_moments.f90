program test_rho_b_collision_moments
    use KIM_kinds_m, only: dp
    use species_m, only: plasma_t
    use constants_m, only: sol
    use config_m, only: artificial_debye_case, resolved_electron_ifunc_conservation_model, &
        IFUNC_MODEL_ENERGY
    use flr2_fourier_kernel_m, only: core_rho_B_sp
    use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
    implicit none
    integer, parameter :: nh = 64
    type(plasma_t) :: background
    complex(dp) :: a(nh, nh), rhs(nh, 2), moments(0:3, 0:3), expected, actual
    complex(dp), parameter :: imaginary = (0.0_dp, 1.0_dp)
    real(dp), parameter :: x1 = 0.73_dp, x2 = -1.21_dp
    integer :: model, j, info, pivots(nh)
    interface
        subroutine getIfunc_model(x1, x2, conservation_model, symbI)
            import dp
            real(dp), intent(in) :: x1, x2
            integer, intent(in) :: conservation_model
            complex(dp), intent(out) :: symbI(0:3, 0:3)
        end subroutine getIfunc_model
    end interface

    allocate(background%spec(0:0))
    allocate(background%spec(0)%lambda_D(1), background%spec(0)%vT(1), &
        background%spec(0)%omega_c(1), background%spec(0)%nu(1), &
        background%spec(0)%A1(1), background%spec(0)%A2(1), &
        background%spec(0)%I01(1, 0:0), background%spec(0)%I21(1, 0:0), &
        background%spec(0)%I03(1, 0:0))
    background%spec(0)%lambda_D = 1.0_dp
    background%spec(0)%vT = 1.0_dp
    background%spec(0)%omega_c = 1.0_dp
    background%spec(0)%nu = 1.0_dp
    background%spec(0)%A1 = 0.0_dp
    background%spec(0)%A2 = 1.0_dp
    artificial_debye_case = 0
    do model = 0, 1
        ! Independent Gaussian-Hermite resolvent of -i*x2+i*x1*u-C.
        ! Number-only C damps Hermite n at rate n; energy restoration removes n=2.
        a = (0.0_dp, 0.0_dp)
        do j = 1, nh
            a(j, j) = real(j-1, dp)-imaginary*x2
        end do
        if (model == 1) a(3, 3) = -imaginary*x2
        do j = 1, nh-1
            a(j, j+1) = imaginary*x1*sqrt(real(j, dp))
            a(j+1, j) = a(j, j+1)
        end do
        ! Maxwellian*u and Maxwellian*u^3 in normalized Hermite coordinates.
        rhs = (0.0_dp, 0.0_dp)
        rhs(2, 1) = 1.0_dp
        rhs(2, 2) = 3.0_dp
        rhs(4, 2) = sqrt(6.0_dp)
        call zgesv(nh, 2, a, nh, pivots, rhs, nh, info)
        if (info /= 0) error stop 'Hermite oracle solve failed'
        if (.not. all(is_finite_complex(rhs))) error stop 'Nonfinite Hermite oracle solution'
        call getIfunc_model(x1, x2, model, moments)
        resolved_electron_ifunc_conservation_model = model
        if (.not. all(is_finite_complex(moments))) error stop 'Nonfinite collision moments'
        if (abs(rhs(1, 1)-moments(0, 1)) > 1.0e-10_dp) error stop 1
        if (abs(rhs(1, 2)-moments(0, 3)) > 1.0e-10_dp) error stop 2
        background%spec(0)%I01(1, 0) = moments(0, 1)
        background%spec(0)%I21(1, 0) = moments(2, 1)
        background%spec(0)%I03(1, 0) = moments(0, 3)
        expected = -(rhs(1, 1)+0.5_dp*rhs(1, 2))/sol
        actual = core_rho_B_sp(background, 0, 0.0_dp, 0.0_dp, 1)
        if (.not. is_finite_complex(expected)) error stop 'Nonfinite density reference'
        if (.not. is_finite_complex(actual)) error stop 'Nonfinite magnetic density response'
        if (abs(expected) == 0.0_dp) error stop 'Density reference must be nonzero'
        if (.not. is_finite_complex(actual/expected)) error stop 'Nonfinite density ratio'
        if (abs(actual/expected-1.0_dp) > 1.0e-10_dp) error stop 3
        if (model == 0) then
            if (abs(moments(0, 3)-moments(2, 1)) < 0.1_dp) &
                error stop 'Oracle must distinguish number-only density moments'
        else
            if (abs(moments(0, 3)-moments(2, 1)) > 1.0e-10_dp) error stop 4
        end if
    end do
    print *, 'PASS: raw rho-B source moments for number-only and conserving collisions'

    ! Energy conservation at x2 = 0 and small x1, near a resonant surface:
    ! particle conservation gives I01 = 1/(i x1), and I03 = I12 = I21 =
    ! (1 + i x2 I02)/(i x1) = 1/(i x1). getIfunc's I03 is inaccurate here, so
    ! the response must come out right only through I21.
    resolved_electron_ifunc_conservation_model = IFUNC_MODEL_ENERGY
    call getIfunc_model(1.0e-4_dp, 0.0_dp, IFUNC_MODEL_ENERGY, moments)
    background%spec(0)%I01(1, 0) = moments(0, 1)
    background%spec(0)%I21(1, 0) = moments(2, 1)
    background%spec(0)%I03(1, 0) = moments(0, 3)
    expected = -(1.5_dp/(imaginary*1.0e-4_dp))/sol
    actual = core_rho_B_sp(background, 0, 0.0_dp, 0.0_dp, 1)
    print '(a,es10.2,a,es10.2)', '  small x1: |I03 - 1/(i x1)|/|1/(i x1)| =', &
        abs(moments(0, 3)*imaginary*1.0e-4_dp - 1.0_dp), &
        '  rho-B response rel. error =', abs(actual/expected - 1.0_dp)
    if (abs(actual/expected - 1.0_dp) > 1.0e-8_dp) &
        error stop 'rho-B response is inaccurate at small x1 with energy conservation'
    print *, 'PASS: rho-B response uses the accurate moment at small x1'

contains

    elemental logical function is_finite_complex(value) result(ok)
        complex(dp), intent(in) :: value
        ok = ieee_is_finite(real(value, dp)) .and. ieee_is_finite(aimag(value))
    end function is_finite_complex

end program test_rho_b_collision_moments
