program test_rho_b_collision_moments
    use KIM_kinds_m, only: dp
    use species_m, only: plasma_t
    use constants_m, only: sol
    use config_m, only: artificial_debye_case
    use flr2_fourier_kernel_m, only: core_rho_B_sp
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
        call getIfunc_model(x1, x2, model, moments)
        if (abs(rhs(1, 1)-moments(0, 1)) > 1.0e-10_dp) error stop 1
        if (abs(rhs(1, 2)-moments(0, 3)) > 1.0e-10_dp) error stop 2
        background%spec(0)%I01(1, 0) = moments(0, 1)
        background%spec(0)%I21(1, 0) = moments(2, 1)
        background%spec(0)%I03(1, 0) = moments(0, 3)
        expected = -(rhs(1, 1)+0.5_dp*rhs(1, 2))/sol
        actual = core_rho_B_sp(background, 0, 0.0_dp, 0.0_dp, 1)
        if (abs(actual/expected-1.0_dp) > 1.0e-10_dp) error stop 3
        if (model == 0) then
            if (abs(moments(0, 3)-moments(2, 1)) < 0.1_dp) &
                error stop 'Oracle must distinguish number-only density moments'
        else
            if (abs(moments(0, 3)-moments(2, 1)) > 1.0e-10_dp) error stop 4
        end if
    end do
    print *, 'PASS: raw rho-B source moments for number-only and conserving collisions'
end program test_rho_b_collision_moments
