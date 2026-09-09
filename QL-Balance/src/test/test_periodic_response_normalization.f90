program test_periodic_response_normalization
    !! Exercise normalization of a complete response, not just its scalar factor.
    use, intrinsic :: ieee_arithmetic, only: ieee_value, ieee_quiet_nan, ieee_positive_inf, &
        ieee_is_finite
    use QLBalance_kinds, only: dp
    use kim_solver_m, only: kim_results_t
    use kim_wave_code_adapter_m, only: kim_normalize_periodic_response
    use control_mod, only: kim_current_floor, kim_current_max_scale, kim_current_relaxation
    implicit none

    integer, parameter :: nr = 4
    real(dp), parameter :: c_light = 2.99792458e10_dp
    real(dp), parameter :: target = 2.0_dp*acos(-1.0_dp)*1.6_dp/c_light
    complex(dp), parameter :: expected_unit = (1.6_dp, 1.6_dp)
    type(kim_results_t) :: initial, response, doubled
    complex(dp) :: unit_current, scale, first_scale
    real(dp) :: invalid_values(2)
    integer :: failures, status, invalid_kind, component

    failures = 0
    kim_current_floor = 1.0e-30_dp
    kim_current_max_scale = 1.0e12_dp
    kim_current_relaxation = 1.0_dp
    call make_response(initial)
    response = initial
    call kim_normalize_periodic_response(response, target, unit_current, scale, status)
    call require(status == 0, 'valid normalization status')
    call require(abs(unit_current-expected_unit) < 1.0e-13_dp, &
        'analytic complex current on clipped trusted interval [1.1,2.1]')
    call require(abs(scale-cmplx(0.5_dp, -0.5_dp, dp)) < 1.0e-13_dp, &
        'complex target amplitude includes c and 2 pi')
    call check_scaled(initial, response, scale)
    call require(abs(2.0_dp*acos(-1.0_dp)*scale*unit_current/c_light-target) &
        < 1.0e-13_dp*target, 'achieved complex current reaches real target')
    first_scale = scale
    doubled = initial
    call kim_normalize_periodic_response(doubled, 2.0_dp*target, unit_current, scale, status)
    call require(status == 0, 'double target status')
    call require(abs(scale-2.0_dp*first_scale) < 1.0e-13_dp, 'target doubling doubles amplitude')
    call require(maxval(abs(linear_fields(doubled)-2.0_dp*linear_fields(response))) &
        < 1.0e-12_dp, 'target doubling doubles every linear response quantity')
    call require(maxval(abs(doubled%D_ion-4.0_dp*response%D_ion)) < 1.0e-12_dp, &
        'target doubling quadruples ion tensor exactly once')

    kim_current_relaxation = 0.25_dp
    response = initial
    call kim_normalize_periodic_response(response, target, unit_current, scale, status)
    call require(status == 0, 'relaxed normalization status')
    call require(abs(scale-(0.75_dp+0.25_dp*first_scale)) < 1.0e-13_dp, &
        'relaxation blends unit amplitude and complex target amplitude')
    call check_scaled(initial, response, scale)
    kim_current_relaxation = 1.0_dp

    response = initial
    call kim_normalize_periodic_response(response, 0.0_dp, unit_current, scale, status)
    call require(all(linear_fields(response) == linear_fields(initial)) &
        .and. all(response%D_ion == initial%D_ion), 'disabled target preserves response')
    call check_background(initial, response)

    response = initial
    response%jpar = (0.0_dp, 0.0_dp)
    response%jpar_e = (0.0_dp, 0.0_dp)
    response%jpar_i = (0.0_dp, 0.0_dp)
    call kim_normalize_periodic_response(response, target, unit_current, scale, status)
    call require(status == 2, 'zero trusted current activates floor guard')
    call check_suppressed(initial, response, scale)

    kim_current_max_scale = 0.1_dp
    response = initial
    call kim_normalize_periodic_response(response, target, unit_current, scale, status)
    call require(status == 3, 'excessive scale activates guard')
    call check_suppressed(initial, response, scale)
    kim_current_max_scale = 1.0e12_dp

    invalid_values = [ieee_value(1.0_dp, ieee_quiet_nan), &
        ieee_value(1.0_dp, ieee_positive_inf)]
    do invalid_kind = 1, size(invalid_values)
        do component = 1, 7
            response = initial
            select case (component)
            case (1)
                response%Es(2) = cmplx(invalid_values(invalid_kind), 0.0_dp, dp)
            case (2)
                response%Phi(2) = cmplx(0.0_dp, invalid_values(invalid_kind), dp)
            case (3)
                response%jrad(2) = cmplx(invalid_values(invalid_kind), 0.0_dp, dp)
            case (4)
                response%D_ion(1,2,2) = invalid_values(invalid_kind)
            case (5)
                response%jpar(2) = cmplx(invalid_values(invalid_kind), 0.0_dp, dp)
            case (6)
                response%jpar_e(2) = cmplx(0.0_dp, invalid_values(invalid_kind), dp)
            case (7)
                response%jpar_i(2) = cmplx(invalid_values(invalid_kind), 0.0_dp, dp)
            end select
            call kim_normalize_periodic_response(response, target, unit_current, scale, status)
            call require(status == 3, 'nonfinite linear quantity or tensor activates guard')
            call check_suppressed(initial, response, scale)
        end do
    end do

    ! Finite input can overflow only after multiplication; validate outputs too.
    response = initial
    response%Es = cmplx(0.75_dp*huge(1.0_dp), 0.0_dp, dp)
    call kim_normalize_periodic_response(response, 4.0_dp*target, unit_current, scale, status)
    call require(status == 3, 'finite field multiplication overflow activates guard')
    call check_suppressed(initial, response, scale)
    response = initial
    response%D_ion = 0.75_dp*huge(1.0_dp)
    call kim_normalize_periodic_response(response, 4.0_dp*target, unit_current, scale, status)
    call require(status == 3, 'finite tensor multiplication overflow activates guard')
    call check_suppressed(initial, response, scale)

    if (failures /= 0) then
        print *, 'failed response normalization assertions: ', failures
        error stop 'periodic response normalization failed'
    end if
    print *, 'complete periodic response normalization passed'

contains

    subroutine require(condition, message)
        logical, intent(in) :: condition
        character(len=*), intent(in) :: message
        if (condition) return
        failures = failures+1
        print *, 'FAIL: ', message
    end subroutine require

    subroutine make_response(res)
        type(kim_results_t), intent(out) :: res
        integer :: i
        res%m = -6
        res%n = 2
        res%r_field = [1.0_dp, 1.4_dp, 1.8_dp, 2.2_dp]
        res%r_resonance = 1.6_dp
        res%dx_asis = 0.5_dp
        res%dx_transition = 0.1_dp
        res%Es = [(cmplx(real(i,dp), 0.3_dp, dp), i=1,nr)]
        res%Ep = 2.0_dp*res%Es
        res%Er = 3.0_dp*res%Es
        res%Etheta = 4.0_dp*res%Es
        res%Ez = 5.0_dp*res%Es
        res%Br = 6.0_dp*res%Es
        res%Bparallel = 7.0_dp*res%Es
        res%jrad = 8.0_dp*res%Es
        res%Phi = 9.0_dp*res%Es
        allocate(res%jpar(nr), res%jpar_e(nr), res%jpar_i(nr), res%D_ion(2,2,nr))
        res%jpar = (1.0_dp, 1.0_dp)
        res%jpar_e = (0.25_dp, 0.25_dp)
        res%jpar_i = (0.75_dp, 0.75_dp)
        do i = 1, nr
            res%D_ion(:,:,i) = real(i,dp)*reshape([2.0_dp, 0.2_dp, 0.2_dp, 3.0_dp], [2,2])
        end do
        res%r_plasma = res%r_field
        res%kp = 0.01_dp*res%r_field
        res%ks = 0.02_dp*res%r_field
        res%om_E = 3.0_dp*res%r_field
        res%nu_e = 4.0_dp*res%r_field
        res%nu_i = 5.0_dp*res%r_field
        res%B0 = 6.0_dp*res%r_field
        res%B0z = 7.0_dp*res%r_field
        res%B0th = 8.0_dp*res%r_field
    end subroutine make_response

    function linear_fields(res) result(fields)
        type(kim_results_t), intent(in) :: res
        complex(dp) :: fields(nr,12)
        fields(:,1) = res%Es
        fields(:,2) = res%Ep
        fields(:,3) = res%Er
        fields(:,4) = res%Etheta
        fields(:,5) = res%Ez
        fields(:,6) = res%Br
        fields(:,7) = res%Bparallel
        fields(:,8) = res%jpar
        fields(:,9) = res%jpar_e
        fields(:,10) = res%jpar_i
        fields(:,11) = res%jrad
        fields(:,12) = res%Phi
    end function linear_fields

    subroutine check_background(before, after)
        type(kim_results_t), intent(in) :: before, after
        call require(all(before%r_field == after%r_field) &
            .and. all(before%r_plasma == after%r_plasma), 'normalization preserves radial grids')
        call require(before%r_resonance == after%r_resonance &
            .and. before%dx_asis == after%dx_asis &
            .and. before%dx_transition == after%dx_transition, &
            'normalization preserves core geometry')
        call require(all(before%kp == after%kp) .and. all(before%ks == after%ks) &
            .and. all(before%om_E == after%om_E), &
            'normalization preserves wave/background frequencies')
        call require(all(before%nu_e == after%nu_e) .and. all(before%nu_i == after%nu_i), &
            'normalization preserves both collision profiles')
        call require(all(before%B0 == after%B0) .and. all(before%B0z == after%B0z) &
            .and. all(before%B0th == after%B0th), 'normalization preserves magnetic background')
        call require(before%m == after%m .and. before%n == after%n, 'normalization preserves mode')
    end subroutine check_background

    subroutine check_scaled(before, after, amplitude)
        type(kim_results_t), intent(in) :: before, after
        complex(dp), intent(in) :: amplitude
        call require(maxval(abs(linear_fields(after)-amplitude*linear_fields(before))) &
            < 1.0e-12_dp, 'all linear fields/current/Phi scale once by common complex amplitude')
        call require(maxval(abs(after%D_ion-abs(amplitude)**2*before%D_ion)) &
            < 1.0e-12_dp, 'ion tensor scales once by squared magnitude')
        call check_background(before, after)
    end subroutine check_scaled

    subroutine check_suppressed(before, after, amplitude)
        type(kim_results_t), intent(in) :: before, after
        complex(dp), intent(in) :: amplitude
        complex(dp) :: fields(nr,12)
        fields = linear_fields(after)
        call require(amplitude == cmplx(0.0_dp, 0.0_dp, dp), 'guard reports zero amplitude')
        call require(all(ieee_is_finite(real(fields,dp))) &
            .and. all(ieee_is_finite(aimag(fields))) &
            .and. all(ieee_is_finite(after%D_ion)), 'guard leaves finite response')
        call require(all(fields == cmplx(0.0_dp, 0.0_dp, dp)) &
            .and. all(after%D_ion == 0.0_dp), 'guard explicitly zeros complete response')
        call check_background(before, after)
    end subroutine check_suppressed
end program test_periodic_response_normalization
