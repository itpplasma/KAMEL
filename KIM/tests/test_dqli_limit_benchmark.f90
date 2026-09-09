program test_dqli_limit_benchmark
    !! Physical-field diagnostic tests; the production FLR tensor is symmetric,
    !! whereas the full legacy tensor can also contain an antisymmetric part.
    use KIM_kinds_m, only: dp
    use constants_m, only: sol
    use config_m, only: resolved_ion_ifunc_conservation_model, IFUNC_MODEL_ENERGY
    use kim_qldiff_m, only: calc_dqli_limit_benchmark, calc_dqli_tensor, calc_dqli11_phi, &
        calc_dqli_benchmark_residuals
    implicit none

    real(dp), parameter :: vt = 2.3e7_dp, nu = 3.7e5_dp, b0 = 1.8e4_dp
    real(dp), parameter :: omega_c = 4.1e7_dp, omega = 3.0e3_dp
    real(dp), parameter :: om_e = omega-0.1_dp*nu, kpar = 0.5_dp*nu/vt
    complex(dp), parameter :: zero = (0.0_dp, 0.0_dp), br = (0.8_dp, 0.3_dp)
    real(dp) :: old(2,2), new(2,2), absolute(2,2), relative(2,2)
    integer :: failures
    character(len=32) :: argument

    failures = 0
    call get_command_argument(1, argument)
    if (trim(argument) == 'reject_zero_ks') then
        call benchmark(0.0_dp, cmplx(1.0e-3_dp, 0.0_dp, dp), br)
        stop 0
    end if
    call test_physical_electric_units()
    call test_real_model_limits()
    call test_amplitude_scaling()
    call test_zero_wave_number()
    call test_residual_boundaries()
    if (failures /= 0) then
        print *, 'failed benchmark assertions: ', failures
        error stop 'drift-kinetic diagnostic benchmark failed'
    end if
    print *, 'physical drift-kinetic diagnostic benchmark passed'

contains

    subroutine benchmark(ks, es, magnetic)
        real(dp), intent(in) :: ks
        complex(dp), intent(in) :: es, magnetic
        call calc_dqli_limit_benchmark(vt, nu, om_e, b0, kpar, ks, omega_c, omega, &
            es, magnetic, old, new, absolute, relative)
    end subroutine benchmark

    subroutine check(condition, message)
        logical, intent(in) :: condition
        character(len=*), intent(in) :: message
        if (condition) return
        failures = failures+1
        print *, 'FAIL: ', message
    end subroutine check

    subroutine test_physical_electric_units()
        complex(dp), parameter :: es = (1.2e-3_dp, -0.7e-3_dp)
        real(dp) :: reference(2,2), electric_d11
        resolved_ion_ifunc_conservation_model = IFUNC_MODEL_ENERGY
        call calc_dqli_tensor(vt, nu, om_e-omega, b0, kpar, es, zero, &
            reference(1,1), reference(1,2), reference(2,1), reference(2,2))
        electric_d11 = calc_dqli11_phi(vt, nu, om_e-omega, b0, kpar, es)
        call benchmark(1.0e-4_dp, es, zero)
        call check(maxval(abs(old-reference))/maxval(abs(reference)) < 2.0e-13_dp, &
            'physical Es and nonzero mode-frequency detuning match established tensor')
        call check(abs(new(1,1)-electric_d11)/abs(electric_d11) < 1.0e-7_dp, &
            'Phi=i Es/ks preserves physical electric D11 including c squared')
    end subroutine test_physical_electric_units

    subroutine test_real_model_limits()
        real(dp), parameter :: ks_values(3) = [1.0e-2_dp, 1.0e-3_dp, 1.0e-4_dp]
        real(dp), parameter :: phases(3) = [0.0_dp, 0.7_dp, 1.5707963267948966_dp]
        real(dp) :: symmetric(2,2), errors(3), gap, scale
        complex(dp) :: es
        integer :: model, phase, ik
        do model = 0, 3
            resolved_ion_ifunc_conservation_model = model
            do phase = 1, size(phases)
                ! c|Es|=vT|Br|, with physical CGS Es, so every channel matters.
                es = vt/sol*br*exp(cmplx(0.0_dp, phases(phase), dp))
                do ik = 1, size(ks_values)
                    call benchmark(ks_values(ik), es, br)
                    symmetric = 0.5_dp*(old+transpose(old))
                    scale = maxval(abs(symmetric))
                    errors(ik) = maxval(abs(new-symmetric))/scale
                    call check(maxval(abs(new-transpose(new)))/scale < 1.0e-13_dp, &
                        'current production tensor remains symmetric')
                    call check(all(absolute >= 0.0_dp) .and. all(relative >= 0.0_dp), &
                        'reported absolute and relative errors are nonnegative')
                end do
                call check(errors(3) < 2.0e-7_dp, 'symmetric physical limit for all models')
                call check(errors(1) > 50.0_dp*errors(2) .and. &
                    errors(2) > 50.0_dp*errors(3), &
                    'finite-FLR symmetric error converges quadratically')
                gap = abs(old(1,2)-old(2,1))/scale
                if ((model == 0 .or. model == 2) .and. phase > 1) then
                    call check(gap > 1.0e-3_dp, &
                        'physical fixture has legacy antisymmetric transport')
                    call check(maxval(relative) > 1.0e-3_dp, &
                        'full-tensor residual honestly retains antisymmetric mismatch')
                else
                    call check(gap < 1.0e-12_dp, 'conserving or in-phase fixture is symmetric')
                    call check(maxval(relative) < 2.0e-7_dp, 'full tensor converges when symmetric')
                end if
            end do
        end do
    end subroutine test_real_model_limits

    subroutine test_amplitude_scaling()
        real(dp), parameter :: factors(3) = [1.0e-12_dp, 1.0_dp, 1.0e12_dp]
        real(dp) :: old_base(2,2), new_base(2,2), abs_base(2,2), rel_base(2,2), factor
        complex(dp) :: es
        integer :: i
        resolved_ion_ifunc_conservation_model = 0
        es = vt/sol*cmplx(1.0_dp, 1.0_dp, dp)
        call benchmark(1.0e-3_dp, es, br)
        old_base = old
        new_base = new
        abs_base = absolute
        rel_base = relative
        do i = 1, size(factors)
            factor = factors(i)
            call benchmark(1.0e-3_dp, factor*es, factor*br)
            call check(maxval(abs(old/factor**2-old_base))/maxval(abs(old_base)) < 1.0e-12_dp, &
                'legacy tensor scales with squared field amplitude')
            call check(maxval(abs(new/factor**2-new_base))/maxval(abs(new_base)) < 1.0e-12_dp, &
                'FLR tensor scales with squared field amplitude')
            call check(maxval(abs(absolute/factor**2-abs_base))/maxval(abs(old_base)) &
                < 1.0e-12_dp, &
                'absolute error scales with squared field amplitude')
            call check(maxval(abs(relative-rel_base)) < 1.0e-12_dp, &
                'fractional error is independent of physical amplitude')
        end do
    end subroutine test_amplitude_scaling

    subroutine test_residual_boundaries()
        real(dp) :: reference(2,2), candidate(2,2), expected(2,2)
        reference = reshape([0.0_dp, 2.0e-200_dp, -2.0_dp, 0.0_dp], [2,2])
        candidate = reshape([3.0e-200_dp, 0.0_dp, 2.0_dp, 0.0_dp], [2,2])
        expected = reshape([1.0_dp, 1.0_dp, 2.0_dp, 0.0_dp], [2,2])
        call calc_dqli_benchmark_residuals(reference, candidate, absolute, relative)
        call check(all(relative == expected), &
            'relative residual handles zero reference, zero candidate, sign reversal and zero/zero')
        call check(all(absolute == abs(candidate-reference)), 'absolute residual uses magnitude')
        reference = reshape([1.0e-200_dp, 3.0_dp, -4.0_dp, 8.0_dp], [2,2])
        candidate = 2.0_dp*reference
        call calc_dqli_benchmark_residuals(reference, candidate, absolute, relative)
        call check(maxval(abs(relative-0.5_dp)) < 1.0e-15_dp, &
            'normalization uses larger coefficient magnitude without a dimensionful floor')
        call calc_dqli_benchmark_residuals(candidate, reference, absolute, relative)
        call check(maxval(abs(relative-0.5_dp)) < 1.0e-15_dp, &
            'fractional residual is symmetric under exchanging the compared tensors')
    end subroutine test_residual_boundaries

    subroutine test_zero_wave_number()
        resolved_ion_ifunc_conservation_model = IFUNC_MODEL_ENERGY
        call benchmark(0.0_dp, zero, br)
        call check(maxval(relative) < 1.0e-12_dp, 'exact zero-ks magnetic tensor limit')
        call benchmark(0.0_dp, zero, zero)
        call check(all(old == 0.0_dp) .and. all(new == 0.0_dp), 'zero fields give zero tensor')
        call check(all(absolute == 0.0_dp) .and. all(relative == 0.0_dp), &
            'zero/zero residual is exactly zero')
    end subroutine test_zero_wave_number
end program test_dqli_limit_benchmark
