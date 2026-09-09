program test_dqli_periodic_limit
    !! Constant-background periodic transport benchmark with physical CGS fields.
    !! Sample an analytic potential, Fourier transform it, and exercise the same
    !! ordered-wave assembly as the periodic solver. This is a constitutive
    !! benchmark for prescribed Phi/Br, not a self-consistent Poisson benchmark.
    use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
    use KIM_kinds_m, only: dp
    use constants_m, only: pi, sol, com_unit
    use config_m, only: resolved_ion_ifunc_conservation_model
    use setup_m, only: mphi_max
    use species_m, only: evaluate_susceptibility
    use rt_electrostatic_periodic_m, only: compute_periodic_ion_tensor_spectrum
    implicit none

    integer, parameter :: nphase = 7
    real(dp), parameter :: vti = 2.3e7_dp, nui = 3.7e5_dp, omega_ci = -4.1e7_dp
    real(dp), parameter :: b0 = 1.8e4_dp, om_e = 1.1e4_dp, omega_mode = 3.0e3_dp
    real(dp), parameter :: kpar = 2.0e-5_dp, electric_scale = 1.0e-3_dp
    real(dp), parameter :: ks_values(4) = [0.1_dp, 0.05_dp, 0.025_dp, 0.0125_dp]
    integer, parameter :: cutoffs(4) = [1, 2, 4, 6], sample_counts(3) = [10, 20, 40]
    complex(dp) :: br, moments(0:3, 0:3)
    real(dp) :: phases(nphase), drift(2, 2, nphase), resolved(2, 2, nphase)
    real(dp) :: actual(2, 2, nphase), sampling_reference(2, 2, nphase)
    real(dp) :: cutoff_error(4), sampling_error(3), flr_error(4), errors(2, 2)
    integer :: i, j

    resolved_ion_ifunc_conservation_model = 1
    mphi_max = 0
    ! c|Es| and vTi|Br| have comparable magnitudes; the electric and mixed
    ! channels cannot disappear underneath a much larger magnetic-only tensor.
    br = cmplx(0.8_dp, -0.4_dp, dp)*sol*electric_scale/vti
    call evaluate_susceptibility(kpar*vti/nui, -(om_e - omega_mode)/nui, &
        resolved_ion_ifunc_conservation_model, moments)
    do j = 1, nphase
        phases(j) = 2.0_dp*pi*(real(j - 1, dp) + 0.37_dp)/real(nphase, dp)
        call drift_ledger(electric_field(phases(j)), br, moments, drift(:, :, j))
    end do

    ! Fourier truncation at fixed fine sampling and fixed finite FLR. The
    ! analytic smooth profile has an infinite Fourier tail, so this is genuine
    ! convergence, rather than padding a finite spectrum with zeros.
    call spectral_tensor(8, 128, ks_values(1), resolved)
    do i = 1, size(cutoffs)
        call spectral_tensor(cutoffs(i), 64, ks_values(1), actual)
        errors = relative_entry_errors(actual, resolved)
        cutoff_error(i) = maxval(errors)
        print '(A,I3,A,4ES13.4)', 'Fourier cutoff M=', cutoffs(i), ' entry errors:', errors
    end do
    call require(cutoff_error(1) > 1.0e-3_dp, 'cutoff fixture has no resolved Fourier tail')
    call require(cutoff_error(2) < 0.5_dp*cutoff_error(1), 'M=2 did not improve truncation')
    call require(cutoff_error(3) < 0.1_dp*cutoff_error(2), 'M=4 did not improve truncation')
    call require(cutoff_error(4) < 1.0e-5_dp, 'M=6 Fourier truncation remains too large')

    ! Sampling/aliasing error is measured separately at fixed M=4, independent
    ! of the physical FLR-to-drift error and the omitted modes above M=4.
    call spectral_tensor(4, 128, ks_values(1), sampling_reference)
    do i = 1, size(sample_counts)
        call spectral_tensor(4, sample_counts(i), ks_values(1), actual)
        errors = relative_entry_errors(actual, sampling_reference)
        sampling_error(i) = maxval(errors)
        print '(A,I3,A,4ES13.4)', 'Sampling N=', sample_counts(i), ' entry errors:', errors
    end do
    call require(sampling_error(1) > 1.0e-7_dp, 'sampling fixture does not detect aliasing')
    call require(sampling_error(2) < 0.01_dp*sampling_error(1), 'sampling did not converge')
    call require(sampling_error(3) < 1.0e-11_dp, 'resolved sampling error exceeds roundoff')

    ! Scale both tangential and radial wavenumbers to zero while retaining
    ! physical Es and Br. Since rho_i is fixed, the leading FLR error is
    ! O((k_perp*rho_i)^2); halving every wavenumber should reduce it fourfold.
    do i = 1, size(ks_values)
        call spectral_tensor(8, 128, ks_values(i), actual)
        errors = relative_entry_errors(actual, drift)
        flr_error(i) = maxval(errors)
        print '(A,ES13.4,A,4ES13.4)', 'k_s*rho_i=', &
            ks_values(i)*vti/abs(omega_ci), ' FLR/drift entry errors:', errors
        if (i > 1) then
            call require(flr_error(i) < 0.35_dp*flr_error(i - 1), &
                'physical FLR error did not decrease quadratically')
            call require(flr_error(i) > 0.15_dp*flr_error(i - 1), &
                'FLR scan no longer resolves the leading quadratic error')
        end if
    end do
    call require(flr_error(4) < 1.0e-3_dp, 'resolved periodic tensor misses drift limit')
    print *, 'PASS: periodic spectral tensor converges to all four drift entries'

contains

    complex(dp) function electric_field(theta) result(es)
        real(dp), intent(in) :: theta
        es = electric_scale*cmplx(exp(0.8_dp*cos(theta)), 0.2_dp*sin(2.0_dp*theta), dp)
    end function electric_field

    subroutine spectral_tensor(cutoff, nsamples, ks, tensor)
        integer, intent(in) :: cutoff, nsamples
        real(dp), intent(in) :: ks
        real(dp), intent(out) :: tensor(2, 2, nphase)
        complex(dp) :: phi_m(2*cutoff + 1), phi_sample
        real(dp) :: theta, period
        integer :: sample, radial_mode, probe

        call require(nsamples > 2*cutoff, 'Fourier sample count is below Nyquist')
        period = 2.0_dp*pi/ks
        phi_m = (0.0_dp, 0.0_dp)
        do sample = 0, nsamples - 1
            theta = 2.0_dp*pi*real(sample, dp)/real(nsamples, dp)
            ! Es=-i*ks*Phi, with Phi in statV and Es in statV/cm.
            phi_sample = com_unit*electric_field(theta)/ks
            do radial_mode = -cutoff, cutoff
                phi_m(radial_mode + cutoff + 1) = phi_m(radial_mode + cutoff + 1) &
                    + phi_sample*exp(-com_unit*real(radial_mode, dp)*theta)/real(nsamples, dp)
            end do
        end do
        do probe = 1, nphase
            call compute_periodic_ion_tensor_spectrum(phi_m, br, period, phases(probe)/ks, &
                ks, kpar, vti, nui, omega_ci, omega_mode, om_e, b0, tensor(:, :, probe))
        end do
        call require(all(ieee_is_finite(tensor)), 'spectral tensor contains non-finite entries')
    end subroutine spectral_tensor

    subroutine drift_ledger(es, magnetic, ifunc, tensor)
        !! Independent zero-FLR Heyn/Markl ledger, using actual susceptibility
        !! moments. The spectral FLR engine exposes only its symmetric part.
        complex(dp), intent(in) :: es, magnetic, ifunc(0:3, 0:3)
        real(dp), intent(out) :: tensor(2, 2)
        real(dp) :: electric2, magnetic2, mixed, factor

        factor = 0.5_dp/(nui*b0**2)
        electric2 = sol**2*abs(es)**2
        magnetic2 = vti**2*abs(magnetic)**2
        mixed = 2.0_dp*sol*vti*real(conjg(es)*magnetic, dp)
        tensor(1, 1) = factor*(electric2*real(ifunc(0, 0), dp) &
            + mixed*real(ifunc(1, 0), dp) + magnetic2*real(ifunc(1, 1), dp))
        tensor(1, 2) = factor*(electric2*real(ifunc(0, 0) + 0.5_dp*ifunc(2, 0), dp) &
            + mixed*real(ifunc(1, 0) + 0.25_dp*(ifunc(3, 0) + ifunc(2, 1)), dp) &
            + magnetic2*real(ifunc(1, 1) + 0.5_dp*ifunc(3, 1), dp))
        tensor(2, 1) = tensor(1, 2)
        tensor(2, 2) = factor*(electric2*real(2.0_dp*ifunc(0, 0) + ifunc(2, 0) &
            + 0.25_dp*ifunc(2, 2), dp) + mixed*real(2.0_dp*ifunc(1, 0) &
            + 0.5_dp*(ifunc(3, 0) + ifunc(2, 1)) + 0.25_dp*ifunc(3, 2), dp) &
            + magnetic2*real(2.0_dp*ifunc(1, 1) + ifunc(3, 1) + 0.25_dp*ifunc(3, 3), dp))
    end subroutine drift_ledger

    function relative_entry_errors(got, expected) result(errors)
        real(dp), intent(in) :: got(2, 2, nphase), expected(2, 2, nphase)
        real(dp) :: errors(2, 2), scale
        integer :: row, col
        do row = 1, 2
            do col = 1, 2
                scale = maxval(abs(expected(row, col, :)))
                call require(scale > tiny(1.0_dp), 'reference tensor entry is not exercised')
                errors(row, col) = maxval(abs(got(row, col, :) - expected(row, col, :)))/scale
            end do
        end do
    end function relative_entry_errors

    subroutine require(condition, message)
        logical, intent(in) :: condition
        character(*), intent(in) :: message
        if (.not. condition) then
            print *, 'FAIL: ', message
            error stop 1
        end if
    end subroutine require
end program test_dqli_periodic_limit
