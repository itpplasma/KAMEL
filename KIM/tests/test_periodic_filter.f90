program test_periodic_filter
    ! Exponential reconstruction filter of the forced-periodicity solver:
    ! weights, invariance of low Fourier modes, suppression of the Gibbs ringing
    ! of a truncated pulse, and an end-to-end periodic solve with and without it.
    use KIM_kinds_m, only: dp

    implicit none

    logical :: all_passed

    all_passed = .true.
    call test_weights(all_passed)
    call test_low_modes_unchanged(all_passed)
    call test_pulse_ringing(all_passed)
    call test_periodic_solve(all_passed)

    if (all_passed) then
        print *, 'All periodic filter tests PASSED'
        stop 0
    else
        print *, 'Periodic filter tests FAILED'
        stop 1
    end if

contains

    subroutine check(name, ok, passed)
        character(*), intent(in) :: name
        logical, intent(in) :: ok
        logical, intent(inout) :: passed

        if (ok) then
            print *, 'PASS: ', name
        else
            print *, 'FAIL: ', name
            passed = .false.
        end if
    end subroutine check

    subroutine test_weights(passed)
        use periodic_solve_m, only: reconstruction_filter
        logical, intent(inout) :: passed
        integer, parameter :: M = 40
        real(dp) :: sigma(2*M + 1)

        sigma = reconstruction_filter(M, 0)
        call check('order 0 leaves every mode unchanged', all(sigma == 1.0_dp), passed)
        sigma = reconstruction_filter(M, 8)
        call check('filter keeps the mean', sigma(M + 1) == 1.0_dp, passed)
        call check('filter removes the top modes', &
            abs(sigma(1) - epsilon(1.0_dp)) < 1.0e-3_dp * epsilon(1.0_dp) .and. &
            abs(sigma(2*M + 1) - epsilon(1.0_dp)) < 1.0e-3_dp * epsilon(1.0_dp), passed)
        call check('filter is symmetric and monotonic', &
            all(sigma(1:M) == sigma(2*M + 1:M + 2:-1)) .and. &
            all(sigma(M + 2:2*M + 1) <= sigma(M + 1:2*M)), passed)
    end subroutine test_weights

    subroutine test_low_modes_unchanged(passed)
        use periodic_solve_m, only: reconstruction_filter, reconstruct_delta_phi
        logical, intent(inout) :: passed
        integer, parameter :: M = 40, n_out = 200
        real(dp), parameter :: L = 3.0_dp
        complex(dp) :: coefficients(2*M + 1)
        complex(dp) :: plain(n_out), filtered(n_out)
        real(dp) :: r(n_out)
        integer :: i, mm

        coefficients = (0.0_dp, 0.0_dp)
        do mm = -M/4, M/4
            coefficients(mm + M + 1) = cmplx(1.0_dp / (1 + abs(mm)), 0.3_dp * mm / M, dp)
        end do
        do i = 1, n_out
            r(i) = 50.0_dp + L * real(i - 1, dp) / real(n_out, dp)
        end do
        plain = reconstruct_delta_phi(coefficients, L, M, r)
        filtered = reconstruct_delta_phi(coefficients, L, M, r, reconstruction_filter(M, 8))
        call check('modes up to M/4 change by less than 1e-3', &
            maxval(abs(filtered - plain)) < 1.0e-3_dp * maxval(abs(plain)), passed)
    end subroutine test_low_modes_unchanged

    subroutine test_pulse_ringing(passed)
        use constants_m, only: pi
        use periodic_solve_m, only: reconstruction_filter, reconstruct_delta_phi
        logical, intent(inout) :: passed

        ! Box of width w on a period L: c_m = (w/L) sinc(m w / L). Away from the
        ! jumps the truncated series rings around 0 or 1; the filter damps it.
        integer, parameter :: M = 60, n_out = 600
        real(dp), parameter :: L = 10.0_dp, w = 2.0_dp
        complex(dp) :: coefficients(2*M + 1)
        complex(dp) :: plain(n_out), filtered(n_out)
        real(dp) :: r(n_out), exact(n_out), x
        logical :: far(n_out)
        integer :: i, mm

        do mm = -M, M
            x = pi * real(mm, dp) * w / L
            if (mm == 0) then
                coefficients(mm + M + 1) = w / L
            else
                coefficients(mm + M + 1) = w / L * sin(x) / x
            end if
        end do
        do i = 1, n_out
            r(i) = -0.5_dp * L + L * real(i - 1, dp) / real(n_out, dp)
            exact(i) = merge(1.0_dp, 0.0_dp, abs(r(i)) < 0.5_dp * w)
            far(i) = abs(abs(r(i)) - 0.5_dp * w) > 1.0_dp
        end do
        plain = reconstruct_delta_phi(coefficients, L, M, r)
        filtered = reconstruct_delta_phi(coefficients, L, M, r, reconstruction_filter(M, 8))
        call check('filter damps ringing 1 away from the jumps by > 10', &
            maxval(abs(filtered - exact), mask=far) < &
            0.1_dp * maxval(abs(plain - exact), mask=far), passed)
        call check('filter keeps the integral of the pulse', &
            abs(sum(real(filtered, dp)) - sum(real(plain, dp))) < 1.0e-10_dp * sum(exact), passed)
    end subroutine test_pulse_ringing

    subroutine write_config(path, filter_order)
        character(*), intent(in) :: path
        integer, intent(in) :: filter_order
        character(len=512) :: line
        integer :: unit_in, unit_out, ios

        open(newunit=unit_in, file='KIM_config_em_small.nml', status='old', action='read')
        open(newunit=unit_out, file=path, status='replace', action='write')
        do
            read(unit_in, '(a)', iostat=ios) line
            if (ios /= 0) exit
            write(unit_out, '(a)') trim(line)
        end do
        close(unit_in)
        write(unit_out, '(a)') '&KIM_PERIODIC'
        write(unit_out, '(a)') '    periodic_calculate_radial_current = .false.'
        write(unit_out, '(a)') '    periodic_calculate_ion_tensor = .false.'
        write(unit_out, '(a,i0)') '    periodic_kmax_scale = 15.0, periodic_n_rg = 512, periodic_filter_order = ', &
            filter_order
        write(unit_out, '(a)') '/'
        close(unit_out)
    end subroutine write_config

    subroutine solve(config, prof, res)
        use kim_solver_m, only: kim_solver_t, kim_results_t, kim_profiles_t, KIM_OK
        character(*), intent(in) :: config
        type(kim_profiles_t), intent(in) :: prof
        type(kim_results_t), intent(out) :: res
        type(kim_solver_t) :: kim
        integer :: ierr

        call kim%init(config, run_type='electrostatic_periodic', profiles=prof, stat=ierr)
        if (ierr == KIM_OK) call kim%solve(m=-6, n=2, stat=ierr)
        if (ierr /= KIM_OK) error stop 'periodic solve failed'
        res = kim%results()
        call kim%finalize()
    end subroutine solve

    subroutine test_periodic_solve(passed)
        use kim_solver_m, only: kim_results_t, kim_profiles_t
        logical, intent(inout) :: passed

        integer, parameter :: npts = 64
        type(kim_profiles_t) :: prof
        type(kim_results_t) :: plain, filtered
        real(dp) :: frac, layer_plain, layer_filtered
        complex(dp) :: integral_plain, integral_filtered
        real(dp), allocatable :: x(:)
        integer :: i

        allocate(prof%r(npts), prof%n(npts), prof%Te(npts), prof%Ti(npts), &
                 prof%q(npts), prof%Er(npts))
        do i = 1, npts
            prof%r(i) = 3.0_dp + 64.0_dp * real(i - 1, dp) / real(npts - 1, dp)
            frac = (prof%r(i) - 3.0_dp) / 64.0_dp
            prof%n(i) = 2.0e13_dp * (1.1_dp - prof%r(i) / 100.0_dp)
            prof%Te(i) = 1.0e3_dp * (1.2_dp - prof%r(i) / 100.0_dp)
            prof%Ti(i) = prof%Te(i)
            prof%q(i) = 1.5_dp + 2.5_dp * frac
            prof%Er(i) = -0.5_dp
        end do

        call write_config('KIM_config_filter_off.nml', 0)
        call write_config('KIM_config_filter_on.nml', 8)
        call solve('KIM_config_filter_off.nml', prof, plain)
        call solve('KIM_config_filter_on.nml', prof, filtered)

        x = plain%r_field - plain%r_resonance
        integral_plain = sum(plain%jpar, mask=abs(x) < 1.0_dp)
        integral_filtered = sum(filtered%jpar, mask=abs(x) < 1.0_dp)
        call check('filter keeps the layer current to 1 %', &
            abs(integral_filtered - integral_plain) < 1.0e-2_dp * abs(integral_plain), passed)
        layer_plain = maxval(abs(plain%jpar))
        layer_filtered = maxval(abs(filtered%jpar))
        call check('filtered current stays finite and of the same size', &
            layer_filtered > 0.3_dp * layer_plain .and. layer_filtered < 1.5_dp * layer_plain, passed)
        call check('filter does not change the window', &
            maxval(abs(filtered%r_field - plain%r_field)) == 0.0_dp, passed)
    end subroutine test_periodic_solve

end program test_periodic_filter
