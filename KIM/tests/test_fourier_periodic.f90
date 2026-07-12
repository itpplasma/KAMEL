program test_fourier_periodic
    ! Behavioral checks of the periodic Fourier assembly and solve on a
    ! uniform background: exact diagonality from discrete orthogonality,
    ! the adiabatic screening diagonal, the screened point-charge solve
    ! against the analytic per-mode solution, and the periodic expansion.

    use KIM_kinds_m, only: dp
    use constants_m, only: pi, com_unit
    use species_m, only: plasma
    use kernels_m, only: kernel_rho_phi_of_kr_krp_rg
    use fourier_periodic_m, only: periodic_k_grid, assemble_periodic_kernel, &
                                  solve_periodic_potential, periodic_field_value
    use kernel_test_background_m, only: setup_uniform_background, require_close

    implicit none

    integer, parameter :: n_modes = 6
    integer, parameter :: n_k = 2*n_modes + 1
    integer, parameter :: n_quad = 64
    real(dp), parameter :: r_mid = 36.5d0
    real(dp), parameter :: period = 30.0d0
    real(dp), parameter :: charge = 3.0d0

    real(dp) :: k_modes(n_k), debye_sum, offdiag, diag_scale
    complex(dp) :: kernel_matrix(n_k, n_k), rhs(n_k), phi_modes(n_k)
    complex(dp) :: analytic, value_a, value_b
    integer :: m, mp, sp

    call setup_uniform_background()
    k_modes = periodic_k_grid(period, n_modes)
    call require_close("mode grid is centered on k = 0", &
                       cmplx(k_modes(n_modes + 1), 0.0d0, dp), &
                       (0.0d0, 0.0d0), 1.0d-15)

    debye_sum = 0.0d0
    do sp = 0, plasma%n_species - 1
        debye_sum = debye_sum + 1.0d0/plasma%spec(sp)%lambda_D(1)**2
    end do

    kernel_matrix = assemble_periodic_kernel(kernel_rho_phi_of_kr_krp_rg, &
                                             k_modes, r_mid, period, n_quad)

    ! Discrete orthogonality: the uniform background is exactly diagonal.
    offdiag = 0.0d0
    diag_scale = 0.0d0
    do mp = 1, n_k
        do m = 1, n_k
            if (m == mp) then
                diag_scale = max(diag_scale, abs(kernel_matrix(m, mp)))
            else
                offdiag = max(offdiag, abs(kernel_matrix(m, mp)))
            end if
        end do
    end do
    call require_close("uniform assembly is exactly diagonal", &
                       cmplx(offdiag/diag_scale, 0.0d0, dp), (0.0d0, 0.0d0), &
                       1.0d-12)

    ! Every diagonal element is the adiabatic screening value.
    do m = 1, n_k
        call require_close("diagonal is -sum(1/lambda_D^2)/(4 pi)", &
                           kernel_matrix(m, m), &
                           cmplx(-debye_sum/(4.0d0*pi), 0.0d0, dp), 1.0d-12)
    end do

    ! Screened point charge: (k^2 + 1/lambda^2) Phi = 4 pi q/L per mode.
    do m = 1, n_k
        rhs(m) = 4.0d0*pi*charge/period*exp(-com_unit*k_modes(m)*r_mid)
    end do
    phi_modes = solve_periodic_potential(k_modes, kernel_matrix, rhs)
    do m = 1, n_k
        analytic = 4.0d0*pi*charge/period*exp(-com_unit*k_modes(m)*r_mid) &
                   /(k_modes(m)**2 + debye_sum)
        call require_close("mode solve matches the screened solution", &
                           phi_modes(m), analytic, 1.0d-12)
    end do

    ! The expansion is L-periodic and real at the charge location.
    value_a = periodic_field_value(k_modes, phi_modes, r_mid + 0.3d0*period)
    value_b = periodic_field_value(k_modes, phi_modes, &
                                   r_mid + 1.3d0*period)
    call require_close("expansion is L-periodic", value_b, value_a, 1.0d-12)
    value_a = periodic_field_value(k_modes, phi_modes, r_mid)
    call require_close("potential is real at the charge", &
                       cmplx(aimag(value_a), 0.0d0, dp), (0.0d0, 0.0d0), &
                       1.0d-12)

    print *, "PASS all fourier periodic checks"

end program test_fourier_periodic
