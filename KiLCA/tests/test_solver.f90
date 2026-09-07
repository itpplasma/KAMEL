module solver_test_rhs_m
    use, intrinsic :: iso_c_binding, only: c_int
    use iso_c_binding
    implicit none
contains
    subroutine diagonal_rhs(t, y, dy, params)
        real(c_double), value :: t
        real(c_double), intent(in) :: y(*)
        real(c_double), intent(out) :: dy(*)
        type(c_ptr), value :: params
        integer :: col
        do col = 0, 1
            dy(4 * col + 1:4 * col + 2) = y(4 * col + 1:4 * col + 2)
            dy(4 * col + 3:4 * col + 4) = -y(4 * col + 3:4 * col + 4)
        end do
    end subroutine

    subroutine reversed_rhs(t, y, dy, params)
        real(c_double), value :: t
        real(c_double), intent(in) :: y(*)
        real(c_double), intent(out) :: dy(*)
        type(c_ptr), value :: params

        call diagonal_rhs(t, y, dy, params)
        dy(1:8) = -dy(1:8)
    end subroutine
end module

program test_solver
    use iso_c_binding
    use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
    use kilca_solver_m, only: integrate_basis_vecs, solver_settings_t
    use solver_test_rhs_m
    implicit none
    type(solver_settings_t), target :: settings
    real(c_double) :: grid(5) = [0d0, 0.25d0, 0.5d0, 0.75d0, 1d0]
    real(c_double), target :: state(8, 5)
    real(c_double) :: expected(8), direction
    integer :: mode, rc, k

    settings%Nort = 10000
    settings%eps_rel = 1d-10
    settings%eps_abs = 1d-12
    settings%debug = 0
    do mode = 1, 3
        settings%norm_fac = 1d6
        if (mode >= 2) settings%norm_fac = 1.01d0
        state = 0d0
        state(1, 1) = 1d0
        state(7, 1) = 1d0
        direction = 1d0
        if (mode == 3) then
            direction = -1d0
            rc = integrate_basis_vecs(reversed_rhs, 2_c_int, 2_c_int, 5_c_int, &
                                      grid, state, c_loc(settings), c_null_ptr)
        else
            rc = integrate_basis_vecs(diagonal_rhs, 2_c_int, 2_c_int, 5_c_int, &
                                      grid, state, c_loc(settings), c_null_ptr)
        end if
        if (rc /= 0) error stop 'Basis integration failed'
        if (.not. all(ieee_is_finite(state))) error stop 'Nonfinite basis'
        if (abs(state(1, 5)) < 0.1d0 .or. abs(state(7, 5)) < 0.1d0) &
            error stop 'Integration lost an independent basis vector'
        if (mode == 1) then
            if (abs(state(1, 5) - exp(1d0)) > 1d-7 .or. &
                abs(state(7, 5) - exp(-1d0)) > 1d-7) &
                error stop 'Unnormalized integration changed the initial amplitudes'
        end if
        do k = 1, 5
            expected = 0d0
            expected(1) = state(1, 5) * exp(direction * (grid(k) - grid(5)))
            expected(7) = state(7, 5) * exp(direction * (grid(5) - grid(k)))
            if (any(abs(state(:, k) - expected) > 1d-7)) &
                error stop 'QR reconstruction changed the differential-equation solution'
        end do
    end do
end program
