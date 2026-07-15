module periodic_solve_m
    !> Fourier discretization of the reduced radial CGS Poisson problem.
    !>
    !> The chapter-15 KIM model uses the radial operator d2/dr2.  In the
    !> exp(+i k_m r) basis its symbol is D_m=-k_m^2, so
    !>
    !>   A = D + 4*pi*K^{rhoPhi},       g = -4*pi*K^{rhoB}*Br.
    !>
    !> The physical solution is represented as deltaPhi=Phi_ref+psi.  Phi_ref
    !> remains a physical-space lifting field; only its kinetic action and
    !> supplied continuous second derivative are projected.  Projecting and
    !> adding Phi_ref in the same finite basis would cancel algebraically and
    !> silently reproduce the old direct-periodic solve.

    use KIM_kinds_m, only: dp

    implicit none
    private

    integer, parameter, public :: PERIODIC_SOLVE_OK = 0
    integer, parameter, public :: PERIODIC_SOLVE_INVALID_INPUT = -100
    integer, parameter, public :: PERIODIC_SOLVE_GAUGE_REQUIRED = -101
    integer, parameter, public :: PERIODIC_SOLVE_INCOMPATIBLE_NULLSPACE = -102
    integer, parameter, public :: PERIODIC_SOLVE_PROJECTION_REJECTED = -103
    integer, parameter, public :: PERIODIC_SOLVE_ILL_CONDITIONED = -104
    integer, parameter, public :: PERIODIC_SOLVE_INACCURATE = -105
    integer, parameter, public :: PERIODIC_SOLVE_NONFINITE = -106
    integer, parameter, public :: PERIODIC_SOLVE_SINGULAR = -107
    real(dp), parameter, public :: PERIODIC_SOLVE_RCOND_TOLERANCE = 1.0e-12_dp
    real(dp), parameter, public :: PERIODIC_SOLVE_RESIDUAL_TOLERANCE = 1.0e-10_dp

    public :: solve_periodic, solve_periodic_deviation
    public :: assemble_periodic_operator, magnetic_drive_rhs
    public :: project_onto_periodic_basis, solve_with_derived_gauge
    public :: reconstruct_delta_phi, reconstruct_delta_phi_from_lift
    public :: build_physical_c2_lift, dense_solve

contains

    subroutine assemble_periodic_operator(Kphi, L, M, A)
        use constants_m, only: pi

        complex(dp), intent(in) :: Kphi(:,:)
        real(dp), intent(in) :: L
        integer, intent(in) :: M
        complex(dp), allocatable, intent(out) :: A(:,:)

        real(dp) :: k_m
        integer :: dim, im, m_row

        dim = 2 * M + 1
        allocate(A(dim, dim))
        A = 4.0_dp * pi * Kphi
        do m_row = -M, M
            im = m_row + M + 1
            k_m = 2.0_dp * pi * real(m_row, dp) / L
            A(im, im) = A(im, im) - k_m * k_m
        end do
    end subroutine assemble_periodic_operator

    subroutine magnetic_drive_rhs(KB, Br_m, rhs)
        use constants_m, only: pi

        complex(dp), intent(in) :: KB(:,:), Br_m(:)
        complex(dp), allocatable, intent(out) :: rhs(:)

        allocate(rhs(size(Br_m)))
        rhs = -4.0_dp * pi * matmul(KB, Br_m)
    end subroutine magnetic_drive_rhs

    !> Compatibility wrapper for the original constant-drive direct solve.
    !> It now uses the same operator, source, and explicit zero-mode policy as
    !> the deviation path.  A vacuum null space without a declared gauge is an
    !> error rather than an arbitrary LAPACK result.
    subroutine solve_periodic(Kphi, KB, L, M, Br_const, Phi_m, info)
        complex(dp), intent(in) :: Kphi(:,:), KB(:,:)
        real(dp), intent(in) :: L
        integer, intent(in) :: M
        complex(dp), intent(in) :: Br_const
        complex(dp), allocatable, intent(out) :: Phi_m(:)
        integer, intent(out) :: info

        complex(dp), allocatable :: A(:,:), Br_m(:), rhs(:)
        integer :: dim

        dim = 2 * M + 1
        if (M < 0 .or. L <= 0.0_dp .or. any(shape(Kphi) /= [dim, dim]) .or. &
            any(shape(KB) /= [dim, dim])) then
            allocate(Phi_m(max(0, dim)))
            Phi_m = (0.0_dp, 0.0_dp)
            info = PERIODIC_SOLVE_INVALID_INPUT
            return
        end if

        allocate(Br_m(dim))
        Br_m = (0.0_dp, 0.0_dp)
        Br_m(M + 1) = Br_const
        call assemble_periodic_operator(Kphi, L, M, A)
        call magnetic_drive_rhs(KB, Br_m, rhs)
        call solve_with_derived_gauge(A, rhs, M, .false., Phi_m, info)
    end subroutine solve_periodic

    !> Solve L psi = g-L Phi_ref.  The caller supplies physical-space samples
    !> of Phi_ref and its continuous radial second derivative on the DFT grid.
    subroutine solve_periodic_deviation(Kphi, KB, L, M, Br_m, r_nodes, &
                                        Phi_ref, d2Phi_ref, projection_tolerance, &
                                        psi_m, projection_residual, info, &
                                        declare_mean_zero_gauge, projection_mask)
        use constants_m, only: pi

        complex(dp), intent(in) :: Kphi(:,:), KB(:,:), Br_m(:)
        real(dp), intent(in) :: L, r_nodes(:), projection_tolerance
        integer, intent(in) :: M
        complex(dp), intent(in) :: Phi_ref(:), d2Phi_ref(:)
        complex(dp), allocatable, intent(out) :: psi_m(:)
        real(dp), intent(out) :: projection_residual
        integer, intent(out) :: info
        logical, intent(in), optional :: declare_mean_zero_gauge
        logical, intent(in), optional :: projection_mask(:)

        complex(dp), allocatable :: A(:,:), rhs(:), ref_m(:), d2ref_m(:)
        complex(dp), allocatable :: projected_ref(:)
        real(dp) :: reference_norm
        integer :: dim
        logical :: use_gauge

        dim = 2 * M + 1
        projection_residual = huge(1.0_dp)
        if (M < 0 .or. L <= 0.0_dp .or. projection_tolerance < 0.0_dp .or. &
            size(r_nodes) == 0 .or. size(Phi_ref) /= size(r_nodes) .or. &
            size(d2Phi_ref) /= size(r_nodes) .or. size(Br_m) /= dim .or. &
            any(shape(Kphi) /= [dim, dim]) .or. any(shape(KB) /= [dim, dim])) then
            allocate(psi_m(max(0, dim)))
            psi_m = (0.0_dp, 0.0_dp)
            info = PERIODIC_SOLVE_INVALID_INPUT
            return
        end if
        if (present(projection_mask)) then
            if (size(projection_mask) /= size(r_nodes) .or. .not. any(projection_mask)) then
                allocate(psi_m(dim))
                psi_m = (0.0_dp, 0.0_dp)
                info = PERIODIC_SOLVE_INVALID_INPUT
                return
            end if
        end if

        call project_onto_periodic_basis(Phi_ref, r_nodes, L, M, ref_m)
        call project_onto_periodic_basis(d2Phi_ref, r_nodes, L, M, d2ref_m)
        projected_ref = reconstruct_delta_phi(ref_m, L, M, r_nodes)
        if (present(projection_mask)) then
            reference_norm = sqrt(sum(abs(Phi_ref)**2, mask=projection_mask))
        else
            reference_norm = sqrt(sum(abs(Phi_ref)**2))
        end if
        if (reference_norm > tiny(1.0_dp)) then
            if (present(projection_mask)) then
                projection_residual = sqrt(sum(abs(projected_ref - Phi_ref)**2, &
                    mask=projection_mask)) / reference_norm
            else
                projection_residual = sqrt(sum(abs(projected_ref - Phi_ref)**2)) / reference_norm
            end if
        else
            if (present(projection_mask)) then
                projection_residual = sqrt(sum(abs(projected_ref - Phi_ref)**2, mask=projection_mask))
            else
                projection_residual = sqrt(sum(abs(projected_ref - Phi_ref)**2))
            end if
        end if
        if (projection_residual > projection_tolerance) then
            allocate(psi_m(dim))
            psi_m = (0.0_dp, 0.0_dp)
            info = PERIODIC_SOLVE_PROJECTION_REJECTED
            return
        end if

        call assemble_periodic_operator(Kphi, L, M, A)
        call magnetic_drive_rhs(KB, Br_m, rhs)
        rhs = rhs - d2ref_m - 4.0_dp * pi * matmul(Kphi, ref_m)
        use_gauge = .false.
        if (present(declare_mean_zero_gauge)) use_gauge = declare_mean_zero_gauge
        call solve_with_derived_gauge(A, rhs, M, use_gauge, psi_m, info)
    end subroutine solve_periodic_deviation

    !> Explicit zero-mode rule.  Screening retains and solves mode zero.  If
    !> the mode is a true decoupled null space, the source must be compatible
    !> and the caller must declare the mean-zero gauge.
    subroutine solve_with_derived_gauge(A_in, rhs_in, M, declare_mean_zero_gauge, &
                                        solution, info)
        complex(dp), intent(in) :: A_in(:,:), rhs_in(:)
        integer, intent(in) :: M
        logical, intent(in) :: declare_mean_zero_gauge
        complex(dp), allocatable, intent(out) :: solution(:)
        integer, intent(out) :: info

        complex(dp), allocatable :: A(:,:), rhs(:)
        real(dp) :: scale, null_tolerance
        integer :: dim, izero
        logical :: zero_is_null

        dim = 2 * M + 1
        allocate(solution(max(0, dim)))
        solution = (0.0_dp, 0.0_dp)
        if (M < 0 .or. any(shape(A_in) /= [dim, dim]) .or. size(rhs_in) /= dim) then
            info = PERIODIC_SOLVE_INVALID_INPUT
            return
        end if

        A = A_in
        rhs = rhs_in
        izero = M + 1
        scale = max(1.0_dp, maxval(abs(A)))
        null_tolerance = 100.0_dp * epsilon(1.0_dp) * scale
        zero_is_null = maxval(abs(A(izero, :))) <= null_tolerance .and. &
                       maxval(abs(A(:, izero))) <= null_tolerance
        if (zero_is_null) then
            if (abs(rhs(izero)) > null_tolerance) then
                info = PERIODIC_SOLVE_INCOMPATIBLE_NULLSPACE
                return
            end if
            if (.not. declare_mean_zero_gauge) then
                info = PERIODIC_SOLVE_GAUGE_REQUIRED
                return
            end if
            A(izero, :) = (0.0_dp, 0.0_dp)
            A(:, izero) = (0.0_dp, 0.0_dp)
            A(izero, izero) = (1.0_dp, 0.0_dp)
            rhs(izero) = (0.0_dp, 0.0_dp)
        end if

        call dense_solve(A, rhs, info)
        if (info == PERIODIC_SOLVE_OK) then
            solution = rhs
        else
            solution = (0.0_dp, 0.0_dp)
        end if
    end subroutine solve_with_derived_gauge

    subroutine project_onto_periodic_basis(values, r_nodes, L, M, coefficients)
        use constants_m, only: pi, com_unit

        complex(dp), intent(in) :: values(:)
        real(dp), intent(in) :: r_nodes(:), L
        integer, intent(in) :: M
        complex(dp), allocatable, intent(out) :: coefficients(:)

        real(dp) :: k_m
        integer :: m_row, im

        allocate(coefficients(2 * M + 1))
        do m_row = -M, M
            im = m_row + M + 1
            k_m = 2.0_dp * pi * real(m_row, dp) / L
            coefficients(im) = sum(values * exp(-com_unit * k_m * r_nodes)) / &
                               real(size(values), dp)
        end do
    end subroutine project_onto_periodic_basis

    !> C2 physical-space lifting field between the two ideal aligned edge
    !> values.  Quintic smootherstep gives zero slope and curvature at both
    !> ends while avoiding the 1/k_parallel singularity at resonance.  Its
    !> unequal endpoint values are intentionally not made periodic: doing so
    !> would restore the old boundary-value problem.
    subroutine build_physical_c2_lift(r_nodes, rm, dx_asis, dx_tr, &
                                      phi_left, phi_right, Phi_ref, d2Phi_ref)
        real(dp), intent(in) :: r_nodes(:), rm, dx_asis, dx_tr
        complex(dp), intent(in) :: phi_left, phi_right
        complex(dp), intent(out) :: Phi_ref(size(r_nodes)), d2Phi_ref(size(r_nodes))

        real(dp) :: x, t, width, s, d2s
        integer :: i

        width = 2.0_dp * (dx_asis + dx_tr)
        do i = 1, size(r_nodes)
            x = r_nodes(i) - rm
            t = max(0.0_dp, min(1.0_dp, (x + 0.5_dp * width) / width))
            call smootherstep_with_second(t, width, s, d2s)
            Phi_ref(i) = phi_left + (phi_right - phi_left) * s
            d2Phi_ref(i) = (phi_right - phi_left) * d2s
        end do
    end subroutine build_physical_c2_lift

    subroutine smootherstep_with_second(t, width, value, second_derivative)
        real(dp), intent(in) :: t, width
        real(dp), intent(out) :: value, second_derivative

        value = 6.0_dp * t**5 - 15.0_dp * t**4 + 10.0_dp * t**3
        second_derivative = (120.0_dp * t**3 - 180.0_dp * t**2 + 60.0_dp * t) / width**2
    end subroutine smootherstep_with_second

    subroutine dense_solve(A, b, info, reciprocal_condition, backward_error)
        use, intrinsic :: ieee_arithmetic, only: ieee_is_finite

        complex(dp), intent(inout) :: A(:,:), b(:)
        integer, intent(out) :: info
        real(dp), intent(out), optional :: reciprocal_condition, backward_error

        integer, allocatable :: ipiv(:)
        complex(dp), allocatable :: A0(:,:), b0(:)
        complex(dp), allocatable :: condition_work(:), residual(:)
        real(dp), allocatable :: condition_real_work(:)
        real(dp) :: anorm, norm_scale, matrix_norm, solution_norm, rhs_norm
        real(dp) :: residual_norm, rcond_local, backward_error_local
        integer :: dim, condition_info

        dim = size(A, 1)
        rcond_local = 0.0_dp
        backward_error_local = huge(1.0_dp)
        if (size(A, 2) /= dim .or. size(b) /= dim .or. dim == 0) then
            info = PERIODIC_SOLVE_INVALID_INPUT
            call export_diagnostics()
            return
        end if
        if (.not. all(ieee_is_finite(real(A, dp))) .or. &
            .not. all(ieee_is_finite(aimag(A))) .or. &
            .not. all(ieee_is_finite(real(b, dp))) .or. &
            .not. all(ieee_is_finite(aimag(b)))) then
            info = PERIODIC_SOLVE_NONFINITE
            call export_diagnostics()
            return
        end if

        A0 = A
        b0 = b
        ! 1-norm of A for ZGECON, taken before ZGESV overwrites A with its LU
        ! factors.  Uniform scaling keeps the column sums from overflowing
        ! without changing the returned value.
        norm_scale = maxval(abs(A0))
        anorm = 0.0_dp
        if (norm_scale > 0.0_dp) then
            anorm = norm_scale * maxval(sum(abs(A0)/norm_scale, dim=1))
        end if
        allocate(ipiv(dim))
        call zgesv(dim, 1, A, dim, ipiv, b, dim, info)
        if (info /= 0) then
            if (info > 0) then
                info = PERIODIC_SOLVE_SINGULAR
            else
                info = PERIODIC_SOLVE_INVALID_INPUT
            end if
            call export_diagnostics()
            return
        end if
        if (.not. all(ieee_is_finite(real(b, dp))) .or. &
            .not. all(ieee_is_finite(aimag(b)))) then
            info = PERIODIC_SOLVE_NONFINITE
            call export_diagnostics()
            return
        end if

        ! ZGECON estimates the reciprocal condition number in the matrix
        ! 1-norm directly from the LU factors ZGESV left in A, so no second
        ! factorization and no matrix copy are needed for the estimate.
        allocate(condition_work(2*dim), condition_real_work(2*dim))
        call zgecon('1', dim, A, dim, anorm, rcond_local, &
                    condition_work, condition_real_work, condition_info)
        if (condition_info /= 0) then
            info = PERIODIC_SOLVE_INVALID_INPUT
            call export_diagnostics()
            return
        end if

        ! Frobenius/L2 norms give a submultiplicative scaled backward residual.
        residual = matmul(A0, b) - b0
        if (.not. all(ieee_is_finite(real(residual, dp))) .or. &
            .not. all(ieee_is_finite(aimag(residual)))) then
            info = PERIODIC_SOLVE_NONFINITE
            call export_diagnostics()
            return
        end if
        matrix_norm = scaled_matrix_norm(A0)
        solution_norm = scaled_vector_norm(b)
        rhs_norm = scaled_vector_norm(b0)
        residual_norm = scaled_vector_norm(residual)
        backward_error_local = scaled_backward_ratio(residual_norm, matrix_norm, &
                                                     solution_norm, rhs_norm)

        if (.not. ieee_is_finite(rcond_local) .or. &
            .not. ieee_is_finite(backward_error_local)) then
            info = PERIODIC_SOLVE_NONFINITE
        else if (rcond_local <= PERIODIC_SOLVE_RCOND_TOLERANCE) then
            info = PERIODIC_SOLVE_ILL_CONDITIONED
        else if (backward_error_local >= PERIODIC_SOLVE_RESIDUAL_TOLERANCE) then
            info = PERIODIC_SOLVE_INACCURATE
        else
            info = PERIODIC_SOLVE_OK
        end if
        call export_diagnostics()

    contains

        subroutine export_diagnostics()
            use, intrinsic :: iso_fortran_env, only: error_unit

            if (present(reciprocal_condition)) reciprocal_condition = rcond_local
            if (present(backward_error)) backward_error = backward_error_local
            ! Healthy solves stay silent: hosts read the diagnostics through
            ! the optional arguments, and printing every solve floods
            ! production and test logs.  Rejections go to error_unit.
            if (info /= PERIODIC_SOLVE_OK) then
                write(error_unit, '(A,ES12.4,A,ES12.4,A,I0)') &
                    ' dense_solve: rcond=', rcond_local, &
                    ' backward_error=', backward_error_local, &
                    ' status=', info
            end if
        end subroutine export_diagnostics

        real(dp) function scaled_vector_norm(values) result(norm)
            complex(dp), intent(in) :: values(:)
            real(dp) :: scale

            scale = maxval(abs(values))
            if (scale > 0.0_dp) then
                norm = scale * sqrt(sum((abs(values)/scale)**2))
            else
                norm = 0.0_dp
            end if
        end function scaled_vector_norm

        real(dp) function scaled_matrix_norm(values) result(norm)
            complex(dp), intent(in) :: values(:,:)
            real(dp) :: scale

            scale = maxval(abs(values))
            if (scale > 0.0_dp) then
                norm = scale * sqrt(sum((abs(values)/scale)**2))
            else
                norm = 0.0_dp
            end if
        end function scaled_matrix_norm

        real(dp) function scaled_backward_ratio(residual_norm, matrix_norm, &
                                                solution_norm, rhs_norm) result(ratio)
            real(dp), intent(in) :: residual_norm, matrix_norm, solution_norm, rhs_norm
            real(dp) :: product_norm, scale

            ! Overflow guard: when ||A||*||x|| would overflow, fall back to
            ! ||r||/(||A||*||x||), dropping the ||b|| term of the regular
            ! denominator.  A smaller denominator only enlarges the reported
            ! backward error, so the reject-gate stays conservative; the
            ! discontinuity at the branch point is accepted for a gate.
            if (matrix_norm > 0.0_dp .and. solution_norm > 0.0_dp .and. &
                matrix_norm > huge(1.0_dp)/solution_norm) then
                ratio = (residual_norm/matrix_norm)/solution_norm
                return
            end if
            product_norm = matrix_norm * solution_norm
            scale = max(product_norm, rhs_norm)
            if (scale > 0.0_dp) then
                ratio = (residual_norm/scale) / &
                        (product_norm/scale + rhs_norm/scale)
            else
                ratio = residual_norm
            end if
        end function scaled_backward_ratio

    end subroutine dense_solve

    function reconstruct_delta_phi(Phi_m, L, M, r_out) result(dPhi)
        use constants_m, only: pi, com_unit

        complex(dp), intent(in) :: Phi_m(:)
        real(dp), intent(in) :: L, r_out(:)
        integer, intent(in) :: M
        complex(dp) :: dPhi(size(r_out))

        real(dp) :: k_m
        integer :: m_row, im

        dPhi = (0.0_dp, 0.0_dp)
        do m_row = -M, M
            im = m_row + M + 1
            k_m = 2.0_dp * pi * real(m_row, dp) / L
            dPhi = dPhi + Phi_m(im) * exp(com_unit * k_m * r_out)
        end do
    end function reconstruct_delta_phi

    function reconstruct_delta_phi_from_lift(psi_m, L, M, r_out, Phi_ref) result(dPhi)
        complex(dp), intent(in) :: psi_m(:), Phi_ref(:)
        real(dp), intent(in) :: L, r_out(:)
        integer, intent(in) :: M
        complex(dp) :: dPhi(size(r_out))

        dPhi = Phi_ref + reconstruct_delta_phi(psi_m, L, M, r_out)
    end function reconstruct_delta_phi_from_lift

end module periodic_solve_m
