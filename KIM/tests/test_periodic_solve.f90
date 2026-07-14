program test_periodic_solve
    ! Validates the derived periodic electrostatic operator, physical-space
    ! deviation lift, explicit gauge policy, and inverse-DFT reconstruction.
    !
    ! System over Fourier modes m = -M..M (k_m = 2*pi*m/L):
    !   [ D_m delta_{m,m'} + 4*pi K^{rhoPhi}_{m,m'} ] Phi_{m'}
    !                                        = -4*pi K^{rhoB}_{m,m'} Br_{m'},
    ! with D_m = -k_m^2 the (diagonal) Fourier-space radial Poisson operator
    ! (matching solve_poisson.f90, whose FEM Laplace operator is the radial
    ! second derivative only), the (2M+1)^2 matrices Kphi/KB from Phase-2.4,
    ! and for type_br_field = 12 a constant Br over the window whose Fourier
    ! coefficient is nonzero only at m' = 0:  Br_{m'} = Br_const*delta_{m',0}.
    ! Hence  b_m = -4*pi*Br_const*KB(m, m'=0) = -4*pi*Br_const*KB(:, M+1).
    ! Solve A Phi = b via LAPACK zgesv, then reconstruct
    !   dPhi(r) = sum_{m=-M}^{M} Phi_m exp(i k_m r).
    !
    ! Assertions:
    !   (a) solver round-trip: a manufactured, diagonally-dominant complex
    !       (2M+1)^2 system solved through the SAME zgesv wiring as
    !       solve_periodic (via the exposed private dense_solve) recovers a
    !       known solution to < 1e-10;
    !   (b) inverse-DFT correctness: reconstruct_delta_phi of a single-mode
    !       coefficient gives exp(i k_m r) to 1e-12, and a two-mode case
    !       confirms superposition;
    !   (c) end-to-end: the REAL assembled (M=6) periodic matrices give a
    !       non-singular solve (info==0), a finite non-zero Phi_m response to a
    !       constant Br drive, and a finite reconstructed dPhi on the window.

    use KIM_kinds_m, only: dp
    use constants_m, only: pi, com_unit
    use periodic_solve_m, only: solve_periodic, reconstruct_delta_phi, dense_solve, &
        assemble_periodic_operator, project_onto_periodic_basis, &
        solve_with_derived_gauge, solve_periodic_deviation, &
        reconstruct_delta_phi_from_lift, build_physical_c2_lift, &
        PERIODIC_SOLVE_OK, PERIODIC_SOLVE_GAUGE_REQUIRED, &
        PERIODIC_SOLVE_INCOMPATIBLE_NULLSPACE, PERIODIC_SOLVE_PROJECTION_REJECTED, &
        PERIODIC_SOLVE_RCOND_TOLERANCE, PERIODIC_SOLVE_RESIDUAL_TOLERANCE

    implicit none

    call test_dense_solve_roundtrip()
    call test_dense_solve_singular_rejection()
    call test_dense_solve_condition_rejection()
    call test_dense_solve_nonfinite_rejection()
    call test_inverse_dft()
    call test_field_operator_oracle()
    call test_zero_mode_policy()
    call test_manufactured_deviation()
    call test_projection_rejection()
    call test_end_to_end()

    print *, 'All tests PASSED'
    stop 0

contains

    !> (a) Manufactured round-trip through the raw zgesv wiring. Isolates the
    !> LAPACK argument pattern shared with solve_periodic: build a well-
    !> conditioned diagonally-dominant complex A0, pick a known Phi, form
    !> b0 = A0*Phi, solve, and assert recovery to < 1e-10.
    subroutine test_dense_solve_roundtrip()
        use, intrinsic :: ieee_arithmetic, only: ieee_is_finite

        integer, parameter :: M = 3
        integer, parameter :: dim = 2 * M + 1
        complex(dp) :: A0(dim, dim), Phi_known(dim), b0(dim), x(dim)
        integer :: i, j, info
        real(dp) :: re, im, rcond, backward_error

        do j = 1, dim
            do i = 1, dim
                ! Deterministic, small off-diagonal entries.
                re = 0.1_dp * real(i - j, dp)
                im = 0.05_dp * real(i + j, dp)
                A0(i, j) = cmplx(re, im, dp)
            end do
            ! Diagonal dominance -> well-conditioned.
            A0(j, j) = cmplx(10.0_dp + real(j, dp), 1.0_dp, dp)
            Phi_known(j) = cmplx(real(j, dp), real(dim - j, dp), dp)
        end do

        b0 = matmul(A0, Phi_known)
        x = b0
        call dense_solve(A0, x, info, rcond, backward_error)

        if (info /= 0) then
            print *, 'FAIL: dense_solve info /= 0:', info
            error stop
        end if
        if (maxval(abs(x - Phi_known)) >= 1.0e-10_dp) then
            print *, 'FAIL: dense_solve round-trip error =', &
                     maxval(abs(x - Phi_known))
            error stop
        end if
        if (.not. ieee_is_finite(rcond) .or. rcond <= PERIODIC_SOLVE_RCOND_TOLERANCE) then
            print *, 'FAIL: dense_solve invalid rcond =', rcond
            error stop
        end if
        if (.not. ieee_is_finite(backward_error) .or. &
            backward_error >= PERIODIC_SOLVE_RESIDUAL_TOLERANCE) then
            print *, 'FAIL: dense_solve backward error =', backward_error
            error stop
        end if
        print *, 'PASS: dense_solve round-trip, max err =', &
                 maxval(abs(x - Phi_known))
        print *, '  rcond =', rcond, ' backward error =', backward_error
    end subroutine test_dense_solve_roundtrip

    subroutine test_dense_solve_singular_rejection()
        use periodic_solve_m, only: PERIODIC_SOLVE_SINGULAR

        complex(dp) :: matrix(2, 2), rhs(2)
        integer :: info

        matrix = (0.0_dp, 0.0_dp)
        matrix(1, 1) = (1.0_dp, 0.0_dp)
        rhs = [(1.0_dp, 0.0_dp), (0.0_dp, 0.0_dp)]
        call dense_solve(matrix, rhs, info)
        if (info /= PERIODIC_SOLVE_SINGULAR) then
            print *, 'FAIL: singular dense solve status =', info
            error stop
        end if
        print *, 'PASS: exactly singular dense solve rejected with host status'
    end subroutine test_dense_solve_singular_rejection

    subroutine test_dense_solve_condition_rejection()
        use periodic_solve_m, only: PERIODIC_SOLVE_ILL_CONDITIONED

        complex(dp) :: matrix(2, 2), rhs(2)
        real(dp) :: rcond, backward_error
        integer :: info

        matrix = (0.0_dp, 0.0_dp)
        matrix(1, 1) = (1.0_dp, 0.0_dp)
        matrix(2, 2) = (1.0e-14_dp, 0.0_dp)
        rhs = [(1.0_dp, 0.0_dp), (1.0e-14_dp, 0.0_dp)]
        call dense_solve(matrix, rhs, info, rcond, backward_error)

        if (info /= PERIODIC_SOLVE_ILL_CONDITIONED) then
            print *, 'FAIL: nearly singular solve status =', info
            print *, '  rcond =', rcond, ' backward error =', backward_error
            error stop
        end if
        if (rcond >= PERIODIC_SOLVE_RCOND_TOLERANCE) then
            print *, 'FAIL: nearly singular fixture rcond =', rcond
            error stop
        end if
        print *, 'PASS: nearly singular dense solve rejected, rcond =', rcond
    end subroutine test_dense_solve_condition_rejection

    subroutine test_dense_solve_nonfinite_rejection()
        use, intrinsic :: ieee_arithmetic, only: ieee_quiet_nan, ieee_value
        use periodic_solve_m, only: PERIODIC_SOLVE_NONFINITE

        complex(dp) :: matrix(1, 1), rhs(1)
        real(dp) :: rcond, backward_error
        integer :: info

        matrix(1, 1) = cmplx(ieee_value(0.0_dp, ieee_quiet_nan), 0.0_dp, dp)
        rhs(1) = (1.0_dp, 0.0_dp)
        call dense_solve(matrix, rhs, info, rcond, backward_error)
        if (info /= PERIODIC_SOLVE_NONFINITE) then
            print *, 'FAIL: non-finite dense solve status =', info
            error stop
        end if
        print *, 'PASS: non-finite dense solve input rejected'
    end subroutine test_dense_solve_nonfinite_rejection

    !> (b) Inverse DFT correctness. A single nonzero coefficient at m=1 must
    !> reconstruct exp(i k_1 r); a two-mode case confirms linear superposition.
    subroutine test_inverse_dft()
        integer, parameter :: M = 3
        integer, parameter :: dim = 2 * M + 1
        integer, parameter :: nr = 5
        real(dp), parameter :: L = 4.0_dp
        complex(dp) :: Phi_m(dim), dPhi(nr), expected(nr)
        real(dp) :: r_out(nr), k1, k2
        integer :: i

        do i = 1, nr
            r_out(i) = -1.5_dp + 0.7_dp * real(i - 1, dp)
        end do

        k1 = 2.0_dp * pi * 1.0_dp / L
        k2 = 2.0_dp * pi * (-2.0_dp) / L

        ! Single-mode: Phi_m = delta_{m,1}.
        Phi_m = (0.0_dp, 0.0_dp)
        Phi_m(1 + M + 1) = (1.0_dp, 0.0_dp)
        dPhi = reconstruct_delta_phi(Phi_m, L, M, r_out)
        do i = 1, nr
            expected(i) = exp(com_unit * k1 * r_out(i))
        end do
        if (maxval(abs(abs(dPhi) - 1.0_dp)) >= 1.0e-12_dp) then
            print *, 'FAIL: single-mode |dPhi| /= 1, err =', &
                     maxval(abs(abs(dPhi) - 1.0_dp))
            error stop
        end if
        if (maxval(abs(dPhi - expected)) >= 1.0e-12_dp) then
            print *, 'FAIL: single-mode dPhi /= exp(i k_1 r), err =', &
                     maxval(abs(dPhi - expected))
            error stop
        end if
        print *, 'PASS: inverse DFT single-mode = exp(i k_1 r), max err =', &
                 maxval(abs(dPhi - expected))

        ! Two-mode superposition: Phi_1 = (1,0), Phi_{-2} = (0.5, -0.25).
        Phi_m = (0.0_dp, 0.0_dp)
        Phi_m(1 + M + 1)    = (1.0_dp, 0.0_dp)
        Phi_m((-2) + M + 1) = (0.5_dp, -0.25_dp)
        dPhi = reconstruct_delta_phi(Phi_m, L, M, r_out)
        do i = 1, nr
            expected(i) = (1.0_dp, 0.0_dp) * exp(com_unit * k1 * r_out(i)) &
                        + (0.5_dp, -0.25_dp) * exp(com_unit * k2 * r_out(i))
        end do
        if (maxval(abs(dPhi - expected)) >= 1.0e-12_dp) then
            print *, 'FAIL: two-mode superposition err =', &
                     maxval(abs(dPhi - expected))
            error stop
        end if
        print *, 'PASS: inverse DFT two-mode superposition, max err =', &
                 maxval(abs(dPhi - expected))
    end subroutine test_inverse_dft

    subroutine test_field_operator_oracle()
        integer, parameter :: M = 2, dim = 2 * M + 1
        complex(dp) :: Kphi(dim, dim), A_expected(dim, dim), psi(dim), rhs(dim)
        complex(dp), allocatable :: A(:,:), solved(:), drive_m(:)
        complex(dp) :: drive_values(4)
        real(dp) :: values(20), nodes(4), L, kappa, hat_nodes(3)
        real(dp) :: hat_mass(3), hat_laplace(3), current_correct, current_printed
        character(len=1024) :: line
        integer :: unit, ios, kind, icase, row, col, idx, seen

        Kphi = (0.0_dp, 0.0_dp)
        A_expected = (0.0_dp, 0.0_dp)
        psi = (0.0_dp, 0.0_dp)
        rhs = (0.0_dp, 0.0_dp)
        seen = 0
        open(newunit=unit, file='field_operator_matrices.dat', status='old', action='read', iostat=ios)
        if (ios /= 0) error stop 'cannot open field_operator_matrices.dat'
        do
            read(unit, '(A)', iostat=ios) line
            if (ios /= 0) exit
            if (len_trim(line) == 0 .or. line(1:1) == '#') cycle
            read(line, *, iostat=ios) kind, icase, row, col, values
            if (ios /= 0) error stop 'invalid field_operator_matrices.dat row'
            seen = seen + 1
            select case (kind)
            case (1)
                idx = row + M + 1
                L = values(1)
                kappa = values(2)
                Kphi(idx, idx) = -kappa * kappa / (4.0_dp * pi)
                A_expected(idx, idx) = cmplx(values(4), values(5), dp)
                psi(idx) = cmplx(values(8), values(9), dp)
                rhs(idx) = cmplx(values(12), values(13), dp)
            case (2)
                hat_nodes(col) = values(2)
                hat_mass(col) = values(3)
                hat_laplace(col) = values(4)
            case (4)
                current_correct = values(11)
                current_printed = values(13)
            end select
        end do
        close(unit)
        if (seen /= 14) error stop 'field/operator oracle row count mismatch'

        call assemble_periodic_operator(Kphi, L, M, A)
        if (maxval(abs(A - A_expected)) > &
            128.0_dp * epsilon(1.0_dp) * max(1.0_dp, maxval(abs(A_expected)))) &
            error stop 'field/operator oracle periodic matrix mismatch'
        call solve_with_derived_gauge(A, rhs, M, .false., solved, ios)
        if (ios /= PERIODIC_SOLVE_OK .or. maxval(abs(solved - psi)) > 1.0e-13_dp) &
            error stop 'field/operator oracle screened solution mismatch'

        nodes = [0.0_dp, 1.0_dp, 2.0_dp, 3.0_dp]
        drive_values = 2.0_dp + exp(com_unit * 2.0_dp * pi * nodes / 4.0_dp)
        call project_onto_periodic_basis(drive_values, nodes, 4.0_dp, M, drive_m)
        if (maxval(abs(drive_m - [cmplx(0.0_dp, 0.0_dp, dp), &
            cmplx(0.0_dp, 0.0_dp, dp), cmplx(2.0_dp, 0.0_dp, dp), &
            cmplx(1.0_dp, 0.0_dp, dp), cmplx(0.0_dp, 0.0_dp, dp)])) > 1.0e-13_dp) &
            error stop 'field/operator oracle drive transform mismatch'
        if (maxval(abs(hat_nodes - [0.0_dp, 1.0_dp, 3.0_dp])) > 1.0e-15_dp .or. &
            maxval(abs(hat_mass - [1.0_dp/6.0_dp, 1.0_dp, 1.0_dp/3.0_dp])) > 1.0e-15_dp .or. &
            maxval(abs(hat_laplace - [1.0_dp, -1.5_dp, 0.5_dp])) > 1.0e-15_dp) &
            error stop 'field/operator oracle hat row mismatch'
        if (abs(current_correct - current_printed) < 1.0e-6_dp) &
            error stop 'field/operator oracle current mutation not distinguished'
        print *, 'PASS: all 14 field/operator oracle rows consumed'
    end subroutine test_field_operator_oracle

    subroutine test_zero_mode_policy()
        integer, parameter :: M = 1, dim = 3
        complex(dp) :: A(dim, dim), rhs(dim)
        complex(dp), allocatable :: solution(:)
        integer :: info

        A = (0.0_dp, 0.0_dp)
        A(1, 1) = (-1.0_dp, 0.0_dp)
        A(3, 3) = (-1.0_dp, 0.0_dp)
        rhs = (0.0_dp, 0.0_dp)
        call solve_with_derived_gauge(A, rhs, M, .false., solution, info)
        if (info /= PERIODIC_SOLVE_GAUGE_REQUIRED) error stop 'missing vacuum gauge was accepted'
        call solve_with_derived_gauge(A, rhs, M, .true., solution, info)
        if (info /= PERIODIC_SOLVE_OK .or. abs(solution(M + 1)) > 1.0e-15_dp) &
            error stop 'declared mean-zero vacuum gauge failed'
        rhs(M + 1) = (1.0_dp, 0.0_dp)
        call solve_with_derived_gauge(A, rhs, M, .true., solution, info)
        if (info /= PERIODIC_SOLVE_INCOMPATIBLE_NULLSPACE) &
            error stop 'incompatible vacuum source was accepted'
        print *, 'PASS: screened/vacuum zero-mode policy'
    end subroutine test_zero_mode_policy

    subroutine test_manufactured_deviation()
        integer, parameter :: M = 12, dim = 25, n = 96
        real(dp), parameter :: L = 30.0_dp, dx_asis = 5.0_dp
        real(dp) :: nodes(n), residual
        logical :: projection_mask(n)
        complex(dp) :: Kphi(dim, dim), KB(dim, dim), Br_m(dim), psi_known(dim)
        complex(dp) :: phi_ref(n), d2phi_ref(n), total(n), expected(n)
        complex(dp), allocatable :: A(:,:), ref_m(:), d2ref_m(:), psi_m(:)
        integer :: i, info

        do i = 1, n
            nodes(i) = -0.5_dp * L + real(i - 1, dp) * L / real(n, dp)
        end do
        call build_physical_c2_lift(nodes, 0.0_dp, dx_asis, 10.0_dp, &
            cmplx(-1.0_dp, 0.4_dp, dp), cmplx(0.7_dp, -0.2_dp, dp), phi_ref, d2phi_ref)
        projection_mask = abs(nodes) <= dx_asis
        Kphi = (0.0_dp, 0.0_dp)
        KB = (0.0_dp, 0.0_dp)
        do i = 1, dim
            Kphi(i, i) = -0.75_dp**2 / (4.0_dp * pi)
            KB(i, i) = (1.0_dp, 0.0_dp)
        end do
        call assemble_periodic_operator(Kphi, L, M, A)
        call project_onto_periodic_basis(phi_ref, nodes, L, M, ref_m)
        call project_onto_periodic_basis(d2phi_ref, nodes, L, M, d2ref_m)
        psi_known = (0.0_dp, 0.0_dp)
        psi_known(M) = (0.2_dp, -0.1_dp)
        psi_known(M + 1) = (-0.3_dp, 0.25_dp)
        psi_known(M + 2) = (0.15_dp, 0.05_dp)
        Br_m = -(d2ref_m + 4.0_dp * pi * matmul(Kphi, ref_m) + &
            matmul(A, psi_known)) / (4.0_dp * pi)
        call solve_periodic_deviation(Kphi, KB, L, M, Br_m, nodes, phi_ref, d2phi_ref, &
            6.0e-2_dp, psi_m, residual, info, projection_mask=projection_mask)
        if (info /= PERIODIC_SOLVE_OK .or. maxval(abs(psi_m - psi_known)) > 2.0e-12_dp) &
            error stop 'manufactured periodic deviation was not recovered'
        total = reconstruct_delta_phi_from_lift(psi_m, L, M, nodes, phi_ref)
        expected = phi_ref + reconstruct_delta_phi(psi_known, L, M, nodes)
        if (maxval(abs(total - expected)) > 2.0e-12_dp) &
            error stop 'physical-space lift reconstruction mismatch'
        print *, 'PASS: manufactured complex deviation, projection residual =', residual
    end subroutine test_manufactured_deviation

    subroutine test_projection_rejection()
        integer, parameter :: M = 2, dim = 5, n = 96
        real(dp), parameter :: L = 30.0_dp
        real(dp) :: nodes(n), residual
        complex(dp) :: Kphi(dim, dim), KB(dim, dim), Br_m(dim)
        complex(dp) :: phi_ref(n), d2phi_ref(n)
        complex(dp), allocatable :: psi_m(:)
        integer :: i, info

        do i = 1, n
            nodes(i) = -0.5_dp * L + real(i - 1, dp) * L / real(n, dp)
        end do
        call build_physical_c2_lift(nodes, 0.0_dp, 5.0_dp, 10.0_dp, &
            cmplx(-1.0_dp, 0.0_dp, dp), cmplx(1.0_dp, 0.0_dp, dp), phi_ref, d2phi_ref)
        Kphi = (0.0_dp, 0.0_dp)
        KB = (0.0_dp, 0.0_dp)
        Br_m = (0.0_dp, 0.0_dp)
        call solve_periodic_deviation(Kphi, KB, L, M, Br_m, nodes, phi_ref, d2phi_ref, &
            1.0e-2_dp, psi_m, residual, info)
        if (info /= PERIODIC_SOLVE_PROJECTION_REJECTED) &
            error stop 'under-resolved physical lift was accepted'
        print *, 'PASS: under-resolved lift rejected, projection residual =', residual
    end subroutine test_projection_rejection

    !> (c) End-to-end sanity with the REAL assembled matrices. Mirrors
    !> test_periodic_assembly's setup: build the (m,n)=(-6,2) periodic plasma,
    !> assemble Kphi/KB with M=6, solve with a constant Br drive, and check the
    !> solve is non-singular, the response is finite and non-zero, and the
    !> reconstructed dPhi on the window grid is finite.
    subroutine test_end_to_end()
        use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
        use config_m, only: profiles_in_memory, nml_config_path
        use species_m, only: plasma, set_profiles_from_arrays
        use kim_base_m, only: kim_t
        use kim_mod_m, only: from_kim_factory_get_kim
        use kim_resonances_m, only: r_res
        use grid_m, only: rg_grid
        use periodic_background_m, only: build_periodic_plasma
        use periodic_assembly_m, only: assemble_periodic_matrices

        integer, parameter :: npts = 201
        integer, parameter :: M = 6

        real(dp) :: r_prof(npts), n_prof(npts), Te_prof(npts)
        real(dp) :: Ti_prof(npts), q_prof(npts), Er_prof(npts)
        class(kim_t), allocatable :: kim_instance

        complex(dp), allocatable :: Kphi(:,:), KB(:,:), Phi_m(:)
        complex(dp), allocatable :: dPhi(:)
        real(dp) :: rm, dx_asis, dx_tr, rho_L_rm, L
        integer :: n_rg, N, dim, info, i

        call make_test_profiles(npts, r_prof, n_prof, Te_prof, Ti_prof, &
                                q_prof, Er_prof)

        call write_test_namelist('./KIM_config_periodic_solve_test.nml')
        nml_config_path = './KIM_config_periodic_solve_test.nml'

        profiles_in_memory = .true.
        call kim_init()
        call set_profiles_from_arrays(r_prof, n_prof, Te_prof, Ti_prof, &
                                      q_prof, Er_prof, npts)

        call from_kim_factory_get_kim('electrostatic', kim_instance)
        call kim_instance%init()

        call prepare_resonances
        if (.not. (r_res > 0.0_dp)) then
            print *, 'FAIL: resonance not found, r_res = ', r_res
            error stop
        end if
        rm = r_res

        rho_L_rm = rho_L_near(rm)
        dx_asis =  5.0_dp * rho_L_rm
        dx_tr   = 10.0_dp * rho_L_rm
        if (dx_asis < 1.0e-3_dp) then
            dx_asis = 1.0_dp
            dx_tr   = 2.0_dp
        end if
        n_rg = 96
        L = 2.0_dp * (dx_asis + dx_tr)

        call build_periodic_plasma(rm, dx_asis, dx_tr, n_rg)
        N = rg_grid%npts_b

        call assemble_periodic_matrices(plasma, L, M, Kphi, KB)

        call solve_periodic(Kphi, KB, L, M, (1.0_dp, 0.0_dp), Phi_m, info)

        dim = 2 * M + 1
        if (info /= 0) then
            print *, 'FAIL: solve_periodic info /= 0:', info
            error stop
        end if
        if (size(Phi_m) /= dim) then
            print *, 'FAIL: Phi_m size /=', dim, ' got ', size(Phi_m)
            error stop
        end if
        do i = 1, dim
            if (.not. ieee_is_finite(real(Phi_m(i), dp)) .or. &
                .not. ieee_is_finite(aimag(Phi_m(i)))) then
                print *, 'FAIL: non-finite Phi_m at', i
                error stop
            end if
        end do
        if (maxval(abs(Phi_m)) <= 0.0_dp) then
            print *, 'FAIL: Phi_m all zero (constant Br should drive response)'
            error stop
        end if
        print *, 'PASS: end-to-end solve info==0, Phi_m finite & non-zero, ', &
                 'max|Phi_m| =', maxval(abs(Phi_m))

        ! Reconstruct on the window grid and assert finite.
        dPhi = reconstruct_delta_phi(Phi_m, L, M, rg_grid%xb)
        do i = 1, size(dPhi)
            if (.not. ieee_is_finite(real(dPhi(i), dp)) .or. &
                .not. ieee_is_finite(aimag(dPhi(i)))) then
                print *, 'FAIL: non-finite dPhi at', i
                error stop
            end if
        end do
        print *, 'PASS: end-to-end reconstruct dPhi on window finite, ', &
                 'max|dPhi| =', maxval(abs(dPhi))
    end subroutine test_end_to_end

    real(dp) function rho_L_near(r) result(val)
        use species_m, only: plasma
        real(dp), intent(in) :: r
        integer :: jn
        jn = minloc(abs(plasma%r_grid - r), dim=1)
        val = plasma%spec(0)%rho_L(jn)
    end function rho_L_near

    subroutine make_test_profiles(npts, r_prof, n_prof, Te_prof, Ti_prof, &
                                  q_prof, Er_prof)
        ! Analytic profiles on r in [2, 60] cm; q crosses 3 inside [10, 50] so
        ! the (m,n) = (-6, 2) mode (|m/n| = 3) has a resonance there.
        integer, intent(in) :: npts
        real(dp), intent(out) :: r_prof(npts), n_prof(npts), Te_prof(npts)
        real(dp), intent(out) :: Ti_prof(npts), q_prof(npts), Er_prof(npts)
        integer :: i

        do i = 1, npts
            r_prof(i) = 2.0_dp + 58.0_dp * real(i - 1, dp) / real(npts - 1, dp)
            n_prof(i) = 2.0e13_dp * (1.1_dp - r_prof(i) / 100.0_dp)
            Te_prof(i) = 1.0e3_dp * (1.2_dp - r_prof(i) / 100.0_dp)
            Ti_prof(i) = Te_prof(i)
            q_prof(i) = 1.5_dp + 2.5_dp * (r_prof(i) - 2.0_dp) / 58.0_dp
            Er_prof(i) = -0.5_dp
        end do
    end subroutine make_test_profiles

    subroutine write_test_namelist(path)
        ! Minimal electrostatic FokkerPlanck configuration; m_mode = -6,
        ! n_mode = 2 makes q resonant at q = 3, type_br_field = 12 (constant Br).
        character(len=*), intent(in) :: path
        integer :: iunit

        open(newunit=iunit, file=path, status='replace', action='write')
        write(iunit, '(A)') '&KIM_CONFIG'
        write(iunit, '(A)') ' number_of_ion_species = 1'
        write(iunit, '(A)') ' artificial_debye_case = 0'
        write(iunit, '(A)') " type_of_run = 'electrostatic'"
        write(iunit, '(A)') " collision_model = 'FokkerPlanck'"
        write(iunit, '(A)') ' read_species_from_namelist = .false.'
        write(iunit, '(A)') ' turn_off_ions = .false.'
        write(iunit, '(A)') ' turn_off_electrons = .false.'
        write(iunit, '(A)') " plasma_type = 'D'"
        write(iunit, '(A)') ' rescale_density = .false.'
        write(iunit, '(A)') ' number_density_rescale = 1.0'
        write(iunit, '(A)') ' ion_flr_scale_factor = 1.0'
        write(iunit, '(A)') '/'
        write(iunit, '(A)') '&WKB_DISPERSION'
        write(iunit, '(A)') '/'
        write(iunit, '(A)') '&KIM_IO'
        write(iunit, '(A)') " profile_location = './'"
        write(iunit, '(A)') " output_path = './out_periodic_solve_test/'"
        write(iunit, '(A)') ' hdf5_input = .false.'
        write(iunit, '(A)') ' hdf5_output = .false.'
        write(iunit, '(A)') ' log_level = 3'
        write(iunit, '(A)') ' data_verbosity = 0'
        write(iunit, '(A)') ' calculate_asymptotics = .false.'
        write(iunit, '(A)') " h5_out_file = ''"
        write(iunit, '(A)') ' write_diagnostics_dat = .false.'
        write(iunit, '(A)') '/'
        write(iunit, '(A)') '&KIM_SETUP'
        write(iunit, '(A)') ' btor = -18000.0'
        write(iunit, '(A)') ' R0 = 165.0'
        write(iunit, '(A)') ' m_mode = -6'
        write(iunit, '(A)') ' n_mode = 2'
        write(iunit, '(A)') ' omega = 0.0'
        write(iunit, '(A)') ' spline_base = 1'
        write(iunit, '(A)') ' type_br_field = 12'
        write(iunit, '(A)') ' collisions_off = .false.'
        write(iunit, '(A)') ' set_profiles_constant = 0'
        write(iunit, '(A)') ' bc_type = 3'
        write(iunit, '(A)') ' mphi_max = 0'
        write(iunit, '(A)') ' Br_boundary_re = 1.0'
        write(iunit, '(A)') ' Br_boundary_im = 0.0'
        write(iunit, '(A)') '/'
        write(iunit, '(A)') '&KIM_GRID'
        write(iunit, '(A)') " grid_spacing_rg = 'equidistant'"
        write(iunit, '(A)') " grid_spacing_xl = 'equidistant'"
        write(iunit, '(A)') ' l_space_dim = 32'
        write(iunit, '(A)') ' rg_space_dim = 32'
        write(iunit, '(A)') " theta_integration = 'GaussLegendre'"
        write(iunit, '(A)') " theta_integration_method = ''"
        write(iunit, '(A)') ' Larmor_skip_factor = 5'
        write(iunit, '(A)') ' gauss_int_nodes_Ntheta = 17'
        write(iunit, '(A)') ' gauss_int_nodes_Nx = 31'
        write(iunit, '(A)') ' gauss_int_nodes_Nxp = 30'
        write(iunit, '(A)') ' r_plas = 50.0'
        write(iunit, '(A)') ' r_min = 10.0'
        write(iunit, '(A)') ' width_res = 0.5'
        write(iunit, '(A)') ' ampl_res = 15.0'
        write(iunit, '(A)') ' hrmax_scaling = 1.0'
        write(iunit, '(A)') ' rkf45_atol = 1.0e-9'
        write(iunit, '(A)') ' rkf45_rtol = 1.0e-6'
        write(iunit, '(A)') ' kernel_taper_skip_threshold = 1.0e-6'
        write(iunit, '(A)') " quadpack_algorithm = 'QAG'"
        write(iunit, '(A)') ' quadpack_key = 6'
        write(iunit, '(A)') ' quadpack_limit = 500'
        write(iunit, '(A)') ' quadpack_epsabs = 1.0e-10'
        write(iunit, '(A)') ' quadpack_epsrel = 1.0e-10'
        write(iunit, '(A)') ' quadpack_use_u_substitution = .true.'
        write(iunit, '(A)') '/'
        write(iunit, '(A)') '&KIM_PROFILES'
        write(iunit, '(A)') '/'
        close(iunit)
    end subroutine write_test_namelist

end program test_periodic_solve
