module periodic_solve_m
    !> Dense solve of the periodic electrostatic system and inverse-DFT
    !> reconstruction of the potential (Phase-2.5).
    !>
    !> Over Fourier modes m = -M..M with k_m = 2*pi*m/L, the electrostatic
    !> forced-periodicity system (design section 3.4) is
    !>
    !>   [ D_m delta_{m,m'} + 4*pi K^{rhoPhi}_{m,m'} ] Phi_{m'}
    !>                                        = -4*pi K^{rhoB}_{m,m'} Br_{m'},
    !>
    !> with D_m = -k_m^2 the (diagonal) Fourier-space radial Poisson operator.
    !> This matches solve_poisson.f90, whose left-hand operator is
    !> A_mat = Laplace + 4*pi K_rho_phi (solve_poisson.f90:49): the FEM Laplace
    !> operator there (prepare_Laplace_matrix) is the RADIAL second derivative
    !> only, which in the Fourier basis maps to -k_m^2. There is NO perpendicular
    !> -k_s^2 term in that operator, so D_m = -k_m^2 as the design specifies.
    !> (Convention question flagged for the P2.7 cross-check: whether a physical
    !> perpendicular k_s^2 contribution belongs on the diagonal.)
    !>
    !> The right-hand side matches create_rhs_vector's sign convention
    !> (solve_poisson.f90:129, rhs = -4*pi * K_rho_B * Br). For type_br_field=12
    !> Br is constant over the window, so its Fourier coefficient is nonzero only
    !> at m' = 0:  Br_{m'} = Br_const * delta_{m',0}. The RHS therefore reduces to
    !>
    !>   b_m = -4*pi * Br_const * KB(m, m'=0) = -4*pi * Br_const * KB(:, M+1).
    !>
    !> The dense complex system is solved with LAPACK zgesv (same argument
    !> pattern as electromagnetic_solver.f90:243), then the potential is
    !> reconstructed on an output radial grid as
    !>
    !>   dPhi(r) = sum_{m=-M}^{M} Phi_m exp(i k_m r).

    use KIM_kinds_m, only: dp

    implicit none
    private

    public :: solve_periodic, reconstruct_delta_phi, reconstruct_delta_phi_derivative
    public :: dense_solve, reconstruct_jpar, reconstruct_jrad
    public :: solve_periodic_em, core_integral_weights

    integer, parameter, public :: PERIODIC_PSI_EXTRAPOLATED = 1
    integer, parameter, public :: PERIODIC_PSI_RESONANCE = 2

contains

    !> Solve the periodic electrostatic system for the Fourier coefficients
    !> Phi_m of the potential, given the Phase-2.4 matrices Kphi/KB, the period
    !> L, the mode cutoff M, and the constant drive amplitude Br_const.
    subroutine solve_periodic(Kphi, KB, L, M, Br_const, Phi_m, info, &
                              KBparallel, Bparallel_const)
        use constants_m, only: pi

        complex(dp), intent(in) :: Kphi(:,:), KB(:,:)
        real(dp), intent(in) :: L
        integer, intent(in) :: M
        complex(dp), intent(in) :: Br_const
        complex(dp), allocatable, intent(out) :: Phi_m(:)
        integer, intent(out) :: info
        complex(dp), intent(in), optional :: KBparallel(:, :)
        complex(dp), intent(in), optional :: Bparallel_const

        complex(dp), allocatable :: A(:,:), b(:)
        real(dp) :: k_m
        integer :: dim, im, m_row

        dim = 2 * M + 1
        allocate(A(dim, dim), b(dim))

        ! A = 4*pi Kphi with the diagonal Fourier-space radial operator
        ! D_m = -k_m^2 (k_m = 2*pi*m/L) added on the diagonal.
        A = 4.0_dp * pi * Kphi
        do m_row = -M, M
            im = m_row + M + 1
            k_m = 2.0_dp * pi * real(m_row, dp) / L
            A(im, im) = A(im, im) - k_m * k_m
        end do

        if (present(KBparallel) .neqv. present(Bparallel_const)) then
            error stop 'solve_periodic requires Bparallel matrix and drive together'
        end if

        ! Constant prescribed drives project onto the m'=0 column only.
        b = -4.0_dp * pi * Br_const * KB(:, M + 1)
        if (present(KBparallel)) then
            if (Bparallel_const /= (0.0_dp, 0.0_dp)) &
                b = b - 4.0_dp * pi * Bparallel_const * KBparallel(:, M + 1)
        end if

        call dense_solve(A, b, info)

        Phi_m = b
    end subroutine solve_periodic

    !> Coupled Poisson-Ampere solve on the period (electromagnetic periodic layer).
    !>
    !> Unknowns: Phi_m and A_m (A = A_parallel), m = -M..M; Br_m = i ks A_m with
    !> the constant ks = ks(r_m) of the resonant surface (window << r).
    !>
    !> Poisson rows (as solve_periodic, Br no longer constant):
    !>   -k_m^2 Phi_m + 4 pi [Kphi Phi + KB Br]_m = -4 pi Bpar KBpar(m, 0).
    !> Ampere rows m /= 0 (nabla_perp^2 A = -(4 pi/c) j, KIM convention):
    !>   -(k_m^2 + ks^2) A_m + (4 pi/c) (j_m - j_0 e_m) = 0,
    !>   e_m = exp(-i k_m (r_m + L/2)).
    !> The net window current L j_0 is returned by a compensating sheet at the
    !> period boundary r_m + L/2 (inside the transition band). Inside the core
    !> the plasma current then produces its local field only (kink convolution,
    !> up to O((ks L)^2)); its far field belongs to the outer region and enters
    !> through the Delta' matching of the k_r = 0 amplitude in QL-Balance.
    !> Row m = 0 fixes the outer amplitude psi = Br_ref/(i ks):
    !>   PERIODIC_PSI_EXTRAPOLATED: even part of A at r_m +- x_a extrapolated
    !>     back to r_m along the kink of the core current,
    !>     [A(r_m+x_a) + A(r_m-x_a)]/2 + (2 pi/c) x_a J_core = psi,
    !>     J_core = int_{r_m-x_a}^{r_m+x_a} j dr (asymptotic matching value,
    !>     requires the current layer to lie inside +-x_a);
    !>   PERIODIC_PSI_RESONANCE: A(r_m) = psi.
    !> With j = 0 both give Br = Br_ref (the electrostatic periodic solution).
    subroutine solve_periodic_em(Kphi, KB, Kjphi, KjB, L, M, r_m, x_a, ks, Br_ref, &
            psi_definition, Phi_m, Br_m, info, KBparallel, KjBparallel, Bparallel_const)
        use constants_m, only: pi, sol, com_unit

        complex(dp), intent(in) :: Kphi(:,:), KB(:,:), Kjphi(:,:), KjB(:,:)
        real(dp), intent(in) :: L, r_m, x_a, ks
        integer, intent(in) :: M, psi_definition
        complex(dp), intent(in) :: Br_ref
        complex(dp), allocatable, intent(out) :: Phi_m(:), Br_m(:)
        integer, intent(out) :: info
        complex(dp), intent(in), optional :: KBparallel(:,:), KjBparallel(:,:)
        complex(dp), intent(in), optional :: Bparallel_const

        complex(dp), allocatable :: A(:,:), b(:), e(:), w_core(:), j_row(:)
        complex(dp) :: alpha, bpar, jbpar0
        real(dp) :: k
        integer :: dim, mm, im, i0, row

        if ((present(KBparallel) .neqv. present(Bparallel_const)) .or. &
                (present(KjBparallel) .neqv. present(Bparallel_const))) then
            error stop 'solve_periodic_em requires Bparallel matrices and drive together'
        end if
        if (abs(ks) <= tiny(1.0_dp)) error stop 'solve_periodic_em requires ks /= 0'
        if (psi_definition /= PERIODIC_PSI_EXTRAPOLATED .and. &
                psi_definition /= PERIODIC_PSI_RESONANCE) then
            error stop 'solve_periodic_em: unknown psi definition'
        end if
        if (psi_definition == PERIODIC_PSI_EXTRAPOLATED .and. &
                (x_a <= 0.0_dp .or. 2.0_dp*x_a >= L)) then
            error stop 'solve_periodic_em: x_a must lie in (0, L/2)'
        end if

        dim = 2*M + 1
        i0 = M + 1
        alpha = com_unit*ks
        bpar = (0.0_dp, 0.0_dp)
        if (present(Bparallel_const)) bpar = Bparallel_const
        allocate(A(2*dim, 2*dim), b(2*dim), e(dim), j_row(2*dim))
        A = (0.0_dp, 0.0_dp)
        b = (0.0_dp, 0.0_dp)
        do mm = -M, M
            e(mm + i0) = exp(-com_unit*2.0_dp*pi*real(mm, dp)/L*(r_m + 0.5_dp*L))
        end do

        ! Poisson rows
        A(1:dim, 1:dim) = 4.0_dp*pi*Kphi
        A(1:dim, dim+1:2*dim) = 4.0_dp*pi*alpha*KB
        do mm = -M, M
            im = mm + i0
            k = 2.0_dp*pi*real(mm, dp)/L
            A(im, im) = A(im, im) - k*k
        end do
        if (bpar /= (0.0_dp, 0.0_dp)) b(1:dim) = -4.0_dp*pi*bpar*KBparallel(:, i0)

        ! Ampere rows m /= 0
        do mm = -M, M
            if (mm == 0) cycle
            im = mm + i0
            row = dim + im
            k = 2.0_dp*pi*real(mm, dp)/L
            A(row, 1:dim) = 4.0_dp*pi/sol*(Kjphi(im, :) - e(im)*Kjphi(i0, :))
            A(row, dim+1:2*dim) = 4.0_dp*pi/sol*alpha*(KjB(im, :) - e(im)*KjB(i0, :))
            A(row, row) = A(row, row) - (k*k + ks*ks)
            if (bpar /= (0.0_dp, 0.0_dp)) b(row) = -4.0_dp*pi/sol*bpar* &
                (KjBparallel(im, i0) - e(im)*KjBparallel(i0, i0))
        end do

        ! Row m = 0: outer amplitude psi = Br_ref/(i ks)
        row = dim + i0
        b(row) = Br_ref/alpha
        select case (psi_definition)
        case (PERIODIC_PSI_RESONANCE)
            do mm = -M, M
                A(row, dim + mm + i0) = exp(com_unit*2.0_dp*pi*real(mm, dp)/L*r_m)
            end do
        case (PERIODIC_PSI_EXTRAPOLATED)
            do mm = -M, M
                k = 2.0_dp*pi*real(mm, dp)/L
                A(row, dim + mm + i0) = 0.5_dp*(exp(com_unit*k*(r_m + x_a)) &
                    + exp(com_unit*k*(r_m - x_a)))
            end do
            ! J_core = w_core . j_m, j_m = Kjphi Phi + alpha KjB A + bpar KjBpar(:,0)
            w_core = core_integral_weights(L, M, r_m, x_a)
            j_row(1:dim) = matmul(w_core, Kjphi)
            j_row(dim+1:2*dim) = alpha*matmul(w_core, KjB)
            A(row, :) = A(row, :) + 2.0_dp*pi/sol*x_a*j_row
            if (bpar /= (0.0_dp, 0.0_dp)) then
                jbpar0 = sum(w_core*KjBparallel(:, i0))*bpar
                b(row) = b(row) - 2.0_dp*pi/sol*x_a*jbpar0
            end if
        end select

        call dense_solve(A, b, info)
        if (info /= 0) return
        Phi_m = b(1:dim)
        Br_m = alpha*b(dim+1:2*dim)
    end subroutine solve_periodic_em

    !> Weights w_m with int_{r_m-x_a}^{r_m+x_a} f dr = sum_m w_m f_m for
    !> f = sum_m f_m exp(i k_m r).
    function core_integral_weights(L, M, r_m, x_a) result(w)
        use constants_m, only: pi, com_unit
        real(dp), intent(in) :: L, r_m, x_a
        integer, intent(in) :: M
        complex(dp) :: w(2*M + 1)
        real(dp) :: k
        integer :: mm

        do mm = -M, M
            if (mm == 0) then
                w(mm + M + 1) = 2.0_dp*x_a
            else
                k = 2.0_dp*pi*real(mm, dp)/L
                w(mm + M + 1) = exp(com_unit*k*r_m)*2.0_dp*sin(k*x_a)/k
            end if
        end do
    end function core_integral_weights

    !> Solve the dense complex system A x = b in place via LAPACK zgesv, using
    !> the same argument pattern as electromagnetic_solver.f90:243. On entry b
    !> is the RHS; on exit b holds the solution and A is LU-overwritten.
    !> info == 0 on success.
    subroutine dense_solve(A, b, info)
        complex(dp), intent(inout) :: A(:,:), b(:)
        integer, intent(out) :: info

        integer, allocatable :: ipiv(:)
        integer :: dim

        dim = size(A, 1)
        allocate(ipiv(dim))
        call zgesv(dim, 1, A, dim, ipiv, b, dim, info)
        deallocate(ipiv)
    end subroutine dense_solve

    !> Reconstruct the parallel current density perturbation j_par(r) on the
    !> output grid r_out, given the solved potential coefficients Phi_m and the
    !> constant drive amplitude Br_const.
    !>
    !> Thesis (11.7) defines the current density directly from the two fields,
    !>
    !>   delta_j_par(x) = K^{jPhi}(x,x') delta_Phi(x') + K^{jB}(x,x') delta_Br(x'),
    !>
    !> which in the Fourier basis over m = -M..M is
    !>
    !>   j_m = sum_{m'} [ Kjphi(m,m') Phi_{m'} + KjB(m,m') Br_{m'} ].
    !>
    !> For type_br_field = 12 the drive is constant over the window, so
    !> Br_{m'} = Br_const * delta_{m',0} and the KjB term collapses to the
    !> m' = 0 column (index M+1) -- the same convention as solve_periodic's RHS.
    !>
    !> hatG_j_phi / hatG_j_B already include the 1/(8*pi^2) Fourier
    !> normalization from equations (14.5)/(14.6). Kjphi/KjB therefore enter
    !> this reconstruction without further scaling.
    !>
    !> j_par is then inverse-DFT'd exactly as the potential is:
    !>   j_par(r) = sum_{m=-M}^{M} j_m exp(i k_m r).
    !> With Br_spectrum (electromagnetic periodic solve) the full Br_m replaces
    !> the constant drive Br_const.
    function reconstruct_jpar(Kjphi, KjB, Phi_m, Br_const, L, M, r_out, &
            KjBparallel, Bparallel_const, Br_spectrum) result(jpar)
        complex(dp), intent(in) :: Kjphi(:,:), KjB(:,:)
        complex(dp), intent(in) :: Phi_m(:)
        complex(dp), intent(in) :: Br_const
        real(dp), intent(in) :: L, r_out(:)
        integer, intent(in) :: M
        complex(dp), intent(in), optional :: KjBparallel(:,:)
        complex(dp), intent(in), optional :: Bparallel_const
        complex(dp), intent(in), optional :: Br_spectrum(:)
        complex(dp) :: jpar(size(r_out))

        complex(dp), allocatable :: j_m(:)

        if (present(KjBparallel) .neqv. present(Bparallel_const)) then
            error stop 'reconstruct_jpar requires Bparallel matrix and drive together'
        end if

        if (present(Br_spectrum)) then
            j_m = matmul(Kjphi, Phi_m) + matmul(KjB, Br_spectrum)
        else
            j_m = matmul(Kjphi, Phi_m) + Br_const * KjB(:, M + 1)
        end if
        if (present(KjBparallel)) then
            if (Bparallel_const /= (0.0_dp, 0.0_dp)) &
                j_m = j_m + Bparallel_const * KjBparallel(:, M + 1)
        end if

        jpar = reconstruct_delta_phi(j_m, L, M, r_out)
    end function reconstruct_jpar

    function reconstruct_delta_phi_derivative(Phi_m, L, M, r_out) result(dphi)
        use constants_m, only: pi, com_unit
        complex(dp), intent(in) :: Phi_m(:)
        real(dp), intent(in) :: L
        integer, intent(in) :: M
        real(dp), intent(in) :: r_out(:)
        complex(dp) :: dphi(size(r_out))
        integer :: i, mm, im
        real(dp) :: k
        dphi = (0.0_dp, 0.0_dp)
        do i = 1, size(r_out)
            do mm = -M, M
                im = mm + M + 1
                k = 2.0_dp*pi*real(mm,dp)/L
                dphi(i) = dphi(i) + com_unit*k*Phi_m(im)*exp(com_unit*k*r_out(i))
            end do
        end do
    end function reconstruct_delta_phi_derivative

    !> Reconstruct j_rad from its potential and magnetic-field response matrices.
    !> Constant Br and optional Bparallel drives occupy only Fourier column m'=0,
    !> as in the Poisson right-hand side and parallel-current reconstruction.
    function reconstruct_jrad(Kjrphi, KjrB, Phi_m, Br_const, L, M, r_out, &
            KjrBparallel, Bparallel_const, Br_spectrum) result(jrad)
        complex(dp), intent(in) :: Kjrphi(:,:), KjrB(:,:)
        complex(dp), intent(in) :: Phi_m(:)
        complex(dp), intent(in) :: Br_const
        real(dp), intent(in) :: L, r_out(:)
        integer, intent(in) :: M
        complex(dp), intent(in), optional :: KjrBparallel(:,:), Bparallel_const
        complex(dp), intent(in), optional :: Br_spectrum(:)
        complex(dp) :: jrad(size(r_out))

        complex(dp), allocatable :: j_m(:)

        if (present(KjrBparallel) .neqv. present(Bparallel_const)) then
            error stop 'KjrBparallel and Bparallel_const must be supplied together'
        end if
        if (present(Br_spectrum)) then
            j_m = matmul(Kjrphi, Phi_m) + matmul(KjrB, Br_spectrum)
        else
            j_m = matmul(Kjrphi, Phi_m) + Br_const * KjrB(:, M + 1)
        end if
        if (present(KjrBparallel)) then
            j_m = j_m + Bparallel_const * KjrBparallel(:, M + 1)
        end if
        jrad = reconstruct_delta_phi(j_m, L, M, r_out)
    end function reconstruct_jrad

    !> Inverse DFT: reconstruct dPhi(r) = sum_{m=-M}^{M} Phi_m exp(i k_m r)
    !> on the output radial grid r_out, with k_m = 2*pi*m/L. Phi_m is indexed
    !> im = m + M + 1 (size 2M+1).
    function reconstruct_delta_phi(Phi_m, L, M, r_out) result(dPhi)
        use constants_m, only: pi, com_unit

        complex(dp), intent(in) :: Phi_m(:)
        real(dp), intent(in) :: L, r_out(:)
        integer, intent(in) :: M
        complex(dp) :: dPhi(size(r_out))

        real(dp) :: k_m
        integer :: i, m_row, im

        dPhi = (0.0_dp, 0.0_dp)
        do m_row = -M, M
            im = m_row + M + 1
            k_m = 2.0_dp * pi * real(m_row, dp) / L
            do i = 1, size(r_out)
                dPhi(i) = dPhi(i) + Phi_m(im) * exp(com_unit * k_m * r_out(i))
            end do
        end do
    end function reconstruct_delta_phi

end module periodic_solve_m
