module periodic_assembly_m
    !> Dense Fourier-matrix assembly for the forced-periodicity solver (Phase-2.4).
    !>
    !> Assembles the complex matrices K^{rhoPhi}_{m,m'} and K^{rhoB}_{m,m'} over
    !> one period from the Phase-2.3 periodic plasma and the Phase-1 fused kernel
    !> (flr2_fourier_kernel_m::hatG_rho_phi / hatG_rho_B).
    !>
    !> Discrete k-grid (memo eq. matelrhoPhi_discr): k_m = 2*pi*m/L for
    !> m = -M..M, mapped to matrix index im = m + M + 1, so the matrices are
    !> (2M+1) x (2M+1).
    !>
    !>   K^{rhoPhi}_{m,m'} = (2*pi/L) * INT_{rm-L/2}^{rm+L/2}
    !>                         G^{rhoPhi}(k_m, k_{m'}, r_g) dr_g.
    !>
    !> The periodic background makes the integrand L-periodic, and the Fourier
    !> phase exp(i(k_m - k_{m'}) r_g) is exactly L-periodic since
    !> (k_m - k_{m'})*L = 2*pi*(m - m'). The periodic trapezoidal (equidistant)
    !> rule is therefore spectrally accurate. The window grid rg_grid%xb is
    !> ENDPOINT-EXCLUSIVE over one period: grid_generate_equidistant sets h = L/N
    !> and xb(j) = r_lo + (j-1)*h (N = rg_grid%npts_b), so the N grid points are the
    !> N DISTINCT samples over [rm-L/2, rm+L/2) and xb(1)+L = rm+L/2 is NOT a node.
    !> The periodic trapezoidal rule therefore sums ALL N samples with equal weight
    !> -- no endpoint half-weights, no dropped sample:
    !>
    !>   K^{rhoPhi}_{m,m'} = (2*pi/N) * sum_{j=1}^{N}
    !>                         hatG_rho_phi(plasma, k_m, k_{m'}, j).
    !>
    !> This equals (2*pi/L)*h*sum with h = L/N; the L cancels. K^{rhoB} is
    !> assembled identically with hatG_rho_B.
    !>
    !> The FLR factors depend on BOTH k_m and k_{m'} (design section 4.2), so K is
    !> in general NOT Toeplitz (not constant along its diagonals).

    use KIM_kinds_m, only: dp
    use species_m, only: plasma_t

    implicit none
    private

    public :: assemble_periodic_matrices, k_of_m

contains

    !> Discrete radial wavenumber for Fourier index m over period L: k_m = 2*pi*m/L.
    pure real(dp) function k_of_m(m, L) result(k)
        use constants_m, only: pi
        integer, intent(in) :: m
        real(dp), intent(in) :: L
        k = 2.0_dp * pi * real(m, dp) / L
    end function k_of_m

    !> Assemble the dense complex Fourier matrices K^{rhoPhi} and K^{rhoB} over
    !> one period, using the periodic trapezoidal rule on the window grid.
    !>
    !> On entry `plasma` must be the periodic window plasma produced by
    !> periodic_background_m::build_periodic_plasma (rg_grid%xb holds the N
    !> equidistant boundary points that hatG_rho_phi / hatG_rho_B index).
    !>
    !> All four matrices are allocated here to (2M+1) x (2M+1), row index
    !> im = m+M+1, column index imp = m'+M+1.
    !>
    !> Kjphi and KjB are the parallel-current analogues, assembled with the same
    !> quadrature from hatG_j_phi / hatG_j_B. They are post-processing only: the
    !> solve uses Kphi/KB, and j_par follows from thesis (11.8) as
    !> Kjphi . Phi_m + KjB . Br_m. The j-kernels already carry their
    !> 1/(8*pi^2) Fourier normalization, so no further scaling is applied here.
    subroutine assemble_periodic_matrices(plasma, L, M, Kphi, KB, Kjphi, KjB)
        use grid_m, only: rg_grid
        use collisionless_fourier_kernel_m, only: configured_hatG_all
        use constants_m, only: pi

        type(plasma_t), intent(in) :: plasma
        real(dp), intent(in) :: L
        integer, intent(in) :: M
        complex(dp), allocatable, intent(out) :: Kphi(:,:), KB(:,:)
        complex(dp), allocatable, intent(out) :: Kjphi(:,:), KjB(:,:)

        integer :: N, dim, m_row, m_col, im, imp, j
        real(dp) :: k_m, k_mp, weight
        complex(dp) :: acc_phi, acc_B, acc_jphi, acc_jB
        complex(dp) :: point_phi, point_B, point_jphi, point_jB

        N = rg_grid%npts_b
        dim = 2 * M + 1
        allocate(Kphi(dim, dim), KB(dim, dim), Kjphi(dim, dim), KjB(dim, dim))

        ! Periodic trapezoidal weight. rg_grid%xb is EQUIDISTANT and ENDPOINT-
        ! EXCLUSIVE over one period (grid_generate_equidistant: h = L/N, xb(N) =
        ! r_hi - h), so the N grid points are the N DISTINCT samples over one
        ! period. Sum ALL N with equal weight (2*pi/L)*h = 2*pi/N -- dropping
        ! sample N and using 2*pi/(N-1) (as if endpoint-inclusive) is an O(1/N)
        ! error that caps r_g-quadrature convergence.
        weight = 2.0_dp * pi / real(N, dp)

        ! Each Fourier-space column is independent.  Parallelize the expensive
        ! radial quadrature while retaining a fixed summation order in j for
        ! bitwise-reproducible entries at a given compiler/architecture.
        !$omp parallel do default(none) schedule(static) &
        !$omp shared(M, L, N, weight, plasma, Kphi, KB, Kjphi, KjB) &
        !$omp private(m_col, k_mp, imp, m_row, k_m, im, j, &
        !$omp         acc_phi, acc_B, acc_jphi, acc_jB, &
        !$omp         point_phi, point_B, point_jphi, point_jB)
        do m_col = -M, M
            k_mp = k_of_m(m_col, L)
            imp = m_col + M + 1
            do m_row = -M, M
                k_m = k_of_m(m_row, L)
                im = m_row + M + 1

                acc_phi  = (0.0_dp, 0.0_dp)
                acc_B    = (0.0_dp, 0.0_dp)
                acc_jphi = (0.0_dp, 0.0_dp)
                acc_jB   = (0.0_dp, 0.0_dp)
                do j = 1, N
                    call configured_hatG_all(plasma, k_m, k_mp, j, &
                        point_phi, point_B, point_jphi, point_jB)
                    acc_phi  = acc_phi  + point_phi
                    acc_B    = acc_B    + point_B
                    acc_jphi = acc_jphi + point_jphi
                    acc_jB   = acc_jB   + point_jB
                end do

                Kphi(im, imp)  = weight * acc_phi
                KB(im, imp)    = weight * acc_B
                Kjphi(im, imp) = weight * acc_jphi
                KjB(im, imp)   = weight * acc_jB
            end do
        end do
        !$omp end parallel do
    end subroutine assemble_periodic_matrices

end module periodic_assembly_m
