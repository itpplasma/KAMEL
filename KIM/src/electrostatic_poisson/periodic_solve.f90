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

    public :: solve_periodic, reconstruct_delta_phi, reconstruct_delta_phi_derivative, dense_solve, reconstruct_jpar

contains

    !> Solve the periodic electrostatic system for the Fourier coefficients
    !> Phi_m of the potential, given the Phase-2.4 matrices Kphi/KB, the period
    !> L, the mode cutoff M, and the constant drive amplitude Br_const.
    subroutine solve_periodic(Kphi, KB, L, M, Br_const, Phi_m, info)
        use constants_m, only: pi

        complex(dp), intent(in) :: Kphi(:,:), KB(:,:)
        real(dp), intent(in) :: L
        integer, intent(in) :: M
        complex(dp), intent(in) :: Br_const
        complex(dp), allocatable, intent(out) :: Phi_m(:)
        integer, intent(out) :: info

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

        ! b = -4*pi Br_const KB(:, m'=0). The constant drive projects onto the
        ! m' = 0 column only (column index M+1).
        b = -4.0_dp * pi * Br_const * KB(:, M + 1)

        call dense_solve(A, b, info)

        Phi_m = b
    end subroutine solve_periodic

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
    function reconstruct_jpar(Kjphi, KjB, Phi_m, Br_const, L, M, r_out) result(jpar)
        complex(dp), intent(in) :: Kjphi(:,:), KjB(:,:)
        complex(dp), intent(in) :: Phi_m(:)
        complex(dp), intent(in) :: Br_const
        real(dp), intent(in) :: L, r_out(:)
        integer, intent(in) :: M
        complex(dp) :: jpar(size(r_out))

        complex(dp), allocatable :: j_m(:)

        j_m = matmul(Kjphi, Phi_m) + Br_const * KjB(:, M + 1)

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
