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

    public :: assemble_periodic_matrices, assemble_periodic_bparallel_matrices
    public :: k_of_m

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
    !> All six aggregate matrices are allocated here to (2M+1) x (2M+1), row index
    !> im = m+M+1, column index imp = m'+M+1.
    !>
    !> Kjphi and KjB are the parallel-current analogues, assembled with the same
    !> quadrature from hatG_j_phi / hatG_j_B. They are post-processing only: the
    !> solve uses Kphi/KB, and j_par follows from thesis (11.8) as
    !> Kjphi . Phi_m + KjB . Br_m. The j-kernels already carry their
    !> 1/(8*pi^2) Fourier normalization, so no further scaling is applied here.
    !> Kjrphi and KjrB similarly reconstruct the radial current as a
    !> post-processing quantity. Prescribed B_parallel columns are assembled
    !> separately so a zero compression drive preserves this established path.
    subroutine assemble_periodic_matrices(plasma, L, M, Kphi, KB, Kjphi, KjB, &
            Kjrphi, KjrB, Kjphi_species, KjB_species, Kphi_species, KB_species)
        use grid_m, only: rg_grid
        use collisionless_fourier_kernel_m, only: configured_hatG_all
        use constants_m, only: pi
        use radial_current_fourier_kernel_m, only: hatG_jrad_phi_br
        use config_m, only: ion_collision_model

        type(plasma_t), intent(in) :: plasma
        real(dp), intent(in) :: L
        integer, intent(in) :: M
        complex(dp), allocatable, intent(out) :: Kphi(:,:), KB(:,:)
        complex(dp), allocatable, intent(out) :: Kjphi(:,:), KjB(:,:)
        complex(dp), allocatable, intent(out), optional :: Kjrphi(:,:), KjrB(:,:)
        complex(dp), allocatable, intent(out), optional :: Kjphi_species(:,:,:)
        complex(dp), allocatable, intent(out), optional :: KjB_species(:,:,:)
        complex(dp), allocatable, intent(out), optional :: Kphi_species(:,:,:)
        complex(dp), allocatable, intent(out), optional :: KB_species(:,:,:)

        integer :: N, dim, m_row, m_col, im, imp, j
        real(dp) :: k_m, k_mp, weight
        complex(dp) :: acc_phi, acc_B, acc_jphi, acc_jB, acc_jrphi, acc_jrB
        complex(dp) :: point_phi, point_B, point_jphi, point_jB, point_jrphi, point_jrB
        complex(dp) :: acc_jphi_species(0:plasma%n_species - 1)
        complex(dp) :: acc_jB_species(0:plasma%n_species - 1)
        complex(dp) :: point_jphi_species(0:plasma%n_species - 1)
        complex(dp) :: point_jB_species(0:plasma%n_species - 1)
        complex(dp) :: acc_phi_species(0:plasma%n_species - 1)
        complex(dp) :: acc_B_species(0:plasma%n_species - 1)
        complex(dp) :: point_phi_species(0:plasma%n_species - 1)
        complex(dp) :: point_B_species(0:plasma%n_species - 1)
        logical :: want_current_species, want_charge_species, want_any_species
        logical :: want_radial_current

        N = rg_grid%npts_b
        dim = 2 * M + 1
        allocate(Kphi(dim, dim), KB(dim, dim), Kjphi(dim, dim), KjB(dim, dim))
        if (present(Kjrphi) .neqv. present(KjrB)) then
            error stop 'assemble_periodic_matrices requires both radial-current matrices'
        end if
        want_radial_current = present(Kjrphi)
        if (want_radial_current) then
            allocate(Kjrphi(dim, dim), KjrB(dim, dim))
        end if
        if (present(Kjphi_species) .neqv. present(KjB_species)) then
            error stop 'assemble_periodic_matrices requires both species-current matrices'
        end if
        if (present(Kphi_species) .neqv. present(KB_species)) then
            error stop 'assemble_periodic_matrices requires both species-charge matrices'
        end if
        want_current_species = present(Kjphi_species)
        want_charge_species = present(Kphi_species)
        want_any_species = want_current_species .or. want_charge_species
        if (want_current_species) then
            allocate(Kjphi_species(dim, dim, 0:plasma%n_species - 1))
            allocate(KjB_species(dim, dim, 0:plasma%n_species - 1))
        end if
        if (want_charge_species) then
            allocate(Kphi_species(dim, dim, 0:plasma%n_species - 1))
            allocate(KB_species(dim, dim, 0:plasma%n_species - 1))
        end if

        ! Periodic trapezoidal weight. rg_grid%xb is EQUIDISTANT and ENDPOINT-
        ! EXCLUSIVE over one period (grid_generate_equidistant: h = L/N, xb(N) =
        ! r_hi - h), so the N grid points are the N DISTINCT samples over one
        ! period. Sum ALL N with equal weight (2*pi/L)*h = 2*pi/N -- dropping
        ! sample N and using 2*pi/(N-1) (as if endpoint-inclusive) is an O(1/N)
        ! error that caps r_g-quadrature convergence.
        weight = 2.0_dp * pi / real(N, dp)

        if (trim(ion_collision_model) == 'FokkerPlanck') then
            call assemble_fp_symmetric(plasma, L, M, N, weight, want_radial_current, &
                want_current_species, want_charge_species, Kphi, KB, Kjphi, KjB, &
                Kjrphi, KjrB, Kjphi_species, KjB_species, Kphi_species, KB_species)
            return
        end if

        ! Each Fourier-space column is independent.  Parallelize the expensive
        ! radial quadrature while retaining a fixed summation order in j for
        ! bitwise-reproducible entries at a given compiler/architecture.
        !$omp parallel do default(none) schedule(static) &
        !$omp shared(M, L, N, weight, plasma, Kphi, KB, Kjphi, KjB, Kjrphi, KjrB, &
        !$omp        want_radial_current, want_current_species, want_charge_species, want_any_species, &
        !$omp        Kjphi_species, KjB_species, Kphi_species, KB_species) &
        !$omp private(m_col, k_mp, imp, m_row, k_m, im, j, &
        !$omp         acc_phi, acc_B, acc_jphi, acc_jB, acc_jrphi, acc_jrB, &
        !$omp         acc_jphi_species, acc_jB_species, acc_phi_species, acc_B_species, &
        !$omp         point_phi, point_B, point_jphi, point_jB, point_jrphi, point_jrB, &
        !$omp         point_jphi_species, point_jB_species, point_phi_species, point_B_species)
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
                if (want_radial_current) then
                    acc_jrphi = (0.0_dp, 0.0_dp)
                    acc_jrB = (0.0_dp, 0.0_dp)
                end if
                if (want_any_species) then
                    acc_phi_species = (0.0_dp, 0.0_dp)
                    acc_B_species = (0.0_dp, 0.0_dp)
                    acc_jphi_species = (0.0_dp, 0.0_dp)
                    acc_jB_species = (0.0_dp, 0.0_dp)
                end if
                do j = 1, N
                    if (want_any_species) then
                        call configured_hatG_all(plasma, k_m, k_mp, j, &
                            point_phi, point_B, point_jphi, point_jB, &
                            point_jphi_species, point_jB_species, &
                            point_phi_species, point_B_species)
                        acc_phi_species = acc_phi_species + point_phi_species
                        acc_B_species = acc_B_species + point_B_species
                        acc_jphi_species = acc_jphi_species + point_jphi_species
                        acc_jB_species = acc_jB_species + point_jB_species
                    else
                        call configured_hatG_all(plasma, k_m, k_mp, j, &
                            point_phi, point_B, point_jphi, point_jB)
                    end if
                    acc_phi  = acc_phi  + point_phi
                    acc_B    = acc_B    + point_B
                    acc_jphi = acc_jphi + point_jphi
                    acc_jB   = acc_jB   + point_jB
                    if (want_radial_current) then
                        call hatG_jrad_phi_br(plasma, k_m, k_mp, j, point_jrphi, point_jrB)
                        acc_jrphi = acc_jrphi + point_jrphi
                        acc_jrB = acc_jrB + point_jrB
                    end if
                end do

                Kphi(im, imp)  = weight * acc_phi
                KB(im, imp)    = weight * acc_B
                Kjphi(im, imp) = weight * acc_jphi
                KjB(im, imp)   = weight * acc_jB
                if (want_radial_current) then
                    Kjrphi(im, imp) = weight * acc_jrphi
                    KjrB(im, imp) = weight * acc_jrB
                end if
                if (want_current_species) then
                    Kjphi_species(im, imp, :) = weight * acc_jphi_species
                    KjB_species(im, imp, :) = weight * acc_jB_species
                end if
                if (want_charge_species) then
                    Kphi_species(im, imp, :) = weight * acc_phi_species
                    KB_species(im, imp, :) = weight * acc_B_species
                end if
            end do
        end do
        !$omp end parallel do
    end subroutine assemble_periodic_matrices

    !> Fokker-Planck assembly. The FP cores depend on k_m and k_m' only through
    !> b_plus and b_cross, which are invariant under sign changes and under
    !> exchange of k_m and k_m'. Each unordered pair of |m| values is therefore
    !> evaluated once per node and added to all of its (up to eight) signed
    !> entries with the phase exp(-i k_m r) exp(i k_m' r). The radial current
    !> depends on the signs, but its Bessel tables do not and are shared the
    !> same way. Every entry keeps its summation order in j; against the
    !> generic loop the values change only by rounding.
    subroutine assemble_fp_symmetric(plasma, L, M, N, weight, want_radial_current, &
            want_current_species, want_charge_species, Kphi, KB, Kjphi, KjB, &
            Kjrphi, KjrB, Kjphi_species, KjB_species, Kphi_species, KB_species)
        use grid_m, only: rg_grid
        use config_m, only: artificial_debye_case, turn_off_electrons, turn_off_ions
        use constants_m, only: pi, com_unit, flr_damping_cutoff
        use flr2_fourier_kernel_m, only: flr_arg_pair_sp, scaled_bessel_pair, &
            core_rho_phi_sp_bessel, core_rho_B_sp_bessel, core_j_phi_sp_bessel, &
            core_j_B_sp_bessel
        use radial_current_fourier_kernel_m, only: hatG_jrad_phi_br, jrad_bessel_tables
        use setup_m, only: mphi_max

        type(plasma_t), intent(in) :: plasma
        real(dp), intent(in) :: L, weight
        integer, intent(in) :: M, N
        logical, intent(in) :: want_radial_current, want_current_species
        logical, intent(in) :: want_charge_species
        complex(dp), intent(inout) :: Kphi(:,:), KB(:,:), Kjphi(:,:), KjB(:,:)
        complex(dp), intent(inout), optional :: Kjrphi(:,:), KjrB(:,:)
        complex(dp), intent(inout), optional :: Kjphi_species(:,:,:), KjB_species(:,:,:)
        complex(dp), intent(inout), optional :: Kphi_species(:,:,:), KB_species(:,:,:)

        integer, parameter :: max_entries = 8
        integer :: mag_a, mag_b, j, sp, e, mode, n_entries, im, imp
        integer :: rows(max_entries), cols(max_entries)
        real(dp) :: k_a, k_b, bplus, bcross, sI0, sIm1
        complex(dp), allocatable :: phase(:,:)
        complex(dp) :: ph, point_jrphi, point_jrB
        complex(dp) :: core(4), core_sp(4, 0:plasma%n_species - 1)
        complex(dp) :: acc(4, max_entries), acc_rad(2, max_entries)
        complex(dp) :: acc_sp(4, 0:plasma%n_species - 1, max_entries)
        real(dp) :: tables(0:mphi_max + 3, 0:plasma%n_species - 1)
        logical :: want_species, fp_active

        want_species = want_current_species .or. want_charge_species
        fp_active = artificial_debye_case == 0 .or. artificial_debye_case == 2

        allocate(phase(-M:M, N))
        do j = 1, N
            do mode = -M, M
                phase(mode, j) = exp(-com_unit * k_of_m(mode, L) * rg_grid%xb(j))
            end do
        end do

        !$omp parallel do schedule(dynamic) default(shared) &
        !$omp private(mag_a, mag_b, j, sp, e, n_entries, im, imp, rows, cols, k_a, k_b, bplus, &
        !$omp         bcross, sI0, sIm1, ph, point_jrphi, point_jrB, core, core_sp, &
        !$omp         acc, acc_rad, acc_sp, tables)
        do mag_a = M, 0, -1
            k_a = k_of_m(mag_a, L)
            do mag_b = 0, mag_a
                k_b = k_of_m(mag_b, L)
                call signed_entries(mag_a, mag_b, rows, cols, n_entries)
                acc = (0.0_dp, 0.0_dp)
                acc_sp = (0.0_dp, 0.0_dp)
                acc_rad = (0.0_dp, 0.0_dp)
                do j = 1, N
                    core = (0.0_dp, 0.0_dp)
                    core_sp = (0.0_dp, 0.0_dp)
                    do sp = 0, plasma%n_species - 1
                        if ((sp == 0 .and. turn_off_electrons) .or. &
                            (sp >= 1 .and. turn_off_ions)) cycle
                        call flr_arg_pair_sp(plasma, sp, k_a, k_b, j, bplus, bcross)
                        sI0 = 0.0_dp
                        sIm1 = 0.0_dp
                        if (fp_active .and. bplus - bcross <= flr_damping_cutoff) &
                            call scaled_bessel_pair(bplus, bcross, sI0, sIm1)
                        core_sp(1, sp) = core_rho_phi_sp_bessel(plasma, sp, bplus, bcross, &
                            j, sI0, sIm1) / (8.0_dp * pi**2)
                        core_sp(2, sp) = core_rho_B_sp_bessel(plasma, sp, bplus, bcross, &
                            j, sI0, sIm1) / (8.0_dp * pi**2)
                        core_sp(3, sp) = core_j_phi_sp_bessel(plasma, sp, bplus, bcross, &
                            j, sI0, sIm1) / (8.0_dp * pi**2)
                        core_sp(4, sp) = core_j_B_sp_bessel(plasma, sp, bplus, bcross, &
                            j, sI0, sIm1) / (8.0_dp * pi**2)
                        core = core + core_sp(:, sp)
                    end do
                    if (want_radial_current) call jrad_bessel_tables(plasma, k_b, k_a, j, tables)
                    do e = 1, n_entries
                        ph = phase(rows(e), j) * conjg(phase(cols(e), j))
                        acc(:, e) = acc(:, e) + ph * core
                        if (want_species) acc_sp(:, :, e) = acc_sp(:, :, e) + ph * core_sp
                        if (want_radial_current) then
                            call hatG_jrad_phi_br(plasma, k_of_m(rows(e), L), &
                                k_of_m(cols(e), L), j, point_jrphi, point_jrB, tables=tables)
                            acc_rad(1, e) = acc_rad(1, e) + point_jrphi
                            acc_rad(2, e) = acc_rad(2, e) + point_jrB
                        end if
                    end do
                end do
                do e = 1, n_entries
                    im = rows(e) + M + 1
                    imp = cols(e) + M + 1
                    Kphi(im, imp) = weight * acc(1, e)
                    KB(im, imp) = weight * acc(2, e)
                    Kjphi(im, imp) = weight * acc(3, e)
                    KjB(im, imp) = weight * acc(4, e)
                    if (want_radial_current) then
                        Kjrphi(im, imp) = weight * acc_rad(1, e)
                        KjrB(im, imp) = weight * acc_rad(2, e)
                    end if
                    ! Assumed-shape species dimension starts at 1 for species 0.
                    if (want_charge_species) then
                        Kphi_species(im, imp, :) = weight * acc_sp(1, :, e)
                        KB_species(im, imp, :) = weight * acc_sp(2, :, e)
                    end if
                    if (want_current_species) then
                        Kjphi_species(im, imp, :) = weight * acc_sp(3, :, e)
                        KjB_species(im, imp, :) = weight * acc_sp(4, :, e)
                    end if
                end do
            end do
        end do
        !$omp end parallel do
    end subroutine assemble_fp_symmetric

    !> Distinct signed (row, column) Fourier indices with magnitudes {a, b}.
    pure subroutine signed_entries(a, b, rows, cols, n_entries)
        integer, intent(in) :: a, b
        integer, intent(out) :: rows(8), cols(8), n_entries
        integer :: candidates(2, 8), c, e
        logical :: seen

        candidates = reshape([a, b, a, -b, -a, b, -a, -b, &
                              b, a, b, -a, -b, a, -b, -a], [2, 8])
        n_entries = 0
        do c = 1, 8
            seen = .false.
            do e = 1, n_entries
                if (rows(e) == candidates(1, c) .and. cols(e) == candidates(2, c)) seen = .true.
            end do
            if (.not. seen) then
                n_entries = n_entries + 1
                rows(n_entries) = candidates(1, c)
                cols(n_entries) = candidates(2, c)
            end if
        end do
    end subroutine signed_entries

    subroutine assemble_periodic_bparallel_matrices(plasma, L, M, &
                                 KBparallel, KjBparallel, KBparallel_species, KjBparallel_species, &
                                                    KjrBparallel)
        ! Assemble only the linear B_parallel charge/current columns. Keeping
        ! this separate leaves the established Phi/Br path and its zero-drive
        ! arithmetic unchanged. The configured kernel enforces the supported
        ! FP, m_phi=0 source-only derivative model before this routine is called.
        use collisionless_fourier_kernel_m, only: configured_hatG_Bparallel_all
        use constants_m, only: pi
        use grid_m, only: rg_grid
        use radial_current_fourier_kernel_m, only: hatG_jrad_bparallel

        type(plasma_t), intent(in) :: plasma
        real(dp), intent(in) :: L
        integer, intent(in) :: M
        complex(dp), allocatable, intent(out) :: KBparallel(:, :), KjBparallel(:, :)
        complex(dp), allocatable, intent(out), optional :: KBparallel_species(:, :, :)
        complex(dp), allocatable, intent(out), optional :: KjBparallel_species(:, :, :)
        complex(dp), allocatable, intent(out), optional :: KjrBparallel(:, :)
        complex(dp) :: point_rho, point_j, point_jrad
        complex(dp) :: acc_rho, acc_j, acc_jrad
        complex(dp) :: point_rho_species(0:plasma%n_species - 1)
        complex(dp) :: point_j_species(0:plasma%n_species - 1)
        complex(dp) :: acc_rho_species(0:plasma%n_species - 1)
        complex(dp) :: acc_j_species(0:plasma%n_species - 1)
        real(dp) :: k_m, k_mp, weight
        integer :: N, dim, m_row, m_col, im, imp, j
        logical :: want_species, want_jrad

        if (present(KBparallel_species) .neqv. present(KjBparallel_species)) then
            error stop 'assemble_periodic_bparallel_matrices requires both species matrices'
        end if
        N = rg_grid%npts_b
        dim = 2 * M + 1
        allocate (KBparallel(dim, dim), KjBparallel(dim, dim))
        want_species = present(KBparallel_species)
        want_jrad = present(KjrBparallel)
        if (want_jrad) allocate (KjrBparallel(dim, dim))
        if (want_species) then
            allocate (KBparallel_species(dim, dim, 0:plasma%n_species - 1))
            allocate (KjBparallel_species(dim, dim, 0:plasma%n_species - 1))
        end if
        weight = 2.0_dp * pi / real(N, dp)

        !$omp parallel do default(none) schedule(static) &
        !$omp shared(M, L, N, weight, plasma, KBparallel, KjBparallel, &
        !$omp        want_species, KBparallel_species, KjBparallel_species, &
        !$omp        want_jrad, KjrBparallel) &
        !$omp private(m_col, k_mp, imp, m_row, k_m, im, j, point_rho, &
        !$omp         point_j, point_jrad, acc_rho, acc_j, acc_jrad, point_rho_species, &
        !$omp         point_j_species, acc_rho_species, acc_j_species)
        do m_col = -M, M
            k_mp = k_of_m(m_col, L)
            imp = m_col + M + 1
            do m_row = -M, M
                k_m = k_of_m(m_row, L)
                im = m_row + M + 1
                acc_rho = (0.0_dp, 0.0_dp)
                acc_j = (0.0_dp, 0.0_dp)
                acc_jrad = (0.0_dp, 0.0_dp)
                if (want_species) then
                    acc_rho_species = (0.0_dp, 0.0_dp)
                    acc_j_species = (0.0_dp, 0.0_dp)
                end if
                do j = 1, N
                    if (want_species) then
                        call configured_hatG_Bparallel_all(plasma, k_m, k_mp, j, &
                                             point_rho, point_j, point_rho_species, point_j_species)
                        acc_rho_species = acc_rho_species + point_rho_species
                        acc_j_species = acc_j_species + point_j_species
                    else
                        call configured_hatG_Bparallel_all(plasma, k_m, k_mp, j, &
                                                           point_rho, point_j)
                    end if
                    acc_rho = acc_rho + point_rho
                    acc_j = acc_j + point_j
                    if (want_jrad) then
                        point_jrad = hatG_jrad_bparallel(plasma, k_m, k_mp, j)
                        acc_jrad = acc_jrad + point_jrad
                    end if
                end do
                KBparallel(im, imp) = weight * acc_rho
                KjBparallel(im, imp) = weight * acc_j
                if (want_jrad) KjrBparallel(im, imp) = weight * acc_jrad
                if (want_species) then
                    KBparallel_species(im, imp, :) = weight * acc_rho_species
                    KjBparallel_species(im, imp, :) = weight * acc_j_species
                end if
            end do
        end do
        !$omp end parallel do
    end subroutine assemble_periodic_bparallel_matrices

end module periodic_assembly_m
