module rhs_balance_m
    !
    ! RHS Balance Module
    !
    ! Implements 1D transport/balance equations for n, Vphi, Te, Ti
    ! [Markl2023 (22)–(25), Heyn2014 (63)–(66)].
    !
    ! Solves the discretized transport equations:
    !   dy/dt = A·y + q
    ! where:
    !   y - vector of plasma parameters at grid points (size 4N)
    !   A - sparse Jacobian matrix
    !   q - source vector
    ! using implicit time stepping.
    ! The implicit time stepping scheme uses:
    !   y' - Δt·A(y)·y' = y + Δt·q(y)
    !
    ! Public API:
    !   - rhs_balance:        Compute sparse Jacobian matrix A
    !   - rhs_balance_source: Compute source vector q
    !
    ! ============================================================================
    ! Frozen-Coefficient Linearization for Jacobian Construction
    ! ============================================================================
    !
    ! The Jacobian matrix A is constructed using frozen-coefficient linearization.
    ! This means certain quantities are evaluated once at the actual plasma state
    ! and held fixed ("frozen") while probing each column of the Jacobian.
    !
    ! FROZEN quantities (evaluated at actual state, held constant during probing):
    !   - Diffusion coefficients D^a, D^ql, D^nc, μ (viscosity)
    !   - Plasma parameters at boundaries: n, Te, Ti, Vphi (from params_b)
    !   - Radial electric field E0r (from Ercov)
    !   - QL particle fluxes Γ_ql_e, Γ_ql_i (from Gamma_ql_*_frozen)
    !
    ! RESPONDING quantities (recomputed for each probe perturbation δy):
    !   - Linearized gradients: δ(∂n/∂r), δ(∂T/∂r), δ(∂V/∂r) from ddr_params_lin
    !   - Linearized electric field δE0r from Ercov_lin
    !   - Linearized fluxes δΓ, δQ computed from above
    !
    ! The Jacobian is built by probing: A[i,j] = δ(dy_i/dt) / δy_j
    !
    ! We use frozen-coefficient linearization where n, T, D are frozen at the
    ! actual state and only gradients respond to the probe perturbation:
    !
    !   ┌─────────────────────────────────────────────────────────────────┐
    !   │  δA1 = (1/n)·δ(∂n/∂r) ± (e/T)·δE0r - (3/2T)·δ(∂T/∂r)            │
    !   │  δA2 = (1/T)·δ(∂T/∂r)                                           │
    !   │                                                                 │
    !   │  δΓ  = -n · D · (δA1, δA2)                                      │
    !   │  δQ  = -n · T · D · (δA1, δA2)                                  │
    !   │                                                                 │
    !   │  where n, T, D are FROZEN (from params_b, not params_b_lin)     │
    !   └─────────────────────────────────────────────────────────────────┘
    !
    ! This works because n largely cancels: Γ ∝ -n·D·(1/n)·∂n/∂r = -D·∂n/∂r
    !
    ! Why this is simpler than finite differences:
    !
    !   Standard finite difference:  J[i,j] = (f(y + δe_j) - f(y)) / δ
    !   → requires computing f twice and subtracting
    !
    !   Our approach: Since everything is linear in gradients with frozen
    !   coefficients, we directly apply the linearized operator to δy:
    !
    !       δy → δ(∂y/∂r) → δA → δΓ, δQ → δ(dy/dt)
    !
    !   Then:  J[i,j] = δ(dy/dt)[i] / δy[j]
    !
    ! References: [Markl2023 (5), (29)], [Heyn2014 (22), (35)]
    !
    use QLBalance_kinds, only: dp
    use baseparam_mod, only: e_charge, p_mass

    implicit none

    private

    !---------------------------------------------------------------------------
    ! Public API
    !---------------------------------------------------------------------------
    public :: initialize_rhs
    public :: rhs_balance
    public :: rhs_balance_source

    !---------------------------------------------------------------------------
    ! Public for testing (pure helper functions)
    !---------------------------------------------------------------------------
    public :: compute_thermodynamic_forces
    public :: compute_particle_fluxes
    public :: compute_total_heat_fluxes
    public :: compute_diffusive_parts
    public :: compute_convective_parts
    public :: compute_source_terms_at_point

    public :: thermodynamic_forces_t
    public :: species_fluxes_t
    public :: transport_fluxes_t

    !---------------------------------------------------------------------------
    ! Derived types
    !---------------------------------------------------------------------------

    !> Thermodynamic forces for a single species
    type :: species_forces_t
        real(dp) :: A1 = 0.0_dp
        real(dp) :: A1_noE = 0.0_dp  ! E0r set to zero
        real(dp) :: A2 = 0.0_dp
    end type species_forces_t

    !> Thermodynamic forces at a single radial point
    type :: thermodynamic_forces_t
        type(species_forces_t) :: e  !< Electron forces
        type(species_forces_t) :: i  !< Ion forces
    end type thermodynamic_forces_t

    !> Particle and heat fluxes for a single species
    type :: species_fluxes_t
        real(dp) :: Gamma_a = 0.0_dp  !< Anomalous particle flux
        real(dp) :: Gamma_ql = 0.0_dp  !< Quasi-linear particle flux
        real(dp) :: Gamma_tot = 0.0_dp  !< Total particle flux
        real(dp) :: Q = 0.0_dp  !< Heat flux
    end type species_fluxes_t

    !> Transport fluxes for electrons and ions
    type :: transport_fluxes_t
        type(species_fluxes_t) :: e  !< Electron fluxes
        type(species_fluxes_t) :: i  !< Ion fluxes
    end type transport_fluxes_t

contains

    !===========================================================================
    ! Main public subroutines
    !===========================================================================

    subroutine initialize_rhs(y, dy)
        !
        ! Initialize the sparse matrix infrastructure.
        !
        ! This subroutine:
        !   1. Calls rhs_balance with isw_rhs=0 to count non-zero elements
        !   2. Allocates sparse matrix arrays (irow, icol, amat, rhsvec)
        !
        ! This is needed because the implicit time stepping requires knowing
        ! the sparsity pattern of the Jacobian matrix before filling it.
        !
        use grid_mod, only: neqset
        use matrix_mod

        implicit none

        real(dp), dimension(neqset), intent(in) :: y
        real(dp), dimension(neqset), intent(out) :: dy
        real(dp) :: dummy

        dummy = 0.0_dp
        isw_rhs = 0

        call rhs_balance(dummy, y, dy)

        isw_rhs = 1
        if (allocated(amat)) deallocate (irow, icol, amat, rhsvec)
        allocate (irow(nz), icol(nz), amat(nz), rhsvec(nsize))

    end subroutine initialize_rhs

    subroutine rhs_balance(x, y, dy)
        !
        ! Compute the sparse Jacobian matrix A in: dy/dt = A·y + q
        !
        ! This subroutine probes each column of the Jacobian by setting
        ! y_lin(iprobe) = 1 one at a time, computing the resulting dy,
        ! and storing non-zero elements in sparse format (irow, icol, amat).
        !
        ! When isw_rhs=0: Only counts non-zeros (for allocation)
        ! When isw_rhs=1: Fills the sparse matrix and calls rhs_balance_source for q
        !
        use grid_mod, only: nbaleqs, neqset, iboutype, npoic, npoib, Sc, Sb, deriv_coef, ipbeg, &
                            ipend, rb, reint_coef, fluxes_dif_lin, fluxes_con_lin, rc, dae11, &
                            dae12, dae22, dai11, dai12, dai22, dni22, visca, gpp_av, dqle11, &
                            dqle12, dqle21, dqle22, dqli11, dqli12, dqli21, dqli22, T_EM_phi_e, &
                            T_EM_phi_i, torque_ntv, sqrt_g_times_B_theta_over_c, Ercov, polforce, &
                            qlheat_e, qlheat_i, Ercov_lin, fluxes_con_nl, Gamma_ql_e_frozen, &
                            Gamma_ql_i_frozen
        use plasma_parameters, only: params, ddr_params_lin, params_lin, params_b_lin, &
                                     params_b, dot_params, ddr_params_nl
        use baseparam_mod, only: Z_i, am
        use matrix_mod, only: isw_rhs, nz, nsize, irow, icol, amat, rhsvec
        use wave_code_data, only: q, Vth

        implicit none

        real(dp), intent(in) :: x
        real(dp), dimension(neqset), intent(in) :: y
        real(dp), dimension(neqset), intent(out) :: dy

        integer :: ipoi, ieq, i, npoi, ibeg, iend, nshift, ibegb, iendb, ibegtot, iendtot, k, iprobe
        real(dp), dimension(neqset) :: y_lin

        ! Local variables for flux computation
        type(thermodynamic_forces_t) :: forces_lin, forces_nl
        real(dp) :: Gamma_tot_e_lin, Gamma_tot_i_lin, Gamma_ql_e_lin, Gamma_ql_i_lin, Qe_lin, Qi_lin
        real(dp) :: Gamma_tot_e_nl, Gamma_tot_i_nl, Gamma_ql_e_nl, Gamma_ql_i_nl, Qe_nl, Qi_nl
        real(dp), dimension(4) :: flux_dif_lin_loc, flux_con_lin_loc
        real(dp), dimension(4) :: flux_dif_nl_loc, flux_con_nl_loc
        real(dp), dimension(4) :: dot_params_loc

        if (iboutype .eq. 1) then
            npoi = npoic - 1
        else
            npoi = npoic
        end if

        ! Initialize arrays for linearized perturbation
        y_lin = 0.0_dp
        params_lin = 0.0_dp
        ddr_params_lin = 0.0_dp
        Ercov_lin = 0.0_dp

        ! Prepare frozen state: evaluate all quantities at actual plasma state
        ! that should be held fixed during Jacobian probing

        ! Copy y to params
        do ipoi = 1, npoi
            do ieq = 1, nbaleqs
                i = nbaleqs * (ipoi - 1) + ieq
                params(ieq, ipoi) = y(i)
            end do
        end do

        ! Interpolate to boundary points
        do ipoi = 1, npoib
            do ieq = 1, nbaleqs
                associate (param => params(ieq, ipbeg(ipoi):ipend(ipoi)))
                    ddr_params_nl(ieq, ipoi) = sum(param * deriv_coef(:, ipoi))
                    params_b(ieq, ipoi) = sum(param * reint_coef(:, ipoi))
                end associate
            end do
        end do

        ! Compute E0r at actual state
        call compute_radial_electric_field(npoib, rb, params_b, ddr_params_nl, &
                                           sqrt_g_times_B_theta_over_c, Vth, q, Z_i, Ercov)

        ! Compute diffusion coefficients at actual state
        call calc_equil_diffusion_coeffs

        ! Compute forces and convective fluxes at ALL boundary points,
        ! storing frozen QL fluxes for complete linearization
        do ipoi = 1, npoib
            call compute_fluxes_at_boundary(ipoi, ddr_params_nl, params_b, Ercov(ipoi), dae11, &
                                            dae12, dae22, dai11, dai12, dai22, dni22, dqle11, &
                                            dqle12, dqle21, dqle22, dqli11, dqli12, dqli21, &
                                            dqli22, visca, gpp_av, Sb, Z_i, forces_nl, &
                                            Gamma_tot_e_nl, Gamma_tot_i_nl, Gamma_ql_e_nl, &
                                            Gamma_ql_i_nl, Qe_nl, Qi_nl, flux_dif_nl_loc, &
                                            flux_con_nl_loc)

            fluxes_con_nl(:, ipoi) = flux_con_nl_loc
            Gamma_ql_e_frozen(ipoi) = Gamma_ql_e_nl
            Gamma_ql_i_frozen(ipoi) = Gamma_ql_i_nl
        end do

        ! Apply boundary condition: zero flux at inner boundary
        fluxes_con_nl(:, 1) = 0.0_dp

        ! Jacobian probing loop
        nshift = 4
        k = 0
        dy = 0.0_dp

        do iprobe = 1, neqset
            ! Determine affected range
            ibeg = iprobe / nbaleqs - nshift
            iend = iprobe / nbaleqs + nshift
            ibegb = ibeg - 1
            iendb = iend + 1
            ibegtot = iprobe - nshift * nbaleqs
            iendtot = iprobe + nshift * nbaleqs
            ibeg = max(1, ibeg)
            iend = min(npoi, iend)
            ibegb = max(1, ibegb)
            iendb = min(npoib, iendb)
            ibegtot = max(1, ibegtot)
            iendtot = min(neqset, iendtot)

            ! Set probe with physically-motivated relative perturbation of actual values
            ! This prevents extreme Jacobian entries from resonance responses
            if (abs(y(iprobe)) > 1.0e-30_dp) then
                y_lin(iprobe) = 1.0e-6_dp * y(iprobe)
            else
                y_lin(iprobe) = 1.0e-6_dp  ! Fallback for zero values
            end if

            ! Copy y_lin to params_lin
            do ipoi = ibeg, iend
                do ieq = 1, nbaleqs
                    i = nbaleqs * (ipoi - 1) + ieq
                    params_lin(ieq, ipoi) = y_lin(i)
                end do
            end do

            ! Interpolate linear perturbation to boundaries
            do ipoi = ibegb, iendb
                do ieq = 1, nbaleqs
                    associate (param => params_lin(ieq, ipbeg(ipoi):ipend(ipoi)))
                        ddr_params_lin(ieq, ipoi) = sum(param * deriv_coef(:, ipoi))
                        params_b_lin(ieq, ipoi) = sum(param * reint_coef(:, ipoi))
                    end associate
                end do
            end do

            ! Compute linearized radial electric field δE0r
            ! From E0r = (Ti/n)·∂n/∂r/ei + (1/ei)·∂Ti/∂r + √g·Bθ/c·(Vphi - q·Vpol/r)
            ! Linearized: δE0r = (Ti/n)·δ(∂n/∂r)/ei + (1/ei)·δ(∂Ti/∂r) + √g·Bθ/c·δVphi
            !
            ! FROZEN (from params_b): Ti, n (in prefactors)
            ! RESPONDING (from ddr_params_lin, params_b_lin): δ(∂n/∂r), δ(∂Ti/∂r), δVphi
            ! Note: Vpol term omitted as poloidal velocity is not evolved (constant)
            Ercov_lin(ibegb:iendb) = sqrt_g_times_B_theta_over_c(ibegb:iendb) * &
                params_b_lin(2, ibegb:iendb) + (params_b(4, ibegb:iendb) * &
                ddr_params_lin(1, ibegb:iendb) / params_b(1, ibegb:iendb) + &
                ddr_params_lin(4, ibegb:iendb)) / (Z_i * e_charge)

            ! Compute linearized fluxes δΓ, δQ at affected boundary points
            ! Using frozen-coefficient linearization:
            !   δA1 = (1/n)·δ(∂n/∂r) ± (e/T)·δE0r - (3/2T)·δ(∂T/∂r)     [n,T frozen]
            !   δA2 = (1/T)·δ(∂T/∂r)                                    [T frozen]
            !   δΓ  = -n·D·(δA1, δA2)                                   [n,D frozen]
            !   δQ  = -n·T·D·(δA1, δA2)                                 [n,T,D frozen]
            !
            ! Arguments to compute_fluxes_at_boundary:
            !   ddr_params_lin: linearized gradients δ(∂/∂r)  [RESPONDING]
            !   params_b:       frozen values n, Te, Ti       [FROZEN - intentional!]
            !   Ercov_lin:      linearized E-field δE0r       [RESPONDING]
            !   D coefficients: frozen diffusion coefficients [FROZEN]
            do ipoi = ibegb, iendb
                call compute_fluxes_at_boundary(ipoi, ddr_params_lin, params_b, Ercov_lin(ipoi), &
                                                dae11, dae12, dae22, dai11, dai12, dai22, dni22, &
                                                dqle11, dqle12, dqle21, dqle22, dqli11, dqli12, &
                                                dqli21, dqli22, visca, gpp_av, Sb, Z_i, &
                                                forces_lin, Gamma_tot_e_lin, Gamma_tot_i_lin, &
                                                Gamma_ql_e_lin, Gamma_ql_i_lin, Qe_lin, Qi_lin, &
                                                flux_dif_lin_loc, flux_con_lin_loc)

                fluxes_dif_lin(:, ipoi) = flux_dif_lin_loc
                fluxes_con_lin(:, ipoi) = flux_con_lin_loc

                ! Compute source terms using linearized QL fluxes and frozen values
                ! Complete linearization: δ(E0r·Γ_ql) = E0r_frozen·δΓ_ql + δE0r·Γ_ql_frozen
                call compute_rmp_induced_sources(Gamma_ql_e_lin, Gamma_ql_i_lin, &
                                                 Gamma_ql_e_frozen(ipoi), Gamma_ql_i_frozen(ipoi), &
                                                 Ercov(ipoi), Ercov_lin(ipoi), &
                                                 sqrt_g_times_B_theta_over_c(ipoi), Z_i, am, &
                                                 torque_ntv(ipoi), polforce(ipoi), qlheat_e(ipoi), &
                                                 qlheat_i(ipoi), T_EM_phi_e(ipoi), T_EM_phi_i(ipoi))
            end do

            ! Apply boundary conditions for linearized fluxes
            fluxes_dif_lin(:, 1) = 0.0_dp
            fluxes_con_lin(:, 1) = 0.0_dp

            ! Compute time derivatives
            do ipoi = ibeg, iend
                call compute_dot_params_at_point(ipoi, npoi, nbaleqs, fluxes_dif_lin, &
                                                 fluxes_con_lin, fluxes_con_nl, params, &
                                                 params_lin, params_b_lin, Sc, rb, rc, gpp_av, &
                                                 polforce, qlheat_e, qlheat_i, Z_i, dot_params_loc)
                dot_params(:, ipoi) = dot_params_loc
            end do

            ! Copy to dy
            do ipoi = ibeg, iend
                do ieq = 1, nbaleqs
                    i = nbaleqs * (ipoi - 1) + ieq
                    dy(i) = dot_params(ieq, ipoi)
                end do
            end do

            ! Store non-zero elements
            ! Divide by perturbation size to get correct Jacobian entry:
            ! J[i,j] = dy[i] / delta_y[j]
            do i = ibegtot, iendtot
                if (abs(dy(i)) > tiny(1.0_dp)) then
                    k = k + 1
                    if (isw_rhs .eq. 1) then
                        irow(k) = i
                        icol(k) = iprobe
                        amat(k) = dy(i) / y_lin(iprobe)
                    end if
                end if
            end do

            ! Clean up for next probe
            y_lin(iprobe) = 0.0_dp
            do ipoi = ibeg, iend
                do ieq = 1, nbaleqs
                    i = nbaleqs * (ipoi - 1) + ieq
                    params_lin(ieq, ipoi) = 0.0_dp
                end do
            end do
            ddr_params_lin(:, ibegb:iendb) = 0.0_dp
            Ercov_lin(ibegb:iendb) = 0.0_dp
            fluxes_dif_lin(:, ibegb:iendb) = 0.0_dp
            fluxes_con_lin(:, ibegb:iendb) = 0.0_dp
            polforce(ibegb:iendb) = 0.0_dp
            qlheat_e(ibegb:iendb) = 0.0_dp
            qlheat_i(ibegb:iendb) = 0.0_dp
            dot_params(:, ibeg:iend) = 0.0_dp
            dy = 0.0_dp
        end do

        if (isw_rhs == 0) then
            nz = k
            nsize = neqset
        else
            call rhs_balance_source(x, y, dy)
            rhsvec = dy
        end if

    end subroutine rhs_balance

    subroutine rhs_balance_source(x, y, dy)
        !
        ! Compute the source vector q in: dy/dt = A·y + q
        !
        ! This computes only the source terms (not flux divergence terms, which
        ! are in matrix A). The sources are:
        !   - polforce: momentum source from RMP-induced radial current
        !   - qlheat_e/i: heat sources from RMP-induced particle transport
        !   - dery_equisource: equilibrium background source
        !
        ! NOTE: This subroutine is only called from rhs_balance, which has already
        ! set up params, params_b, ddr_params_nl, Ercov, and diffusion coefficients.
        !
        use grid_mod, only: nbaleqs, neqset, iboutype, npoic, npoib, gpp_av, dery_equisource, &
                            T_EM_phi_e, T_EM_phi_i, torque_ntv, sqrt_g_times_B_theta_over_c, &
                            Ercov, polforce, qlheat_e, qlheat_i, Gamma_ql_e_frozen, &
                            Gamma_ql_i_frozen
        use plasma_parameters, only: params, dot_params
        use baseparam_mod, only: Z_i, am

        implicit none

        real(dp), intent(in) :: x
        real(dp), dimension(neqset), intent(in) :: y
        real(dp), dimension(neqset), intent(out) :: dy

        integer :: ipoi, ieq, i, npoi
        real(dp), dimension(4) :: dot_params_loc

        if (iboutype .eq. 1) then
            npoi = npoic - 1
        else
            npoi = npoic
        end if

        ! Compute RMP-induced sources using frozen QL fluxes (already computed in rhs_balance)
        do ipoi = 1, npoib
            ! For actual state computation: no frozen contribution (zeros)
            ! Result: qlheat = E0r * Gamma_ql (full product, not linearized)
            ! => set linearized quantities to 0
            call compute_rmp_induced_sources(Gamma_ql_e_frozen(ipoi), Gamma_ql_i_frozen(ipoi), &
                                             0.0_dp, 0.0_dp, Ercov(ipoi), 0.0_dp, &
                                             sqrt_g_times_B_theta_over_c(ipoi), Z_i, am, &
                                             torque_ntv(ipoi), polforce(ipoi), qlheat_e(ipoi), &
                                             qlheat_i(ipoi), T_EM_phi_e(ipoi), T_EM_phi_i(ipoi))
        end do

        ! Compute source terms only (flux divergence is in matrix A)
        do ipoi = 1, npoi
            call compute_source_terms_at_point(ipoi, params, gpp_av, polforce, qlheat_e, qlheat_i, &
                                               Z_i, dot_params_loc)
            dot_params(:, ipoi) = dot_params_loc
        end do

        ! Copy to dy
        do ipoi = 1, npoi
            do ieq = 1, nbaleqs
                i = nbaleqs * (ipoi - 1) + ieq
                dy(i) = dot_params(ieq, ipoi)
            end do
        end do

        ! Add equilibrium source (background source term)
        dy = dy + dery_equisource

    end subroutine rhs_balance_source

    !===========================================================================
    ! Core computation routines (shared by rhs_balance and rhs_balance_source)
    !===========================================================================

    pure subroutine compute_radial_electric_field(npoib, rb, params_b, ddr_params, &
                                                  sqrt_g_Bth_over_c, V_pol_arr, q_arr, Z, E0r)
        !
        ! Compute equilibrium radial electric field at all boundary points.
        !   E0r = 1/(ei ni) ∂(ni Ti)/∂r + sqrt(g) B^θ_0 / c (V^φ - q V^θ)
        !     [Markl2023 (26), Heyn2014 (71)]
        !   V^θ = V_pol / r
        !     [Markl2023 directly above (40)]
        !

        implicit none

        integer, intent(in) :: npoib
        real(dp), intent(in) :: rb(:)
        real(dp), intent(in) :: params_b(:, :)  ! (4, npoib)
        real(dp), intent(in) :: ddr_params(:, :)  ! (4, npoib)
        real(dp), intent(in) :: sqrt_g_Bth_over_c(:)
        real(dp), intent(in) :: V_pol_arr(:), q_arr(:)
        real(dp), intent(in) :: Z
        real(dp), intent(out) :: E0r(:)

        associate ( &
            Vphi => params_b(2, 1:npoib), &
            Ti => params_b(4, 1:npoib), &
            dn_dr => ddr_params(1, 1:npoib), &
            n => params_b(1, 1:npoib), &
            dTi_dr => ddr_params(4, 1:npoib), &
            q => q_arr(1:npoib), &
            V_pol => V_pol_arr(1:npoib), &
            r => rb(1:npoib) &
        )
            E0r(1:npoib) = &
                (Ti / n * dn_dr + dTi_dr) / (Z * e_charge) + &
                sqrt_g_Bth_over_c(1:npoib) * (Vphi - q * V_pol / r)
        end associate

    end subroutine compute_radial_electric_field

    pure subroutine compute_fluxes_at_boundary(ipoi, ddr_params, params_b, E0r, Dae11, Dae12, &
                                               Dae22, Dai11, Dai12, Dai22, Dni22, Dqle11, Dqle12, &
                                               Dqle21, Dqle22, Dqli11, Dqli12, Dqli21, Dqli22, &
                                               visca, g_phi_phi, S, Z, forces, Gamma_tot_e, &
                                               Gamma_tot_i, Gamma_ql_e, Gamma_ql_i, Qe, Qi, &
                                               flux_diffusion, flux_convection)
        !
        ! Compute all fluxes at a single boundary point.
        !
        ! This subroutine is used in two contexts:
        !
        ! 1. ACTUAL STATE (rhs_balance frozen state preparation):
        !    - ddr_params = actual gradients (ddr_params_nl)
        !    - params_b   = actual values
        !    - E0r        = actual electric field (Ercov)
        !    Result: fluxes at actual plasma state, stored in Gamma_ql_*_frozen
        !
        ! 2. LINEARIZED STATE (rhs_balance Jacobian probing):
        !    - ddr_params = linearized gradients δ(∂/∂r) from probe (ddr_params_lin)
        !    - params_b   = ACTUAL values n, Te, Ti (NOT params_b_lin!)
        !    - E0r        = linearized electric field δE0r (Ercov_lin)
        !    Result: linearized flux response δΓ, δQ
        !
        ! In the linearized case, params_b contains frozen values (actual state)
        ! because the flux formulas have n, T appearing in two roles:
        !   Γ = -n · D · [(1/n)·∂n/∂r - ...]
        !
        ! For frozen-coefficient linearization, both the prefactor n and the 1/n
        ! in the force are frozen at actual state. Only the gradient ∂n/∂r responds.
        ! This avoids numerical instability from dividing by tiny perturbation values.
        !

        implicit none

        integer, intent(in) :: ipoi
        real(dp), intent(in) :: ddr_params(:, :)  !< Gradients (actual or linearized)
        real(dp), intent(in) :: params_b(:, :)  !< Values n,Vphi,Te,Ti (ALWAYS actual state)
        real(dp), intent(in) :: E0r
        real(dp), intent(in) :: Dae11(:), Dae12(:), Dae22(:)
        real(dp), intent(in) :: Dai11(:), Dai12(:), Dai22(:), Dni22(:)
        real(dp), intent(in) :: Dqle11(:), Dqle12(:), Dqle21(:), Dqle22(:)
        real(dp), intent(in) :: Dqli11(:), Dqli12(:), Dqli21(:), Dqli22(:)
        real(dp), intent(in) :: visca(:), g_phi_phi(:), S(:)
        real(dp), intent(in) :: Z

        type(thermodynamic_forces_t), intent(out) :: forces
        real(dp), intent(out) :: Gamma_tot_e, Gamma_tot_i, Gamma_ql_e, Gamma_ql_i
        real(dp), intent(out) :: Qe, Qi
        real(dp), dimension(4), intent(out) :: flux_diffusion, flux_convection

        type(transport_fluxes_t) :: fluxes
        real(dp) :: n, Te, Ti, dn_dr, dVphi_dr, dTe_dr, dTi_dr

        ! Extract local values
        n = params_b(1, ipoi)
        Te = params_b(3, ipoi)
        Ti = params_b(4, ipoi)
        dn_dr = ddr_params(1, ipoi)
        dVphi_dr = ddr_params(2, ipoi)
        dTe_dr = ddr_params(3, ipoi)
        dTi_dr = ddr_params(4, ipoi)

        call compute_thermodynamic_forces(dn_dr, dTe_dr, dTi_dr, n, Te, Ti, E0r, Z, forces)

        call compute_particle_fluxes(forces, n, Z, Dae11(ipoi), Dae12(ipoi), Dqle11(ipoi), &
                                     Dqle12(ipoi), Dai11(ipoi), Dai12(ipoi), Dqli11(ipoi), &
                                     Dqli12(ipoi), fluxes)

        ! Extract individual flux values for output and further use
        Gamma_tot_e = fluxes%e%Gamma_tot
        Gamma_tot_i = fluxes%i%Gamma_tot
        Gamma_ql_e = fluxes%e%Gamma_ql
        Gamma_ql_i = fluxes%i%Gamma_ql

        call compute_total_heat_fluxes(forces, n, Te, Ti, Z, Dae12(ipoi), Dqle21(ipoi), &
                                       Dae22(ipoi), Dqle22(ipoi), Dai12(ipoi), Dqli21(ipoi), &
                                       Dai22(ipoi), Dni22(ipoi), Dqli22(ipoi), Qe, Qi)

        ! extract contributions from total heat fluxes
        call compute_diffusive_parts(dn_dr, dVphi_dr, dTe_dr, dTi_dr, n, Te, Ti, Z, S(ipoi), &
                                     Dae11(ipoi), Dqle11(ipoi), Dae22(ipoi), Dqle22(ipoi), &
                                     Dai22(ipoi), Dni22(ipoi), Dqli22(ipoi), Dqli21(ipoi), &
                                     visca(ipoi), g_phi_phi(ipoi), flux_diffusion)

        call compute_convective_parts(Gamma_tot_e, Qe, Qi, n, Te, Ti, S(ipoi), flux_diffusion, &
                                      flux_convection)

    end subroutine compute_fluxes_at_boundary

    pure subroutine compute_rmp_induced_sources(Gamma_ql_e, Gamma_ql_i, Gamma_ql_e_froz, &
                                                Gamma_ql_i_froz, E0r, E0r_lin, sqrt_g_Bth_over_c, &
                                                Z, am, torque_ntv, polforce, qlheat_e, &
                                                qlheat_i, torque_e, torque_i)
        !
        ! Compute RMP-induced source terms of the four balance equations.
        !
        ! For actual state computation (rhs_balance_source):
        !   Pass Gamma_ql_*_froz = 0 and E0r_lin = 0
        !   Result: qlheat = E0r * Gamma_ql
        !
        ! For Jacobian probing (rhs_balance):
        !   Pass all quantities
        !   Result: δ(qlheat) = E0r_frozen * δΓ_ql + δE0r * Γ_ql_frozen
        !
        ! The RHS for the four equations (with external sources set to zero) are as follows:
        !   1.) = 0
        !         [Markl2023 (22), Heyn2014 (63)]
        !   2.) = Σ_{e,i} -sqrt(g) B^θ_0 / c e Γ^EM = T^EM_φ,e + T^EM_φ,i
        !         [Markl2023 (23), Heyn2014 (64), (67)]
        !   3.) = -E0r ee Γ^EM_e = +E0r e Γ^EM_e
        !         [Markl2023 (24), Heyn2014 (65)]
        !   4.) = -E0r ei Γ^EM_i = -E0r Z e Γ^EM_i
        !         [Markl2023 (24), Heyn2014 (65)]
        !

        implicit none

        real(dp), intent(in) :: Gamma_ql_e, Gamma_ql_i
        real(dp), intent(in) :: Gamma_ql_e_froz, Gamma_ql_i_froz
        real(dp), intent(in) :: E0r, E0r_lin
        real(dp), intent(in) :: sqrt_g_Bth_over_c
        real(dp), intent(in) :: Z, am
        real(dp), intent(in) :: torque_ntv
        real(dp), intent(out) :: polforce, qlheat_e, qlheat_i
        real(dp), intent(out) :: torque_e, torque_i

        ! Eq. 1: Particle density
        ! no internal sources

        ! Eq. 2: Toroidal rotation frequency
        ! Divide by the factor mi to get velocity instead of momentum.
        ! This code evolves velocity directly.
        torque_e = +sqrt_g_Bth_over_c * e_charge * Gamma_ql_e
        torque_i = -sqrt_g_Bth_over_c * Z * e_charge * Gamma_ql_i
        ! Net radial current drive entering the momentum evolution as a source term.
        ! torque_{e,i} has units [dyn/cm²] = [g·cm⁻¹·s⁻²] (force per area, same as pressure).
        ! polforce = torque / mass has units [cm⁻¹·s⁻²].
        ! In compute_dot_params_at_point, polforce is multiplied by Z/(n·g_φφ) to give
        ! angular acceleration [s⁻²] for the toroidal rotation frequency evolution.
        polforce = (torque_e + torque_i + torque_ntv) / (am * p_mass)

        ! Eq. 3: Electron temperature
        ! Complete linearization: δ(E0r·Γ_ql) = E0r·δΓ_ql + δE0r·Γ_ql_frozen
        qlheat_e = +e_charge * (E0r * Gamma_ql_e + E0r_lin * Gamma_ql_e_froz)

        ! Eq. 4: Ion temperature
        ! Complete linearization: δ(E0r·Γ_ql) = E0r·δΓ_ql + δE0r·Γ_ql_frozen
        qlheat_i = -Z * e_charge * (E0r * Gamma_ql_i + E0r_lin * Gamma_ql_i_froz)
    end subroutine compute_rmp_induced_sources

    pure subroutine compute_dot_params_at_point(ipoi, npoi, nbaleqs, fluxes_dif, fluxes_con, &
                                                fluxes_con_nl, params, params_lin, params_b_lin, &
                                                Sc, rb, rc, gpp_av, polforce, qlheat_e, qlheat_i, &
                                                Z_i, dot_params_out)
        !
        ! Compute time derivatives at a single grid point.
        !
        ! This combines:
        !  - Flux divergence
        !  - Upstream convection
        !  - Sources
        !  - Quantity conversions (momentum -> omega, nT -> T)
        !
        ! The splitting into diffusive and convective flux is to improve numerical stability.
        ! The upstream convection has two branches s.t. the resulting directional derivative
        ! always points in the same direction.
        !   if convection_velocity > 0: (p(i+1) - p(i)) / (r(i+1) - r(i))
        !   else:                       (p(i-1) - p(i)) / (r(i-1) - r(i))
        !

        implicit none

        integer, intent(in) :: ipoi, npoi, nbaleqs
        real(dp), intent(in) :: fluxes_dif(:, :), fluxes_con(:, :), fluxes_con_nl(:, :)
        real(dp), intent(in) :: params(:, :), params_lin(:, :), params_b_lin(:, :)
        real(dp), intent(in) :: Sc(:), rb(:), rc(:), gpp_av(:)
        real(dp), intent(in) :: polforce(:), qlheat_e(:), qlheat_i(:)
        real(dp), intent(in) :: Z_i
        real(dp), dimension(4), intent(out) :: dot_params_out

        integer :: ieq
        real(dp) :: convec_velocity, dr

        dr = rb(ipoi + 1) - rb(ipoi)

        ! Flux divergence for all 4 equations:
        do ieq = 1, nbaleqs
            dot_params_out(ieq) = -(fluxes_dif(ieq, ipoi + 1) - fluxes_dif(ieq, ipoi)) / &
                (Sc(ipoi) * dr) - (fluxes_con(ieq, ipoi + 1) - fluxes_con(ieq, ipoi)) / &
                (Sc(ipoi) * dr) * params(ieq, ipoi)

            ! Upstream convection:
            convec_velocity = 0.5_dp * (fluxes_con_nl(ieq, ipoi + 1) + fluxes_con_nl(ieq, ipoi)) / &
                              Sc(ipoi)
            if (convec_velocity .gt. 0.0_dp) then
                dot_params_out(ieq) = dot_params_out(ieq) - convec_velocity * &
                                      (params_lin(ieq, ipoi + 1) - params_lin(ieq, ipoi)) / &
                                      (rc(ipoi + 1) - rc(ipoi))
            else
                if (ipoi .gt. 1) then
                    dot_params_out(ieq) = dot_params_out(ieq) - convec_velocity * &
                                          (params_lin(ieq, ipoi - 1) - params_lin(ieq, ipoi)) / &
                                          (rc(ipoi - 1) - rc(ipoi))
                else  ! no center point before first grid point
                    dot_params_out(ieq) = dot_params_out(ieq) - convec_velocity * &
                        (params_b_lin(ieq, 1) - params_lin(ieq, 1)) / (rb(1) - rc(1))
                end if
            end if
        end do

        ! Add internal sources:
        ! Momentum:
        dot_params_out(2) = dot_params_out(2) + 0.5_dp * (polforce(ipoi) + polforce(ipoi + 1))

        ! Heat into electrons:
        dot_params_out(3) = dot_params_out(3) + 0.5_dp * (qlheat_e(ipoi) + qlheat_e(ipoi + 1))

        ! Heat into ions:
        dot_params_out(4) = dot_params_out(4) + 0.5_dp * (qlheat_i(ipoi) + qlheat_i(ipoi + 1))

        ! Convert momentum time derivative to rotation frequency:
        dot_params_out(2) = dot_params_out(2) * Z_i / params(1, ipoi) * 2.0_dp / (gpp_av(ipoi + 1) &
                                                                                  + gpp_av(ipoi))

        ! Convert d(nT)/dt to dT/dt:
        dot_params_out(3) = (-params(3, ipoi) * dot_params_out(1) + dot_params_out(3) / 1.5_dp) / &
                            params(1, ipoi)
        dot_params_out(4) = (-params(4, ipoi) * dot_params_out(1) + dot_params_out(4) / 1.5_dp) / &
                            params(1, ipoi)

    end subroutine compute_dot_params_at_point

    pure subroutine compute_source_terms_at_point(ipoi, params, gpp_av, polforce, qlheat_e, &
                                                  qlheat_i, Z_i, dot_params_out)
        !
        ! Compute source contributions only (no flux divergence).
        !
        ! This subroutine is used by rhs_balance_source where flux terms belong
        ! in the matrix A, not in the source vector q. It computes:
        !   - polforce: momentum source from radial current (divided by mi)
        !   - qlheat_e/i: heat sources from QL fluxes
        !
        ! And applies the necessary conversions:
        !   - Momentum → rotation frequency (divide by n, gpp)
        !   - d(nT)/dt → dT/dt (using density evolution)
        !
        ! Note: Unlike compute_dot_params_at_point, this does NOT include:
        !   - Flux divergence terms
        !   - Upstream convection terms
        ! Those belong in matrix A for the source computation.
        !

        implicit none

        integer, intent(in) :: ipoi
        real(dp), intent(in) :: params(:, :)
        real(dp), intent(in) :: gpp_av(:)
        real(dp), intent(in) :: polforce(:), qlheat_e(:), qlheat_i(:)
        real(dp), intent(in) :: Z_i
        real(dp), dimension(4), intent(out) :: dot_params_out

        ! Initialize - no flux contributions here
        dot_params_out = 0.0_dp

        ! Eq 1: Particle density - no internal sources
        dot_params_out(1) = 0.0_dp

        ! Eq 2: Momentum source (averaged to cell center)
        dot_params_out(2) = 0.5_dp * (polforce(ipoi) + polforce(ipoi + 1))

        ! Eq 3: Electron heat source (averaged to cell center)
        dot_params_out(3) = 0.5_dp * (qlheat_e(ipoi) + qlheat_e(ipoi + 1))

        ! Eq 4: Ion heat source (averaged to cell center)
        dot_params_out(4) = 0.5_dp * (qlheat_i(ipoi) + qlheat_i(ipoi + 1))

        ! Convert momentum time derivative to rotation frequency:
        ! d(Omega)/dt = d(p/mi)/dt * Z / (n * gpp)
        dot_params_out(2) = dot_params_out(2) * Z_i / params(1, ipoi) * &
                            2.0_dp / (gpp_av(ipoi + 1) + gpp_av(ipoi))

        ! Convert d(nT)/dt to dT/dt:
        ! Since dot_params_out(1) = 0 for source-only computation,
        ! d(nT)/dt = n*dT/dt + T*dn/dt simplifies to:
        ! dT/dt = (1/n) * [source / 1.5]
        dot_params_out(3) = dot_params_out(3) / (1.5_dp * params(1, ipoi))
        dot_params_out(4) = dot_params_out(4) / (1.5_dp * params(1, ipoi))

    end subroutine compute_source_terms_at_point

    !===========================================================================
    ! Testing helper functions (pure, for unit testing)
    !===========================================================================

    pure subroutine compute_thermodynamic_forces(dn_dr, dTe_dr, dTi_dr, n, Te, Ti, E0r, Z, forces)
        implicit none

        real(dp), intent(in) :: dn_dr, dTe_dr, dTi_dr
        real(dp), intent(in) :: n, Te, Ti
        real(dp), intent(in) :: E0r, Z
        type(thermodynamic_forces_t), intent(out) :: forces

        !
        ! Compute the thermodynamic forces.
        !   A1 = 1/n ∂n/∂r - e/T E0r - 3/(2T) ∂T/∂r
        !   A2 = 1/T ∂T/∂r
        !     [Markl2023 (5), Heyn2014 (22)]
        !
        ! A*_noE denotes the same expressions with the electric field term dropped.
        ! This is the form used for the anomalous E0r=0-frame model [Heyn2014 below (68)].
        !

        forces%e%A1_noE = dn_dr / n - 1.5_dp * dTe_dr / Te
        forces%e%A1 = forces%e%A1_noE + e_charge * E0r / Te
        forces%e%A2 = dTe_dr / Te

        forces%i%A1_noE = dn_dr / n - 1.5_dp * dTi_dr / Ti
        forces%i%A1 = forces%i%A1_noE - Z * e_charge * E0r / Ti
        forces%i%A2 = dTi_dr / Ti
    end subroutine compute_thermodynamic_forces

    pure subroutine compute_particle_fluxes(forces, n, Z, Dae11, Dae12, Dqle11, Dqle12, Dai11, &
                                            Dai12, Dqli11, Dqli12, fluxes)
        !
        ! Compute particle fluxes for electrons and ions.
        !   Γ = -n (D11 A1 + D12 A2)
        !     [Markl2023 (29), Heyn2014 (35)]
        ! Same equations for RMP (EM) and turbulence (A) due to the E0r=0 frame.
        ! EM equation uses D^ql coefficients, A uses D^a.
        !

        implicit none

        type(thermodynamic_forces_t), intent(in) :: forces
        real(dp), intent(in) :: n, Z
        real(dp), intent(in) :: Dae11, Dae12, Dqle11, Dqle12
        real(dp), intent(in) :: Dai11, Dai12, Dqli11, Dqli12
        type(transport_fluxes_t), intent(out) :: fluxes

        ! Electron particle fluxes
        fluxes%e%Gamma_a = -(Dae11 * forces%e%A1_noE + Dae12 * forces%e%A2) * n
        fluxes%e%Gamma_ql = -(Dqle11 * forces%e%A1 + Dqle12 * forces%e%A2) * n
        fluxes%e%Gamma_tot = fluxes%e%Gamma_a + fluxes%e%Gamma_ql

        ! Ion particle fluxes
        fluxes%i%Gamma_a = -(Dai11 * forces%i%A1_noE + Dai12 * forces%i%A2) * n / Z
        fluxes%i%Gamma_ql = -(Dqli11 * forces%i%A1 + Dqli12 * forces%i%A2) * n / Z
        fluxes%i%Gamma_tot = fluxes%i%Gamma_a + fluxes%i%Gamma_ql
    end subroutine compute_particle_fluxes

    pure subroutine compute_total_heat_fluxes(forces, n, Te, Ti, Z, Dae21, Dqle21, Dae22, Dqle22, &
                                              Dai21, Dqli21, Dai22, Dni22, Dqli22, Qe, Qi)
        !
        ! Compute heat flux densities for electrons and ions.
        !
        ! Flux-force relation for heat flux:
        !   Q = -n T (D21 A1 + D22 A2)
        !     [Markl2023 (29), Heyn2014 (35)]
        ! Same equations for RMP (EM) and turbulence (A) due to the E0r=0 frame.
        ! EM equation uses D^ql coefficients, A uses D^a.
        ! Relations between D^a's:
        !   D^a_12 = D^a_21 = 3/2 D^a_11
        !   D^a_22 = 15/4 D^a_11
        !     [Heyn2014 below (68)]
        ! The total heat flux is then
        !   Q_e = Q^EM_e + Q^A_e
        !   Q_i = Q^EM_i + Q^A_i + Q^NEO_i
        !     [Markl2023 (24)-(25), Heyn2014 (65)-(66)]
        ! The neoclassical ion heat flux contribution is included via Dni22.
        !
        ! Note: Rearranging the computation breaks the golden record test.
        ! Needs investigation.
        !

        implicit none

        type(thermodynamic_forces_t), intent(in) :: forces
        real(dp), intent(in) :: n, Te, Ti, Z
        real(dp), intent(in) :: Dae21, Dqle21, Dae22, Dqle22
        real(dp), intent(in) :: Dai21, Dqli21, Dai22, Dni22, Dqli22
        real(dp), intent(out) :: Qe, Qi

        real(dp) :: De22, Di22

        De22 = Dae22 + Dqle22
        Di22 = Dai22 + Dni22 + Dqli22
        Qe = -(Dae21 * forces%e%A1_noE + Dqle21 * forces%e%A1 + De22 * forces%e%A2) * n * Te
        Qi = -(Dai21 * forces%i%A1_noE + Dqli21 * forces%i%A1 + Di22 * forces%i%A2) * n / Z * Ti
    end subroutine compute_total_heat_fluxes

    pure subroutine compute_diffusive_parts(dn_dr, dVphi_dr, dTe_dr, dTi_dr, n, Te, Ti, Z, S, &
                                            Dae11, Dqle11, Dae22, Dqle22, Dai22, Dni22, Dqli22, &
                                            Dqli21, mu_perp, g_phi_phi, flux_diffusion)
        !
        ! Compute diffusive fluxes for all four balance equations.
        ! These are only the operands of the divergence operator in each
        ! equation [Markl2023 (22)-(25), Heyn2014 (63)-(66)].
        !
        ! These expressions are algebraic rearrangements of the flux-force relations.
        ! From RMP:
        !   Γ^EM = -n (D^ql_11 A1 + D^ql_12 A2)
        !   Q^EM = -n T (D^ql_21 A1 + D^ql_22 A2)
        !     [Markl2023 (29), Heyn2014 (35)]
        ! From turbulence:
        !   Γ^A = -D^a_11 ∂n/∂r
        !   Q^A = -3/2 D^a_11 ∂(ni Ti)/∂r
        !     [Markl2023 (27), Heyn2014 (68)]
        ! In the E0r=0-frame used for turbulence [Markl2023 (29), Heyn2014 (35)] reduce
        ! to [Markl2023 (27), Heyn2014 (68)] as well:
        !   Γ^A = -n (D^a_11 A1 + D^a_12 A2)
        !   Q^A = -n T (D^a_21 A1 + D^a_22 A2)
        !     [Markl2023 (29), Heyn2014 (35)]
        ! Relevant formulas to compute the rearrangements by hand:
        !   A1 = 1/n ∂n/∂r - e/T E0r - 3/(2T) ∂T/∂r
        !   A2 = 1/T ∂T/∂r
        !     [Markl2023 (5), Heyn2014 (22)]
        !   E0r = 1/(ei ni) ∂(ni Ti)/∂r + sqrt(g) B^θ_0 / c (V^φ - q V^θ)
        !     [Markl2023 (26), Heyn2014 (71)]
        ! and for completeness:
        !   ni = ne / Z = n / Z
        !   ei = -ee Z = e Z
        !
        ! Note: Rearranging the computation breaks the golden record test.
        ! Needs investigation.
        !

        implicit none

        real(dp), intent(in) :: dn_dr, dVphi_dr, dTe_dr, dTi_dr
        real(dp), intent(in) :: n, Te, Ti, Z, S
        real(dp), intent(in) :: Dae11, Dqle11, Dae22, Dqle22
        real(dp), intent(in) :: Dai22, Dni22, Dqli22, Dqli21
        real(dp), intent(in) :: mu_perp, g_phi_phi
        real(dp), dimension(4), intent(out) :: flux_diffusion

        real(dp) :: dfluxvphi

        ! Eq. 1: Particle flux (diffusive part)
        ! Insert A1e, E0r and A2e into Γe
        ! Drop all terms without ∂ne/∂r = Z ∂ni/∂r
        flux_diffusion(1) = -S * dn_dr * (Dae11 + Dqle11 * (1.0_dp + Ti / (Z * Te)))

        ! Eq. 2: Toroidal momentum diffusion (viscous term)
        ! Divide the equation by the factor mi to get velocity instead of momentum.
        ! This code evolves velocity directly.
        dfluxvphi = -mu_perp * dVphi_dr * n / Z * g_phi_phi
        flux_diffusion(2) = S * dfluxvphi

        ! Eq. 3: Electron heat flux (diffusive part)
        ! TODO: try: -S * n * dTe_dr * (-1.5_dp * Dqle21 + Dqle22 + 0.4_dp * Dae22)
        flux_diffusion(3) = -S * (Dae22 + Dqle22) * n * dTe_dr

        ! Eq. 4: Ion heat flux (diffusive part)
        ! Dni22 is due to Markl2023 (28)
        ! TODO: try: -S * n / Z * dTi_dr * (-2.5_dp * Dqli21 + Dqli22 + 0.4_dp * Dai22 + Dni22)
        flux_diffusion(4) = -S * (Dai22 + Dni22 + Dqli22 - 2.5_dp * Dqli21) * n / Z * dTi_dr
    end subroutine compute_diffusive_parts

    pure subroutine compute_convective_parts(Gamma_tot_e, Qe, Qi, n, Te, Ti, S, flux_diffusion, &
                                             flux_convection)
        !
        ! Compute convective fluxes for all four balance equations.
        !   convective flux = (total flux - diffusive flux) / parameter value
        !
        ! This decomposition is for numerical stability.
        !

        implicit none

        real(dp), intent(in) :: Gamma_tot_e, Qe, Qi
        real(dp), intent(in) :: n, Te, Ti, S
        real(dp), dimension(4), intent(in) :: flux_diffusion
        real(dp), dimension(4), intent(out) :: flux_convection

        ! Eq 1: Particle convection
        flux_convection(1) = (S * Gamma_tot_e - flux_diffusion(1)) / n

        ! Eq 2: Momentum convection
        ! Zero -- only viscous diffusion
        flux_convection(2) = 0.0_dp

        ! Eq 3: Electron heat convection
        flux_convection(3) = (S * Qe - flux_diffusion(3)) / Te

        ! Eq 4: Ion heat convection
        flux_convection(4) = (S * Qi - flux_diffusion(4)) / Ti
    end subroutine compute_convective_parts

end module rhs_balance_m
