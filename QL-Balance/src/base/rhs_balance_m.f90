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
    !
    ! This module provides:
    !   - initialize_rhs: Count non-zeros and allocate sparse matrix
    !   - rhs_balance: Compute sparse Jacobian matrix A
    !   - rhs_balance_source: Compute source vector q
    !
    ! The implicit time stepping scheme uses:
    !   y' - Δt·A(y)·y' = y + Δt·q(y)
    !
    use QLBalance_kinds, only: dp
    use baseparam_mod, only: e_charge

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
    public :: apply_boundary_conditions
    public :: compute_total_fluxes_at_point
    public :: compute_time_derivatives

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
                            ipend, rb, reint_coef, fluxes_dif_lin, fluxes_con_lin, rc, dae11, dae12, &
                            dae22, dai11, dai12, dai22, dni22, visca, gpp_av, dqle11, dqle12, &
                            dqle21, dqle22, dqli11, dqli12, dqli21, dqli22, T_EM_phi_e, &
                            T_EM_phi_i, sqrt_g_times_B_theta_over_c, Ercov, polforce, &
                            qlheat_e, qlheat_i, Ercov_lin, fluxes_con_nl
        use plasma_parameters, only: params, ddr_params_lin, params_lin, ddr_params_nl, params_b_lin, &
                                     params_b, dot_params
        use baseparam_mod, only: Z_i, am, p_mass, c
        use wave_code_data, only: q, Vth
        use matrix_mod, only: isw_rhs, nz, nsize, irow, icol, amat, rhsvec

        implicit none

        real(dp), intent(in) :: x
        real(dp), dimension(neqset), intent(in) :: y
        real(dp), dimension(neqset), intent(out) :: dy

        integer :: ipoi, ieq, i, npoi, ibeg, iend, nshift, ibegb, iendb, ibegtot, iendtot, k, iprobe
        real(dp), dimension(neqset) :: y_lin

        ! Local variables for flux computation
        ! _lin suffix: evaluated at linearized state (for Jacobian construction)
        ! _nl suffix: evaluated at actual plasma state
        type(thermodynamic_forces_t) :: forces_lin, forces_nl
        real(dp) :: Gamma_e_lin, Gamma_i_lin, Gamma_ql_e_lin, Gamma_ql_i_lin, Qe_lin, Qi_lin
        real(dp) :: Gamma_e_nl, Gamma_i_nl, Gamma_ql_e_nl, Gamma_ql_i_nl, Qe_nl, Qi_nl
        real(dp), dimension(4) :: flux_dif_lin_loc, flux_con_lin_loc
        real(dp), dimension(4) :: flux_dif_nl_loc, flux_con_nl_loc
        real(dp), dimension(4) :: dot_params_loc

        if (iboutype .eq. 1) then
            npoi = npoic - 1
        else
            npoi = npoic
        end if

        ! Initialize arrays
        y_lin = 0.0_dp
        params_lin = 0.0_dp
        ddr_params_lin = 0.0_dp
        Ercov_lin = 0.0_dp

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

        call compute_radial_electric_field(npoib, rb, params_b, ddr_params_nl, &
                                           sqrt_g_times_B_theta_over_c, Vth, q, Z_i, &
                                           Ercov)

        call calc_equil_diffusion_coeffs

        ! Compute nonlinear fluxes at all boundary points (needed for fluxes_con_nl)
        do ipoi = 1, npoib
            call compute_fluxes_at_boundary(ipoi, ddr_params_nl, params_b, Ercov(ipoi), dae11, &
                                            dae12, dae22, dai11, dai12, dai22, dni22, dqle11, &
                                            dqle12, dqle21, dqle22, dqli11, dqli12, dqli21, &
                                            dqli22, visca, gpp_av, Sb, Z_i, &
                                            forces_nl, Gamma_e_nl, Gamma_i_nl, Gamma_ql_e_nl, &
                                            Gamma_ql_i_nl, Qe_nl, Qi_nl, flux_dif_nl_loc, &
                                            flux_con_nl_loc)

            call compute_nonlinear_convective_flux(ipoi, Gamma_e_nl, Qe_nl, Qi_nl, &
                                                   ddr_params_nl, params_b, Sb, dae11, dqle11, &
                                                   dae22, dqle22, dai22, dni22, dqli22, &
                                                   dqli21, Z_i, flux_con_nl_loc)

            fluxes_con_nl(:, ipoi) = flux_con_nl_loc
        end do

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

            ! Set probe
            y_lin(iprobe) = 1.0_dp

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

            ! Compute linearized Ercov
            ! TODO: what is this? It's different from the formula in the subroutine
            ! here: Ercov_lin = sqrt(g) B^θ/c * V^φ + (Ti/n * dn/dr + dTi/dr)/(Ze)
            Ercov_lin(ibegb:iendb) = sqrt_g_times_B_theta_over_c(ibegb:iendb) * &
                params_b_lin(2, ibegb:iendb) + (params_b(4, ibegb:iendb) * &
                ddr_params_lin(1, ibegb:iendb) / params_b(1, ibegb:iendb) + &
                ddr_params_lin(4, ibegb:iendb)) / (Z_i * e_charge)

            ! Compute linearized fluxes at affected boundary points
            do ipoi = ibegb, iendb
                call compute_fluxes_at_boundary(ipoi, ddr_params_lin, params_b, Ercov_lin(ipoi), &
                                                dae11, dae12, dae22, dai11, dai12, dai22, dni22, &
                                                dqle11, dqle12, dqle21, dqle22, dqli11, dqli12, &
                                                dqli21, dqli22, visca, gpp_av, Sb, Z_i, &
                                                forces_lin, Gamma_e_lin, Gamma_i_lin, Gamma_ql_e_lin, &
                                                Gamma_ql_i_lin, Qe_lin, Qi_lin, flux_dif_lin_loc, &
                                                flux_con_lin_loc)

                fluxes_dif_lin(:, ipoi) = flux_dif_lin_loc
                fluxes_con_lin(:, ipoi) = flux_con_lin_loc

                ! Recompute nonlinear ql fluxes for T_EM_phi using:
                ! - stale forces_nl (from pre-loop's last ipoi=npoib)
                ! - current diffusion coefficients and density (ipoi)
                ! This matches the original code's behavior
                Gamma_ql_e_nl = -(dqle11(ipoi) * forces_nl%e%A1 + dqle12(ipoi) * &
                                  forces_nl%e%A2) * params_b(1, ipoi)
                Gamma_ql_i_nl = -(dqli11(ipoi) * forces_nl%i%A1 + dqli12(ipoi) * &
                                  forces_nl%i%A2) * params_b(1, ipoi) / Z_i

                ! Compute source terms
                ! Only internal sources (due to the RMP fields are handled).
                ! External forces (heating, fueling, NBI, ...) are set to 0.
                call compute_internal_sources(Gamma_e_lin, Gamma_i_lin, Gamma_ql_e_lin, &
                                              Gamma_ql_i_lin, Gamma_ql_e_nl, Gamma_ql_i_nl, &
                                              Ercov(ipoi), sqrt_g_times_B_theta_over_c(ipoi), Z_i, &
                                              am, p_mass, polforce(ipoi), qlheat_e(ipoi), &
                                              qlheat_i(ipoi), T_EM_phi_e(ipoi), T_EM_phi_i(ipoi))
            end do

            ! Apply boundary conditions
            fluxes_dif_lin(:, 1) = 0.0_dp
            fluxes_con_lin(:, 1) = 0.0_dp
            fluxes_con_nl(:, 1) = 0.0_dp

            ! Compute time derivatives
            do ipoi = ibeg, iend
                call compute_dot_params_at_point(ipoi, npoi, nbaleqs, fluxes_dif_lin, fluxes_con_lin, &
                                                 fluxes_con_nl, params, params_lin, params_b_lin, &
                                                 Sc, rb, rc, gpp_av, polforce, qlheat_e, &
                                                 qlheat_i, Z_i, dot_params_loc)
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
            do i = ibegtot, iendtot
                if (dy(i) /= 0.0_dp) then
                    k = k + 1
                    if (isw_rhs .eq. 1) then
                        irow(k) = i
                        icol(k) = iprobe
                        amat(k) = dy(i)
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
        ! This computes the full RHS using the nonlinear (actual) values.
        ! By setting params_lin = params (not zero), the linear flux contributions
        ! are computed with the actual solution, giving the source vector.
        !
        use grid_mod, only: nbaleqs, neqset, iboutype, npoic, npoib, Sc, Sb, deriv_coef, ipbeg, &
                            ipend, rb, reint_coef, fluxes_dif_lin, fluxes_con_lin, rc, dae11, dae12, &
                            dae22, dai11, dai12, dai22, dni22, visca, gpp_av, dery_equisource, &
                            dqle11, dqle12, dqle21, dqle22, dqli11, dqli12, dqli21, dqli22, &
                            T_EM_phi_e_source, T_EM_phi_i_source, sqrt_g_times_B_theta_over_c, &
                            Ercov, polforce, polforce_ql, qlheat_e, qlheat_i, Ercov_lin, &
                            fluxes_con_nl
        use plasma_parameters, only: params, ddr_params_lin, params_b, params_lin, params_b_lin, &
                                     ddr_params_nl, dot_params
        use baseparam_mod, only: Z_i, am, p_mass, c
        use wave_code_data, only: q, Vth

        implicit none

        real(dp), intent(in) :: x
        real(dp), dimension(neqset), intent(in) :: y
        real(dp), dimension(neqset), intent(out) :: dy

        integer :: ipoi, ieq, i, npoi
        real(dp), dimension(neqset) :: y_lin

        ! Local variables for flux computation
        ! _lin suffix: evaluated at linearized state (y_lin=0 for source computation)
        ! _nl suffix: evaluated at actual plasma state
        type(thermodynamic_forces_t) :: forces_lin, forces_nl
        real(dp) :: Gamma_e_lin, Gamma_i_lin, Gamma_ql_e_lin, Gamma_ql_i_lin, Q_e_lin, Q_i_lin
        real(dp) :: Gamma_e_nl, Gamma_i_nl, Gamma_ql_e_nl, Gamma_ql_i_nl, Q_e_nl, Q_i_nl
        real(dp), dimension(4) :: flux_dif_lin_loc, flux_con_lin_loc
        real(dp), dimension(4) :: flux_dif_nl_loc, flux_con_nl_loc
        real(dp), dimension(4) :: dot_params_loc

        if (iboutype .eq. 1) then
            npoi = npoic - 1
        else
            npoi = npoic
        end if

        ! For source computation: y_lin = 0 (no linear perturbation)
        ! but params_lin = params (actual values)
        y_lin = 0.0_dp
        params_lin = params

        ! Copy y to params
        do ipoi = 1, npoi
            do ieq = 1, nbaleqs
                i = nbaleqs * (ipoi - 1) + ieq
                params(ieq, ipoi) = y(i)
                params_lin(ieq, ipoi) = y_lin(i)
            end do
        end do

        ! Interpolate to boundary points
        do ipoi = 1, npoib
            do ieq = 1, nbaleqs
                ddr_params_nl(ieq, ipoi) = sum(params(ieq, ipbeg(ipoi):ipend(ipoi)) * &
                                               deriv_coef(:, ipoi))
                ddr_params_lin(ieq, ipoi) = sum(params_lin(ieq, ipbeg(ipoi):ipend(ipoi)) * &
                                            deriv_coef(:, ipoi))
                params_b(ieq, ipoi) = sum(params(ieq, ipbeg(ipoi):ipend(ipoi)) * &
                                          reint_coef(:, ipoi))
                params_b_lin(ieq, ipoi) = sum(params_lin(ieq, ipbeg(ipoi):ipend(ipoi)) * &
                                              reint_coef(:, ipoi))
            end do
        end do

        ! Compute radial electric fields
        call compute_radial_electric_field(npoib, rb, params_b, ddr_params_nl, &
                                           sqrt_g_times_B_theta_over_c, Vth, q, Z_i, &
                                           Ercov)

        ! Linear Ercov (with y_lin=0, this simplifies)
        Ercov_lin(1:npoib) = sqrt_g_times_B_theta_over_c(1:npoib) * params_b_lin(2, 1:npoib) + &
            (params_b(4, 1:npoib) * ddr_params_lin(1, 1:npoib) / params_b(1, 1:npoib) + &
            ddr_params_lin(4, 1:npoib)) / (Z_i * e_charge)

        call calc_equil_diffusion_coeffs

        ! Compute fluxes at all boundary points
        do ipoi = 1, npoib
            ! Linearized fluxes (with ddr_params_lin from y_lin=0)
            call compute_fluxes_at_boundary(ipoi, ddr_params_lin, params_b, Ercov_lin(ipoi), dae11, &
                                            dae12, dae22, dai11, dai12, dai22, dni22, dqle11, &
                                            dqle12, dqle21, dqle22, dqli11, dqli12, dqli21, &
                                            dqli22, visca, gpp_av, Sb, Z_i, forces_lin, &
                                            Gamma_e_lin, Gamma_i_lin, Gamma_ql_e_lin, Gamma_ql_i_lin, &
                                            Q_e_lin, Q_i_lin, flux_dif_lin_loc, flux_con_lin_loc)

            fluxes_dif_lin(:, ipoi) = flux_dif_lin_loc
            fluxes_con_lin(:, ipoi) = flux_con_lin_loc

            ! Nonlinear fluxes (with ddr_params_nl from actual y)
            call compute_fluxes_at_boundary(ipoi, ddr_params_nl, params_b, Ercov(ipoi), dae11, &
                                            dae12, dae22, dai11, dai12, dai22, dni22, dqle11, &
                                            dqle12, dqle21, dqle22, dqli11, dqli12, dqli21, &
                                            dqli22, visca, gpp_av, Sb, Z_i, forces_nl, &
                                            Gamma_e_nl, Gamma_i_nl, Gamma_ql_e_nl, Gamma_ql_i_nl, &
                                            Q_e_nl, Q_i_nl, flux_dif_nl_loc, flux_con_nl_loc)

            call compute_nonlinear_convective_flux(ipoi, Gamma_e_nl, Q_e_nl, Q_i_nl, &
                                                   ddr_params_nl, params_b, Sb, dae11, dqle11, &
                                                   dae22, dqle22, dai22, dni22, dqli22, dqli21, &
                                                   Z_i, flux_con_nl_loc)

            fluxes_con_nl(:, ipoi) = flux_con_nl_loc

            ! Source terms
            call compute_internal_sources(Gamma_e_lin, Gamma_i_lin, Gamma_ql_e_lin, &
                                          Gamma_ql_i_lin, Gamma_ql_e_nl, Gamma_ql_i_nl, &
                                          Ercov(ipoi), sqrt_g_times_B_theta_over_c(ipoi), Z_i, am, &
                                          p_mass, polforce(ipoi), qlheat_e(ipoi), &
                                          qlheat_i(ipoi), T_EM_phi_e_source(ipoi), &
                                          T_EM_phi_i_source(ipoi))

            ! Additional source terms specific to rhs_balance_source
            polforce_ql(ipoi) = (T_EM_phi_i_source(ipoi) - T_EM_phi_e_source(ipoi)) / (am * p_mass)
        end do

        ! Apply boundary conditions
        fluxes_dif_lin(:, 1) = 0.0_dp
        fluxes_con_lin(:, 1) = 0.0_dp
        fluxes_con_nl(:, 1) = 0.0_dp

        ! Compute time derivatives
        do ipoi = 1, npoi
            call compute_dot_params_at_point(ipoi, npoi, nbaleqs, fluxes_dif_lin, fluxes_con_lin, &
                                             fluxes_con_nl, params, params_lin, params_b_lin, Sc, &
                                             rb, rc, gpp_av, polforce, qlheat_e, qlheat_i, Z_i, &
                                             dot_params_loc)
            dot_params(:, ipoi) = dot_params_loc
        end do

        ! Copy to dy
        do ipoi = 1, npoi
            do ieq = 1, nbaleqs
                i = nbaleqs * (ipoi - 1) + ieq
                dy(i) = dot_params(ieq, ipoi)
            end do
        end do

        ! Add equilibrium source
        dy = dy + dery_equisource

    end subroutine rhs_balance_source

    !===========================================================================
    ! Core computation routines (shared by rhs_balance and rhs_balance_source)
    !===========================================================================

    subroutine compute_radial_electric_field(npoib, rb, params_b, ddr_params, sqrt_g_Bth_over_c, &
                                             V_pol_arr, q_arr, Z, E0r)
        !
        ! Compute equilibrium radial electric field at all boundary points.
        !   E0r = 1/(ei ni) ∂(ni Ti)/∂r + sqrt(g) B^θ_0 / c (V^φ - q V^θ)
        !     [Markl2023 (26), Heyn2014 (71)]
        !   V^θ = V_pol / r
        !     [Markl2023 direclty above (40)]
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

    subroutine compute_fluxes_at_boundary( &
        ! inputs:
        ipoi, ddr_params, params_b, E0r, Dae11, Dae12, Dae22, Dai11, Dai12, Dai22, Dni22, &
        Dqle11, Dqle12, Dqle21, Dqle22, Dqli11, Dqli12, Dqli21, Dqli22, visca, g_phi_phi, S, Z, &
        ! outputs
        forces, Gamma_e, Gamma_i, Gamma_ql_e, Gamma_ql_i, Qe, Qi, flux_diffusion, flux_convection)
        !
        ! Compute all fluxes at a single boundary point.
        !

        implicit none

        integer, intent(in) :: ipoi
        real(dp), intent(in) :: ddr_params(:, :)  ! (4, npoib)
        real(dp), intent(in) :: params_b(:, :)  ! (4, npoib)
        real(dp), intent(in) :: E0r
        real(dp), intent(in) :: Dae11(:), Dae12(:), Dae22(:)
        real(dp), intent(in) :: Dai11(:), Dai12(:), Dai22(:), Dni22(:)
        real(dp), intent(in) :: Dqle11(:), Dqle12(:), Dqle21(:), Dqle22(:)
        real(dp), intent(in) :: Dqli11(:), Dqli12(:), Dqli21(:), Dqli22(:)
        real(dp), intent(in) :: visca(:), g_phi_phi(:), S(:)
        real(dp), intent(in) :: Z

        type(thermodynamic_forces_t), intent(out) :: forces
        real(dp), intent(out) :: Gamma_e, Gamma_i, Gamma_ql_e, Gamma_ql_i
        real(dp), intent(out) :: Qe, Qi
        real(dp), dimension(4), intent(out) :: flux_diffusion, flux_convection

        type(transport_fluxes_t) :: fluxes
        real(dp) :: n, Te, Ti, dn_dr, dVphi_dr, dTe_dr, dTi_dr
        real(dp) :: De22, Di22

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
        Gamma_e = fluxes%e%Gamma_tot
        Gamma_i = fluxes%i%Gamma_tot
        Gamma_ql_e = fluxes%e%Gamma_ql
        Gamma_ql_i = fluxes%i%Gamma_ql

        De22 = Dae22(ipoi) + Dqle22(ipoi)
        Di22 = Dai22(ipoi) + Dni22(ipoi) + Dqli22(ipoi)
        call compute_total_heat_fluxes(forces, n, Te, Ti, Z, Dae12(ipoi), Dqle21(ipoi), De22, &
                                       Dai12(ipoi), Dqli21(ipoi), Di22, Qe, Qi)

        ! extract contributions from total heat fluxes
        call compute_diffusive_parts(dn_dr, dVphi_dr, dTe_dr, dTi_dr, n, Te, Ti, Z, S(ipoi), &
                                     Dae11(ipoi), Dqle11(ipoi), Dae22(ipoi), Dqle22(ipoi), &
                                     Dai22(ipoi), Dni22(ipoi), Dqli22(ipoi), Dqli21(ipoi), &
                                     visca(ipoi), g_phi_phi(ipoi), flux_diffusion)

        call compute_convective_parts(Gamma_e, Qe, Qi, n, Te, Ti, S(ipoi), flux_diffusion, &
                                      flux_convection)

    end subroutine compute_fluxes_at_boundary

    subroutine compute_nonlinear_convective_flux(ipoi, Gamma_e_nl, Qe_nl, Qi_nl, ddr_params_nl, &
                                                 params_b, S, Dae11, Dqle11, Dae22, Dqle22, Dai22, &
                                                 Dni22, Dqli22, Dqli21, Z, flux_convection_nl)
        !
        ! Compute nonlinear convective flux at a single boundary point.
        ! This is used for the upstream convection scheme.
        !

        implicit none

        integer, intent(in) :: ipoi
        real(dp), intent(in) :: Gamma_e_nl, Qe_nl, Qi_nl
        real(dp), intent(in) :: ddr_params_nl(:, :)
        real(dp), intent(in) :: params_b(:, :)
        real(dp), intent(in) :: S(:)
        real(dp), intent(in) :: Dae11(:), Dqle11(:)
        real(dp), intent(in) :: Dae22(:), Dqle22(:)
        real(dp), intent(in) :: Dai22(:), Dni22(:), Dqli22(:), Dqli21(:)
        real(dp), intent(in) :: Z
        real(dp), dimension(4), intent(out) :: flux_convection_nl

        real(dp) :: n, Te, Ti, dn_dr, dTe_dr, dTi_dr

        ! extract local values
        n = params_b(1, ipoi)
        Te = params_b(3, ipoi)
        Ti = params_b(4, ipoi)
        dn_dr = ddr_params_nl(1, ipoi)
        dTe_dr = ddr_params_nl(3, ipoi)
        dTi_dr = ddr_params_nl(4, ipoi)

        ! Particle flux
        ! TODO: the subtrahend is just the diffusive part calculated in compute_diffusive_parts
        ! maybe we can reuse it here without computing it anew
        flux_convection_nl(1) = (S(ipoi) * Gamma_e_nl - (-S(ipoi) * dn_dr * &
            (Dae11(ipoi) + Dqle11(ipoi) * (1.0_dp + Ti / (Z * Te))))) / n

        ! Momentum flux
        flux_convection_nl(2) = 0.0_dp

        ! Electron heat flux
        ! TODO: same thing here as for particle flux
        flux_convection_nl(3) = (S(ipoi) * Qe_nl - (-S(ipoi) * (Dae22(ipoi) + Dqle22(ipoi)) * n * &
            dTe_dr)) / Te

        ! Ion heat flux
        !colli - 2.5d0*dqli12 was changed to dqli21
        ! TODO: and again, same thing
        flux_convection_nl(4) = (S(ipoi) * Qi_nl - (-S(ipoi) * (Dai22(ipoi) + Dni22(ipoi) + &
            Dqli22(ipoi) - 2.5_dp * Dqli21(ipoi)) * n / Z * dTi_dr)) / Ti
    end subroutine compute_nonlinear_convective_flux

    pure subroutine compute_internal_sources(Gamma_e_lin, Gamma_i_lin, Gamma_ql_e_lin, &
                                             Gamma_ql_i_lin, Gamma_ql_e_nl, Gamma_ql_i_nl, E0r, &
                                             sqrt_g_Bth_over_c, Z, am, p_mass, polforce, &
                                             qlheat_e, qlheat_i, torque_e_nl, torque_i_nl)
        !
        ! Compute internal source terms of the four balance equations.
        ! The RHS for the four equations (with external sources set to zero) are as follows:
        !   1.) = 0
        !         [Markl2023 (22), Heyn2014 (63)]
        !   2.) = Σ_{e,i} -sqrt(g) B^θ_0 / c e Γ^EM = T^EM_φ,e + T^EM_φ,i
        !         [Markl2023 (23), Heyn2014 (64), (67)]
        !   3.) = -E0r ee Γ^EM_e = +E0r e Γ^EM_e
        !         [Markl2023 (24), Heyn2014 (65)]
        !   4.) = -E0r ei Γ^EM_i = E0r Z e Γ^EM_i
        !         [Markl2023 (24), Heyn2014 (65)]
        !

        implicit none

        real(dp), intent(in) :: Gamma_e_lin, Gamma_i_lin
        real(dp), intent(in) :: Gamma_ql_e_lin, Gamma_ql_i_lin
        real(dp), intent(in) :: Gamma_ql_e_nl, Gamma_ql_i_nl
        real(dp), intent(in) :: E0r
        real(dp), intent(in) :: sqrt_g_Bth_over_c
        real(dp), intent(in) :: Z, am, p_mass
        real(dp), intent(out) :: polforce, qlheat_e, qlheat_i
        real(dp), intent(out) :: torque_e_nl, torque_i_nl

        ! Eq. 1: Particle density
        ! no internal sources

        ! Eq. 2: Toroidal rotation frequency
        ! Divide by the factor mi to get velocity instead of momentum.
        ! This code evolves velocity directly.
        torque_e_nl = +sqrt_g_Bth_over_c * e_charge * Gamma_ql_e_nl
        torque_i_nl = -sqrt_g_Bth_over_c * Z * e_charge * Gamma_ql_i_nl

        ! Eq. 3: Electron temperature
        qlheat_e = -e_charge * E0r * Gamma_ql_e_lin  !+????

        ! Eq. 4: Ion temperature
        qlheat_i = Z * e_charge * E0r * Gamma_ql_i_lin

        ! TODO: Why do we need this?
        ! Net radial current drive (Γe - Z Γi) entering the momentum evolution as a source term.
        ! This is just a difference of torque, but with other Γ
        polforce = (Gamma_e_lin - Z * Gamma_i_lin) * e_charge * sqrt_g_Bth_over_c / (am * p_mass)
    end subroutine compute_internal_sources

    subroutine compute_dot_params_at_point(ipoi, npoi, nbaleqs, fluxes_dif, fluxes_con, &
                                           fluxes_con_nl, params, params_lin, params_b_lin, Sc, &
                                           rb, rc, gpp_av, polforce, qlheat_e, qlheat_i, Z_i, &
                                           dot_params_out)
        !
        ! Compute time derivatives at a single grid point.
        !
        ! This combines:
        !   1. Flux divergence
        !   2. Upstream convection
        !   3. Internal sources
        !   4. Unit conversions (momentum -> omega, nT -> T)
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
        real(dp) :: convel, dr

        dr = rb(ipoi + 1) - rb(ipoi)

        ! Flux divergence for all 4 equations:
        do ieq = 1, nbaleqs
            dot_params_out(ieq) = -(fluxes_dif(ieq, ipoi + 1) - fluxes_dif(ieq, ipoi)) / &
                (Sc(ipoi) * dr) - (fluxes_con(ieq, ipoi + 1) - fluxes_con(ieq, ipoi)) / &
                (Sc(ipoi) * dr) * params(ieq, ipoi)

            ! Upstream convection:
            convel = 0.5_dp * (fluxes_con_nl(ieq, ipoi + 1) + fluxes_con_nl(ieq, ipoi)) / Sc(ipoi)
            if (convel .gt. 0.0_dp) then
                dot_params_out(ieq) = dot_params_out(ieq) - convel * &
                    (params_lin(ieq, ipoi + 1) - params_lin(ieq, ipoi)) / (rc(ipoi + 1) - rc(ipoi))
            else
                if (ipoi .gt. 1) then
                    dot_params_out(ieq) = dot_params_out(ieq) - convel * &
                        (params_lin(ieq, ipoi - 1) - params_lin(ieq, ipoi)) / &
                        (rc(ipoi - 1) - rc(ipoi))
                else
                    dot_params_out(ieq) = dot_params_out(ieq) - convel * &
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

    !===========================================================================
    ! Testing helper functions (pure, for unit testing)
    !===========================================================================

    pure subroutine compute_thermodynamic_forces(dn_dr, dTe_dr, dTi_dr, n, Te, Ti, E0r, Z, &
                                                 forces)
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
        ! A_noE_* denote the same expressions with the electric-field term dropped.
        ! This is the form used for the anomalous E0r=0-frame model [Heyn2014 below (68)].
        !

        forces%e%A1_noE = dn_dr / n - 1.5_dp * dTe_dr / Te
        forces%e%A1 = forces%e%A1_noE + e_charge * E0r / Te
        forces%e%A2 = dTe_dr / Te

        forces%i%A1_noE = dn_dr / n - 1.5_dp * dTi_dr / Ti
        forces%i%A1 = forces%i%A1_noE - Z * e_charge * E0r / Ti
        forces%i%A2 = dTi_dr / Ti
    end subroutine compute_thermodynamic_forces

    pure subroutine compute_particle_fluxes(forces, n, Z, Dae11, Dae12, Dqle11, Dqle12, &
                                            Dai11, Dai12, Dqli11, Dqli12, fluxes)
        !
        ! Compute particle fluxes for electrons and ions.
        !   Γ^EM = -n (D^ql_11 A1 + D^ql_12 A2)
        !     [Markl2023 (29), Heyn2014 (35)]
        !
        ! We split Γ = Γ(A) + Γ(EM) into:
        ! Gamma_a  = -n (D11_a * A_noE_1 + D12_a * A_noE_2)
        ! Gamma_ql = -n (D11_ql * A_1 + D12_ql * A_noE_2)
        ! TODO: I don't get this.
        !
        ! The anomalous part follows the simple E0r=0-frame model which justifies the same
        ! equation as for the RMP-induced flux [Heyn2014 (35), referenced in text below (68)].
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

    pure subroutine compute_total_heat_fluxes(forces, n, Te, Ti, Z, Dae12, Dqle21, De22, Dai12, &
                                              Dqli21, Di22, Qe, Qi)
        !
        ! Compute heat flux densities for electrons and ions.
        !
        ! Flux-force relation for heat flux:
        !   Q^EM = -n T (D^ql_21 A1 + D^ql_22 A2)
        !     [Markl2023 (29), Heyn2014 (35)]
        !
        ! We split anomalous vs quasilinear contributions by coefficient sets:
        ! TODO: I don't see where this formula comes from. Why the dependence on D12? Why sometimes E0r=0?
        ! Q_e = -(D12_a_e*A_noE_1e + D21_ql_e*A_1e + D22_e*A_noE_2e) * n * Te
        ! Q_i = -(D12_a_i*A_noE_1i + D21_ql_i*A_1i + D22_i*A_noE_2i) * n/Z * Ti
        !
        ! The neoclassical ion heat flux contribution is included via dni22 in D22_i
        ! (Heyn2014 Eq. (69); Markl2023 Eq. (25)).
        !
        ! Note: D12 (anomalous) couples to A_noE_1 (no E-field)
        !       D21 (quasi-linear) couples to A_1 (with E-field)
        !       D22 couples to temperature gradient force A_noE_2
        !

        implicit none

        type(thermodynamic_forces_t), intent(in) :: forces
        real(dp), intent(in) :: n, Te, Ti, Z
        real(dp), intent(in) :: Dae12, Dqle21, De22
        real(dp), intent(in) :: Dai12, Dqli21, Di22
        real(dp), intent(out) :: Qe, Qi

        ! Note: Here too. We're outside of machine precision already. Rearranging the computation
        ! significantly changes the results. We need to consider if quad precision is
        ! viable here or if we can normalize some variables to have a "nicer" decimal range.
        Qe = -(Dae12 * forces%e%A1_noE + Dqle21 * forces%e%A1 + De22 * forces%e%A2) * n * Te
        Qi = -(Dai12 * forces%i%A1_noE + Dqli21 * forces%i%A1 + Di22 * forces%i%A2) * n / Z * Ti
    end subroutine compute_total_heat_fluxes

    pure subroutine compute_diffusive_parts(dn_dr, dVphi_dr, dTe_dr, dTi_dr, n, Te, Ti, Z, S, &
                                            Dae11, Dqle11, Dae22, Dqle22, Dai22, Dni22, Dqli22, &
                                            Dqli21, mu_perb, g_phi_phi, flux_diffusion)
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
        ! Relevant formulas to compute this by hand:
        !   A1 = 1/n ∂n/∂r - e/T E0r - 3/(2T) ∂T/∂r
        !   A2 = 1/T ∂T/∂r
        !     [Markl2023 (5), Heyn2014 (22)]
        !   E0r = 1/(ei ni) ∂(ni Ti)/∂r + sqrt(g) B^θ_0 / c (V^φ - q V^θ)
        !     [Markl2023 (26), Heyn2014 (71)]
        ! and for completeness:
        !   ni = ne / Z = n / Z
        !   ei = -ee Z = -e Z
        !
        ! Note: We're outside of machine precision already. Rearranging the computation
        ! significantly changes the results. We need to consider if quad precision is
        ! viable here or if we can normalize some variables to have a "nicer" decimal range.
        !

        implicit none

        real(dp), intent(in) :: dn_dr, dVphi_dr, dTe_dr, dTi_dr
        real(dp), intent(in) :: n, Te, Ti, Z, S
        real(dp), intent(in) :: Dae11, Dqle11, Dae22, Dqle22
        real(dp), intent(in) :: Dai22, Dni22, Dqli22, Dqli21
        real(dp), intent(in) :: mu_perb, g_phi_phi
        real(dp), dimension(4), intent(out) :: flux_diffusion

        real(dp) :: dfluxvphi

        ! Eq. 1: Particle flux (diffusive part)
        ! Insert A1e, E0r and A2e into Γe
        ! Drop all terms without ∂ne/∂r = Z ∂ni/∂r
        flux_diffusion(1) = -S * dn_dr * (Dae11 + Dqle11 * (1.0_dp + Ti / (Z * Te)))

        ! Eq. 2: Toroidal momentum diffusion (viscous term)
        ! Divide by the factor mi to get velocity instead of momentum.
        ! This code evolves velocity directly.
        dfluxvphi = -mu_perb * dVphi_dr * n / Z * g_phi_phi
        flux_diffusion(2) = S * dfluxvphi

        ! Eq. 3: Electron heat flux (diffusive part)
        ! TODO: don't get it. I would arrive at:
        ! S n dTe/dr (3/2 Dqle21 - Dqle22 - 3/2 Dae22)
        flux_diffusion(3) = -S * (Dae22 + Dqle22) * n * dTe_dr

        ! Eq. 4: Ion heat flux (diffusive part)
        ! TODO: don't get it. How do Q^NEO and Dni22 relate?
        flux_diffusion(4) = -S * (Dai22 + Dni22 + Dqli22 - 2.5_dp * Dqli21) * &
                            n / Z * dTi_dr
    end subroutine compute_diffusive_parts

    pure subroutine compute_convective_parts(Gamma_e, Qe, Qi, n, Te, Ti, S, flux_diffusion, &
                                             flux_convection)
        !
        ! Compute convective fluxes for all 4 balance equations.
        ! Convective flux = (total flux - diffusive flux) / parameter value
        !
        ! This is a numerical decomposition used to build an affine RHS:
        !   flux ≈ flux_dif + flux_con * (state variable)
        ! so that dy/dt can be written as A·y + q for implicit stepping.
        !

        implicit none

        real(dp), intent(in) :: Gamma_e, Qe, Qi
        real(dp), intent(in) :: n, Te, Ti, S
        real(dp), dimension(4), intent(in) :: flux_diffusion
        real(dp), dimension(4), intent(out) :: flux_convection

        ! Eq 1: Particle convection
        flux_convection(1) = (S * Gamma_e - flux_diffusion(1)) / n

        ! Eq 2: Momentum convection (zero - only viscous diffusion)
        flux_convection(2) = 0.0_dp

        ! Eq 3: Electron heat convection
        flux_convection(3) = (S * Qe - flux_diffusion(3)) / Te

        ! Eq 4: Ion heat convection
        flux_convection(4) = (S * Qi - flux_diffusion(4)) / Ti
    end subroutine compute_convective_parts

    pure subroutine apply_boundary_conditions(fluxes_dif, fluxes_con, fluxes_con_nl, nbaleqs)
        implicit none

        integer, intent(in) :: nbaleqs
        real(dp), dimension(nbaleqs, *), intent(inout) :: fluxes_dif
        real(dp), dimension(nbaleqs, *), intent(inout) :: fluxes_con
        real(dp), dimension(nbaleqs, *), intent(inout) :: fluxes_con_nl

        fluxes_dif(:, 1) = 0.0_dp
        fluxes_con(:, 1) = 0.0_dp
        fluxes_con_nl(:, 1) = 0.0_dp
    end subroutine apply_boundary_conditions

    pure subroutine compute_total_fluxes_at_point(Gamma_e, Gamma_e_nl, Gamma_i, Q_e, Q_e_nl, Q_i, &
                                                  Q_i_nl, ddr_n, ddr_n_nl, ddr_Te, ddr_Te_nl, &
                                                  ddr_Ti, ddr_Ti_nl, ddr_vphi, n_b, Te_b, Ti_b, &
                                                  Z_i_val, Sb_val, dae11_val, dqle11_val, &
                                                  dae22_val, dqle22_val, dai22_val, dni22_val, &
                                                  dqli22_val, dqli21_val, visca_val, gpp_av_val, &
                                                  flux_dif, flux_con, flux_con_nl)
        implicit none

        real(dp), intent(in) :: Gamma_e, Gamma_e_nl, Gamma_i
        real(dp), intent(in) :: Q_e, Q_e_nl, Q_i, Q_i_nl
        real(dp), intent(in) :: ddr_n, ddr_n_nl, ddr_Te, ddr_Te_nl, ddr_Ti, ddr_Ti_nl, ddr_vphi
        real(dp), intent(in) :: n_b, Te_b, Ti_b, Z_i_val, Sb_val
        real(dp), intent(in) :: dae11_val, dqle11_val, dae22_val, dqle22_val
        real(dp), intent(in) :: dai22_val, dni22_val, dqli22_val, dqli21_val
        real(dp), intent(in) :: visca_val, gpp_av_val
        real(dp), dimension(4), intent(out) :: flux_dif, flux_con, flux_con_nl
        real(dp) :: dfluxvphi

        flux_dif(1) = -Sb_val * ddr_n * (dae11_val + dqle11_val * (1.0_dp + Ti_b / Te_b / Z_i_val))
        flux_con(1) = (Sb_val * Gamma_e - flux_dif(1)) / n_b
        flux_con_nl(1) = (Sb_val * Gamma_e_nl - (-Sb_val * ddr_n_nl * &
                                                 (dae11_val + dqle11_val * &
                                                  (1.0_dp + Ti_b / Te_b / Z_i_val)))) / n_b

        dfluxvphi = -visca_val * ddr_vphi * n_b / Z_i_val * gpp_av_val
        flux_dif(2) = Sb_val * dfluxvphi
        flux_con(2) = 0.0_dp
        flux_con_nl(2) = 0.0_dp

        flux_dif(3) = -Sb_val * (dae22_val + dqle22_val) * n_b * ddr_Te
        flux_con(3) = (Sb_val * Q_e - flux_dif(3)) / Te_b
        flux_con_nl(3) = (Sb_val * Q_e_nl - (-Sb_val * (dae22_val + dqle22_val) * n_b * &
                                             ddr_Te_nl)) / Te_b

        flux_dif(4) = -Sb_val * (dai22_val + dni22_val + dqli22_val - 2.5_dp * dqli21_val) * n_b / &
                      Z_i_val * ddr_Ti
        flux_con(4) = (Sb_val * Q_i - flux_dif(4)) / Ti_b
        flux_con_nl(4) = (Sb_val * Q_i_nl - (-Sb_val * &
                                             (dai22_val + dni22_val + dqli22_val - 2.5_dp * &
                                              dqli21_val) * n_b / Z_i_val * ddr_Ti_nl)) / Ti_b
    end subroutine compute_total_fluxes_at_point

    pure subroutine compute_time_derivatives(dot_params_in, n_val, Te_val, Ti_val, Z_i_val, &
                                             gpp_av_avg, dot_params_out)
        implicit none

        real(dp), dimension(4), intent(in) :: dot_params_in
        real(dp), intent(in) :: n_val, Te_val, Ti_val, Z_i_val, gpp_av_avg
        real(dp), dimension(4), intent(out) :: dot_params_out

        dot_params_out(1) = dot_params_in(1)
        dot_params_out(2) = dot_params_in(2) * Z_i_val / n_val * 2.0_dp / gpp_av_avg
        dot_params_out(3) = (-Te_val * dot_params_in(1) + dot_params_in(3) / 1.5_dp) / n_val
        dot_params_out(4) = (-Ti_val * dot_params_in(1) + dot_params_in(4) / 1.5_dp) / n_val
    end subroutine compute_time_derivatives

end module rhs_balance_m
