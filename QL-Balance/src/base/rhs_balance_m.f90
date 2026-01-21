module rhs_balance_m
    !
    ! RHS Balance Module
    !
    ! Solves the discretized transport equations:
    !   dy/dt = A·y + q
    !
    ! where:
    !   y   - vector of plasma parameters at grid points (size 4N)
    !   A   - sparse Jacobian matrix
    !   q   - source vector
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
    public :: compute_heat_fluxes
    public :: compute_diffusive_fluxes
    public :: compute_convective_fluxes
    public :: compute_source_terms
    public :: apply_boundary_conditions
    public :: compute_total_fluxes_at_point
    public :: compute_time_derivatives

    public :: thermodynamic_forces_t
    public :: species_fluxes_t
    public :: transport_fluxes_t

    !---------------------------------------------------------------------------
    ! Derived types
    !---------------------------------------------------------------------------

    !> Thermodynamic forces at a single radial point
    type :: thermodynamic_forces_t
        real(dp) :: A_noE_1e = 0.0_dp  !< Electron force without E-field
        real(dp) :: A_noE_2e = 0.0_dp  !< Electron temperature gradient force
        real(dp) :: A_noE_1i = 0.0_dp  !< Ion force without E-field
        real(dp) :: A_noE_2i = 0.0_dp  !< Ion temperature gradient force
        real(dp) :: A_1e = 0.0_dp  !< Electron force with E-field
        real(dp) :: A_1i = 0.0_dp  !< Ion force with E-field
    end type thermodynamic_forces_t

    !> Particle and heat fluxes for a single species
    type :: species_fluxes_t
        real(dp) :: gamma_a = 0.0_dp  !< Anomalous particle flux
        real(dp) :: gamma_ql = 0.0_dp  !< Quasi-linear particle flux
        real(dp) :: gamma = 0.0_dp  !< Total particle flux
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

        real(dp) :: x
        real(dp), dimension(neqset) :: y, dy

        x = 0.0_dp
        isw_rhs = 0

        call rhs_balance(x, y, dy)

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
                            T_EM_phi_i, sqrt_g_times_B_theta_over_c, Ercov, polforce, &
                            qlheat_e, qlheat_i, Ercov_lin, fluxes_con
        use plasma_parameters, only: params, ddr_params_lin, params_lin, ddr_params, params_b_lin, &
                                     params_b, dot_params
        use baseparam_mod, only: Z_i, am, p_mass, c
        use wave_code_data, only: q, Vth
        use matrix_mod, only: isw_rhs, nz, nsize, irow, icol, amat, rhsvec

        implicit none

        integer :: ipoi, ieq, i, npoi, ibeg, iend, nshift, ibegb, iendb, ibegtot, iendtot, k, iprobe
        real(dp) :: x
        real(dp), dimension(neqset) :: y, dy, y_lin

        ! Local variables for flux computation
        ! _lin suffix: evaluated at linearized/probed state (for Jacobian construction)
        ! no suffix: evaluated at actual plasma state
        type(thermodynamic_forces_t) :: forces_lin, forces
        real(dp) :: gamma_e_lin, gamma_i_lin, gamma_ql_e_lin, gamma_ql_i_lin, Q_e_lin, Q_i_lin
        real(dp) :: gamma_e, gamma_i, gamma_ql_e, gamma_ql_i, Q_e, Q_i
        real(dp) :: flux_dif_lin_loc(4), flux_con_lin_loc(4), flux_dif_loc(4), flux_con_loc(4)
        real(dp) :: dot_params_loc(4)

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
                ddr_params(ieq, ipoi) = sum(params(ieq, ipbeg(ipoi):ipend(ipoi)) * &
                                            deriv_coef(:, ipoi))
                params_b(ieq, ipoi) = sum(params(ieq, ipbeg(ipoi):ipend(ipoi)) * &
                                          reint_coef(:, ipoi))
            end do
        end do

        call compute_radial_electric_field(npoib, rb, params_b, ddr_params, &
                                           sqrt_g_times_B_theta_over_c, Vth, q, Z_i, Ercov)

        call calc_equil_diffusion_coeffs

        ! Compute fluxes at actual state for all boundary points (needed for fluxes_con)
        do ipoi = 1, npoib
            call compute_fluxes_at_boundary(ipoi, ddr_params, params_b, Ercov(ipoi), dae11, &
                                            dae12, dae22, dai11, dai12, dai22, dni22, dqle11, &
                                            dqle12, dqle21, dqle22, dqli11, dqli12, dqli21, &
                                            dqli22, visca, gpp_av, Sb, Z_i, forces, &
                                            gamma_e, gamma_i, gamma_ql_e, &
                                            gamma_ql_i, Q_e, Q_i, flux_dif_loc, &
                                            flux_con_loc)

            call compute_nonlinear_convective_flux(ipoi, gamma_e, Q_e, Q_i, &
                                                   ddr_params, params_b, Sb, dae11, dqle11, &
                                                   dae22, dqle22, dai22, dni22, dqli22, &
                                                   dqli21, Z_i, flux_con_loc)

            fluxes_con(:, ipoi) = flux_con_loc
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
                    ddr_params_lin(ieq, ipoi) = sum(params_lin(ieq, ipbeg(ipoi):ipend(ipoi)) * &
                                                    deriv_coef(:, ipoi))
                    params_b_lin(ieq, ipoi) = sum(params_lin(ieq, ipbeg(ipoi):ipend(ipoi)) * &
                                                  reint_coef(:, ipoi))
                end do
            end do

            ! Compute linearized Ercov
            Ercov_lin(ibegb:iendb) = sqrt_g_times_B_theta_over_c(ibegb:iendb) * &
                params_b_lin(2, ibegb:iendb) + (params_b(4, ibegb:iendb) * &
                ddr_params_lin(1, ibegb:iendb) / params_b(1, ibegb:iendb) + &
                ddr_params_lin(4, ibegb:iendb)) / (Z_i * e_charge)

            ! Compute linearized fluxes at affected boundary points
            do ipoi = ibegb, iendb
                call compute_fluxes_at_boundary(ipoi, ddr_params_lin, params_b, Ercov_lin(ipoi), &
                                                dae11, dae12, dae22, dai11, dai12, dai22, dni22, &
                                                dqle11, dqle12, dqle21, dqle22, dqli11, dqli12, &
                                                dqli21, dqli22, visca, gpp_av, Sb, Z_i, forces_lin, &
                                                gamma_e_lin, gamma_i_lin, gamma_ql_e_lin, &
                                                gamma_ql_i_lin, Q_e_lin, Q_i_lin, flux_dif_lin_loc, &
                                                flux_con_lin_loc)

                fluxes_dif_lin(:, ipoi) = flux_dif_lin_loc
                fluxes_con_lin(:, ipoi) = flux_con_lin_loc

                ! Recompute ql fluxes for T_EM_phi using:
                ! - stale forces (from pre-loop's last ipoi=npoib)
                ! - current diffusion coefficients and density (ipoi)
                ! This matches the original code's behavior
                gamma_ql_e = -(dqle11(ipoi) * forces%A_1e + dqle12(ipoi) * &
                               forces%A_noE_2e) * params_b(1, ipoi)
                gamma_ql_i = -(dqli11(ipoi) * forces%A_1i + dqli12(ipoi) * &
                               forces%A_noE_2i) * params_b(1, ipoi) / Z_i

                ! Compute source terms (using actual ql fluxes for T_EM_phi)
                call compute_internal_sources(ipoi, gamma_e_lin, gamma_i_lin, gamma_ql_e_lin, &
                                              gamma_ql_i_lin, gamma_ql_e, gamma_ql_i, Ercov(ipoi), &
                                              sqrt_g_times_B_theta_over_c, Z_i, am, p_mass, &
                                              polforce(ipoi), qlheat_e(ipoi), qlheat_i(ipoi), &
                                              T_EM_phi_e(ipoi), T_EM_phi_i(ipoi))
            end do

            ! Apply boundary conditions
            fluxes_dif_lin(:, 1) = 0.0_dp
            fluxes_con_lin(:, 1) = 0.0_dp
            fluxes_con(:, 1) = 0.0_dp

            ! Compute time derivatives
            do ipoi = ibeg, iend
                call compute_dot_params_at_point(ipoi, npoi, nbaleqs, fluxes_dif_lin, &
                                                 fluxes_con_lin, fluxes_con, params, params_lin, &
                                                 params_b_lin, Sc, rb, rc, gpp_av, polforce, &
                                                 qlheat_e, qlheat_i, Z_i, dot_params_loc)
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
                            ipend, rb, reint_coef, fluxes_dif_lin, fluxes_con_lin, rc, dae11, &
                            dae12, dae22, dai11, dai12, dai22, dni22, visca, gpp_av, &
                            dery_equisource, dqle11, dqle12, dqle21, dqle22, dqli11, dqli12, &
                            dqli21, dqli22, T_EM_phi_e_source, T_EM_phi_i_source, &
                            sqrt_g_times_B_theta_over_c, Ercov, polforce, polforce_ql, qlheat_e, &
                            qlheat_i, Ercov_lin, fluxes_con
        use plasma_parameters, only: params, ddr_params_lin, params_b, params_lin, params_b_lin, &
                                     ddr_params, dot_params
        use baseparam_mod, only: Z_i, am, p_mass, c
        use wave_code_data, only: q, Vth

        implicit none

        integer :: ipoi, ieq, i, npoi
        real(dp) :: x
        real(dp), dimension(neqset) :: y, dy, y_lin

        ! Local variables for flux computation
        ! _lin suffix: evaluated at linearized state (y_lin=0 for source computation)
        ! no suffix: evaluated at actual plasma state
        type(thermodynamic_forces_t) :: forces_lin, forces
        real(dp) :: gamma_e_lin, gamma_i_lin, gamma_ql_e_lin, gamma_ql_i_lin, Q_e_lin, Q_i_lin
        real(dp) :: gamma_e, gamma_i, gamma_ql_e, gamma_ql_i, Q_e, Q_i
        real(dp) :: flux_dif_lin_loc(4), flux_con_lin_loc(4), flux_dif_loc(4), flux_con_loc(4)
        real(dp) :: dot_params_loc(4)

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
                ddr_params(ieq, ipoi) = sum(params(ieq, ipbeg(ipoi):ipend(ipoi)) * &
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
        call compute_radial_electric_field(npoib, rb, params_b, ddr_params, &
                                           sqrt_g_times_B_theta_over_c, Vth, q, Z_i, Ercov)

        ! Linearized Ercov (with y_lin=0, this simplifies)
        Ercov_lin(1:npoib) = sqrt_g_times_B_theta_over_c(1:npoib) * params_b_lin(2, 1:npoib) + &
            (params_b(4, 1:npoib) * ddr_params_lin(1, 1:npoib) / params_b(1, 1:npoib) + &
            ddr_params_lin(4, 1:npoib)) / (Z_i * e_charge)

        call calc_equil_diffusion_coeffs

        ! Compute fluxes at all boundary points
        do ipoi = 1, npoib
            ! Linearized fluxes (with ddr_params_lin from y_lin=0)
            call compute_fluxes_at_boundary(ipoi, ddr_params_lin, params_b, Ercov_lin(ipoi), &
                                            dae11, dae12, dae22, dai11, dai12, dai22, dni22, &
                                            dqle11, dqle12, dqle21, dqle22, dqli11, dqli12, &
                                            dqli21, dqli22, visca, gpp_av, Sb, Z_i, forces_lin, &
                                            gamma_e_lin, gamma_i_lin, gamma_ql_e_lin, &
                                            gamma_ql_i_lin, Q_e_lin, Q_i_lin, flux_dif_lin_loc, &
                                            flux_con_lin_loc)

            fluxes_dif_lin(:, ipoi) = flux_dif_lin_loc
            fluxes_con_lin(:, ipoi) = flux_con_lin_loc

            ! Fluxes at actual state (with ddr_params from actual y)
            call compute_fluxes_at_boundary(ipoi, ddr_params, params_b, Ercov(ipoi), dae11, &
                                            dae12, dae22, dai11, dai12, dai22, dni22, dqle11, &
                                            dqle12, dqle21, dqle22, dqli11, dqli12, dqli21, &
                                            dqli22, visca, gpp_av, Sb, Z_i, forces, &
                                            gamma_e, gamma_i, gamma_ql_e, &
                                            gamma_ql_i, Q_e, Q_i, flux_dif_loc, &
                                            flux_con_loc)

            call compute_nonlinear_convective_flux(ipoi, gamma_e, Q_e, Q_i, &
                                                   ddr_params, params_b, Sb, dae11, dqle11, &
                                                   dae22, dqle22, dai22, dni22, dqli22, dqli21, &
                                                   Z_i, flux_con_loc)

            fluxes_con(:, ipoi) = flux_con_loc

            ! Source terms
            call compute_internal_sources(ipoi, gamma_e_lin, gamma_i_lin, gamma_ql_e_lin, &
                                          gamma_ql_i_lin, gamma_ql_e, gamma_ql_i, Ercov(ipoi), &
                                          sqrt_g_times_B_theta_over_c, Z_i, am, p_mass, &
                                          polforce(ipoi), qlheat_e(ipoi), qlheat_i(ipoi), &
                                          T_EM_phi_e_source(ipoi), T_EM_phi_i_source(ipoi))

            ! Additional source terms specific to rhs_balance_source
            polforce_ql(ipoi) = (T_EM_phi_i_source(ipoi) - T_EM_phi_e_source(ipoi)) / (am * p_mass)
        end do

        ! Apply boundary conditions
        fluxes_dif_lin(:, 1) = 0.0_dp
        fluxes_con_lin(:, 1) = 0.0_dp
        fluxes_con(:, 1) = 0.0_dp

        ! Compute time derivatives
        do ipoi = 1, npoi
            call compute_dot_params_at_point(ipoi, npoi, nbaleqs, fluxes_dif_lin, fluxes_con_lin, &
                                             fluxes_con, params, params_lin, params_b_lin, Sc, &
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
                                             Vth_arr, q_arr, Z_i, Ercov_out)
        !
        ! Compute equilibrium radial electric field at all boundary points
        ! E0r = -∂r Φ0 = (sqrt(g) B0θ / c) (Viφ - q Viθ) + ∂r(ni Ti) / (ei ni)
        ! (Heyn2014 (71), Markl2023 (26))
        !

        implicit none

        integer, intent(in) :: npoib
        real(dp), intent(in) :: rb(:)
        real(dp), intent(in) :: params_b(:, :)  ! (4, npoib)
        real(dp), intent(in) :: ddr_params(:, :)  ! (4, npoib)
        real(dp), intent(in) :: sqrt_g_Bth_over_c(:)
        real(dp), intent(in) :: Vth_arr(:), q_arr(:)
        real(dp), intent(in) :: Z_i
        real(dp), intent(out) :: Ercov_out(:)

        associate ( &
            Vphi => params_b(2, 1:npoib), &
            Ti => params_b(4, 1:npoib), &
            dni_dr => ddr_params(1, 1:npoib), &
            ni => params_b(1, 1:npoib), &
            dTi_dr => ddr_params(4, 1:npoib), &
            q => q_arr(1:npoib), &
            rbp => rb(1:npoib) &
        )
            Ercov_out(1:npoib) = &
                (Ti / ni * dni_dr + dTi_dr) / (Z_i * e_charge) + &
                sqrt_g_Bth_over_c(1:npoib) * (Vphi - q * Vth_arr(1:npoib) / rbp)
        end associate

    end subroutine compute_radial_electric_field

    subroutine compute_fluxes_at_boundary( &
        ! inputs:
        ipoi, ddr_params, params_b, Ercov_val, dae11, dae12, dae22, dai11, dai12, dai22, dni22, &
        dqle11, dqle12, dqle21, dqle22, dqli11, dqli12, dqli21, dqli22, visca, gpp_av, Sb, &
        Z_i, &
        ! outputs
        forces, gamma_e, gamma_i, gamma_ql_e, gamma_ql_i, Q_e, Q_i, flux_dif, flux_con)
        !
        ! Compute all fluxes at a single boundary point.
        ! This is the core shared computation used by both rhs_balance and rhs_balance_source.
        !
        ! Delegates to pure helper functions for each computation step.
        !

        implicit none

        integer, intent(in) :: ipoi
        real(dp), intent(in) :: ddr_params(:, :)  ! (4, npoib)
        real(dp), intent(in) :: params_b(:, :)  ! (4, npoib)
        real(dp), intent(in) :: Ercov_val
        real(dp), intent(in) :: dae11(:), dae12(:), dae22(:)
        real(dp), intent(in) :: dai11(:), dai12(:), dai22(:), dni22(:)
        real(dp), intent(in) :: dqle11(:), dqle12(:), dqle21(:), dqle22(:)
        real(dp), intent(in) :: dqli11(:), dqli12(:), dqli21(:), dqli22(:)
        real(dp), intent(in) :: visca(:), gpp_av(:), Sb(:)
        real(dp), intent(in) :: Z_i

        type(thermodynamic_forces_t), intent(out) :: forces
        real(dp), intent(out) :: gamma_e, gamma_i, gamma_ql_e, gamma_ql_i
        real(dp), intent(out) :: Q_e, Q_i
        real(dp), intent(out) :: flux_dif(4), flux_con(4)

        type(transport_fluxes_t) :: fluxes
        real(dp) :: n_b, Te_b, Ti_b

        ! Extract local values
        n_b = params_b(1, ipoi)
        Te_b = params_b(3, ipoi)
        Ti_b = params_b(4, ipoi)

        call compute_thermodynamic_forces(ddr_params(1, ipoi), ddr_params(3, ipoi), ddr_params(4, &
                                           ipoi), n_b, Te_b, Ti_b, Ercov_val, Z_i, forces)

        call compute_particle_fluxes(forces, n_b, Z_i, dae11(ipoi), dae12(ipoi), dqle11(ipoi), &
                                     dqle12(ipoi), dai11(ipoi), dai12(ipoi), dqli11(ipoi), &
                                     dqli12(ipoi), fluxes)

        ! Extract individual flux values for output and further use
        gamma_e = fluxes%e%gamma
        gamma_i = fluxes%i%gamma
        gamma_ql_e = fluxes%e%gamma_ql
        gamma_ql_i = fluxes%i%gamma_ql

        call compute_heat_fluxes(forces, n_b, Te_b, Ti_b, Z_i, dae12(ipoi), dqle21(ipoi), &
                                 dae22(ipoi) + dqle22(ipoi), dai12(ipoi), dqli21(ipoi), &
                                 dai22(ipoi) + dni22(ipoi) + dqli22(ipoi), Q_e, Q_i)

        call compute_diffusive_fluxes(ddr_params(1, ipoi), ddr_params(2, ipoi), &
                                      ddr_params(3, ipoi), ddr_params(4, ipoi), n_b, Te_b, Ti_b, &
                                      Z_i, Sb(ipoi), dae11(ipoi), dqle11(ipoi), dae22(ipoi), &
                                      dqle22(ipoi), dai22(ipoi), dni22(ipoi), dqli22(ipoi), &
                                      dqli21(ipoi), visca(ipoi), gpp_av(ipoi), flux_dif)

        call compute_convective_fluxes(gamma_e, Q_e, Q_i, n_b, Te_b, Ti_b, Sb(ipoi), &
                                       flux_dif, flux_con)

    end subroutine compute_fluxes_at_boundary

    subroutine compute_nonlinear_convective_flux(ipoi, gamma_e, Q_e, Q_i, ddr_params, &
                                                 params_b, Sb, dae11, dqle11, dae22, dqle22, &
                                                 dai22, dni22, dqli22, dqli21, Z_i, flux_con)
        !
        ! Compute convective flux at actual state for a single boundary point.
        ! This is used for the upstream convection scheme.
        !

        implicit none

        integer, intent(in) :: ipoi
        real(dp), intent(in) :: gamma_e, Q_e, Q_i
        real(dp), intent(in) :: ddr_params(:, :)
        real(dp), intent(in) :: params_b(:, :)
        real(dp), intent(in) :: Sb(:)
        real(dp), intent(in) :: dae11(:), dqle11(:)
        real(dp), intent(in) :: dae22(:), dqle22(:)
        real(dp), intent(in) :: dai22(:), dni22(:), dqli22(:), dqli21(:)
        real(dp), intent(in) :: Z_i
        real(dp), intent(out) :: flux_con(4)

        real(dp) :: n_b, Te_b, Ti_b

        n_b = params_b(1, ipoi)
        Te_b = params_b(3, ipoi)
        Ti_b = params_b(4, ipoi)

        ! Particle flux
        flux_con(1) = (Sb(ipoi) * gamma_e - (-Sb(ipoi) * ddr_params(1, ipoi) * &
            (dae11(ipoi) + dqle11(ipoi) * (1.0_dp + Ti_b / Te_b / Z_i)))) / n_b

        ! Momentum flux
        flux_con(2) = 0.0_dp

        ! Electron heat flux
        flux_con(3) = (Sb(ipoi) * Q_e - (-Sb(ipoi) * (dae22(ipoi) + dqle22(ipoi)) * n_b * &
            ddr_params(3, ipoi))) / Te_b

        ! Ion heat flux
        flux_con(4) = (Sb(ipoi) * Q_i - (-Sb(ipoi) * &
            (dai22(ipoi) + dni22(ipoi) + dqli22(ipoi) - 2.5_dp * dqli21(ipoi)) * n_b / Z_i * &
            ddr_params(4, ipoi))) / Ti_b

    end subroutine compute_nonlinear_convective_flux

    subroutine compute_internal_sources(ipoi, gamma_e, gamma_i, gamma_ql_e, gamma_ql_i, &
                                        gamma_ql_e_nl, gamma_ql_i_nl, Ercov_val, sqrt_g_B_theta_c, &
                                        Z_i, am, p_mass, polforce_out, qlheat_e_out, &
                                        qlheat_i_out, T_EM_phi_e_out, T_EM_phi_i_out)
        !
        ! Compute internal source terms at a single boundary point.
        !
        ! Note: T_EM_phi uses nonlinear ql fluxes, while polforce/qlheat use linear fluxes.
        !       This is handled by calling the pure helper twice with different flux values.
        !

        implicit none

        integer, intent(in) :: ipoi
        real(dp), intent(in) :: gamma_e, gamma_i
        real(dp), intent(in) :: gamma_ql_e, gamma_ql_i
        real(dp), intent(in) :: gamma_ql_e_nl, gamma_ql_i_nl
        real(dp), intent(in) :: Ercov_val
        real(dp), intent(in) :: sqrt_g_B_theta_c(:)
        real(dp), intent(in) :: Z_i, am, p_mass
        real(dp), intent(out) :: polforce_out, qlheat_e_out, qlheat_i_out
        real(dp), intent(out) :: T_EM_phi_e_out, T_EM_phi_i_out

        real(dp) :: T_EM_phi_e_dummy, T_EM_phi_i_dummy
        real(dp) :: polforce_dummy, qlheat_e_dummy, qlheat_i_dummy

        ! Get T_EM_phi from nonlinear ql fluxes
        call compute_source_terms(gamma_e, gamma_i, gamma_ql_e_nl, gamma_ql_i_nl, Ercov_val, &
                                  sqrt_g_B_theta_c(ipoi), Z_i, am, p_mass, &
                                  polforce_dummy, qlheat_e_dummy, qlheat_i_dummy, T_EM_phi_e_out, &
                                  T_EM_phi_i_out)

        ! Get polforce and qlheat from linear fluxes
        call compute_source_terms(gamma_e, gamma_i, gamma_ql_e, gamma_ql_i, Ercov_val, &
                                  sqrt_g_B_theta_c(ipoi), Z_i, am, p_mass, polforce_out, &
                                  qlheat_e_out, qlheat_i_out, T_EM_phi_e_dummy, T_EM_phi_i_dummy)
    end subroutine compute_internal_sources

    subroutine compute_dot_params_at_point(ipoi, npoi, nbaleqs, fluxes_dif_lin, fluxes_con_lin, &
                                           fluxes_con, params, params_lin, params_b_lin, Sc, &
                                           rb, rc, gpp_av, polforce, qlheat_e, qlheat_i, Z_i, &
                                           dot_params_out)
        !
        ! Compute time derivatives at a single grid point.
        !
        ! This combines:
        !   1. Flux divergence (using linearized fluxes)
        !   2. Upstream convection (using actual state fluxes for velocity)
        !   3. Internal sources
        !   4. Unit conversions (momentum -> omega, nT -> T)
        !

        implicit none

        integer, intent(in) :: ipoi, npoi, nbaleqs
        real(dp), intent(in) :: fluxes_dif_lin(:, :), fluxes_con_lin(:, :), fluxes_con(:, :)
        real(dp), intent(in) :: params(:, :), params_lin(:, :), params_b_lin(:, :)
        real(dp), intent(in) :: Sc(:), rb(:), rc(:), gpp_av(:)
        real(dp), intent(in) :: polforce(:), qlheat_e(:), qlheat_i(:)
        real(dp), intent(in) :: Z_i
        real(dp), intent(out) :: dot_params_out(4)

        integer :: ieq
        real(dp) :: convel, dr

        dr = rb(ipoi + 1) - rb(ipoi)

        ! Flux divergence for all 4 equations (using linearized fluxes):
        do ieq = 1, nbaleqs
            dot_params_out(ieq) = -(fluxes_dif_lin(ieq, ipoi + 1) - fluxes_dif_lin(ieq, ipoi)) / &
                (Sc(ipoi) * dr) - (fluxes_con_lin(ieq, ipoi + 1) - fluxes_con_lin(ieq, ipoi)) / &
                (Sc(ipoi) * dr) * params(ieq, ipoi)

            ! Upstream convection (using actual state fluxes for velocity):
            convel = 0.5_dp * (fluxes_con(ieq, ipoi + 1) + fluxes_con(ieq, ipoi)) / Sc(ipoi)
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

    pure subroutine compute_thermodynamic_forces(dn_dr, dTe_dr, dTi_dr, n, Te, Ti, E0r, Z_i, &
                                                 forces)
        implicit none

        real(dp), intent(in) :: dn_dr, dTe_dr, dTi_dr
        real(dp), intent(in) :: n, Te, Ti
        real(dp), intent(in) :: E0r, Z_i
        type(thermodynamic_forces_t), intent(out) :: forces

        ! Thermodynamic forces A1, A2:
        !   Markl2023 (5), Heyn2014 (22)
        !   A1α = (1/nα) ∂r nα  - (eα E0r)/Tα - (3/2)(1/Tα) ∂r Tα
        !   A2α = (1/Tα) ∂r Tα
        !
        ! A_noE_* denote the same expressions with the electric-field term dropped.
        ! This is the form used for the anomalous E0r=0-frame model (Heyn2014 (72), (68), (70);
        ! Markl2023 (27))
        forces%A_noE_1e = dn_dr / n - 1.5_dp * dTe_dr / Te
        forces%A_noE_2e = dTe_dr / Te
        forces%A_noE_1i = dn_dr / n - 1.5_dp * dTi_dr / Ti
        forces%A_noE_2i = dTi_dr / Ti
        forces%A_1e = forces%A_noE_1e + e_charge * E0r / Te
        forces%A_1i = forces%A_noE_1i - Z_i * e_charge * E0r / Ti
    end subroutine compute_thermodynamic_forces

    pure subroutine compute_particle_fluxes(forces, n_b, Z_i_val, D11_a_e, D12_a_e, D11_ql_e, &
                                            D12_ql_e, D11_a_i, D12_a_i, D11_ql_i, D12_ql_i, fluxes)
        !
        ! Compute particle fluxes for electrons and ions.
        !
        ! gamma_a  = -(D11_a*A_noE_1 + D12_a*A_noE_2) * n
        ! gamma_ql = -(D11_ql*A_1 + D12_ql*A_noE_2) * n
        !

        implicit none

        type(thermodynamic_forces_t), intent(in) :: forces
        real(dp), intent(in) :: n_b, Z_i_val
        real(dp), intent(in) :: D11_a_e, D12_a_e, D11_ql_e, D12_ql_e
        real(dp), intent(in) :: D11_a_i, D12_a_i, D11_ql_i, D12_ql_i
        type(transport_fluxes_t), intent(out) :: fluxes

        ! Electron particle fluxes
        fluxes%e%gamma_a = -(D11_a_e * forces%A_noE_1e + D12_a_e * forces%A_noE_2e) * n_b
        fluxes%e%gamma_ql = -(D11_ql_e * forces%A_1e + D12_ql_e * forces%A_noE_2e) * n_b
        fluxes%e%gamma = fluxes%e%gamma_a + fluxes%e%gamma_ql

        ! Ion particle fluxes
        fluxes%i%gamma_a = -(D11_a_i * forces%A_noE_1i + D12_a_i * forces%A_noE_2i) * n_b / Z_i_val
        fluxes%i%gamma_ql = -(D11_ql_i * forces%A_1i + D12_ql_i * forces%A_noE_2i) * n_b / Z_i_val
        fluxes%i%gamma = fluxes%i%gamma_a + fluxes%i%gamma_ql
    end subroutine compute_particle_fluxes

    pure subroutine compute_heat_fluxes(forces, n_b, Te_b, Ti_b, Z_i_val, D12_a_e, D21_ql_e, &
                                        D22_e, D12_a_i, D21_ql_i, D22_i, Q_e, Q_i)
        !
        ! Compute heat flux densities for electrons and ions.
        !
        ! Q_e = -(D12_a_e*A_noE_1e + D21_ql_e*A_1e + D22_e*A_noE_2e) * n * Te
        ! Q_i = -(D12_a_i*A_noE_1i + D21_ql_i*A_1i + D22_i*A_noE_2i) * n/Z * Ti
        !
        ! Note: D12 (anomalous) couples to A_noE_1 (no E-field)
        !       D21 (quasi-linear) couples to A_1 (with E-field)
        !       D22 couples to temperature gradient force A_noE_2
        !

        implicit none

        type(thermodynamic_forces_t), intent(in) :: forces
        real(dp), intent(in) :: n_b, Te_b, Ti_b, Z_i_val
        real(dp), intent(in) :: D12_a_e, D21_ql_e, D22_e
        real(dp), intent(in) :: D12_a_i, D21_ql_i, D22_i
        real(dp), intent(out) :: Q_e, Q_i

        Q_e = -(D12_a_e * forces%A_noE_1e + D21_ql_e * forces%A_1e + D22_e * forces%A_noE_2e) &
              * n_b * Te_b
        Q_i = -(D12_a_i * forces%A_noE_1i + D21_ql_i * forces%A_1i + D22_i * forces%A_noE_2i) * &
              n_b / Z_i_val * Ti_b
    end subroutine compute_heat_fluxes

    pure subroutine compute_diffusive_fluxes(ddr_n, ddr_vphi, ddr_Te, ddr_Ti, n_b, Te_b, Ti_b, &
                                             Z_i_val, Sb_val, dae11_val, dqle11_val, dae22_val, &
                                             dqle22_val, dai22_val, dni22_val, dqli22_val, &
                                             dqli21_val, visca_val, gpp_av_val, flux_dif)
        !
        ! Compute diffusive fluxes for all 4 balance equations.
        ! Diffusive flux depends only on gradients and transport coefficients.
        !

        implicit none

        real(dp), intent(in) :: ddr_n, ddr_vphi, ddr_Te, ddr_Ti
        real(dp), intent(in) :: n_b, Te_b, Ti_b, Z_i_val, Sb_val
        real(dp), intent(in) :: dae11_val, dqle11_val, dae22_val, dqle22_val
        real(dp), intent(in) :: dai22_val, dni22_val, dqli22_val, dqli21_val
        real(dp), intent(in) :: visca_val, gpp_av_val
        real(dp), intent(out) :: flux_dif(4)

        real(dp) :: dfluxvphi

        ! Eq 1: Particle flux (diffusive part)
        flux_dif(1) = -Sb_val * ddr_n * (dae11_val + dqle11_val * (1.0_dp + Ti_b / Te_b / Z_i_val))

        ! Eq 2: Momentum flux (viscous)
        dfluxvphi = -visca_val * ddr_vphi * n_b / Z_i_val * gpp_av_val
        flux_dif(2) = Sb_val * dfluxvphi

        ! Eq 3: Electron heat flux (diffusive part)
        flux_dif(3) = -Sb_val * (dae22_val + dqle22_val) * n_b * ddr_Te

        ! Eq 4: Ion heat flux (diffusive part)
        flux_dif(4) = -Sb_val * (dai22_val + dni22_val + dqli22_val - 2.5_dp * dqli21_val) * n_b / &
                      Z_i_val * ddr_Ti
    end subroutine compute_diffusive_fluxes

    pure subroutine compute_convective_fluxes(gamma_e, Q_e, Q_i, n_b, Te_b, Ti_b, Sb_val, &
                                              flux_dif, flux_con)
        !
        ! Compute convective fluxes for all 4 balance equations.
        ! Convective flux = (total flux - diffusive flux) / parameter value
        !

        implicit none

        real(dp), intent(in) :: gamma_e, Q_e, Q_i
        real(dp), intent(in) :: n_b, Te_b, Ti_b, Sb_val
        real(dp), intent(in) :: flux_dif(4)
        real(dp), intent(out) :: flux_con(4)

        ! Eq 1: Particle convection
        flux_con(1) = (Sb_val * gamma_e - flux_dif(1)) / n_b

        ! Eq 2: Momentum convection (zero - only viscous diffusion)
        flux_con(2) = 0.0_dp

        ! Eq 3: Electron heat convection
        flux_con(3) = (Sb_val * Q_e - flux_dif(3)) / Te_b

        ! Eq 4: Ion heat convection
        flux_con(4) = (Sb_val * Q_i - flux_dif(4)) / Ti_b
    end subroutine compute_convective_fluxes

    pure subroutine compute_source_terms(gamma_e, gamma_i, gamma_ql_e, gamma_ql_i, Ercov_val, &
                                         sqrt_g_B_theta_c, Z_i_val, am_val, p_mass_val, &
                                         polforce_out, qlheat_e_out, qlheat_i_out, &
                                         T_EM_phi_e_out, T_EM_phi_i_out)
        implicit none

        real(dp), intent(in) :: gamma_e, gamma_i, gamma_ql_e, gamma_ql_i
        real(dp), intent(in) :: Ercov_val, sqrt_g_B_theta_c
        real(dp), intent(in) :: Z_i_val, am_val, p_mass_val
        real(dp), intent(out) :: polforce_out, qlheat_e_out, qlheat_i_out
        real(dp), intent(out) :: T_EM_phi_e_out, T_EM_phi_i_out

        T_EM_phi_e_out = +e_charge * sqrt_g_B_theta_c * gamma_ql_e
        T_EM_phi_i_out = -Z_i_val * e_charge * sqrt_g_B_theta_c * gamma_ql_i
        polforce_out = (gamma_e - Z_i_val * gamma_i) * e_charge * sqrt_g_B_theta_c / &
                       (am_val * p_mass_val)
        qlheat_e_out = -Ercov_val * gamma_ql_e * e_charge
        qlheat_i_out = Z_i_val * Ercov_val * gamma_ql_i * e_charge
    end subroutine compute_source_terms

    pure subroutine apply_boundary_conditions(fluxes_dif_lin, fluxes_con_lin, fluxes_con, nbaleqs)
        implicit none

        integer, intent(in) :: nbaleqs
        real(dp), intent(inout) :: fluxes_dif_lin(nbaleqs, *)
        real(dp), intent(inout) :: fluxes_con_lin(nbaleqs, *)
        real(dp), intent(inout) :: fluxes_con(nbaleqs, *)

        fluxes_dif_lin(:, 1) = 0.0_dp
        fluxes_con_lin(:, 1) = 0.0_dp
        fluxes_con(:, 1) = 0.0_dp
    end subroutine apply_boundary_conditions

    pure subroutine compute_total_fluxes_at_point(gamma_e_lin, gamma_e, gamma_i_lin, Q_e_lin, &
                                                  Q_e, Q_i_lin, Q_i, ddr_n_lin, ddr_n, ddr_Te_lin, &
                                                  ddr_Te, ddr_Ti_lin, ddr_Ti, ddr_vphi_lin, n_b, &
                                                  Te_b, Ti_b, Z_i_val, Sb_val, dae11_val, &
                                                  dqle11_val, dae22_val, dqle22_val, dai22_val, &
                                                  dni22_val, dqli22_val, dqli21_val, visca_val, &
                                                  gpp_av_val, flux_dif_lin, flux_con_lin, flux_con)
        implicit none

        real(dp), intent(in) :: gamma_e_lin, gamma_e, gamma_i_lin
        real(dp), intent(in) :: Q_e_lin, Q_e, Q_i_lin, Q_i
        real(dp), intent(in) :: ddr_n_lin, ddr_n, ddr_Te_lin, ddr_Te, ddr_Ti_lin, ddr_Ti
        real(dp), intent(in) :: ddr_vphi_lin
        real(dp), intent(in) :: n_b, Te_b, Ti_b, Z_i_val, Sb_val
        real(dp), intent(in) :: dae11_val, dqle11_val, dae22_val, dqle22_val
        real(dp), intent(in) :: dai22_val, dni22_val, dqli22_val, dqli21_val
        real(dp), intent(in) :: visca_val, gpp_av_val
        real(dp), intent(out) :: flux_dif_lin(4), flux_con_lin(4), flux_con(4)
        real(dp) :: dfluxvphi

        flux_dif_lin(1) = -Sb_val * ddr_n_lin * (dae11_val + dqle11_val * (1.0_dp + Ti_b / Te_b / &
                          Z_i_val))
        flux_con_lin(1) = (Sb_val * gamma_e_lin - flux_dif_lin(1)) / n_b
        flux_con(1) = (Sb_val * gamma_e - (-Sb_val * ddr_n * &
                       (dae11_val + dqle11_val * (1.0_dp + Ti_b / Te_b / Z_i_val)))) / n_b

        dfluxvphi = -visca_val * ddr_vphi_lin * n_b / Z_i_val * gpp_av_val
        flux_dif_lin(2) = Sb_val * dfluxvphi
        flux_con_lin(2) = 0.0_dp
        flux_con(2) = 0.0_dp

        flux_dif_lin(3) = -Sb_val * (dae22_val + dqle22_val) * n_b * ddr_Te_lin
        flux_con_lin(3) = (Sb_val * Q_e_lin - flux_dif_lin(3)) / Te_b
        flux_con(3) = (Sb_val * Q_e - (-Sb_val * (dae22_val + dqle22_val) * n_b * ddr_Te)) / Te_b

        flux_dif_lin(4) = -Sb_val * (dai22_val + dni22_val + dqli22_val - 2.5_dp * dqli21_val) * &
                          n_b / Z_i_val * ddr_Ti_lin
        flux_con_lin(4) = (Sb_val * Q_i_lin - flux_dif_lin(4)) / Ti_b
        flux_con(4) = (Sb_val * Q_i - (-Sb_val * (dai22_val + dni22_val + dqli22_val - 2.5_dp * &
                       dqli21_val) * n_b / Z_i_val * ddr_Ti)) / Ti_b
    end subroutine compute_total_fluxes_at_point

    pure subroutine compute_time_derivatives(dot_params_in, n_val, Te_val, Ti_val, Z_i_val, &
                                             gpp_av_avg, dot_params_out)
        implicit none

        real(dp), intent(in) :: dot_params_in(4)
        real(dp), intent(in) :: n_val, Te_val, Ti_val, Z_i_val, gpp_av_avg
        real(dp), intent(out) :: dot_params_out(4)

        dot_params_out(1) = dot_params_in(1)
        dot_params_out(2) = dot_params_in(2) * Z_i_val / n_val * 2.0_dp / gpp_av_avg
        dot_params_out(3) = (-Te_val * dot_params_in(1) + dot_params_in(3) / 1.5_dp) / n_val
        dot_params_out(4) = (-Ti_val * dot_params_in(1) + dot_params_in(4) / 1.5_dp) / n_val
    end subroutine compute_time_derivatives

end module rhs_balance_m
