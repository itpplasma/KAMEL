subroutine rhs_balance_stell(x, y, dy)

    ! version of rhs_balance that includes the 1/nu transport used to study stellarators
    ! Mostly, it is a copy of rhs_balance

    use grid_mod, only: nbaleqs, neqset, iboutype, npoic, npoib &
                        , Sc, Sb, deriv_coef &
                        , ipbeg, ipend, rb, reint_coef &
                        , fluxes_dif_lin, fluxes_con_lin, rc &
                        , dae11, dae12, dae22, dai11, dai12, dai22 &
                        , dni22, visca, gpp_av &
                        , dqle11, dqle12, dqle21, dqle22 &
                        , dqli11, dqli12, dqli21, dqli22 &
                        , sqrt_g_times_B_theta_over_c, Ercov, polforce, qlheat_e, qlheat_i &
                        , Ercov_lin, fluxes_con, Donue11, Donue12, Donue21, Donue22 &
                        , Donui11, Donui12, Donui21, Donui22, cneo

    use plasma_parameters, only: params, ddr_params_lin, params_lin, ddr_params &
                        , params_b_lin, params_b, dot_params
    use baseparam_mod, only: Z_i, e_charge, am, p_mass, c
    use wave_code_data, only: q, Vth
    use matrix_mod, only: isw_rhs, nz, nsize, irow, icol, amat, rhsvec
    use KAMEL_hdf5_tools
    use QLBalance_kinds, only: dp
    use time_evolution_stellarator, only: set_Q_neo_to_zero, turn_off_heat_sources

    implicit none

    integer :: ipoi, ieq, i, npoi, ibeg, iend, nshift, ibegb, iendb, ibegtot, iendtot, k, iprobe
    real(dp) :: x, A_noE_1e_lin, A_noE_2e_lin, A_noE_1i_lin, A_noE_2i_lin, convel
    real(dp) :: A_noE_1e, A_noE_2e, A_noE_1i, A_noE_2i
    real(dp) :: gamma_e_lin, gamma_i_lin, dfluxvphi, Q_e_lin, Q_i_lin, A_1e_lin, A_1i_lin
    real(dp) :: gamma_e, Q_e, Q_i, A_1e, A_1i
    real(dp) :: gamma_ql_e_lin, gamma_ql_i_lin, gamma_onu_e_lin, gamma_onu_i_lin
    real(dp), dimension(neqset) :: y, dy, y_lin
    real(dp) :: De11, De12, De21, De22, Di11, Di12, Di21, Di22

    if (iboutype .eq. 1) then
        npoi = npoic - 1
    else
        npoi = npoic
    end if

    y_lin = 0.0d0
    params_lin = 0.0d0
    ddr_params_lin = 0.0d0
    Ercov_lin = 0.0d0

    do ipoi = 1, npoi
        do ieq = 1, nbaleqs
            i = nbaleqs*(ipoi - 1) + ieq
            params(ieq, ipoi) = y(i)
        end do
    end do

    do ipoi = 1, npoib
        do ieq = 1, nbaleqs
            ! radial derivatives of equilibrium parameters at cell boundaries:
            ddr_params(ieq, ipoi) &
                = sum(params(ieq, ipbeg(ipoi):ipend(ipoi))*deriv_coef(:, ipoi))
            ! equilibrium parameters at cell boundaries:
            params_b(ieq, ipoi) &
                = sum(params(ieq, ipbeg(ipoi):ipend(ipoi))*reint_coef(:, ipoi))
        end do
    end do


    Ercov = sqrt_g_times_B_theta_over_c*(params_b(2, :) - Vth*q/rb) &
            + (params_b(4, :)*ddr_params(1, :)/params_b(1, :) + ddr_params(4, :)) &
            /(Z_i*e_charge)

    if (set_Q_neo_to_zero) cneo = 0.0d0 ! sets Dni22 to zero
    call calc_equil_diffusion_coeffs

    do ipoi = 1, npoib
        De11 = dqle11(ipoi) + Donue11(ipoi)
        De12 = dqle12(ipoi) + Donue12(ipoi)
        De21 = dqle21(ipoi) + Donue21(ipoi)
        De22 = dqle22(ipoi) + Donue22(ipoi)

        Di11 = dqli11(ipoi) + Donui11(ipoi)
        Di12 = dqli12(ipoi) + Donui12(ipoi)
        Di21 = dqli21(ipoi) + Donui21(ipoi)
        Di22 = dqli22(ipoi) + Donui22(ipoi)

        ! Thermodynamic forces for zero radial electric field (actual state):
        A_noE_1e = ddr_params(1, ipoi)/params_b(1, ipoi) &
                    - 1.5d0*ddr_params(3, ipoi)/params_b(3, ipoi)
        A_noE_2e = ddr_params(3, ipoi)/params_b(3, ipoi)
        A_noE_1i = ddr_params(1, ipoi)/params_b(1, ipoi) &
                    - 1.5d0*ddr_params(4, ipoi)/params_b(4, ipoi)
        A_noE_2i = ddr_params(4, ipoi)/params_b(4, ipoi)

        ! Thermodynamic forces for finite radial electric field (actual state):
        A_1e = A_noE_1e + Ercov(ipoi)*e_charge/params_b(3, ipoi)
        A_1i = A_noE_1i - Ercov(ipoi)*e_charge*Z_i/params_b(4, ipoi)

        ! particle flux densities (actual state):
        gamma_e = -(dae11(ipoi)*A_noE_1e + dae12(ipoi)*A_noE_2e &
                    + De11*A_1e + De12*A_noE_2e &
                    )*params_b(1, ipoi)

        ! total particle flux (actual state):
        fluxes_con(1, ipoi) = (Sb(ipoi)*gamma_e - &
                               (-Sb(ipoi)*ddr_params(1, ipoi)*(dae11(ipoi) &
                            + De11*(1.d0 + params_b(4, ipoi)/params_b(3, ipoi)/Z_i))))&
                            /params_b(1, ipoi)

        ! toroidal moment flux (actual state):
        fluxes_con(2, ipoi) = 0.d0

        ! electron heat flux density (actual state):
        Q_e = -(dae12(ipoi)*A_noE_1e + De21*A_1e &
                + (dae22(ipoi) + De22)*A_noE_2e) &
               *params_b(1, ipoi)*params_b(3, ipoi)

        ! ion heat flux density (actual state):
        Q_i = -(dai12(ipoi)*A_noE_1i + Di21*A_1i &
                + (dai22(ipoi) + dni22(ipoi) + Di22)*A_noE_2i) &
               *params_b(1, ipoi)/Z_i*params_b(4, ipoi)

        ! total heat fluxes (actual state):
        fluxes_con(3, ipoi) = (Sb(ipoi)*Q_e - &
                              (-Sb(ipoi)*(dae22(ipoi) + De22)*params_b(1, ipoi) &
                              * ddr_params(3, ipoi))) /params_b(3, ipoi)
        fluxes_con(4, ipoi) = (Sb(ipoi)*Q_i - &
                              (-Sb(ipoi)*(dai22(ipoi) + dni22(ipoi) + Di22 - 2.5d0*Di21) &
                                 *params_b(1, ipoi)/Z_i*ddr_params(4, ipoi))) &
                              /params_b(4, ipoi)
    end do

    nshift = 4
    k = 0
    dy = 0.d0

    do iprobe = 1, neqset
        ibeg = iprobe/nbaleqs - nshift
        iend = iprobe/nbaleqs + nshift
        ibegb = ibeg - 1
        iendb = iend + 1
        ibegtot = iprobe - nshift*nbaleqs
        iendtot = iprobe + nshift*nbaleqs
        ibeg = max(1, ibeg)
        iend = min(npoi, iend)
        ibegb = max(1, ibegb)
        iendb = min(npoib, iendb)
        ibegtot = max(1, ibegtot)
        iendtot = min(neqset, iendtot)

        y_lin(iprobe) = 1.d0

        do ipoi = ibeg, iend
            do ieq = 1, nbaleqs
                i = nbaleqs*(ipoi - 1) + ieq
                params_lin(ieq, ipoi) = y_lin(i)
            end do
        end do

        ! Compute fluxes and internal sources:

        do ipoi = ibegb, iendb
            do ieq = 1, nbaleqs
                ! radial derivatives of linearized parameters at cell boundaries:
                ddr_params_lin(ieq, ipoi) &
                    = sum(params_lin(ieq, ipbeg(ipoi):ipend(ipoi))*deriv_coef(:, ipoi))
                params_b_lin(ieq, ipoi) &
                    = sum(params_lin(ieq, ipbeg(ipoi):ipend(ipoi))*reint_coef(:, ipoi))
            end do
        end do

        Ercov_lin(ibegb:iendb) &
            = sqrt_g_times_B_theta_over_c(ibegb:iendb)*params_b_lin(2, ibegb:iendb) &
                + (params_b(4, ibegb:iendb)*ddr_params_lin(1, ibegb:iendb) &
                /params_b(1, ibegb:iendb) + ddr_params_lin(4, ibegb:iendb)) &
                /(Z_i*e_charge)


        do ipoi = ibegb, iendb
            De11 = dqle11(ipoi) + Donue11(ipoi)
            De12 = dqle12(ipoi) + Donue12(ipoi)
            De21 = dqle21(ipoi) + Donue21(ipoi)
            De22 = dqle22(ipoi) + Donue22(ipoi)

            Di11 = dqli11(ipoi) + Donui11(ipoi)
            Di12 = dqli12(ipoi) + Donui12(ipoi)
            Di21 = dqli21(ipoi) + Donui21(ipoi)
            Di22 = dqli22(ipoi) + Donui22(ipoi)

            ! Thermodynamic forces for zero radial electric field (linearized):
            A_noE_1e_lin = ddr_params_lin(1, ipoi)/params_b(1, ipoi) &
                        - 1.5d0*ddr_params_lin(3, ipoi)/params_b(3, ipoi)
            A_noE_2e_lin = ddr_params_lin(3, ipoi)/params_b(3, ipoi)
            A_noE_1i_lin = ddr_params_lin(1, ipoi)/params_b(1, ipoi) &
                        - 1.5d0*ddr_params_lin(4, ipoi)/params_b(4, ipoi)
            A_noE_2i_lin = ddr_params_lin(4, ipoi)/params_b(4, ipoi)

            ! Thermodynamic forces for finite radial electric field (linearized):
            A_1e_lin = A_noE_1e_lin + Ercov_lin(ipoi)*e_charge/params_b(3, ipoi)
            A_1i_lin = A_noE_1i_lin - Ercov_lin(ipoi)*e_charge*Z_i/params_b(4, ipoi)

            ! particle flux densities (linearized):
            gamma_e_lin = -(dae11(ipoi)*A_noE_1e_lin + dae12(ipoi)*A_noE_2e_lin)*params_b(1, ipoi)
            gamma_ql_e_lin = -(dqle11(ipoi)*A_1e_lin + dqle12(ipoi)*A_noE_2e_lin)*params_b(1, ipoi)
            gamma_onu_e_lin = -(Donue11(ipoi)*A_1e_lin + Donue12(ipoi)*A_noE_2e_lin)*params_b(1, ipoi)
            gamma_e_lin = gamma_e_lin + gamma_ql_e_lin + gamma_onu_e_lin
            gamma_i_lin = -(dai11(ipoi)*A_noE_1i_lin + dai12(ipoi)*A_noE_2i_lin)*params_b(1, ipoi)/Z_i
            gamma_ql_i_lin = -(dqli11(ipoi)*A_1i_lin + dqli12(ipoi)*A_noE_2i_lin)*params_b(1, ipoi)/Z_i
            gamma_onu_i_lin = -(Donui11(ipoi)*A_1i_lin + Donui12(ipoi)*A_noE_2i_lin)*params_b(1, ipoi) / Z_i
            gamma_i_lin = gamma_i_lin + gamma_ql_i_lin + gamma_onu_i_lin

            ! total particle flux (linearized):
            fluxes_dif_lin(1, ipoi) = -Sb(ipoi)*ddr_params_lin(1, ipoi)*(dae11(ipoi) &
                                    + De11*(1.d0 + params_b(4, ipoi)/params_b(3, ipoi)/Z_i))
            fluxes_con_lin(1, ipoi) = (Sb(ipoi)*gamma_e_lin - fluxes_dif_lin(1, ipoi))/params_b(1, ipoi)

            ! toroidal moment flux density divided by mass (linearized):
            dfluxvphi = -visca(ipoi)*ddr_params_lin(2, ipoi)*params_b(1, ipoi)/Z_i*gpp_av(ipoi)

            ! total toroidal moment flux (linearized):
            fluxes_dif_lin(2, ipoi) = Sb(ipoi)*dfluxvphi
            fluxes_con_lin(2, ipoi) = 0.d0

            ! electron heat flux density (linearized):
            Q_e_lin = -(dae12(ipoi)*A_noE_1e_lin + De21*A_1e_lin &
                    + (dae22(ipoi) + De22)*A_noE_2e_lin) &
                  *params_b(1, ipoi)*params_b(3, ipoi)

            ! ion heat flux density (linearized):
            Q_i_lin = -(dai12(ipoi)*A_noE_1i_lin + Di21*A_1i_lin &
                    + (dai22(ipoi) + dni22(ipoi) + Di22)*A_noE_2i_lin) &
                  *params_b(1, ipoi)/Z_i*params_b(4, ipoi)

            ! total heat fluxes (linearized):
            fluxes_dif_lin(3, ipoi) = -Sb(ipoi)*(dae22(ipoi) + De22) &
                                  *params_b(1, ipoi)*ddr_params_lin(3, ipoi)
            fluxes_con_lin(3, ipoi) = (Sb(ipoi)*Q_e_lin - fluxes_dif_lin(3, ipoi))/params_b(3, ipoi)

            fluxes_dif_lin(4, ipoi) = -Sb(ipoi)*(dai22(ipoi) + dni22(ipoi) + Di22 &
                                - 2.5d0*Di21)*params_b(1, ipoi)/Z_i*ddr_params_lin(4, ipoi)
            fluxes_con_lin(4, ipoi) = (Sb(ipoi)*Q_i_lin - fluxes_dif_lin(4, ipoi))/params_b(4, ipoi)

            ! Momentum source due to the polarization current:
            polforce(ipoi) = (gamma_e_lin - Z_i*gamma_i_lin)*e_charge*sqrt_g_times_B_theta_over_c(ipoi) &
                            /(am*p_mass)

            ! Heat sources due to the radial QL drift in the equilibrium electric field:
            qlheat_e(ipoi) = -Ercov(ipoi)*(gamma_ql_e_lin + gamma_onu_e_lin)*e_charge
            qlheat_i(ipoi) = Z_i*Ercov(ipoi)*(gamma_ql_i_lin + gamma_onu_i_lin)*e_charge
        end do

        if (turn_off_heat_sources) then
            qlheat_e = 0.d0
            qlheat_i = 0.d0
        end if

        ! Condition of zero flux at the inner boundary:
        fluxes_dif_lin(:, 1) = 0.d0
        fluxes_con_lin(:, 1) = 0.d0
        fluxes_con(:, 1) = 0.d0

        ! Partial time derivatives of equilibrium parameters:
        do ipoi = ibeg, iend
        ! Flux divergence:
            do ieq = 1, nbaleqs
                dot_params(ieq, ipoi) = -(fluxes_dif_lin(ieq, ipoi + 1) - fluxes_dif_lin(ieq, ipoi)) &
                                        /(Sc(ipoi)*(rb(ipoi + 1) - rb(ipoi))) &
                                        - (fluxes_con_lin(ieq, ipoi + 1) - fluxes_con_lin(ieq, ipoi)) &
                                        /(Sc(ipoi)*(rb(ipoi + 1) - rb(ipoi)))*params(ieq, ipoi)
                convel = 0.5d0*(fluxes_con(ieq, ipoi + 1) + fluxes_con(ieq, ipoi))/Sc(ipoi)
        ! upstream convection:
                if (convel .gt. 0.d0) then
                !if(convel.lt.0.d0) then
                    dot_params(ieq, ipoi) = dot_params(ieq, ipoi) &
                                            - convel*(params_lin(ieq, ipoi + 1) - &
                                            params_lin(ieq, ipoi))/(rc(ipoi + 1) - rc(ipoi))
                else
                    if (ipoi .gt. 1) then
                        dot_params(ieq, ipoi) = dot_params(ieq, ipoi) &
                                            - convel*(params_lin(ieq, ipoi - 1) - &
                                            params_lin(ieq, ipoi))/(rc(ipoi - 1) - rc(ipoi))
                    else
                        dot_params(ieq, ipoi) = dot_params(ieq, ipoi) &
                                                - convel*(params_b_lin(ieq, 1) - params_lin(ieq, 1))/(rb(1) - rc(1))
                    end if
                end if
            end do

            ! Add internal sources:
            ! Momentum:
            dot_params(2, ipoi) = dot_params(2, ipoi) &
                                + 0.5d0*(polforce(ipoi) + polforce(ipoi + 1))

            ! Heat into electrons:
            dot_params(3, ipoi) = dot_params(3, ipoi) &
                                + 0.5d0*(qlheat_e(ipoi) + qlheat_e(ipoi + 1))

            ! Heat into ions:
            dot_params(4, ipoi) = dot_params(4, ipoi) &
                                + 0.5d0*(qlheat_i(ipoi) + qlheat_i(ipoi + 1))

            ! Covert momentum time derivative to time derivative of the rotation frequency:
            dot_params(2, ipoi) = dot_params(2, ipoi)*Z_i/params(1, ipoi) &
                                  *2.d0/(gpp_av(ipoi + 1) + gpp_av(ipoi))

            ! Convert dot_params from d(nT_{e,i})/dt to d(T_{e,i})/dt:
            dot_params(3, ipoi) = (-params(3, ipoi)*dot_params(1, ipoi) &
                                + dot_params(3, ipoi)/1.5d0)/params(1, ipoi)

            dot_params(4, ipoi) = (-params(4, ipoi)*dot_params(1, ipoi) &
                                + dot_params(4, ipoi)/1.5d0)/params(1, ipoi)

        end do

        ! RHS vector of ODE system
        do ipoi = ibeg, iend
            do ieq = 1, nbaleqs
                i = nbaleqs*(ipoi - 1) + ieq
                dy(i) = dot_params(ieq, ipoi)
            end do
        end do

        do i = ibegtot, iendtot
            if (dy(i) .ne. 0.d0) then
                k = k + 1
                if (isw_rhs .eq. 1) then
                    irow(k) = i
                    icol(k) = iprobe
                    amat(k) = dy(i)
                end if
            end if
        end do

        ! Clean:

        y_lin(iprobe) = 0.d0
        do ipoi = ibeg, iend
            do ieq = 1, nbaleqs
                i = nbaleqs*(ipoi - 1) + ieq
                params_lin(ieq, ipoi) = y_lin(i)
            end do
        end do
        do ipoi = ibegb, iendb
            do ieq = 1, nbaleqs
                ddr_params_lin(ieq, ipoi) = 0.d0
            end do
        end do
        Ercov_lin(ibegb:iendb) = 0.d0
        fluxes_dif_lin(:, ibegb:iendb) = 0.d0
        fluxes_con_lin(:, ibegb:iendb) = 0.d0
        polforce(ibegb:iendb) = 0.d0
        qlheat_e(ibegb:iendb) = 0.d0
        qlheat_i(ibegb:iendb) = 0.d0
        dot_params(:, ibeg:iend) = 0.d0
        do ipoi = ibeg, iend
            do ieq = 1, nbaleqs
                i = nbaleqs*(ipoi - 1) + ieq
                dy(i) = dot_params(ieq, ipoi)
            end do
        end do
    end do

    if (isw_rhs .eq. 0) then
        !  print *,'Number of non-zero elements = ',k,' out of: ',neqset**2
        nz = k
        nsize = neqset
    else
        call rhs_balance_source_stell(x, y, dy)
        rhsvec = dy
    end if

end subroutine rhs_balance_stell


subroutine initialize_rhs_stell(y, dy)

    use grid_mod, only: neqset
    use matrix_mod
    use QLBalance_kinds, only: dp

    implicit none
    real(dp) :: x
    real(dp), dimension(neqset) :: y, dy

    x = 0.d0
    isw_rhs = 0

    call rhs_balance_stell(x, y, dy)

    isw_rhs = 1
    if (allocated(amat)) deallocate (irow, icol, amat, rhsvec)
    allocate (irow(nz), icol(nz), amat(nz), rhsvec(nsize))

end subroutine initialize_rhs_stell


subroutine rhs_balance_source_stell(x, y, dy)

    use grid_mod, only: nbaleqs, neqset, iboutype, npoic, npoib &
                        , Sc, Sb, deriv_coef &
                        , ipbeg, ipend, rb, reint_coef &
                        , fluxes_dif_lin, fluxes_con_lin, rc &
                        , dae11, dae12, dae22, dai11, dai12, dai22 &
                        , dni22, visca, gpp_av, dery_equisource &
                        , dqle11, dqle12, dqle21, dqle22 &
                        , dqli11, dqli12, dqli21, dqli22 &
                        , sqrt_g_times_B_theta_over_c, Ercov, polforce, qlheat_e, qlheat_i &
                        , Ercov_lin, fluxes_con, cneo, Donue11, Donue12, Donue21, Donue22 &
                        , Donui11, Donui12, Donui21, Donui22

    use plasma_parameters, only: params, ddr_params_lin, params_b, params_lin &
                        , params_b_lin, ddr_params, dot_params
    use baseparam_mod, only: Z_i, e_charge, am, p_mass, c
    use wave_code_data, only: q, Vth
    use time_evolution_stellarator, only: set_Q_neo_to_zero, turn_off_heat_sources
    use QLBalance_kinds, only: dp

    implicit none

    integer :: ipoi, ieq, i, npoi
    real(dp) :: x, A_noE_1e_lin, A_noE_2e_lin, A_noE_1i_lin, A_noE_2i_lin, convel
    real(dp) :: A_noE_1e, A_noE_2e, A_noE_1i, A_noE_2i
    real(dp) :: gamma_e_lin, gamma_i_lin, dfluxvphi, Q_e_lin, Q_i_lin, A_1e_lin, A_1i_lin
    real(dp) :: gamma_e, Q_e, Q_i, A_1e, A_1i
    real(dp) :: gamma_ql_e_lin, gamma_ql_i_lin
    real(dp), dimension(neqset) :: y, dy, y_lin
    real(dp) :: De11, De12, De21, De22, Di11, Di12, Di21, Di22

    if (iboutype .eq. 1) then
        npoi = npoic - 1
    else
        npoi = npoic
    end if
    
    y_lin = 0.0d0
    params_lin = params

    do ipoi = 1, npoi
        do ieq = 1, nbaleqs
            i = nbaleqs*(ipoi - 1) + ieq
            params(ieq, ipoi) = y(i)
            params_lin(ieq, ipoi) = y_lin(i)
        end do
    end do

    !
    ! Interpolation:
    !
    do ipoi = 1, npoib
        do ieq = 1, nbaleqs
            ! radial derivatives of equilibrium parameters at cell boundaries:
            ddr_params(ieq, ipoi) &
                = sum(params(ieq, ipbeg(ipoi):ipend(ipoi))*deriv_coef(:, ipoi))
            ddr_params_lin(ieq, ipoi) &
                = sum(params_lin(ieq, ipbeg(ipoi):ipend(ipoi))*deriv_coef(:, ipoi))
            ! equilibrium parameters at cell boundaries:
            params_b(ieq, ipoi) &
                = sum(params(ieq, ipbeg(ipoi):ipend(ipoi))*reint_coef(:, ipoi))
            params_b_lin(ieq, ipoi) &
                = sum(params_lin(ieq, ipbeg(ipoi):ipend(ipoi))*reint_coef(:, ipoi))
        end do
    end do

    !
    ! Compute radial electric field:
    !
    Ercov = sqrt_g_times_B_theta_over_c*(params_b(2, :) - Vth*q/rb) &
            + (params_b(4, :)*ddr_params(1, :)/params_b(1, :) + ddr_params(4, :)) &
            /(Z_i*e_charge)

    Ercov_lin = sqrt_g_times_B_theta_over_c*params_b_lin(2, :) &
                + (params_b(4, :)*ddr_params_lin(1, :)/params_b(1, :) + ddr_params_lin(4, :)) &
                /(Z_i*e_charge)

    ! Compute diffusion coefficient matrices:
    if (set_Q_neo_to_zero) cneo = 0.0d0
    call calc_equil_diffusion_coeffs
    !
    ! Compute fluxes and internal sources:
    !

    do ipoi = 1, npoib
        De11 = dqle11(ipoi) + Donue11(ipoi)
        De12 = dqle12(ipoi) + Donue12(ipoi)
        De21 = dqle21(ipoi) + Donue21(ipoi)
        De22 = dqle22(ipoi) + Donue22(ipoi)

        Di11 = dqli11(ipoi) + Donui11(ipoi)
        Di12 = dqli12(ipoi) + Donui12(ipoi)
        Di21 = dqli21(ipoi) + Donui21(ipoi)
        Di22 = dqli22(ipoi) + Donui22(ipoi)

        ! Thermodynamic forces for zero radial electric field (linearized):
        A_noE_1e_lin = ddr_params_lin(1, ipoi)/params_b(1, ipoi) &
                    - 1.5d0*ddr_params_lin(3, ipoi)/params_b(3, ipoi)
        A_noE_2e_lin = ddr_params_lin(3, ipoi)/params_b(3, ipoi)
        A_noE_1i_lin = ddr_params_lin(1, ipoi)/params_b(1, ipoi) &
                    - 1.5d0*ddr_params_lin(4, ipoi)/params_b(4, ipoi)
        A_noE_2i_lin = ddr_params_lin(4, ipoi)/params_b(4, ipoi)

        ! Thermodynamic forces for zero radial electric field (actual state):
        A_noE_1e = ddr_params(1, ipoi)/params_b(1, ipoi) &
                    - 1.5d0*ddr_params(3, ipoi)/params_b(3, ipoi)
        A_noE_2e = ddr_params(3, ipoi)/params_b(3, ipoi)
        A_noE_1i = ddr_params(1, ipoi)/params_b(1, ipoi) &
                    - 1.5d0*ddr_params(4, ipoi)/params_b(4, ipoi)
        A_noE_2i = ddr_params(4, ipoi)/params_b(4, ipoi)

        ! Thermodynamic forces for finite radial electric field (linearized):
        A_1e_lin = A_noE_1e_lin + Ercov_lin(ipoi)*e_charge/params_b(3, ipoi)
        A_1i_lin = A_noE_1i_lin - Ercov_lin(ipoi)*e_charge*Z_i/params_b(4, ipoi)

        ! Thermodynamic forces for finite radial electric field (actual state):
        A_1e = A_noE_1e + Ercov(ipoi)*e_charge/params_b(3, ipoi)
        A_1i = A_noE_1i - Ercov(ipoi)*e_charge*Z_i/params_b(4, ipoi)


        ! particle flux densities (linearized):
        gamma_e_lin = -(dae11(ipoi)*A_noE_1e_lin + dae12(ipoi)*A_noE_2e_lin)*params_b(1, ipoi)
        gamma_ql_e_lin = -(De11*A_1e_lin + De12*A_noE_2e_lin)*params_b(1, ipoi)
        gamma_e_lin = gamma_e_lin + gamma_ql_e_lin
        gamma_i_lin = -(dai11(ipoi)*A_noE_1i_lin + dai12(ipoi)*A_noE_2i_lin)*params_b(1, ipoi)/Z_i
        gamma_ql_i_lin = -(Di11*A_1i_lin + Di12*A_noE_2i_lin)*params_b(1, ipoi)/Z_i
        gamma_i_lin = gamma_i_lin + gamma_ql_i_lin

        ! particle flux density (actual state):
        gamma_e = -(dae11(ipoi)*A_noE_1e + dae12(ipoi)*A_noE_2e &
                        + De11*A_1e + De12*A_noE_2e &
                        )*params_b(1, ipoi)

        ! total particle flux (linearized):
        fluxes_dif_lin(1, ipoi) = -Sb(ipoi)*ddr_params_lin(1, ipoi)*(dae11(ipoi) &
                            + De11*(1.d0 + params_b(4, ipoi)/params_b(3, ipoi)/Z_i))
        fluxes_con_lin(1, ipoi) = (Sb(ipoi)*gamma_e_lin - fluxes_dif_lin(1, ipoi))/params_b(1, ipoi)
        ! total particle flux (actual state):
        fluxes_con(1, ipoi) = (Sb(ipoi)*gamma_e - &
                                (-Sb(ipoi)*ddr_params(1, ipoi)*(dae11(ipoi) &
                            + De11*(1.d0 + params_b(4, ipoi)/params_b(3, ipoi)/Z_i))))/params_b(1, ipoi)

        ! toroidal moment flux density divided by mass (linearized):
        dfluxvphi = -visca(ipoi)*ddr_params_lin(2, ipoi)*params_b(1, ipoi)/Z_i*gpp_av(ipoi)
        ! total toroidal moment flux:
        fluxes_dif_lin(2, ipoi) = Sb(ipoi)*dfluxvphi
        fluxes_con_lin(2, ipoi) = 0.d0
        fluxes_con(2, ipoi) = 0.d0

        ! electron heat flux density (linearized):
        Q_e_lin = -(dae12(ipoi)*A_noE_1e_lin + De21*A_1e_lin &
                + (dae22(ipoi) + De22)*A_noE_2e_lin) &
              *params_b(1, ipoi)*params_b(3, ipoi)

        ! electron heat flux density (actual state):
        Q_e = -(dae12(ipoi)*A_noE_1e + De21*A_1e &
                + (dae22(ipoi) + De22)*A_noE_2e) &
                *params_b(1, ipoi)*params_b(3, ipoi)

        ! ion heat flux density (linearized):
        Q_i_lin = -(dai12(ipoi)*A_noE_1i_lin + Di21*A_1i_lin &
                + (dai22(ipoi) + dni22(ipoi) + Di22)*A_noE_2i_lin) &
              *params_b(1, ipoi)/Z_i*params_b(4, ipoi)

        ! ion heat flux density (actual state):
        Q_i = -(dai12(ipoi)*A_noE_1i + Di21*A_1i &
                + (dai22(ipoi) + dni22(ipoi) + Di22)*A_noE_2i) &
                 *params_b(1, ipoi)/Z_i*params_b(4, ipoi)

        ! total heat fluxes (linearized):
        fluxes_dif_lin(3, ipoi) = -Sb(ipoi)*(dae22(ipoi) + De22) &
                              *params_b(1, ipoi)*ddr_params_lin(3, ipoi)
        fluxes_con_lin(3, ipoi) = (Sb(ipoi)*Q_e_lin - fluxes_dif_lin(3, ipoi))/params_b(3, ipoi)
        ! total heat flux (actual state):
        fluxes_con(3, ipoi) = (Sb(ipoi)*Q_e - &
                            (-Sb(ipoi)*(dae22(ipoi) + De22)*params_b(1, ipoi)*ddr_params(3, ipoi))) &
                            /params_b(3, ipoi)

        fluxes_dif_lin(4, ipoi) = -Sb(ipoi)*(dai22(ipoi) + dni22(ipoi) + Di22 &
                            - 2.5d0*Di21)*params_b(1, ipoi)/Z_i*ddr_params_lin(4, ipoi)
        fluxes_con_lin(4, ipoi) = (Sb(ipoi)*Q_i_lin - fluxes_dif_lin(4, ipoi))/params_b(4, ipoi)
        fluxes_con(4, ipoi) = (Sb(ipoi)*Q_i - &
                            (-Sb(ipoi)*(dai22(ipoi) + dni22(ipoi) + Di22 - 2.5d0*Di21) &
                            *params_b(1, ipoi)/Z_i*ddr_params(4, ipoi))) &
                            /params_b(4, ipoi)

        ! Momentum source due to the polarization current:
        polforce(ipoi) = (gamma_e_lin - Z_i*gamma_i_lin)*e_charge*sqrt_g_times_B_theta_over_c(ipoi) &
                        /(am*p_mass)

        ! Heat sources due to the radial QL drift in the equilibrium electric field:
        qlheat_e(ipoi) = -Ercov(ipoi)*gamma_ql_e_lin*e_charge
        qlheat_i(ipoi) = Z_i*Ercov(ipoi)*gamma_ql_i_lin*e_charge
    end do
    ! Condition of zero flux at the inner boundary:
    fluxes_dif_lin(:, 1) = 0.d0
    fluxes_con_lin(:, 1) = 0.d0
    fluxes_con(:, 1) = 0.d0

    if (turn_off_heat_sources) then
        qlheat_e = 0.d0
        qlheat_i = 0.d0
    end if

    ! Partial time derivatives of equilibrium parameters:
    do ipoi = 1, npoi
        ! Flux divergence:
        do ieq = 1, nbaleqs
            dot_params(ieq, ipoi) = -(fluxes_dif_lin(ieq, ipoi + 1) - fluxes_dif_lin(ieq, ipoi)) &
                                    /(Sc(ipoi)*(rb(ipoi + 1) - rb(ipoi))) &
                                    - (fluxes_con_lin(ieq, ipoi + 1) - fluxes_con_lin(ieq, ipoi)) &
                                    /(Sc(ipoi)*(rb(ipoi + 1) - rb(ipoi)))*params(ieq, ipoi)
            convel = 0.5d0*(fluxes_con(ieq, ipoi + 1) + fluxes_con(ieq, ipoi))/Sc(ipoi)
            ! upstream convection:
            if (convel .gt. 0.d0) then
                dot_params(ieq, ipoi) = dot_params(ieq, ipoi) &
                                        - convel*(params_lin(ieq, ipoi + 1) - &
                                        params_lin(ieq, ipoi))/(rc(ipoi + 1) - rc(ipoi))
            else
                if (ipoi .gt. 1) then
                    dot_params(ieq, ipoi) = dot_params(ieq, ipoi) &
                                            - convel*(params_lin(ieq, ipoi - 1) - &
                                            params_lin(ieq, ipoi))/(rc(ipoi - 1) - rc(ipoi))
                else
                    dot_params(ieq, ipoi) = dot_params(ieq, ipoi)                    !    &
                end if
            end if
        end do


        ! Add internal sources:
        ! Momentum:
        dot_params(2, ipoi) = dot_params(2, ipoi) &
                            + 0.5d0*(polforce(ipoi) + polforce(ipoi + 1))

        ! Heat into electrons:
        dot_params(3, ipoi) = dot_params(3, ipoi) &
                            + 0.5d0*(qlheat_e(ipoi) + qlheat_e(ipoi + 1))

        ! Heat into ions:
        dot_params(4, ipoi) = dot_params(4, ipoi) &
                            + 0.5d0*(qlheat_i(ipoi) + qlheat_i(ipoi + 1))

        ! Convert momentum time derivative to time derivative of the rotation frequency:
        dot_params(2, ipoi) = dot_params(2, ipoi)*Z_i/params(1, ipoi) &
                              *2.d0/(gpp_av(ipoi + 1) + gpp_av(ipoi))

        ! Convert dot_params from d(nT_{e,i})/dt to d(T_{e,i})/dt:
        dot_params(3, ipoi) = (-params(3, ipoi)*dot_params(1, ipoi) &
                            + dot_params(3, ipoi)/1.5d0)/params(1, ipoi)

        dot_params(4, ipoi) = (-params(4, ipoi)*dot_params(1, ipoi) &
                            + dot_params(4, ipoi)/1.5d0)/params(1, ipoi)

    end do


    ! RHS vector of ODE system
    do ipoi = 1, npoi
        do ieq = 1, nbaleqs
            i = nbaleqs*(ipoi - 1) + ieq
            dy(i) = dot_params(ieq, ipoi)
        end do
    end do

    dy = dy + dery_equisource

end subroutine rhs_balance_source_stell

