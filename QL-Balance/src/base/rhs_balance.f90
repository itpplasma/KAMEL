
subroutine rhs_balance(x, y, dy)

    use grid_mod, only: nbaleqs, neqset, iboutype, npoic, npoib &
                        , Sc, Sb, deriv_coef &
                        , ipbeg, ipend, rb, reint_coef &
                        , fluxes_dif, fluxes_con, rc &
                        , dae11, dae12, dae22, dai11, dai12, dai22 &
                        , dni22, visca, gpp_av &
                        , dqle11, dqle12, dqle21, dqle22 &
                        , dqli11, dqli12, dqli21, dqli22 &
                        , sqg_bthet_overc, Ercov, polforce, qlheat_e, qlheat_i &
                        , Ercov_lin, fluxes_con_nl
                        
    use plasma_parameters, only: params, ddr_params, params_lin, ddr_params_nl &
                        , params_b_lin, params_b, dot_params
    use baseparam_mod, only: Z_i, e_charge, am, p_mass, c
    use control_mod, only: iwrite, irf
    use wave_code_data, only: q, Vth
    use matrix_mod, only: isw_rhs, nz, nsize, irow, icol, amat, rhsvec
    use QLBalance_hdf5_tools

    implicit none

    integer :: ipoi, ieq, i, npoi, ibeg, iend, nshift, ibegb, iendb, ibegtot, iendtot, k, iprobe
    double precision :: x, A_noE_1e, A_noE_2e, A_noE_1i, A_noE_2i, convel
    double precision :: A_noE_1e_nl, A_noE_2e_nl, A_noE_1i_nl, A_noE_2i_nl
    double precision :: gamma_e, gamma_i, dfluxvphi, Q_e, Q_i, A_1e, A_1i
    double precision :: gamma_e_nl, Q_e_nl, Q_i_nl, A_1e_nl, A_1i_nl
    double precision :: gamma_ql_e, gamma_ql_i
    double precision, dimension(neqset) :: y, dy, y_lin ! y is a profile vector, holding the data of all profiles

    if (iboutype .eq. 1) then
        npoi = npoic - 1
    else
        npoi = npoic
    end if

    iwrite = 1

    ! isw_rhs is the switch for initializing the RHS vector, if isw_rhs=0, then the RHS vector is initialized

    y_lin = 0.0d0
    params_lin = 0.0d0
    ddr_params = 0.0d0
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
            ddr_params_nl(ieq, ipoi) &
                = sum(params(ieq, ipbeg(ipoi):ipend(ipoi))*deriv_coef(:, ipoi))
            ! equilibrium parameters at cell boundaries:
            params_b(ieq, ipoi) &
                = sum(params(ieq, ipbeg(ipoi):ipend(ipoi))*reint_coef(:, ipoi))
        end do
    end do

    
    Ercov = sqg_bthet_overc*(params_b(2, :) - Vth*q/rb) &
            + (params_b(4, :)*ddr_params_nl(1, :)/params_b(1, :) + ddr_params_nl(4, :)) &
            /(Z_i*e_charge)

    call calc_equil_diffusion_coeffs

    do ipoi = 1, npoib
        ! Thermodynamic forces for zero radial electric field:
        A_noE_1e_nl = ddr_params_nl(1, ipoi)/params_b(1, ipoi) &
                    - 1.5d0*ddr_params_nl(3, ipoi)/params_b(3, ipoi)
        A_noE_2e_nl = ddr_params_nl(3, ipoi)/params_b(3, ipoi)
        A_noE_1i_nl = ddr_params_nl(1, ipoi)/params_b(1, ipoi) &
                    - 1.5d0*ddr_params_nl(4, ipoi)/params_b(4, ipoi)
        A_noE_2i_nl = ddr_params_nl(4, ipoi)/params_b(4, ipoi)
        
        ! Thermodynamic forces for finite radial electric field:
        A_1e_nl = A_noE_1e_nl + Ercov(ipoi)*e_charge/params_b(3, ipoi)
        A_1i_nl = A_noE_1i_nl - Ercov(ipoi)*e_charge*Z_i/params_b(4, ipoi) !<-FIXED

        ! particle flux densities:
        gamma_e_nl = -(dae11(ipoi)*A_noE_1e_nl + dae12(ipoi)*A_noE_2e_nl &
                        + dqle11(ipoi)*A_1e_nl + dqle12(ipoi)*A_noE_2e_nl &
                        )*params_b(1, ipoi)

        ! total particle flux:
        fluxes_con_nl(1, ipoi) = (Sb(ipoi)*gamma_e_nl - &
                                    (-Sb(ipoi)*ddr_params_nl(1, ipoi)*(dae11(ipoi) &
                                + dqle11(ipoi)*(1.d0 + params_b(4, ipoi)/params_b(3, ipoi)/Z_i))))&
                                /params_b(1, ipoi)
                                                
        ! toroidal moment flux density divided by mass:
        ! total toroidal moment flux:
        fluxes_con_nl(2, ipoi) = 0.d0

        ! electron heat flux density:
        !colli    Q_e_nl=-(dae12(ipoi)*A_noE_1e_nl+dqle12(ipoi)*A_1e_nl    &
        Q_e_nl = -(dae12(ipoi)*A_noE_1e_nl + dqle21(ipoi)*A_1e_nl &
                + (dae22(ipoi) + dqle22(ipoi))*A_noE_2e_nl) &
                 *params_b(1, ipoi)*params_b(3, ipoi)

        ! ion heat flux density:
        !colli    Q_i_nl=-(dai12(ipoi)*A_noE_1i_nl+dqli12(ipoi)*A_1i_nl    &
        Q_i_nl = -(dai12(ipoi)*A_noE_1i_nl + dqli21(ipoi)*A_1i_nl &
                + (dai22(ipoi) + dni22(ipoi) + dqli22(ipoi))*A_noE_2i_nl) &
                 *params_b(1, ipoi)/Z_i*params_b(4, ipoi)

        ! total heat fluxes:
        fluxes_con_nl(3, ipoi) = (Sb(ipoi)*Q_e_nl - &
                                (-Sb(ipoi)*(dae22(ipoi) + dqle22(ipoi))*params_b(1, ipoi) &
                                * ddr_params_nl(3, ipoi))) /params_b(3, ipoi)
        fluxes_con_nl(4, ipoi) = (Sb(ipoi)*Q_i_nl - &
                                !colli    (-Sb(ipoi)*(dai22(ipoi)+dni22(ipoi)+dqli22(ipoi)-2.5d0*dqli12(ipoi)) &
                                (-Sb(ipoi)*(dai22(ipoi) + dni22(ipoi) + dqli22(ipoi) - 2.5d0*dqli21(ipoi)) &
                                   *params_b(1, ipoi)/Z_i*ddr_params_nl(4, ipoi))) &
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
                ! radial derivatives of equilibrium parameters at cell boundaries:
                ddr_params(ieq, ipoi) &
                    = sum(params_lin(ieq, ipbeg(ipoi):ipend(ipoi))*deriv_coef(:, ipoi))
                params_b_lin(ieq, ipoi) &
                    = sum(params_lin(ieq, ipbeg(ipoi):ipend(ipoi))*reint_coef(:, ipoi))
            end do
        end do

        Ercov_lin(ibegb:iendb) &
            = sqg_bthet_overc(ibegb:iendb)*params_b_lin(2, ibegb:iendb) &
                + (params_b(4, ibegb:iendb)*ddr_params(1, ibegb:iendb) &
                /params_b(1, ibegb:iendb) + ddr_params(4, ibegb:iendb)) &
                /(Z_i*e_charge)


        do ipoi = ibegb, iendb
        ! Thermodynamic forces for zero radial electric field:
            A_noE_1e = ddr_params(1, ipoi)/params_b(1, ipoi) &
                        - 1.5d0*ddr_params(3, ipoi)/params_b(3, ipoi)
            A_noE_2e = ddr_params(3, ipoi)/params_b(3, ipoi)
            A_noE_1i = ddr_params(1, ipoi)/params_b(1, ipoi) &
                        - 1.5d0*ddr_params(4, ipoi)/params_b(4, ipoi)
            A_noE_2i = ddr_params(4, ipoi)/params_b(4, ipoi)
        ! Thermodynamic forces for finite radial electric field:
            A_1e = A_noE_1e + Ercov_lin(ipoi)*e_charge/params_b(3, ipoi)
        !ERROR    A_1i=A_noE_1e-Ercov_lin(ipoi)*e_charge*Z_i/params_b(4,ipoi)
            A_1i = A_noE_1i - Ercov_lin(ipoi)*e_charge*Z_i/params_b(4, ipoi) !<-FIXED
        ! particle flux densities:
            gamma_e = -(dae11(ipoi)*A_noE_1e + dae12(ipoi)*A_noE_2e)*params_b(1, ipoi)
            gamma_ql_e = -(dqle11(ipoi)*A_1e + dqle12(ipoi)*A_noE_2e)*params_b(1, ipoi)
            gamma_e = gamma_e + gamma_ql_e
            gamma_i = -(dai11(ipoi)*A_noE_1i + dai12(ipoi)*A_noE_2i)*params_b(1, ipoi)/Z_i
            gamma_ql_i = -(dqli11(ipoi)*A_1i + dqli12(ipoi)*A_noE_2i)*params_b(1, ipoi)/Z_i
            gamma_i = gamma_i + gamma_ql_i
        ! total particle flux:
            fluxes_dif(1, ipoi) = -Sb(ipoi)*ddr_params(1, ipoi)*(dae11(ipoi) &
                                    + dqle11(ipoi)*(1.d0 + params_b(4, ipoi)/params_b(3, ipoi)/Z_i))
            fluxes_con(1, ipoi) = (Sb(ipoi)*gamma_e - fluxes_dif(1, ipoi))/params_b(1, ipoi)
        ! toroidal moment flux density divided by mass:
            dfluxvphi = -visca(ipoi)*ddr_params(2, ipoi)*params_b(1, ipoi)/Z_i*gpp_av(ipoi)
        ! total toroidal moment flux:
            fluxes_dif(2, ipoi) = Sb(ipoi)*dfluxvphi
            fluxes_con(2, ipoi) = 0.d0
        ! electron heat flux density:
        !colli    Q_e=-(dae12(ipoi)*A_noE_1e+dqle12(ipoi)*A_1e             &
            Q_e = -(dae12(ipoi)*A_noE_1e + dqle21(ipoi)*A_1e &
                    + (dae22(ipoi) + dqle22(ipoi))*A_noE_2e) &
                  *params_b(1, ipoi)*params_b(3, ipoi)
        ! ion heat flux density:
        !colli    Q_i=-(dai12(ipoi)*A_noE_1i+dqli12(ipoi)*A_1i             &
            Q_i = -(dai12(ipoi)*A_noE_1i + dqli21(ipoi)*A_1i &
                    + (dai22(ipoi) + dni22(ipoi) + dqli22(ipoi))*A_noE_2i) &
                  *params_b(1, ipoi)/Z_i*params_b(4, ipoi)
        ! total heat fluxes:
            fluxes_dif(3, ipoi) = -Sb(ipoi)*(dae22(ipoi) + dqle22(ipoi)) &
                                  *params_b(1, ipoi)*ddr_params(3, ipoi)
            fluxes_con(3, ipoi) = (Sb(ipoi)*Q_e - fluxes_dif(3, ipoi))/params_b(3, ipoi)

            fluxes_dif(4, ipoi) = -Sb(ipoi)*(dai22(ipoi) + dni22(ipoi) + dqli22(ipoi) &
                                !colli             - 2.5d0*dqli12(ipoi))*params_b(1,ipoi)/Z_i*ddr_params(4,ipoi)
                                - 2.5d0*dqli21(ipoi))*params_b(1, ipoi)/Z_i*ddr_params(4, ipoi)
            fluxes_con(4, ipoi) = (Sb(ipoi)*Q_i - fluxes_dif(4, ipoi))/params_b(4, ipoi)
        ! Momentum source due to the polarization current:
            polforce(ipoi) = (gamma_e - Z_i*gamma_i)*e_charge*sqg_bthet_overc(ipoi) &
                            /(am*p_mass)

        ! Heat sources due to the radial QL drift in the equilibrium electric field:
            qlheat_e(ipoi) = -Ercov(ipoi)*gamma_ql_e*e_charge
            qlheat_i(ipoi) = Z_i*Ercov(ipoi)*gamma_ql_i*e_charge
        end do

        ! Condition of zero flux at the inner boundary:
        fluxes_dif(:, 1) = 0.d0
        fluxes_con(:, 1) = 0.d0
        fluxes_con_nl(:, 1) = 0.d0

    if (irf .eq. 100) then
        if (iwrite .eq. 1) then
            open (15, file = 'fluxes_dif.dat')
            open (16, file = 'fluxes_con.dat')
            open (17, file = 'fluxes_con_nl.dat')
            open (18, file = 'polforce.dat')
            open (19, file = 'qlheat_e.dat')
            open (20, file = 'qlheat_i.dat')
            open (22, file = 'dqle22.dat')
            open (24, file = 'dae22.dat')
            open (23, file = 'params_b.dat')
            open(21, file = 'Sb.dat')
            open(25, file = 'ddr_params.dat')
            open(26, file = 'params_lin.dat')
            do ipoi = 1, npoi
                write (15, *) rb(ipoi), fluxes_dif(1, ipoi), fluxes_dif(2, ipoi), fluxes_dif(3, ipoi), fluxes_dif(4, ipoi)
                write (16, *) rb(ipoi), fluxes_con(1, ipoi), fluxes_con(2, ipoi), fluxes_con(3, ipoi), fluxes_con(4, ipoi)
                write (17, *) rb(ipoi), fluxes_con_nl(1, ipoi), fluxes_con_nl(2, ipoi), fluxes_con_nl(3, ipoi), &
                fluxes_con_nl(4, ipoi)
                write (18, *) rb(ipoi), polforce(ipoi)
                write (19, *) rb(ipoi), qlheat_e(ipoi)
                write (20, *) rb(ipoi), qlheat_i(ipoi)
                write (21, *) rb(ipoi), Sb(ipoi)
                write (22, *) rb(ipoi), dqle22(ipoi)
                write (24, *) rb(ipoi), dae22(ipoi)
                write (23, *) rb(ipoi), params_b(:, ipoi)
                write (25, *) rb(ipoi), ddr_params(:, ipoi)
                write (26, *) rb(ipoi), params_lin(:, ipoi)
                write (27, *) rb(ipoi), deriv_coef(:, ipoi)
            end do
            close(15)
            close(16)
            close(17)
            close(18)
            close(19)
            close(20)
            close(21)
            close(22)
            close(23)
            close(24)
            close(25)
            close(26)
            close(27)


        end if
        print *, 'After writing fluxes'
    end if

        ! Partial time derivatives of equilibrium parameters:
        do ipoi = ibeg, iend
        ! Flux divergence:
            do ieq = 1, nbaleqs
                dot_params(ieq, ipoi) = -(fluxes_dif(ieq, ipoi + 1) - fluxes_dif(ieq, ipoi)) &
                                        /(Sc(ipoi)*(rb(ipoi + 1) - rb(ipoi))) &
                                        - (fluxes_con(ieq, ipoi + 1) - fluxes_con(ieq, ipoi)) &
                                        /(Sc(ipoi)*(rb(ipoi + 1) - rb(ipoi)))*params(ieq, ipoi)
                convel = 0.5d0*(fluxes_con_nl(ieq, ipoi + 1) + fluxes_con_nl(ieq, ipoi))/Sc(ipoi)
        ! upstream convection:
                if (convel .gt. 0.d0) then
        !      if(convel.lt.0.d0) then
                    dot_params(ieq, ipoi) = dot_params(ieq, ipoi) &
                                            - convel*(params_lin(ieq, ipoi + 1) - params_lin(ieq, ipoi))/(rc(ipoi + 1) - rc(ipoi))
                else
                    if (ipoi .gt. 1) then
                        dot_params(ieq, ipoi) = dot_params(ieq, ipoi) &
                                            - convel*(params_lin(ieq, ipoi - 1) - params_lin(ieq, ipoi))/(rc(ipoi - 1) - rc(ipoi))
                    else
                        dot_params(ieq, ipoi) = dot_params(ieq, ipoi) &
                                                - convel*(params_b_lin(ieq, 1) - params_lin(ieq, 1))/(rb(1) - rc(1))
                    end if
                end if
            end do

            !
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
                ddr_params(ieq, ipoi) = 0.d0
            end do
        end do
        Ercov_lin(ibegb:iendb) = 0.d0
        fluxes_dif(:, ibegb:iendb) = 0.d0
        fluxes_con(:, ibegb:iendb) = 0.d0
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
        call rhs_balance_source(x, y, dy)
        rhsvec = dy
    end if

end subroutine rhs_balance


subroutine initialize_rhs(y, dy)

    use grid_mod, only: neqset
    use matrix_mod

    implicit none
    double precision :: x
    double precision, dimension(neqset) :: y, dy

    x = 0.d0
    isw_rhs = 0

    call rhs_balance(x, y, dy)

    isw_rhs = 1
    if (allocated(amat)) deallocate (irow, icol, amat, rhsvec)
    allocate (irow(nz), icol(nz), amat(nz), rhsvec(nsize))

end subroutine initialize_rhs


subroutine rhs_balance_source(x, y, dy)

    use grid_mod, only: nbaleqs, neqset, iboutype, npoic, npoib &
                        , Sc, Sb, deriv_coef &
                        , ipbeg, ipend, rb, reint_coef &
                        , fluxes_dif, fluxes_con, rc &
                        , dae11, dae12, dae22, dai11, dai12, dai22 &
                        , dni22, visca, gpp_av, dery_equisource &
                        , dqle11, dqle12, dqle21, dqle22 &
                        , dqli11, dqli12, dqli21, dqli22 &
                        , sqg_bthet_overc, Ercov, polforce, qlheat_e, qlheat_i &
                        , Ercov_lin, fluxes_con_nl 
                        
    use plasma_parameters, only: params, ddr_params, params_b, params_lin &
                        , params_b_lin, ddr_params_nl, dot_params
    use baseparam_mod, only: Z_i, e_charge, am, p_mass, c
    use wave_code_data, only: q, Vth

    implicit none

    integer :: ipoi, ieq, i, npoi
    double precision :: x, A_noE_1e, A_noE_2e, A_noE_1i, A_noE_2i, convel
    double precision :: A_noE_1e_nl, A_noE_2e_nl, A_noE_1i_nl, A_noE_2i_nl
    double precision :: gamma_e, gamma_i, dfluxvphi, Q_e, Q_i, A_1e, A_1i
    double precision :: gamma_e_nl, Q_e_nl, Q_i_nl, A_1e_nl, A_1i_nl
    double precision :: gamma_ql_e, gamma_ql_i
    double precision, dimension(neqset) :: y, dy, y_lin

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
            ddr_params_nl(ieq, ipoi) &
                = sum(params(ieq, ipbeg(ipoi):ipend(ipoi))*deriv_coef(:, ipoi))
            ddr_params(ieq, ipoi) &
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
    Ercov = sqg_bthet_overc*(params_b(2, :) - Vth*q/rb) &
            + (params_b(4, :)*ddr_params_nl(1, :)/params_b(1, :) + ddr_params_nl(4, :)) &
            /(Z_i*e_charge)

    Ercov_lin = sqg_bthet_overc*params_b_lin(2, :) &
                + (params_b(4, :)*ddr_params(1, :)/params_b(1, :) + ddr_params(4, :)) &
                /(Z_i*e_charge)

    ! Compute diffusion coefficient matrices:
    call calc_equil_diffusion_coeffs
    !
    ! Compute fluxes and internal sources:
    !

    do ipoi = 1, npoib
        ! Thermodynamic forces for zero radial electric field:
        A_noE_1e = ddr_params(1, ipoi)/params_b(1, ipoi) &
                    - 1.5d0*ddr_params(3, ipoi)/params_b(3, ipoi)
        A_noE_2e = ddr_params(3, ipoi)/params_b(3, ipoi)
        A_noE_1i = ddr_params(1, ipoi)/params_b(1, ipoi) &
                    - 1.5d0*ddr_params(4, ipoi)/params_b(4, ipoi)
        A_noE_2i = ddr_params(4, ipoi)/params_b(4, ipoi)

        A_noE_1e_nl = ddr_params_nl(1, ipoi)/params_b(1, ipoi) &
                    - 1.5d0*ddr_params_nl(3, ipoi)/params_b(3, ipoi)
        A_noE_2e_nl = ddr_params_nl(3, ipoi)/params_b(3, ipoi)
        A_noE_1i_nl = ddr_params_nl(1, ipoi)/params_b(1, ipoi) &
                    - 1.5d0*ddr_params_nl(4, ipoi)/params_b(4, ipoi)
        A_noE_2i_nl = ddr_params_nl(4, ipoi)/params_b(4, ipoi)

        ! Thermodynamic forces for finite radial electric field:
        A_1e = A_noE_1e + Ercov_lin(ipoi)*e_charge/params_b(3, ipoi)
        A_1i = A_noE_1i - Ercov_lin(ipoi)*e_charge*Z_i/params_b(4, ipoi) !<-FIXED

        A_1e_nl = A_noE_1e_nl + Ercov(ipoi)*e_charge/params_b(3, ipoi)
        A_1i_nl = A_noE_1i_nl - Ercov(ipoi)*e_charge*Z_i/params_b(4, ipoi) !<-FIXED


        ! particle flux densities:
        gamma_e = -(dae11(ipoi)*A_noE_1e + dae12(ipoi)*A_noE_2e)*params_b(1, ipoi)
        gamma_ql_e = -(dqle11(ipoi)*A_1e + dqle12(ipoi)*A_noE_2e)*params_b(1, ipoi)
        gamma_e = gamma_e + gamma_ql_e
        gamma_i = -(dai11(ipoi)*A_noE_1i + dai12(ipoi)*A_noE_2i)*params_b(1, ipoi)/Z_i
        gamma_ql_i = -(dqli11(ipoi)*A_1i + dqli12(ipoi)*A_noE_2i)*params_b(1, ipoi)/Z_i
        gamma_i = gamma_i + gamma_ql_i

        gamma_e_nl = -(dae11(ipoi)*A_noE_1e_nl + dae12(ipoi)*A_noE_2e_nl &
                        + dqle11(ipoi)*A_1e_nl + dqle12(ipoi)*A_noE_2e_nl &
                        )*params_b(1, ipoi)

        ! total particle flux:
        fluxes_dif(1, ipoi) = -Sb(ipoi)*ddr_params(1, ipoi)*(dae11(ipoi) &
                            + dqle11(ipoi)*(1.d0 + params_b(4, ipoi)/params_b(3, ipoi)/Z_i))
        fluxes_con(1, ipoi) = (Sb(ipoi)*gamma_e - fluxes_dif(1, ipoi))/params_b(1, ipoi)
        fluxes_con_nl(1, ipoi) = (Sb(ipoi)*gamma_e_nl - &
                                (-Sb(ipoi)*ddr_params_nl(1, ipoi)*(dae11(ipoi) &
                            + dqle11(ipoi)*(1.d0 + params_b(4, ipoi)/params_b(3, ipoi)/Z_i))))/params_b(1, ipoi)

        ! toroidal moment flux density divided by mass:
        dfluxvphi = -visca(ipoi)*ddr_params(2, ipoi)*params_b(1, ipoi)/Z_i*gpp_av(ipoi)
        ! total toroidal moment flux:
        fluxes_dif(2, ipoi) = Sb(ipoi)*dfluxvphi
        fluxes_con(2, ipoi) = 0.d0
        fluxes_con_nl(2, ipoi) = 0.d0

        ! electron heat flux density:
        Q_e = -(dae12(ipoi)*A_noE_1e + dqle21(ipoi)*A_1e &
                + (dae22(ipoi) + dqle22(ipoi))*A_noE_2e) &
              *params_b(1, ipoi)*params_b(3, ipoi)

        Q_e_nl = -(dae12(ipoi)*A_noE_1e_nl + dqle21(ipoi)*A_1e_nl &
                + (dae22(ipoi) + dqle22(ipoi))*A_noE_2e_nl) &
                *params_b(1, ipoi)*params_b(3, ipoi)

        ! ion heat flux density:
        Q_i = -(dai12(ipoi)*A_noE_1i + dqli21(ipoi)*A_1i &
                + (dai22(ipoi) + dni22(ipoi) + dqli22(ipoi))*A_noE_2i) &
              *params_b(1, ipoi)/Z_i*params_b(4, ipoi)

        Q_i_nl = -(dai12(ipoi)*A_noE_1i_nl + dqli21(ipoi)*A_1i_nl &
                + (dai22(ipoi) + dni22(ipoi) + dqli22(ipoi))*A_noE_2i_nl) &
                 *params_b(1, ipoi)/Z_i*params_b(4, ipoi)

        ! total heat fluxes:
        fluxes_dif(3, ipoi) = -Sb(ipoi)*(dae22(ipoi) + dqle22(ipoi)) &
                              *params_b(1, ipoi)*ddr_params(3, ipoi)
        fluxes_con(3, ipoi) = (Sb(ipoi)*Q_e - fluxes_dif(3, ipoi))/params_b(3, ipoi)
        fluxes_con_nl(3, ipoi) = (Sb(ipoi)*Q_e_nl - &
                            (-Sb(ipoi)*(dae22(ipoi) + dqle22(ipoi))*params_b(1, ipoi)*ddr_params_nl(3, ipoi))) &
                            /params_b(3, ipoi)

        fluxes_dif(4, ipoi) = -Sb(ipoi)*(dai22(ipoi) + dni22(ipoi) + dqli22(ipoi) &
                        !colli             - 2.5d0*dqli12(ipoi))*params_b(1,ipoi)/Z_i*ddr_params(4,ipoi)
                            - 2.5d0*dqli21(ipoi))*params_b(1, ipoi)/Z_i*ddr_params(4, ipoi)
        fluxes_con(4, ipoi) = (Sb(ipoi)*Q_i - fluxes_dif(4, ipoi))/params_b(4, ipoi)
        fluxes_con_nl(4, ipoi) = (Sb(ipoi)*Q_i_nl - &
                            !colli    (-Sb(ipoi)*(dai22(ipoi)+dni22(ipoi)+dqli22(ipoi)-2.5d0*dqli12(ipoi)) &
                            (-Sb(ipoi)*(dai22(ipoi) + dni22(ipoi) + dqli22(ipoi) - 2.5d0*dqli21(ipoi)) &
                            *params_b(1, ipoi)/Z_i*ddr_params_nl(4, ipoi))) &
                            /params_b(4, ipoi)

        ! Momentum source due to the polarization current:
        polforce(ipoi) = (gamma_e - Z_i*gamma_i)*e_charge*sqg_bthet_overc(ipoi) &
                        /(am*p_mass)

        ! Heat sources due to the radial QL drift in the equilibrium electric field:
        qlheat_e(ipoi) = -Ercov(ipoi)*gamma_ql_e*e_charge
        qlheat_i(ipoi) = Z_i*Ercov(ipoi)*gamma_ql_i*e_charge
    end do
    ! Condition of zero flux at the inner boundary:
    fluxes_dif(:, 1) = 0.d0
    fluxes_con(:, 1) = 0.d0
    fluxes_con_nl(:, 1) = 0.d0

    ! Partial time derivatives of equilibrium parameters:
    do ipoi = 1, npoi
        ! Flux divergence:
        do ieq = 1, nbaleqs
            dot_params(ieq, ipoi) = -(fluxes_dif(ieq, ipoi + 1) - fluxes_dif(ieq, ipoi)) &
                                    /(Sc(ipoi)*(rb(ipoi + 1) - rb(ipoi))) &
                                    - (fluxes_con(ieq, ipoi + 1) - fluxes_con(ieq, ipoi)) &
                                    /(Sc(ipoi)*(rb(ipoi + 1) - rb(ipoi)))*params(ieq, ipoi)
            convel = 0.5d0*(fluxes_con_nl(ieq, ipoi + 1) + fluxes_con_nl(ieq, ipoi))/Sc(ipoi)
            ! upstream convection:
            if (convel .gt. 0.d0) then
                dot_params(ieq, ipoi) = dot_params(ieq, ipoi) &
                                        - convel*(params_lin(ieq, ipoi + 1) - params_lin(ieq, ipoi))/(rc(ipoi + 1) - rc(ipoi))
            else
                if (ipoi .gt. 1) then
                    dot_params(ieq, ipoi) = dot_params(ieq, ipoi) &
                                            - convel*(params_lin(ieq, ipoi - 1) - params_lin(ieq, ipoi))/(rc(ipoi - 1) - rc(ipoi))
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

end subroutine rhs_balance_source