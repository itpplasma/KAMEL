
subroutine rhs_balance(x, y, dy)

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
                        
    use plasma_parameters, only: params, ddr_params, params_lin, ddr_params_nl &
                        , init_params, params_b_lin, params_b, dot_params
    use baseparam_mod, only: Z_i, e_charge, am, p_mass, c, btor
    use control_mod, only: iwrite, irf
    use wave_code_data, only: q, Vth
    use matrix_mod, only: isw_rhs, nz, nsize, irow, icol, amat, rhsvec
    use hdf5_tools

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
        !ERROR    A_1i_nl=A_noE_1e_nl-Ercov(ipoi)*e_charge*Z_i/params_b(4,ipoi)
        A_1i_nl = A_noE_1i_nl - Ercov(ipoi)*e_charge*Z_i/params_b(4, ipoi) !<-FIXED

        ! particle flux densities:
        gamma_e_nl = -(dae11(ipoi)*A_noE_1e_nl + dae12(ipoi)*A_noE_2e_nl &
                       + dqle11(ipoi)*A_1e_nl + dqle12(ipoi)*A_noE_2e_nl &
                       )*params_b(1, ipoi)

        ! total particle flux:
        fluxes_con_nl(1, ipoi) = (Sb(ipoi)*gamma_e_nl - &
                                  (-Sb(ipoi)*ddr_params_nl(1, ipoi)*(dae11(ipoi) &
                                                + dqle11(ipoi)*(1.d0 + params_b(4, ipoi)/params_b(3, ipoi)/Z_i))))/params_b(1, ipoi)
                                                
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
                                  (-Sb(ipoi)*(dae22(ipoi) + dqle22(ipoi))*params_b(1, ipoi)*ddr_params_nl(3, ipoi))) &
                                 /params_b(3, ipoi)
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
!
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!
subroutine initialize_rhs(y, dy)

    use grid_mod, only: neqset
    use matrix_mod

    implicit none
    double precision :: x
    double precision, dimension(neqset) :: y, dy

    x = 0.d0
    isw_rhs = 0
!
    call rhs_balance(x, y, dy)
!
    isw_rhs = 1
    if (allocated(amat)) deallocate (irow, icol, amat, rhsvec)
    allocate (irow(nz), icol(nz), amat(nz), rhsvec(nsize))
!
end subroutine initialize_rhs
!
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!
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
                        , init_params, params_b_lin, ddr_params_nl, dot_params
    use baseparam_mod, only: Z_i, e_charge, am, p_mass, c, btor
    use control_mod, only: iwrite
    use wave_code_data, only: q, Vth

    implicit none

    integer :: ipoi, ieq, i, npoi
    double precision :: x, A_noE_1e, A_noE_2e, A_noE_1i, A_noE_2i, convel
    double precision :: A_noE_1e_nl, A_noE_2e_nl, A_noE_1i_nl, A_noE_2i_nl
    double precision :: gamma_e, gamma_i, dfluxvphi, Q_e, Q_i, A_1e, A_1i
    double precision :: gamma_e_nl, Q_e_nl, Q_i_nl, A_1e_nl, A_1i_nl
    double precision :: gamma_ql_e, gamma_ql_i
    double precision, dimension(neqset) :: y, dy, y_lin
!
    if (iboutype .eq. 1) then
        npoi = npoic - 1
    else
        npoi = npoic
    end if
!npoi=npoic
!
! equilibrium parameters:
!
    y_lin = 0.0d0
    params_lin = 0.0d0
    !params_lin = 0.0d0

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
!   Ercov=sqg_bthet_overc*params_b(2,:)                                 & !OLD
    Ercov = sqg_bthet_overc*(params_b(2, :) - Vth*q/rb) &
            + (params_b(4, :)*ddr_params_nl(1, :)/params_b(1, :) + ddr_params_nl(4, :)) &
            /(Z_i*e_charge)

    Ercov_lin = sqg_bthet_overc*params_b_lin(2, :) &
                + (params_b(4, :)*ddr_params(1, :)/params_b(1, :) + ddr_params(4, :)) &
                /(Z_i*e_charge)

! Compute diffusion coefficient matrices:
!
    call calc_equil_diffusion_coeffs
!
! Compute fluxes and internal sources:
!
!open(8765,file='torqueden.dat')
    do ipoi = 1, npoib
!
! Thermodynamic forces for zero radial electric field:
        A_noE_1e = ddr_params(1, ipoi)/params_b(1, ipoi) &
                   - 1.5d0*ddr_params(3, ipoi)/params_b(3, ipoi)
        A_noE_2e = ddr_params(3, ipoi)/params_b(3, ipoi)
        A_noE_1i = ddr_params(1, ipoi)/params_b(1, ipoi) &
                   - 1.5d0*ddr_params(4, ipoi)/params_b(4, ipoi)
        A_noE_2i = ddr_params(4, ipoi)/params_b(4, ipoi)
!
        A_noE_1e_nl = ddr_params_nl(1, ipoi)/params_b(1, ipoi) &
                      - 1.5d0*ddr_params_nl(3, ipoi)/params_b(3, ipoi)
        A_noE_2e_nl = ddr_params_nl(3, ipoi)/params_b(3, ipoi)
        A_noE_1i_nl = ddr_params_nl(1, ipoi)/params_b(1, ipoi) &
                      - 1.5d0*ddr_params_nl(4, ipoi)/params_b(4, ipoi)
        A_noE_2i_nl = ddr_params_nl(4, ipoi)/params_b(4, ipoi)
!
! Thermodynamic forces for finite radial electric field:
        A_1e = A_noE_1e + Ercov_lin(ipoi)*e_charge/params_b(3, ipoi)
!ERROR    A_1i=A_noE_1e-Ercov_lin(ipoi)*e_charge*Z_i/params_b(4,ipoi)
        A_1i = A_noE_1i - Ercov_lin(ipoi)*e_charge*Z_i/params_b(4, ipoi) !<-FIXED
!
        A_1e_nl = A_noE_1e_nl + Ercov(ipoi)*e_charge/params_b(3, ipoi)
!ERROR    A_1i_nl=A_noE_1e_nl-Ercov(ipoi)*e_charge*Z_i/params_b(4,ipoi)
        A_1i_nl = A_noE_1i_nl - Ercov(ipoi)*e_charge*Z_i/params_b(4, ipoi) !<-FIXED

!
! particle flux densities:
        gamma_e = -(dae11(ipoi)*A_noE_1e + dae12(ipoi)*A_noE_2e)*params_b(1, ipoi)
        gamma_ql_e = -(dqle11(ipoi)*A_1e + dqle12(ipoi)*A_noE_2e)*params_b(1, ipoi)
        gamma_e = gamma_e + gamma_ql_e
        gamma_i = -(dai11(ipoi)*A_noE_1i + dai12(ipoi)*A_noE_2i)*params_b(1, ipoi)/Z_i
        gamma_ql_i = -(dqli11(ipoi)*A_1i + dqli12(ipoi)*A_noE_2i)*params_b(1, ipoi)/Z_i
        gamma_i = gamma_i + gamma_ql_i
!
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
!
! toroidal moment flux density divided by mass:
        dfluxvphi = -visca(ipoi)*ddr_params(2, ipoi)*params_b(1, ipoi)/Z_i*gpp_av(ipoi)
! total toroidal moment flux:
        fluxes_dif(2, ipoi) = Sb(ipoi)*dfluxvphi
        fluxes_con(2, ipoi) = 0.d0
        fluxes_con_nl(2, ipoi) = 0.d0
!
! electron heat flux density:
!colli    Q_e=-(dae12(ipoi)*A_noE_1e+dqle12(ipoi)*A_1e             &
        Q_e = -(dae12(ipoi)*A_noE_1e + dqle21(ipoi)*A_1e &
                + (dae22(ipoi) + dqle22(ipoi))*A_noE_2e) &
              *params_b(1, ipoi)*params_b(3, ipoi)
!
!colli    Q_e_nl=-(dae12(ipoi)*A_noE_1e_nl+dqle12(ipoi)*A_1e_nl    &
        Q_e_nl = -(dae12(ipoi)*A_noE_1e_nl + dqle21(ipoi)*A_1e_nl &
                   + (dae22(ipoi) + dqle22(ipoi))*A_noE_2e_nl) &
                 *params_b(1, ipoi)*params_b(3, ipoi)
!
! ion heat flux density:
!colli    Q_i=-(dai12(ipoi)*A_noE_1i+dqli12(ipoi)*A_1i             &
        Q_i = -(dai12(ipoi)*A_noE_1i + dqli21(ipoi)*A_1i &
                + (dai22(ipoi) + dni22(ipoi) + dqli22(ipoi))*A_noE_2i) &
              *params_b(1, ipoi)/Z_i*params_b(4, ipoi)
!
!colli    Q_i_nl=-(dai12(ipoi)*A_noE_1i_nl+dqli12(ipoi)*A_1i_nl    &
        Q_i_nl = -(dai12(ipoi)*A_noE_1i_nl + dqli21(ipoi)*A_1i_nl &
                   + (dai22(ipoi) + dni22(ipoi) + dqli22(ipoi))*A_noE_2i_nl) &
                 *params_b(1, ipoi)/Z_i*params_b(4, ipoi)
!
! total heat fluxes:
!
        fluxes_dif(3, ipoi) = -Sb(ipoi)*(dae22(ipoi) + dqle22(ipoi)) &
                              *params_b(1, ipoi)*ddr_params(3, ipoi)
        fluxes_con(3, ipoi) = (Sb(ipoi)*Q_e - fluxes_dif(3, ipoi))/params_b(3, ipoi)
        fluxes_con_nl(3, ipoi) = (Sb(ipoi)*Q_e_nl - &
                                  (-Sb(ipoi)*(dae22(ipoi) + dqle22(ipoi))*params_b(1, ipoi)*ddr_params_nl(3, ipoi))) &
                                 /params_b(3, ipoi)
!
        fluxes_dif(4, ipoi) = -Sb(ipoi)*(dai22(ipoi) + dni22(ipoi) + dqli22(ipoi) &
                                         !colli             - 2.5d0*dqli12(ipoi))*params_b(1,ipoi)/Z_i*ddr_params(4,ipoi)
                                         - 2.5d0*dqli21(ipoi))*params_b(1, ipoi)/Z_i*ddr_params(4, ipoi)
        fluxes_con(4, ipoi) = (Sb(ipoi)*Q_i - fluxes_dif(4, ipoi))/params_b(4, ipoi)
        fluxes_con_nl(4, ipoi) = (Sb(ipoi)*Q_i_nl - &
                                  !colli    (-Sb(ipoi)*(dai22(ipoi)+dni22(ipoi)+dqli22(ipoi)-2.5d0*dqli12(ipoi)) &
                                  (-Sb(ipoi)*(dai22(ipoi) + dni22(ipoi) + dqli22(ipoi) - 2.5d0*dqli21(ipoi)) &
                                   *params_b(1, ipoi)/Z_i*ddr_params_nl(4, ipoi))) &
                                 /params_b(4, ipoi)
!
! Momentum source due to the polarization current:
        polforce(ipoi) = (gamma_e - Z_i*gamma_i)*e_charge*sqg_bthet_overc(ipoi) &
                         /(am*p_mass)
!
! Heat sources due to the radial QL drift in the equilibrium electric field:
        qlheat_e(ipoi) = -Ercov(ipoi)*gamma_ql_e*e_charge
        qlheat_i(ipoi) = Z_i*Ercov(ipoi)*gamma_ql_i*e_charge
!
!write(8765,*) rb(ipoi), -e_charge*sqg_bthet_overc(ipoi)            &
!*(dqle11(ipoi)*A_1e_nl+dqle12(ipoi)*A_noE_2e_nl)*params_b(1,ipoi), &
!Z_i*e_charge*sqg_bthet_overc(ipoi)                                 &
!*(dqli11(ipoi)*A_1i_nl+dqli12(ipoi)*A_noE_2i_nl)*params_b(1,ipoi), &
!A_1e_nl,A_noE_2e_nl,dqle11(ipoi),dqle12(ipoi),sqg_bthet_overc(ipoi),&
!A_1i_nl,A_noE_2i_nl,dqli11(ipoi),dqli12(ipoi)
    end do
!close(8765)
!
! Condition of zero flux at the inner boundary:
    fluxes_dif(:, 1) = 0.d0
    fluxes_con(:, 1) = 0.d0
    fluxes_con_nl(:, 1) = 0.d0

!
!
! Partial time derivatives of equilibrium parameters:
!
    do ipoi = 1, npoi
!
! Flux divergence:
        do ieq = 1, nbaleqs
            dot_params(ieq, ipoi) = -(fluxes_dif(ieq, ipoi + 1) - fluxes_dif(ieq, ipoi)) &
                                    /(Sc(ipoi)*(rb(ipoi + 1) - rb(ipoi))) &
                                    - (fluxes_con(ieq, ipoi + 1) - fluxes_con(ieq, ipoi)) &
                                    /(Sc(ipoi)*(rb(ipoi + 1) - rb(ipoi)))*params(ieq, ipoi)
            convel = 0.5d0*(fluxes_con_nl(ieq, ipoi + 1) + fluxes_con_nl(ieq, ipoi))/Sc(ipoi)
! upstream convection:
!      if(convel.lt.0.d0) then
            if (convel .gt. 0.d0) then
                dot_params(ieq, ipoi) = dot_params(ieq, ipoi) &
                                        - convel*(params_lin(ieq, ipoi + 1) - params_lin(ieq, ipoi))/(rc(ipoi + 1) - rc(ipoi))
            else
                if (ipoi .gt. 1) then
                    dot_params(ieq, ipoi) = dot_params(ieq, ipoi) &
                                            - convel*(params_lin(ieq, ipoi - 1) - params_lin(ieq, ipoi))/(rc(ipoi - 1) - rc(ipoi))
                else
                    dot_params(ieq, ipoi) = dot_params(ieq, ipoi)                    !    &
!               -convel*(params_b(ieq,1))/(rb(1)-rc(1))
                end if
            end if
        end do
        if (iwrite .eq. 1) write (13, *) rc(ipoi), dot_params(2, ipoi) &
            , dot_params(3, ipoi), dot_params(4, ipoi) &
            , 0.5d0*(polforce(ipoi) + polforce(ipoi + 1)) &
            , 0.5d0*(qlheat_e(ipoi) + qlheat_e(ipoi + 1)) &
            , 0.5d0*(qlheat_i(ipoi) + qlheat_i(ipoi + 1))

        !
        ! Add internal sources:
        ! Momentum:
        dot_params(2, ipoi) = dot_params(2, ipoi) &
                              + 0.5d0*(polforce(ipoi) + polforce(ipoi + 1))
        if (iwrite .eq. 1) write (11, *) rc(ipoi), dot_params(2, ipoi)

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

        if (iwrite .eq. 1) write (10, *) rc(ipoi), dot_params(2, ipoi) &
            , dot_params(3, ipoi), dot_params(4, ipoi) &
            , 0.5d0*(polforce(ipoi) + polforce(ipoi + 1)) &
            , 0.5d0*(qlheat_e(ipoi) + qlheat_e(ipoi + 1)) &
            , 0.5d0*(qlheat_i(ipoi) + qlheat_i(ipoi + 1))
    end do

    if (iwrite .eq. 1) then
        close (10)
        close (13)
        close (21)
    end if

!
! RHS vector of ODE system
!
    do ipoi = 1, npoi
        do ieq = 1, nbaleqs
            i = nbaleqs*(ipoi - 1) + ieq
            dy(i) = dot_params(ieq, ipoi)
        end do
    end do
!
    dy = dy + dery_equisource
!
end subroutine rhs_balance_source
!
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!
subroutine calc_equil_diffusion_coeffs
    !
    use grid_mod, only: npoib, dae11, dae12, dae22, dai11, dai12, dai22 &
                        , dni22, visca, rb, cneo
    use plasma_parameters, only: params_b
    use baseparam_mod, only: dperp, rsepar
    ! added by Markus Markl, 12.04.2021
    use h5mod
    use control_mod, only: ihdf5IO, debug_mode
    use paramscan_mod, only: viscosity_factor
    use PolyLagrangeInterpolation
    !
    implicit none
    !
    integer :: ipoi
    double precision :: rnorm, weight

    character(1024) :: fname;
    integer :: nr, i, iunit_res
    double precision :: r
    double precision, dimension(:), allocatable :: r_raw
    double precision, dimension(:), allocatable :: Da_raw
    integer :: lb, ub

    !This subroutine was changed to include estimated Da from outside of this code
    !file Da.dat which is located in profiles is read and the other Da.. calculated based on this
    !Changed by Philipp Ulbl 18.05.2020
    ! Added the option to read from hdf5 file by Markus Markl 12.04.2021

    if (ihdf5IO .eq. 1) then
        ! read Da data from hdf5 input file
        fname = "/da_estimation/" ! fname is used for the group name in the hdf5 version
        CALL h5_init()
        CALL h5_open_rw(path2inp, h5_id)
        CALL h5_get_bounds_1(h5_id, trim(fname)//'Da', lb, ub)
        allocate (r_raw(ub), Da_raw(ub))
        CALL h5_get_double_1(h5_id, trim(fname)//'Da', Da_raw)
        CALL h5_get_double_1(h5_id, trim(fname)//'r', r_raw)
        CALL h5_close(h5_id)
        CALL h5_deinit()
        nr = ub
    else
        !code from gengrid to read Da file
        fname = 'profiles/Da.dat';
        iunit_res = 157
        nr = 0
        open (iunit_res, file=fname)
        do
            read (iunit_res, *, end=1)
            nr = nr + 1
        end do
1       continue
        close (iunit_res)
        allocate (r_raw(nr), Da_raw(nr))
        !
        open (iunit_res, file=fname)
        do i = 1, nr
            read (iunit_res, *) r_raw(i), Da_raw(i)
        end do
        close (iunit_res)
    end if
    !code from amn_of_r to interpolate (+do loop)

    !interpolate on balance grid
    if (.not. allocated(coef)) allocate (coef(0:nder, nlagr))

    !open(77, FILE='dae12.dat')
    do ipoi = 1, npoib

        r = rb(ipoi);
        if (r .gt. r_raw(nr)) then
            r = r_raw(nr)
        end if

        !binsearch
        call binsrc(r_raw, 1, nr, r, i)
        call getIndicesForLagrangeInterp(i)
        !lagrange interpolation with order 4 only for function (0)

        call plag_coeff(nlagr, nder, r, r_raw(indBeginInterp:indEndInterp), coef)
        !
        !Da estimated is dae12 -> see notes on conversion
        dae12(ipoi) = sum(coef(0, :)*Da_raw(indBeginInterp:indEndInterp))
        !call localizer(-1.d0, rsepar, rsepar + 0.5d0, rb(ipoi), weight)
        !dae12(ipoi) = dae12(ipoi)*(1.d0 - weight) + weight*1d6
    !    write(77, *) r, dae12(ipoi)
    end do
    !close(77)
        call localizer(-1.d0, rsepar, rsepar + 0.5d0, rb(ipoi), weight)
        dae12(ipoi) = dae12(ipoi)*(1.d0 - weight) + weight*1d6


    !get other da
    dae11 = dae12/1.499999d0 !previously used instead of 1.5d0, no idea why
    dae22 = 3.75d0*dae11
    dai11 = dae11
    dai12 = dae12
    dai22 = dae22
    visca = dae11 * viscosity_factor
    
    dni22 = cneo*params_b(1, :)/sqrt(abs(params_b(4, :)))
    !
end subroutine calc_equil_diffusion_coeffs
!
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!
!
!> @brief subroutine get_dql(timeStep). Calculates quasilinear diffusion coefficients.
subroutine get_dql

    use grid_mod, only: nbaleqs, neqset, iboutype, npoic, npoib &
                        , deriv_coef &
                        , ipbeg, ipend, rb, reint_coef &
                        , rc, sqg_bthet_overc, Ercov &
                        , y, mwind &
                        , dqle11, dqle12, dqle21, dqle22 &
                        , dqli11, dqli12, dqli21, dqli22 &
                        , de11, de12, de21, de22, di11, di12, di21, di22 &
                        , rb_cut_in, re_cut_in, rb_cut_out, re_cut_out, rb &
                        , r_resonant, rmax, d11_misalign, Es_pert_flux
    use plasma_parameters
    use baseparam_mod, only: Z_i, e_charge, am, p_mass, c, btor, e_mass, ev, rtor, pi, rsepar
    use control_mod, only: irf, write_formfactors, ihdf5IO, &
                           diagnostics_output, suppression_mode, &
                           debug_mode, misalign_diffusion
    use time_evolution, only: save_prof_time_step, timeIndex
    use h5mod
    use wave_code_data
    use parallelTools
    use diag_mod, only: write_diag, iunit_diag, write_diag_b, iunit_diag_b, i_mn_loop
    use PolyLagrangeInterpolation    

    implicit none
    !logical :: suppression_mode = .true.

    integer :: modpernode, imin, imax;
    integer :: ipoi, ieq, i, npoi, i_mn, ierr, mwind_save
    double precision, dimension(:), allocatable :: dummy
    double complex, dimension(npoib) :: amn_psi, amn_theta, amn_theta_cyl

    double precision, dimension(npoib) :: spec_weight
    double precision :: weight
    double precision, dimension(npoib) :: vT_e, vT_i, nu_e, nu_i

    DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE :: dqle11_loc
    DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE :: dqle12_loc
    DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE :: dqle21_loc
    DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE :: dqle22_loc
    DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE :: dqli11_loc
    DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE :: dqli12_loc
    DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE :: dqli21_loc
    DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE :: dqli22_loc
    DOUBLE COMPLEX, DIMENSION(:), ALLOCATABLE :: Es_pert_flux_temp
    !DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE :: d11_misalign ! diffusion due to misalignment of equipotentials and flux surfaces
    double complex, dimension(:), allocatable :: formfactor

    ! added variables for interpolation of Brvac
    integer :: ibrabsres
    double precision :: brvac_interp
    double precision :: MI_width

    CHARACTER(LEN=1024) :: tempch

    allocate (dqle11_loc(npoib))
    allocate (dqle12_loc(npoib))
    allocate (dqle21_loc(npoib))
    allocate (dqle22_loc(npoib))
    allocate (dqli11_loc(npoib))
    allocate (dqli12_loc(npoib))
    allocate (dqli21_loc(npoib))
    allocate (dqli22_loc(npoib))
    allocate (formfactor(npoib))  
    if (.not. allocated(d11_misalign)) allocate (d11_misalign(npoib))
    if (.not. allocated(Es_pert_flux)) allocate (Es_pert_flux(npoib))
    if (.not. allocated(Es_pert_flux_temp)) allocate (Es_pert_flux_temp(npoib))

    dqle11_loc = 0.0d0
    dqle12_loc = 0.0d0
    dqle21_loc = 0.0d0
    dqle22_loc = 0.0d0
    dqli11_loc = 0.0d0
    dqli12_loc = 0.0d0
    dqli21_loc = 0.0d0
    dqli22_loc = 0.0d0

    if (irf .eq. 0) then
        return
    elseif (irf .eq. 2) then
        dqle11 = 0.0d0
        dqle12 = 0.0d0
        dqle21 = 0.0d0
        dqle22 = 0.0d0
        dqli11 = 0.0d0
        dqli12 = 0.0d0
        dqli21 = 0.0d0
        dqli22 = 0.0e0
        return
    end if

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

    ! Smooth input for KILCA
    if (.false.) then
        allocate (dummy(npoib))
        do ieq = 1, nbaleqs
            call smooth_array_gauss(npoib, mwind, ddr_params_nl(ieq, :), dummy)
            ddr_params_nl(ieq, :) = dummy
            call smooth_array_gauss(npoib, mwind, params_b(ieq, :), dummy)
            params_b(ieq, :) = dummy
        end do
        deallocate (dummy)
    end if

    ! Compute radial electric field:
    Ercov = sqg_bthet_overc*(params_b(2, :) - Vth*q/rb) &
            + (params_b(4, :)*ddr_params_nl(1, :)/params_b(1, :) + ddr_params_nl(4, :)) &
            /(Z_i*e_charge)

    call MPI_Comm_rank(MPI_COMM_WORLD, irank, ierror);
    if (irank .eq. 0) then
        if (write_diag_b) then !write_diag_b
            write(*,*) "writing params in ", iunit_diag_b
            open(iunit_diag_b)
            do ipoi = 1, npoib
                write (iunit_diag_b, *) rb(ipoi), params_b(1:4, ipoi)
            end do
            close(iunit_diag_b)
        end if
    end if


    ! Compute diffusion coefficient matrices:

    call MPI_Comm_size(MPI_COMM_WORLD, np_num, ierror);
    call MPI_Comm_rank(MPI_COMM_WORLD, irank, ierror);
    !sum over modes:
    if (np_num .gt. dim_mn) then
        print *, ' '
        print *, 'Number of processes', np_num, 'is larger than number of modes', dim_mn
        call MPI_finalize(ierror)
        stop
    end if

    modpernode = ceiling(float(dim_mn)/float(np_num));
    imin = modpernode*irank + 1;
    imax = min(dim_mn, modpernode*(irank + 1));

    if (irf .eq. 1) call update_background_files(path2profs);
    if (irf .eq. 1) call get_wave_code_data(imin, imax);
    if (irf .eq. 1) call get_background_magnetic_fields_from_wave_code(flre_cd_ptr(imin), dim_r, r, B0t, B0z, B0);
    if (irf .eq. 1) call get_collision_frequences_from_wave_code(flre_cd_ptr(imin), dim_r, r, nui, nue);

    !  nu_e=15.4d-6*params_b(1,:)/sqrt(params_b(3,:)/ev)**3            &
    !      *(23.d0-0.5d0*log(params_b(1,:)/(params_b(3,:)/ev)**3))
    !  nu_i=1.d-7*params_b(1,:)/sqrt(params_b(4,:)/ev)**3              &
    !      *(23.d0-0.5d0*log(2.d0*params_b(1,:)/(params_b(4,:)/ev)**3))

    nu_e = nue
    nu_i = nui

    !initialization before summing up over modes:
    dqle11 = 0.0d0
    dqle12 = 0.0d0
    dqle21 = 0.0d0
    dqle22 = 0.0d0
    dqli11 = 0.0d0
    dqli12 = 0.0d0
    dqli21 = 0.0d0
    dqli22 = 0.0e0
    Es_pert_flux = 0.0d0

    !sum over modes:
    do i_mn = imin, imax
        call get_wave_vectors_from_wave_code(flre_cd_ptr(i_mn), dim_r, r, &
                                             m_vals(i_mn), n_vals(i_mn), ks, kp)
        call get_wave_fields_from_wave_code(flre_cd_ptr(i_mn), dim_r, r, &
                                            m_vals(i_mn), n_vals(i_mn), Er, Es, Ep, Et, Ez, Br, Bs, Bp, Bt, Bz)
        om_E = ks*c*dPhi0/B0
        vT_e = sqrt(params_b(3, :)/e_mass)
        vT_i = sqrt(params_b(4, :)/p_mass/am)

        i_mn_loop = i_mn

        ! TODO: add switch to choose calculation of collisionless transport coefficients.
        if (.false.) then
            call calc_transport_coeffs_collisionless(npoib, vT_e, de11, de12, de22)
            de21 = de12
            call calc_transport_coeffs_collisionless(npoib, vT_i, di11, di12, di22)
            di21 = di12
        else
        if (.true.) then
                call calc_transport_coeffs_ornuhl(npoib, vT_e, nu_e, de11, de12, de21, de22)
                call calc_transport_coeffs_ornuhl(npoib, vT_i, nu_i, di11, di12, di21, di22)
            else
                call calc_transport_coeffs_ornuhl_drift(1, npoib, de11, de12, de21, de22)
                call calc_transport_coeffs_ornuhl_drift(2, npoib, di11, di12, di21, di22)
            end if
        end if

        if (misalign_diffusion .eqv. .true.) then
            call get_wave_fields_from_wave_code(flre_cd_ptr(i_mn), dim_r, r, &
                                            m_vals(i_mn), n_vals(i_mn), Er, Es, Ep, Et, Ez, Br, Bs, Bp, Bt, Bz)


            ! caluclate part of perpendicular electric field perturbation that comes from
            if (.not. allocated(coef)) allocate(coef(0:nder,nlagr))
			if (debug_mode) write(*,*) "at r_resonant(i_mn) = ", r_resonant(i_mn)
            call binsrc(rb, 1, npoib, r_resonant(i_mn), ibrabsres)
            if (debug_mode) write(*,*) "binary search found ibrabsres = ", ibrabsres

            call getIndicesForLagrangeInterp(ibrabsres)

            call plag_coeff(nlagr, nder, r_resonant(i_mn), rb(indBeginInterp:indEndInterp), coef)

            CALL magnetic_island_width(coef, nder, nlagr, indBeginInterp, indEndInterp, m_vals(i_mn), MI_width)
 
            ! the perturbed flux surfaces
            !Es_pert_flux_temp = (-dPhi0) * Br * (m_vals(i_mn) * rtor**2d0 - n_vals(i_mn) * r**2d0 / qsaf) &
            !/ (B0 * r * rtor * (n_vals(i_mn) + (m_vals(i_mn)) / qsaf))
            Es_pert_flux_temp = (-dPhi0) * Br * ks / (B0 * kp) 

            ! cut magnetic island from diffusion 
            !do ipoi = 1, npoi
            !    if (r(ipoi) .gt. r_resonant(i_mn) - MI_width/2d0 .and. &
            !    r(ipoi) .lt. r_resonant(i_mn) + MI_width/2d0) then
            !        Es_pert_flux_temp(ipoi) = 0d0
            !    end if
            !end do
            Es_pert_flux = Es_pert_flux + Es_pert_flux_temp
        end if

        call get_wave_fields_from_wave_code(vac_cd_ptr(i_mn), dim_r, r, &
                                            m_vals(i_mn), n_vals(i_mn), Bz, Bz, Bz, Bz, Bz, Br, Bz, Bz, Bz, Bz)

        formfactor = (1.d0, 0.d0)/Br

        ! spec_weight was wrongly set to 2.0. In case that the tmhd code uses double sided Fourier series, it
        ! must be set to 4.0d0, since the factor 2.0d0 should occur in the fields.
        spec_weight = 1.0d0

        do ipoi = 1, npoib
            call localizer(1.d0, rb_cut_out, re_cut_out, r(ipoi), weight)
            spec_weight(ipoi) = spec_weight(ipoi)*weight
            call localizer(-1.d0, rb_cut_in, rb_cut_in, r(ipoi), weight)
            spec_weight(ipoi) = spec_weight(ipoi)*weight
        end do

        if (irank .eq. 0) then
            if (timeIndex .le. 1) then
                if (ihdf5IO .eq. 1) then
                    if (.not. allocated(coef)) allocate(coef(0:nder,nlagr))
                    call binsrc(rb, 1, npoib, r_resonant(1), ibrabsres)

                    call getIndicesForLagrangeInterp(ibrabsres)

                    call plag_coeff(nlagr, nder, r_resonant(1), rb(indBeginInterp:indEndInterp), coef)

                    brvac_interp = sum(coef(0,:)*abs(Br(indBeginInterp:indEndInterp)))

                    if (debug_mode) write(*,*) "Debug: writing Brvac interpolation"
                    if (debug_mode) write(*,*) "Debug: Brvac = ", brvac_interp, " at r_res = ", r_resonant(1)
                    CALL h5_init()
                    CALL h5_open_rw(path2out, h5_id)
                    tempch = "/"//trim(h5_mode_groupname)//"/Brvac_res"
                    CALL h5_add_double_0(h5_id, trim(tempch), brvac_interp)
                    CALL h5_close(h5_id)
                    CALL h5_deinit()

                    if (diagnostics_output) then
                      write (*, *) "writing Brvac.dat"
                      CALL h5_init()
                      CALL h5_open_rw(path2out, h5_id)
                      tempch = "/"//trim(h5_mode_groupname)//"/Brvac.dat"

                      CALL h5_obj_exists(h5_id, trim(tempch), h5_exists_log)
                      if (h5_exists_log) then
                        CALL h5_delete(h5_id, trim(tempch))
                      end if

                      CALL h5_define_unlimited_matrix(h5_id, trim(tempch), &
                                                    H5T_NATIVE_DOUBLE, (/-1, 2/), dataset_id)
                      CALL h5_append_double_1(dataset_id, r, 1)
                      CALL h5_append_double_1(dataset_id, abs(Br), 2)
                      CALL h5_close(h5_id)
                      CALL h5_deinit()
                    end if

                else
                    open (7000, file='Brvac.dat')
                    do ipoi = 1, npoib
                        write (7000, *) r(ipoi), abs(Br(ipoi))
                    end do
                    close (7000)
                end if
            end if
        end if

        call get_wave_fields_from_wave_code(flre_cd_ptr(i_mn), dim_r, r, &
                                            m_vals(i_mn), n_vals(i_mn), Bz, Bz, Bz, Bz, Bz, Br, Bz, Bz, Bz, Bz)
        formfactor = Br*formfactor
        !write(*,*) "sum(abs(formfactor))/size(formfactor) = ", sum(abs(formfactor))/size(formfactor)
        if (write_formfactors) then
            do ipoi = 1, npoib
                write (10000 + n_vals(i_mn)*1000 + m_vals(i_mn), *) r(ipoi), abs(formfactor(ipoi))
            end do
            close (10000 + n_vals(i_mn)*1000 + m_vals(i_mn))
        end if
        !  spec_weight=1.0d0  ! for DIII-D

        dqle11_loc = dqle11_loc + de11*spec_weight
        dqle12_loc = dqle12_loc + de12*spec_weight
        dqle21_loc = dqle21_loc + de21*spec_weight
        dqle22_loc = dqle22_loc + de22*spec_weight
        dqli11_loc = dqli11_loc + di11*spec_weight
        dqli12_loc = dqli12_loc + di12*spec_weight
        dqli21_loc = dqli21_loc + di21*spec_weight
        dqli22_loc = dqli22_loc + di22*spec_weight

        call get_current_densities_from_wave_code(flre_cd_ptr(i_mn), dim_r, r, &
                                                  m_vals(i_mn), n_vals(i_mn), Jri, Jsi, Jpi, Jre, Jse, Jpe)

    end do

    
    ! calculate diffusion due to misalignment of equipotentials and flux surfaces
    if (misalign_diffusion .eqv. .true.) then
        ! rsepar/rtor is the inverse aspect ratio
        d11_misalign = (16.0d0*sqrt(2.0d0) / (9.0d0 * pi**1.5d0)) * (c * abs(Es_pert_flux + Es) / B0)**2.0d0 &
         * (rsepar / rtor)**(1.5d0) / (0.5d0 * nu_e)
        !write(*,*) B0
        !write(*,*) "r         Es     abs(Es)     nu_e"
        !do ipoi = 1, npoib
        !    write(*,*) r(ipoi), Es(ipoi), abs(Es(ipoi)), nu_e(ipoi)
        !end do
         !dqle11_loc = dqle11_loc + d11_misalign
         !dqle12_loc = dqle12_loc + 3 * d11_misalign
         !dqle21_loc = dqle21_loc + 3 * d11_misalign
         !dqle22_loc = dqle22_loc + 12 * d11_misalign
   end if ! misalign_diffusion .eqv. .true.


    call MPI_Allreduce(dqle11_loc, dqle11, npoib, MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_WORLD, ierror);
    call MPI_Allreduce(dqle12_loc, dqle12, npoib, MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_WORLD, ierror);
    call MPI_Allreduce(dqle21_loc, dqle21, npoib, MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_WORLD, ierror);
    call MPI_Allreduce(dqle22_loc, dqle22, npoib, MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_WORLD, ierror);
    call MPI_Allreduce(dqli11_loc, dqli11, npoib, MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_WORLD, ierror);
    call MPI_Allreduce(dqli12_loc, dqli12, npoib, MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_WORLD, ierror);
    call MPI_Allreduce(dqli21_loc, dqli21, npoib, MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_WORLD, ierror);
    call MPI_Allreduce(dqli22_loc, dqli22, npoib, MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_WORLD, ierror);
    call MPI_Barrier(MPI_COMM_WORLD, ierror);
    deallocate (dqle11_loc);
    deallocate (dqle12_loc);
    deallocate (dqle21_loc);
    deallocate (dqle22_loc);
    deallocate (dqli11_loc);
    deallocate (dqli12_loc);
    deallocate (dqli21_loc);
    deallocate (dqli22_loc);
    deallocate (formfactor)
    

    call calc_parallel_current_directly
    call calc_ion_parallel_current_directly


    if (.true.) then
        mwind_save = mwind
        mwind = 30
        allocate (dummy(npoib))
        call smooth_array_gauss(npoib, mwind, dqle11, dummy)
        dqle11 = dummy
        call smooth_array_gauss(npoib, mwind, dqle12, dummy)
        dqle12 = dummy
        call smooth_array_gauss(npoib, mwind, dqle21, dummy)
        dqle21 = dummy
        call smooth_array_gauss(npoib, mwind, dqle22, dummy)
        dqle22 = dummy
        mwind = 30
        call smooth_array_gauss(npoib, mwind, dqli11, dummy)
        dqli11 = dummy
        call smooth_array_gauss(npoib, mwind, dqli12, dummy)
        dqli12 = dummy
        call smooth_array_gauss(npoib, mwind, dqli21, dummy)
        dqli21 = dummy
        call smooth_array_gauss(npoib, mwind, dqli22, dummy)
        dqli22 = dummy
        mwind = mwind_save
        deallocate (dummy)
    end if

    ! set ion particle flux coefficients to zero
    if (.false.) then
        mwind_save = mwind
        mwind = 30
        allocate (dummy(npoib))
        !call smooth_array_gauss(npoib, mwind, dqle11, dummy)
        dqle11 = 0.d0!dummy
        !call smooth_array_gauss(npoib, mwind, dqle12, dummy)
        dqle12 = 0.d0!dummy
        call smooth_array_gauss(npoib, mwind, dqle21, dummy)
        dqle21 = dummy
        call smooth_array_gauss(npoib, mwind, dqle22, dummy)
        dqle22 = dummy
        mwind = 30
        call smooth_array_gauss(npoib, mwind, dqli12, dummy)
        dqli12 = dummy
        call smooth_array_gauss(npoib, mwind, dqli21, dummy)
        dqli21 = dummy
        call smooth_array_gauss(npoib, mwind, dqli21, dummy)
        dqli21 = dummy
        call smooth_array_gauss(npoib, mwind, dqli22, dummy)
        dqli22 = dummy
        mwind = mwind_save
        deallocate (dummy)
    end if


    if (debug_mode) print *, "Debug: writeFieldsCurrentsAndTranspCoeffsToH5"

    if (irank .eq. 0) then
        if (modulo(timeIndex, save_prof_time_step) .eq. 0) then
            if (suppression_mode .eqv. .false.) then
                CALL writeFieldsCurrentsAndTranspCoeffsToH5
            end if
        end if
    end if

    if (debug_mode) write(*,*) "Debug: going out of get_dql"

end subroutine get_dql

subroutine initialize_get_dql

    use control_mod, only: irf

    implicit none

    irf = 2
    call get_dql
    irf = 1

end subroutine

!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

subroutine calc_transport_coeffs_collisionless(dim, vT, D_11, D_12, D_22)

    use baseparam_mod
!use wave_code_data, only: om_E, kp, ks, Ep, Br
    use wave_code_data, only: om_E, kp, ks, Ep, Br, Bs, Bp, r, B0, Te

    integer, intent(in) :: dim
    real(8), dimension(dim), intent(in) :: vT
    real(8), dimension(dim), intent(out) :: D_11, D_12, D_22
    real(8), dimension(dim) :: Z
!complex(8), dimension(dim) :: field_fac
    double precision, dimension(dim) :: field_fac

    Z = om_E/sqrt(2.0)/kp/vT
!Z = om_E/sqrt(2.0)/sqrt(kp**2+1.2e-5**2)/vT

!field_fac = c*ks*Ep - om_E*Br
    field_fac = abs(c*ks*Ep - om_E*Br)**2

!D_11 = sqrt(pi)*vT**2/btor**2*(abs(Z/om_E))**3*exp(-Z**2)*(abs(field_fac))**2
    D_11 = sqrt(pi)*vT**2/btor**2*(abs(Z/om_E))**3*exp(-Z**2)*field_fac
    D_12 = (1.0 + Z**2)*D_11
    D_22 = (1.0 + (1.0 + Z**2)**2)*D_11

!do i=1,dim
!write(333,*) r(i),sqrt(field_fac(i))/abs(c*ks(i)), &
!abs(kp(i)*Bp(i)*Te(i)*ev/e_charge/B0(i))
!write(333,*) r(i),abs(Br(i)),abs(kp(i)*Bs(i))
!enddo
!stop
!do i=1,dim
!write(335,*) r(i),D_11(i),D_12(i),D_12(i),D_22(i)
!enddo
!stop

end subroutine

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine calc_transport_coeffs_ornuhl(dim, vT, nu, D_11, D_12, D_21, D_22)

    use baseparam_mod
    use wave_code_data, only: om_E, kp, Es, Br, B0, r, ks, Ep
    use diag_mod, only: i_mn_loop
    use grid_mod, only: r_resonant, gg_width, rb

    integer, parameter :: mnmax = 3
    integer, intent(in) :: dim
    real(8), dimension(dim), intent(in) :: vT, nu
    real(8), dimension(dim), intent(out) :: D_11, D_12, D_21, D_22
    double precision, dimension(:), allocatable :: x1, x2, comfac, d_12a
    double precision, dimension(:), allocatable :: epm2, brm2, epbr_re, epbr_im
    double complex, dimension(:, :, :), allocatable :: symbI

    allocate (comfac(dim), d_12a(dim), epm2(dim), brm2(dim), epbr_re(dim), epbr_im(dim))
    allocate (x1(dim), x2(dim), symbI(0:mnmax, 0:mnmax, dim))

    symbI = 0.d0
!
!    if  Br=c*kp*Es/om_E diffusion tensor iz zero

    comfac = 0.5d0/(nu*B0**2)
    epm2 = c**2*abs(Es)**2
    brm2 = vT**2*abs(Br)**2
    epbr_re = 2.d0*c*vT*real(conjg(Es)*Br)
    epbr_im = 2.d0*c*vT*dimag(conjg(Es)*Br)
!epm2=0.0d0 !c**2*abs(Es)**2
!brm2=1.0d0 !vT**2*abs(Br)**2
!epbr_re=0.0d0 !2.d0*c*vT*real(conjg(Es)*Br)
!epbr_im=0.0d0 !2.d0*c*vT*dimag(conjg(Es)*Br)

    x1 = kp*vT/nu
    x2 = -om_E/nu
    
    do i = 1, dim
        if (rb(i) .lt. r_resonant(i_mn_loop) - 2.d0*gg_width) cycle
        if (rb(i) .gt. r_resonant(i_mn_loop) + 2.d0*gg_width) cycle
        call getIfunc(x1(i), x2(i), symbI(:, :, i))
    end do

    D_11 = comfac*(epm2*real(symbI(0, 0, :)) &
                   + epbr_re*real(symbI(1, 0, :)) &
                   + brm2*real(symbI(1, 1, :)))
    D_12 = comfac*(epm2*real(symbI(0, 0, :) + 0.5d0*symbI(2, 0, :)) &
                   + epbr_re*real(symbI(1, 0, :) + 0.25d0*(symbI(3, 0, :) + symbI(2, 1, :))) &
                   + brm2*real(symbI(1, 1, :) + 0.5d0*symbI(3, 1, :)))
    D_21 = D_12
    D_22 = comfac*(epm2*real(2.d0*symbI(0, 0, :) + symbI(2, 0, :) &
                             + 0.25d0*symbI(2, 2, :)) &
                   + epbr_re*real(2.d0*symbI(1, 0, :) &
                                  + 0.5d0*(symbI(3, 0, :) + symbI(2, 1, :)) &
                                  + 0.25d0*symbI(3, 2, :)) &
                   + brm2*real(2.d0*symbI(1, 1, :) + symbI(3, 1, :) &
                               + 0.25d0*symbI(3, 3, :)))

    D_12a = comfac*epbr_im*0.25d0*dimag(symbI(2, 1, :) - symbI(3, 0, :))

    D_12 = D_12 + D_12a
    D_21 = D_21 - D_12a

!D_11=comfac*epm2*real(symbI(0,0,:))
!D_12=comfac*epbr_re*real(symbI(1,0,:))
!D_22=comfac*brm2*real(symbI(1,1,:))
!do i=1,dim
!write(7000,*) r(i),D_11(i),D_12(i),D_22(i),D_11(i)+D_12(i)+D_22(i)
!enddo
!D_11=comfac*real(symbI(0,0,:))
!D_12=comfac*real(symbI(1,0,:))
!D_22=comfac*real(symbI(1,1,:))
!do i=1,dim
!write(7001,*) r(i),D_11(i),D_12(i),D_22(i),x1(i),x2(i)
!enddo
!stop

    deallocate (x1, x2, symbI)
    deallocate (comfac, d_12a, epm2, brm2, epbr_re, epbr_im)

end subroutine

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
subroutine getIfunc(x1, x2, symbI)
    integer, parameter :: mnmax = 3
    integer :: m, n
    double precision, intent(in) :: x1, x2
    double precision :: z
    double complex :: denom
    double complex, dimension(0:mnmax, 0:mnmax), intent(out) :: symbI
    double complex, dimension(0:mnmax, 0:mnmax) :: Imn
!
!  if(.true.) then
    if (.false.) then
! collisionless case:
        symbI = (0.d0, 0.d0)
        z = x2/(sqrt(2.d0)*x1)
!
        symbI(0, 0) = sqrt(2.d0)*exp(-z**2)/abs(x1)
        symbI(1, 0) = symbI(0, 0)*x2/x1
        symbI(1, 1) = symbI(1, 0)*x2/x1
!
        symbI(2, 0) = symbI(0, 0)*(x2/x1)**2
        symbI(2, 1) = symbI(1, 0)*(x2/x1)**2
        symbI(3, 0) = symbI(2, 1)
        symbI(3, 1) = symbI(1, 1)*(x2/x1)**2
!
        symbI(2, 2) = symbI(2, 0)*(x2/x1)**2
        symbI(3, 2) = symbI(2, 1)*(x2/x1)**2
        symbI(3, 3) = symbI(3, 1)*(x2/x1)**2
    else
! collisional case:
!
        call W2_arr(x1, x2, Imn)
!
        if (.true.) then
!    if(.false.) then
! energy conservation:
            denom = (1.d0, 0.d0) - Imn(0, 0) + (2.d0, 0.d0)*Imn(2, 0) - Imn(2, 2)
            do m = 0, 3
                do n = 0, 3
                    symbI(m, n) = Imn(m, n) + (Imn(m, 0) - Imn(m, 2))*(Imn(n, 0) - Imn(n, 2))/denom
                end do
            end do
        else
            symbI = Imn
        end if
!
    end if
!
end subroutine getIfunc
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!

subroutine smooth_array_gauss(dimx, ngauss, y, ys)

    implicit none

    integer :: dimx, ngauss, i
    real(8) :: sgauss, dummy
    real(8), dimension(dimx) :: y, ys
    double precision, dimension(:), allocatable :: wgauss
    integer :: k, j;
    real(8) :: sum

    sgauss = dfloat(ngauss)/5.d0
    allocate (wgauss(-ngauss:ngauss))
    do i = -ngauss, ngauss
        wgauss(i) = exp(-(dfloat(i)/sgauss)**2)
    end do
    dummy = sum(wgauss)
    wgauss = wgauss/dummy

    do i = ngauss + 1, dimx - ngauss
        ys(i) = sum(wgauss*y(i - ngauss:i + ngauss))
    end do
    do i = 1, ngauss
        ys(i) = sum(wgauss(1 - i:ngauss)*y(1:i + ngauss))/sum(wgauss(1 - i:ngauss))
    end do
    do i = dimx - ngauss + 1, dimx
        ys(i) = sum(wgauss(-ngauss:dimx - i)*y(i - ngauss:dimx))/sum(wgauss(-ngauss:dimx - i))
    end do

end subroutine
!
subroutine equipotentials
!
    use baseparam_mod, only: c, rtor
    use wave_code_data, only: dim_r, m_vals, n_vals, r, om_E, ks, kp, Es, Br, B0 &
                              , vac_cd_ptr
    use diag_mod, only: iunit_diag

    use mpi
!
    implicit none
!
    integer :: ierror, irank;
    integer :: nr, i, m, n, ipoi, ierr
    double complex :: amnp, amnt, ampl
    double precision, dimension(:), allocatable :: psi0, phi0
    double complex, dimension(:), allocatable :: psi1, phi1, dum, Brv
    integer :: ind = 1
!
    nr = dim_r
    m = m_vals(ind)
    n = n_vals(ind)
!
    allocate (psi0(nr), phi0(nr), psi1(nr), phi1(nr), dum(nr), Brv(nr))
!
    psi0(1) = 0.d0
    phi0(1) = 0.d0
    ipoi = 1
!
    do i = 2, nr
        psi0(i) = psi0(i - 1) + 0.5d0*(r(i) - r(i - 1))*(B0(i)*kp(i) + B0(i - 1)*kp(i - 1))
        phi0(i) = phi0(i - 1) + 0.5d0*(r(i) - r(i - 1))/c &
                  *(om_E(i)*B0(i)/ks(i) + om_E(i - 1)*B0(i - 1)/ks(i - 1))
        if (kp(i)*kp(i - 1) .lt. 0) ipoi = i
    end do
!
    !Not used anymore
    !Changed by Philipp Ulbl 13.05.2020
    !call amn_of_r(-m,n,r(ipoi),amnp,amnt,ierr)
    call get_wave_fields_from_wave_code(vac_cd_ptr(ind), nr, r, m, n, &
                                        dum, dum, dum, dum, dum, Brv, dum, dum, dum, dum)

!
    !amn not used anymore, use value for DIII-D
    !Changed by Philipp Ulbl 13.05.2020
    !ampl = 2.d0*amnt*dble(n)/(r(ipoi)*rtor*Brv(ipoi))
    ampl = 2.d0*dble(n)/(r(ipoi)*rtor*Brv(ipoi))  ! DIII-D
    print *, 'equipotentials: AMN ARE NOT SUPPORTED ANYMORE, THIS SUBROUTINE GIVE WRONG VALUES!!!'
!
    psi1 = (0.d0, 1.d0)*Br*ampl
    phi1 = (0.d0, 1.d0)*Es/ks*ampl
!
    call MPI_Comm_rank(MPI_COMM_WORLD, irank, ierror);
    if (irank .eq. 0) then
        open (iunit_diag, file='equipotentials.dat')
        do i = 1, nr
            write (iunit_diag, *) r(i), psi0(i), phi0(i), real(psi1(i)), dimag(psi1(i)) &
                , real(phi1(i)), dimag(phi1(i)), abs(Br(i)) &
                , abs(Br(i) - c*kp(i)*Es(i)/om_E(i))
        end do
        close (iunit_diag)
    end if

    deallocate (psi0, phi0, psi1, phi1, dum, Brv)
!
end subroutine equipotentials
!
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!
subroutine calc_parallel_current_directly
    ! this subroutine calculates the electron parallel current (eq. (60) in Heyn et. al 2014)

    use grid_mod, only: npoib, rb, Ercov
    use plasma_parameters, only: params_b, ddr_params_nl
    use baseparam_mod, only: Z_i, e_charge, am, p_mass, c, btor, e_mass, ev, rtor
    use control_mod, only: ihdf5IO, diagnostics_output, write_gyro_current, &
        gyro_current_study
    use h5mod
    use wave_code_data
    use mpi
!
    implicit none
!
    integer :: ierror, irank;
    integer :: ipoi, i, iunit, mnmax
    double precision, dimension(:), allocatable :: dummy, x1, x2, vT, A1, A2
    double complex, dimension(:), allocatable :: curr_e_par
    double complex, dimension(:, :, :), allocatable :: symbI
    double precision, dimension(1) :: coll_fac = (/1/)!(/10.0, 100.0, 1.0e3, 1.0e4, 1.0e5/)
    integer :: study_i_omE
    integer :: study_j_nue
!
    character(len=1024) :: tempch
!
    iunit = 731
    mnmax = 3
    allocate (x1(npoib), x2(npoib), A1(npoib), A2(npoib), symbI(0:mnmax, 0:mnmax, npoib))
    allocate (curr_e_par(npoib), vT(npoib))
!
    vT = sqrt(params_b(3, :)/e_mass)
!

    if (gyro_current_study .eq. 0) then
        x1 = kp*vT/nue
        x2 = -om_E/nue
!
        do i = 1, npoib
            call getIfunc(x1(i), x2(i), symbI(:, :, i))
        end do


        A2 = ddr_params_nl(3, :)/params_b(3, :)
        A1 = ddr_params_nl(1, :)/params_b(1, :) + e_charge*Ercov/params_b(3, :) - 1.5d0*A2
!
        curr_e_par = e_charge*params_b(1, :)*vT/(nue*B0) &
                 *(c*Es*((A1 + A2)*symbI(1, 0, :) + 0.5d0*A2*symbI(2, 1, :)) &
                   + vT*Br*((A1 + A2)*symbI(1, 1, :) + 0.5d0*A2*symbI(3, 1, :)))
!
        call MPI_Comm_rank(MPI_COMM_WORLD, irank, ierror);
        if (irank .eq. 0) then
            CALL h5_init()
            CALL h5_open_rw(path2out, h5_id)
            tempch = "/"//trim(h5_mode_groupname)//"/par_current_e/"
            write(*,*) "In group: "//trim(tempch)

            CALL h5_obj_exists(h5_id, trim(tempch), h5_exists_log)
            if (.not. h5_exists_log) then
                CALL h5_define_group(h5_id, trim(tempch), group_id_1)
                CALL h5_close_group(group_id_1)
            end if
            

            CALL h5_add_double_1(h5_id, trim(tempch)//"rb", &
                    rb, lbound(rb), ubound(rb))
            CALL h5_add_double_1(h5_id, trim(tempch)//"x2", &
                    x2, lbound(x2), ubound(x2)) 

            if (write_gyro_current) then
                write(*,*) "writing par_current_e.dat"
                write(*,*) " - "
                ! Write out gyro current which is different to KiLCA current Jpe.
            ! The gyro current is calculated from (60) in Heyn et. al 2014
               !CALL h5_add_double_1(h5_id, trim(tempch)//"par_current_e_real", &
                !   real(curr_e_par), lbound(real(curr_e_par)), ubound(real(curr_e_par)))
                !ALL h5_add_double_1(h5_id, trim(tempch)//"par_current_e_imag", &
                !   dimag(curr_e_par), lbound(dimag(curr_e_par)), ubound(dimag(curr_e_par))) 
                !ALL h5_add_double_1(h5_id, trim(tempch)//"Jpe_real", &
                !   real(Jpe), lbound(real(Jpe)), ubound(real(Jpe)))
                !ALL h5_add_double_1(h5_id, trim(tempch)//"Jpe_imag", &
                !   dimag(Jpe), lbound(dimag(Jpe)), ubound(dimag(Jpe))) 
                !ALL h5_add_double_1(h5_id, trim(tempch)//"Jse_real", &
                !   real(Jse), lbound(real(Jse)), ubound(real(Jse)))
                !ALL h5_add_double_1(h5_id, trim(tempch)//"Jse_imag", &
                !   dimag(Jse), lbound(dimag(Jse)), ubound(dimag(Jse))) 
                !ALL h5_add_double_1(h5_id, trim(tempch)//"Jre_real", &
                !   real(Jpe), lbound(real(Jpe)), ubound(real(Jpe)))
                !ALL h5_add_double_1(h5_id, trim(tempch)//"Jre_imag", &
                !   dimag(Jre), lbound(dimag(Jre)), ubound(dimag(Jre))) 
                CALL h5_add_double_1(h5_id, trim(tempch)//"kp", &
                    kp, lbound(kp), ubound(kp))
                CALL h5_add_double_1(h5_id, trim(tempch)//"ks", &
                    ks, lbound(ks), ubound(ks)) 
                CALL h5_add_double_1(h5_id, trim(tempch)//"x1", &
                    x1, lbound(x1), ubound(x1))
                CALL h5_add_double_1(h5_id, trim(tempch)//"nue", &
                    nue, lbound(nue), ubound(nue)) 
                CALL h5_add_double_1(h5_id, trim(tempch)//"om_E", &
                    om_E, lbound(om_E), ubound(om_E)) 
                CALL h5_add_double_1(h5_id, trim(tempch)//"vT", &
                    vT, lbound(vT), ubound(vT)) 

                CALL h5_add_double_1(h5_id, trim(tempch)//"I00_re", &
                    real(symbI(0,0,:)), lbound(symbI(0,0,:)), ubound(symbI(0,0,:))) 
                CALL h5_add_double_1(h5_id, trim(tempch)//"I00_im", &
                    dimag(symbI(0,0,:)), lbound(symbI(0,0,:)), ubound(symbI(0,0,:))) 
 
                CALL h5_add_double_1(h5_id, trim(tempch)//"I20_re", &
                    real(symbI(2,0,:)), lbound(symbI(2,0,:)), ubound(symbI(2,0,:))) 
                CALL h5_add_double_1(h5_id, trim(tempch)//"I20_im", &
                    dimag(symbI(2,0,:)), lbound(symbI(2,0,:)), ubound(symbI(2,0,:))) 
 
                CALL h5_add_double_1(h5_id, trim(tempch)//"I22_re", &
                    real(symbI(2,2,:)), lbound(symbI(2,2,:)), ubound(symbI(2,2,:))) 
                CALL h5_add_double_1(h5_id, trim(tempch)//"I22_im", &
                    dimag(symbI(2,2,:)), lbound(symbI(2,2,:)), ubound(symbI(2,2,:))) 

            end if ! write_gyro_current
            CALL h5_close(h5_id)
            CALL h5_deinit()



            if (diagnostics_output) then
                if (ihdf5IO .eq. 1) then
                    write(*,*) "writing par_current_e.dat"
                    write(*,*) " - "
                    CALL h5_init()
                    CALL h5_open_rw(path2out, h5_id)

                    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
                    ! par_current_e data
                    tempch = "/"//trim(h5_mode_groupname)//"/par_current_e.dat"

                    CALL h5_obj_exists(h5_id, trim(tempch), h5_exists_log)
                    if (h5_exists_log) then
                        CALL h5_delete(h5_id, trim(tempch))
                    end if

                    CALL h5_define_unlimited_matrix(h5_id, trim(tempch), &
                                                H5T_NATIVE_DOUBLE, (/-1, 5/), dataset_id)
                    CALL h5_append_double_1(dataset_id, rb, 1)
                    CALL h5_append_double_1(dataset_id, real(curr_e_par), 2)
                    CALL h5_append_double_1(dataset_id, dimag(curr_e_par), 3)
                    CALL h5_append_double_1(dataset_id, real(Jpe), 4)
                    CALL h5_append_double_1(dataset_id, dimag(Jpe), 5)

                    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
                    ! cond_e data
                    tempch = "/"//trim(h5_mode_groupname)//"/cond_e.dat"

                    CALL h5_obj_exists(h5_id, trim(tempch), h5_exists_log)
                    if (h5_exists_log) then
                        CALL h5_delete(h5_id, trim(tempch))
                    end if

                    CALL h5_define_unlimited_matrix(h5_id, trim(tempch), &
                                                H5T_NATIVE_DOUBLE, (/-1, 9/), dataset_id)
                    CALL h5_append_double_1(dataset_id, rb, 1)
                    CALL h5_append_double_1(dataset_id, real(symbI(1, 0, :)), 2)
                    CALL h5_append_double_1(dataset_id, dimag(symbI(1, 0, :)), 3)
                    CALL h5_append_double_1(dataset_id, real(symbI(1, 1, :)), 4)
                    CALL h5_append_double_1(dataset_id, dimag(symbI(1, 1, :)), 5)
                    CALL h5_append_double_1(dataset_id, real(symbI(2, 1, :)), 6)
                    CALL h5_append_double_1(dataset_id, dimag(symbI(2, 1, :)), 7)
                    CALL h5_append_double_1(dataset_id, real(symbI(3, 1, :)), 8)
                    CALL h5_append_double_1(dataset_id, dimag(symbI(3, 1, :)), 9)

                    CALL h5_close(h5_id)
                    CALL h5_deinit()
                else ! ihdf5IO .eq. 1
                    open (iunit, file='par_current_e.dat')
                    open (10000, file='cond_e.dat')
                    do ipoi = 1, npoib
                        write (iunit, *) rb(ipoi) &
                            , real(curr_e_par(ipoi)), dimag(curr_e_par(ipoi)) &
                            , real(Jpe(ipoi)), dimag(Jpe(ipoi))
                        write (10000, *) rb(ipoi) &
                            , real(symbI(1, 0, ipoi)), dimag(symbI(1, 0, ipoi)) &
                            , real(symbI(1, 1, ipoi)), dimag(symbI(1, 1, ipoi)) &
                            , real(symbI(2, 1, ipoi)), dimag(symbI(2, 1, ipoi)) &
                            , real(symbI(3, 1, ipoi)), dimag(symbI(3, 1, ipoi))

                    end do
                    close (10000)
                    close (iunit)
                end if
            end if
        end if ! irank .eq. 0

    elseif (gyro_current_study .eq. 1) then ! used to scan over nue and omega_E, deprecated
        write(*,*) " - - - - - - - - - "
        write(*,*) "Gyro current study"

        CALL h5_init()
        CALL h5_open_rw(path2out, h5_id)
        tempch = "/"//trim(h5_mode_groupname)//"/gyro_current_study/"
        write(*,*) "In group: "//trim(tempch)

        CALL h5_define_group(h5_id, trim(tempch), group_id_1)
        CALL h5_close_group(group_id_1)

        ! Quantities that are not affected by parameter study:
        CALL h5_add_double_1(h5_id, trim(tempch)//"rb", &
            rb, lbound(rb), ubound(rb))
        CALL h5_add_double_1(h5_id, trim(tempch)//"nue0", &
            nue, lbound(nue), ubound(nue))
        CALL h5_add_double_1(h5_id, trim(tempch)//"om_E", &
            om_E, lbound(om_E), ubound(om_E))

        CALL h5_add_double_1(h5_id, trim(tempch)//"B0", &
            B0, lbound(B0), ubound(B0))

        CALL h5_add_double_1(h5_id, trim(tempch)//"Br_real", &
            real(Br), lbound(real(Br)), ubound(real(Br)))
        CALL h5_add_double_1(h5_id, trim(tempch)//"Es_real", &
            real(Es), lbound(real(Es)), ubound(real(Es)))

        CALL h5_add_double_1(h5_id, trim(tempch)//"Br_imag", &
            dimag(Br), lbound(dimag(Br)), ubound(dimag(Br)))
        CALL h5_add_double_1(h5_id, trim(tempch)//"Es_imag", &
            dimag(Es), lbound(dimag(Es)), ubound(dimag(Es)))
        ! - - - - - - - - - - - - - - - - - - - - - - - - -

        A2 = ddr_params_nl(3, :)/params_b(3, :)
        A1 = ddr_params_nl(1, :)/params_b(1, :) + e_charge*Ercov/params_b(3, :) - 1.5d0*A2
!
        CALL h5_add_double_1(h5_id, trim(tempch)//"A1", &
            A1, lbound(A1), ubound(A1))
        CALL h5_add_double_1(h5_id, trim(tempch)//"A2", &
            A2, lbound(A2), ubound(A2))

        ! vary om_E/nue
        do study_i_omE = 1, 1
            tempch = "/"//trim(h5_mode_groupname)//"/gyro_current_study/"
            write(tempch, "(A,i4.4,A)") trim(tempch), study_i_omE, "/"

            CALL h5_define_group(h5_id, trim(tempch), group_id_1)
            CALL h5_close_group(group_id_1)

            do study_j_nue = 1, size(coll_fac)

                om_E = om_E * study_i_omE
                nue = nue * coll_fac(study_j_nue)
                x2 = -om_E/(nue)! * study_i * 0.5
                x1 = kp*vT/(nue)! * study_i * 0.5
!
                do i = 1, npoib
                    call getIfunc(x1(i), x2(i), symbI(:, :, i))
                end do
                tempch = "/"//trim(h5_mode_groupname)//"/gyro_current_study/"
                write(tempch, "(A,i4.4,A,i4.4,A)") trim(tempch), study_i_omE, "/", &
                    study_j_nue, "/"
                write(*,*) "In group: ", trim(tempch)


                CALL h5_define_group(h5_id, trim(tempch), group_id_1)
                CALL h5_close_group(group_id_1)
                ! write I functions 
                ! real part
                !CALL h5_add_double_1(h5_id, trim(tempch)//"I10_real", &
                !    real(symbI(1,0,:)), lbound(real(symbI(1,0,:))), &
                !    ubound(real(symbI(1,0,:))))
                !CALL h5_add_double_1(h5_id, trim(tempch)//"I11_real", &
                !    real(symbI(1,1,:)), lbound(real(symbI(1,1,:))), &
                !    ubound(real(symbI(1,1,:))))
                !CALL h5_add_double_1(h5_id, trim(tempch)//"I21_real", &
                !    real(symbI(2,1,:)), lbound(real(symbI(2,1,:))), &
                !    ubound(real(symbI(2,1,:))))
                !CALL h5_add_double_1(h5_id, trim(tempch)//"I31_real", &
                !    real(symbI(3,1,:)), lbound(real(symbI(3,1,:))), &
                !    ubound(real(symbI(3,1,:))))
                ! imag part

                !CALL h5_add_double_1(h5_id, trim(tempch)//"I10_imag", &
                !    dimag(symbI(1,0,:)), lbound(dimag(symbI(1,0,:))),&
                !    ubound(dimag(symbI(1,0,:))))
                !CALL h5_add_double_1(h5_id, trim(tempch)//"I11_imag", &
                !    dimag(symbI(1,1,:)), lbound(dimag(symbI(1,1,:))), &
                !    ubound(dimag(symbI(1,1,:))))
                !CALL h5_add_double_1(h5_id, trim(tempch)//"I21_imag", &
                !    dimag(symbI(2,1,:)), lbound(dimag(symbI(2,1,:))), &
                !    ubound(dimag(symbI(2,1,:))))
                !CALL h5_add_double_1(h5_id, trim(tempch)//"I31_imag", &
                !    dimag(symbI(3,1,:)), lbound(dimag(symbI(3,1,:))), &
                !    ubound(dimag(symbI(3,1,:))))

                curr_e_par = e_charge*params_b(1, :)*vT/(nue*B0) &
                 *(c*Es*((A1 + A2)*symbI(1, 0, :) + 0.5d0*A2*symbI(2, 1, :)) &
                   + vT*Br*((A1 + A2)*symbI(1, 1, :) + 0.5d0*A2*symbI(3, 1, :)))

                call get_current_densities_from_wave_code(flre_cd_ptr(1), dim_r, r, &
                                                  m_vals(1), n_vals(1), Jri, Jsi, Jpi, Jre, Jse, Jpe)
!
                CALL h5_add_double_1(h5_id, trim(tempch)//"je_real", &
                    real(curr_e_par), lbound(real(curr_e_par)), ubound(real(curr_e_par)))
                CALL h5_add_double_1(h5_id, trim(tempch)//"je_imag", &
                    dimag(curr_e_par), lbound(dimag(curr_e_par)), ubound(dimag(curr_e_par))) 

                CALL h5_add_double_1(h5_id, trim(tempch)//"Je_real", &
                    real(Jpe), lbound(real(Jpe)), ubound(real(Jpe)))
                CALL h5_add_double_1(h5_id, trim(tempch)//"Je_imag", &
                    dimag(Jpe), lbound(dimag(Jpe)), ubound(dimag(Jpe))) 
            end do ! study_j_nue
        end do ! study_i_omE

        CALL h5_close(h5_id)
        CALL h5_deinit()

        write(*,*) " - - - - - - - - - "

    elseif (gyro_current_study .eq. 2) then
        write(*,*) " - - - - - - - - - "
        write(*,*) "Writing currents"
        x1 = kp*vT/nue
        x2 = -om_E/nue


        CALL h5_init()
        CALL h5_open_rw(path2out, h5_id)
        tempch = "/"//trim(h5_mode_groupname)//"/currents/"
        write(*,*) "In group: "//trim(tempch)

        CALL h5_define_group(h5_id, trim(tempch), group_id_1)
        CALL h5_close_group(group_id_1)

        ! Quantities that are not affected by parameter study:
        CALL h5_add_double_1(h5_id, trim(tempch)//"rb", &
            rb, lbound(rb), ubound(rb))
        CALL h5_add_double_1(h5_id, trim(tempch)//"nue0", &
            nue, lbound(nue), ubound(nue))
        CALL h5_add_double_1(h5_id, trim(tempch)//"om_E", &
            om_E, lbound(om_E), ubound(om_E))

        CALL h5_add_double_1(h5_id, trim(tempch)//"B0", &
            B0, lbound(B0), ubound(B0))

        CALL h5_add_double_1(h5_id, trim(tempch)//"Br_real", &
            real(Br), lbound(real(Br)), ubound(real(Br)))
        CALL h5_add_double_1(h5_id, trim(tempch)//"Es_real", &
            real(Es), lbound(real(Es)), ubound(real(Es)))

        CALL h5_add_double_1(h5_id, trim(tempch)//"Br_imag", &
            dimag(Br), lbound(dimag(Br)), ubound(dimag(Br)))
        CALL h5_add_double_1(h5_id, trim(tempch)//"Es_imag", &
            dimag(Es), lbound(dimag(Es)), ubound(dimag(Es)))
        ! - - - - - - - - - - - - - - - - - - - - - - - - -

        A2 = ddr_params_nl(3, :)/params_b(3, :)
        A1 = ddr_params_nl(1, :)/params_b(1, :) + e_charge*Ercov/params_b(3, :) - 1.5d0*A2
!
        CALL h5_add_double_1(h5_id, trim(tempch)//"A1", &
            A1, lbound(A1), ubound(A1))
        CALL h5_add_double_1(h5_id, trim(tempch)//"A2", &
            A2, lbound(A2), ubound(A2))

        do i = 1, npoib
            call getIfunc(x1(i), x2(i), symbI(:, :, i))
        end do
        curr_e_par = e_charge*params_b(1, :)*vT/(nue*B0) &
                    *(c*Es*((A1 + A2)*symbI(1, 0, :) + 0.5d0*A2*symbI(2, 1, :)) &
                   + vT*Br*((A1 + A2)*symbI(1, 1, :) + 0.5d0*A2*symbI(3, 1, :)))

        call get_current_densities_from_wave_code(flre_cd_ptr(1), dim_r, r, &
                                    m_vals(1), n_vals(1), Jri, Jsi, Jpi, Jre, Jse, Jpe)
!
        CALL h5_add_double_1(h5_id, trim(tempch)//"je_real", & ! drift kinetic current (gyrocurrent)
            real(curr_e_par), lbound(real(curr_e_par)), ubound(real(curr_e_par)))
        CALL h5_add_double_1(h5_id, trim(tempch)//"je_imag", &
            dimag(curr_e_par), lbound(dimag(curr_e_par)), ubound(dimag(curr_e_par))) 

        CALL h5_add_double_1(h5_id, trim(tempch)//"Je_real", & ! kinetic current from KiLCA
            real(Jpe), lbound(real(Jpe)), ubound(real(Jpe)))
        CALL h5_add_double_1(h5_id, trim(tempch)//"Je_imag", &
            dimag(Jpe), lbound(dimag(Jpe)), ubound(dimag(Jpe))) 

        CALL h5_close(h5_id)
        CALL h5_deinit()

        write(*,*) " - - - - - - - - - "
 
    end if

!
    deallocate (x1, x2, symbI, curr_e_par, vT)
!
end subroutine calc_parallel_current_directly
!
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!
subroutine calc_ion_parallel_current_directly

    use grid_mod, only: npoib, rb, Ercov
    use plasma_parameters, only: params_b, ddr_params_nl
    use baseparam_mod, only: Z_i, e_charge, am, p_mass, c, btor, e_mass, ev, rtor
    use wave_code_data
    use control_mod, only: ihdf5IO, diagnostics_output, write_gyro_current
    use h5mod
    use mpi
!
    implicit none
!
    integer :: ierror, irank;
    integer :: ipoi, i, iunit, mnmax
    double precision :: ei_charge
    double precision, dimension(:), allocatable :: dummy, x1, x2, vT
    double complex, dimension(:), allocatable :: curr_i_par
    double complex, dimension(:, :, :), allocatable :: symbI
!
    character(len=1024) :: tempch
!
    iunit = 731
    mnmax = 3
    allocate (x1(npoib), x2(npoib), symbI(0:mnmax, 0:mnmax, npoib))
    allocate (curr_i_par(npoib), vT(npoib))
!
    ei_charge = -Z_i*e_charge
!
!
    vT = sqrt(params_b(4, :)/p_mass)
!
    x1 = kp*vT/nui
    x2 = -om_E/nui
!
    do i = 1, npoib
        call getIfunc(x1(i), x2(i), symbI(:, :, i))
    end do
!
! Here x1 and x2 are used for A_1 and A_2:
    x2 = ddr_params_nl(4, :)/params_b(4, :)
    x1 = ddr_params_nl(1, :)/params_b(1, :) + ei_charge*Ercov/params_b(4, :) - 1.5d0*x2
!
    curr_i_par = ei_charge*params_b(1, :)*vT/(nui*B0) &
                 *(c*Es*((x1 + x2)*symbI(1, 0, :) + 0.5d0*x2*symbI(2, 1, :)) &
                   + vT*Br*((x1 + x2)*symbI(1, 1, :) + 0.5d0*x2*symbI(3, 1, :)))
!
    call MPI_Comm_rank(MPI_COMM_WORLD, irank, ierror);
    if (irank .eq. 0) then
        if (write_gyro_current) then
            write(*,*) "writing par_current_i.dat"
            write(*,*) " - "
            ! Write out gyro current which is different to KiLCA current Jpe.
            ! The gyro current is calculated from (60) in Heyn et. al 2014
            CALL h5_init()
            CALL h5_open_rw(path2out, h5_id)
            tempch = "/"//trim(h5_mode_groupname)//"/par_current_i/"
            write(*,*) "In group: "//trim(tempch)

            CALL h5_define_group(h5_id, trim(tempch), group_id_1)
            CALL h5_close_group(group_id_1)

                !ALL h5_add_double_1(h5_id, trim(tempch)//"rb", &
                !   rb, lbound(rb), ubound(rb))
                !ALL h5_add_double_1(h5_id, trim(tempch)//"par_current_i_real", &
                !   real(curr_i_par), lbound(real(curr_i_par)), ubound(real(curr_i_par)))
                !ALL h5_add_double_1(h5_id, trim(tempch)//"par_current_i_imag", &
                !   dimag(curr_i_par), lbound(dimag(curr_i_par)), ubound(dimag(curr_i_par))) 
                !ALL h5_add_double_1(h5_id, trim(tempch)//"Jpi_real", &
                !   real(Jpi), lbound(real(Jpi)), ubound(real(Jpi)))
                !ALL h5_add_double_1(h5_id, trim(tempch)//"Jpi_imag", &
                !   dimag(Jpi), lbound(dimag(Jpi)), ubound(dimag(Jpi))) 
                !ALL h5_add_double_1(h5_id, trim(tempch)//"Jsi_real", &
                !   real(Jsi), lbound(real(Jsi)), ubound(real(Jsi)))
                !ALL h5_add_double_1(h5_id, trim(tempch)//"Jsi_imag", &
                !   dimag(Jsi), lbound(dimag(Jsi)), ubound(dimag(Jsi))) 
                !ALL h5_add_double_1(h5_id, trim(tempch)//"Jri_real", &
                !   real(Jpi), lbound(real(Jpi)), ubound(real(Jpi)))
                !ALL h5_add_double_1(h5_id, trim(tempch)//"Jri_imag", &
                !   dimag(Jri), lbound(dimag(Jri)), ubound(dimag(Jri))) 
                CALL h5_add_double_1(h5_id, trim(tempch)//"kp", &
                   kp, lbound(kp), ubound(kp))
                CALL h5_add_double_1(h5_id, trim(tempch)//"ks", &
                    ks, lbound(ks), ubound(ks)) 
                CALL h5_add_double_1(h5_id, trim(tempch)//"x1", &
                    x1, lbound(x1), ubound(x1))
                CALL h5_add_double_1(h5_id, trim(tempch)//"x2", &
                    x2, lbound(x2), ubound(x2)) 
                CALL h5_add_double_1(h5_id, trim(tempch)//"nui", &
                    nui, lbound(nui), ubound(nui)) 




                CALL h5_close(h5_id)
                CALL h5_deinit()

            end if ! write_gyro_current



        if (diagnostics_output) then
            if (ihdf5IO .eq. 1) then
                print *, "writing par_current_i.dat"
                CALL h5_init()
                CALL h5_open_rw(path2out, h5_id)
                tempch = "/"//trim(h5_mode_groupname)//"/par_current_i.dat"

                CALL h5_obj_exists(h5_id, trim(tempch), h5_exists_log)
                if (h5_exists_log) then
                    CALL h5_delete(h5_id, trim(tempch))
                end if

                CALL h5_define_unlimited_matrix(h5_id, trim(tempch), &
                                                H5T_NATIVE_DOUBLE, (/-1, 5/), dataset_id)
                CALL h5_append_double_1(dataset_id, rb, 1)
                CALL h5_append_double_1(dataset_id, real(curr_i_par), 2)
                CALL h5_append_double_1(dataset_id, dimag(curr_i_par), 3)
                CALL h5_append_double_1(dataset_id, real(Jpi), 4)
                CALL h5_append_double_1(dataset_id, dimag(Jpi), 5)

                CALL h5_close(h5_id)
                CALL h5_deinit()

            else
                open (iunit, file='par_current_i.dat')
                do ipoi = 1, npoib
                    write (iunit, *) rb(ipoi), real(curr_i_par(ipoi)), dimag(curr_i_par(ipoi)) &
                        , real(Jpi(ipoi)), dimag(Jpi(ipoi))
                end do
                close (iunit)
            end if
        end if
    end if
!
    deallocate (x1, x2, symbI, curr_i_par, vT)
!
end subroutine calc_ion_parallel_current_directly

