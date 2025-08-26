module flr2_asymptotics_m

    implicit none


    contains

    subroutine calc_flr2_asymptotic_Phi_MA(plasma_in, EBdat)

        use KIM_kinds_m, only: dp
        use species_m, only: plasma_t
        use fields_m, only: EBdat_t
        use constants_m, only: pi, com_unit, sol, e_charge, ev
        use grid_m, only: xl_grid
        use equilibrium_m, only: B0

        implicit none

        type(plasma_t), intent(in) :: plasma_in
        type(EBdat_t), intent(inout) :: EBdat

        integer :: sp, j
        complex(dp), allocatable :: F0(:), F2(:), H(:)
        real(dp), allocatable :: inv_lambda_tot_squared(:)

        integer :: nlagr = 4
        integer :: nder = 0
        integer :: ibeg, iend, ir
        real(dp), dimension(:,:), allocatable :: coef

        if (.not. allocated(coef)) allocate(coef(0:nder, nlagr))

        allocate(F0(size(EBdat%Br)))
        allocate(F2(size(EBdat%Br)))
        allocate(H(size(EBdat%Br)))
        allocate(inv_lambda_tot_squared(size(EBdat%Br)))
        if (.not. allocated(EBdat%Phi_MA_asymptotic)) allocate(EBdat%Phi_MA_asymptotic(size(EBdat%Br)))

        F0 = 0.0d0
        H = 0.0d0
        inv_lambda_tot_squared = 0.0d0

        do sp = 0, plasma_in%n_species-1
            do j = 1, size(EBdat%r_grid)
                call binsrc(plasma_in%r_grid, 1, size(plasma_in%r_grid), EBdat%r_grid(j), ir) 
                ibeg = max(1, ir - nlagr/2)
                iend = ibeg + nlagr - 1
                if (iend .gt. size(plasma_in%r_grid)) then
                    iend = size(plasma_in%r_grid)
                    ibeg = iend -nlagr + 1
                end if

                call plag_coeff(nlagr, nder, EBdat%r_grid(j), plasma_in%r_grid(ibeg:iend), coef)

                inv_lambda_tot_squared(j) = 1.0d0 / sum(coef(0,:) * plasma_in%spec(sp)%lambda_D(ibeg:iend))**2.0d0
            end do
        end do

        ! this uses expressions from FLR2 for the asymptotics far away from the resonant surface
        do j = 1, size(EBdat%Br)
            do sp = 0, plasma_in%n_species-1
                F0(j) = F0(j) + inv_lambda_tot_squared(j)
                H(j) = H(j) + plasma_in%spec(sp)%rho_L(j)**2.0d0 / (inv_lambda_tot_squared(j)) &
                    * (1.0d0 - ev * (plasma_in%spec(sp)%dndr(j) * plasma_in%spec(sp)%T(j)&
                        + plasma_in%spec(sp)%n(j) * plasma_in%spec(sp)%dTdr(j)) / &
                        (plasma_in%spec(sp)%Zspec * e_charge * plasma_in%spec(sp)%n(j) &
                        * plasma_in%Er(j)) )
            end do
            H(j) = H(j) / (8.0d0 * pi)
            F2(j) = - com_unit * plasma_in%om_E(j) * plasma_in%Er(j) / plasma_in%kp(j)**2.0d0 * H(j)
            F0(j) = F0(j) * (- com_unit) * (sol * plasma_in%om_E(j) * plasma_in%Er(j) / (4.0d0 * pi *plasma_in%kp(j)**2.0d0))
        end do

        EBdat%Phi_MA_asymptotic = 0.0d0
        !do j = 1, xl_grid%npts_c
            !print *, xl_grid%ipbeg(j), xl_grid%ipend(j)
            !EBdat%Phi_MA_asymptotic(j) = plasma_in%Er(j) * F2(j) / (com_unit * plasma_in%kp(j) * F0(j)) &
                !* sum(EBdat%Br(xl_grid%ipbeg(j):xl_grid%ipend(j)) / B0(xl_grid%ipbeg(j):xl_grid%ipend(j)) &
                    !* xl_grid%deriv2_coef(:,j)) &
                !+ 1.0d0 / F0(j) * ( &
                        !com_unit * F2(j) - plasma_in%om_E(j) * plasma_in%Er(j) * H(j)/ plasma_in%kp(j)**2.0d0 &
                    !) &
                    !* sum(EBdat%Br(xl_grid%ipbeg(j):xl_grid%ipend(j)) * plasma_in%Er(xl_grid%ipbeg(j):xl_grid%ipend(j)) &
                        !/ (plasma_in%kp(xl_grid%ipbeg(j):xl_grid%ipend(j)) * B0(xl_grid%ipbeg(j):xl_grid%ipend(j)))&
                        !* xl_grid%reint_coef(:,j)) &
                !- plasma_in%om_E(j) * plasma_in%Er(j) / (plasma_in%kp(j)**2.0d0 * F0(j)) &
                    !* sum(EBdat%Br(xl_grid%ipbeg(j):xl_grid%ipend(j)) * plasma_in%Er(xl_grid%ipbeg(j):xl_grid%ipend(j)) &
                        !* H(xl_grid%ipbeg(j):xl_grid%ipend(j))  &
                        !/ (B0(xl_grid%ipbeg(j):xl_grid%ipend(j)) * plasma_in%kp(xl_grid%ipbeg(j):xl_grid%ipend(j))) &
                        !* xl_grid%deriv_coef(:,j))
        !end do
        !EBdat%Phi_MA_asymptotic = EBdat%Phi_MA_asymptotic * sol

        do j = 1, xl_grid%npts_c
            EBdat%Phi_MA_asymptotic(j) = - 4.0d0 * pi * com_unit / inv_lambda_tot_squared(j) &
                * (&
                    H(j) * plasma_in%Er(j) / plasma_in%kp(j) &
                    * sum(EBdat%Br(xl_grid%ipbeg(j):xl_grid%ipend(j)) / B0(xl_grid%ipbeg(j):xl_grid%ipend(j)) &
                        * xl_grid%deriv2_coef(:,j)) &
                    + sum(EBdat%Br(xl_grid%ipbeg(j):xl_grid%ipend(j)) * plasma_in%Er(xl_grid%ipbeg(j):xl_grid%ipend(j)) &
                        * H(xl_grid%ipbeg(j):xl_grid%ipend(j))  &
                        / (plasma_in%kp(xl_grid%ipbeg(j):xl_grid%ipend(j)) * B0(xl_grid%ipbeg(j):xl_grid%ipend(j)))&
                        * xl_grid%reint_coef(:,j)) &
                )
        end do

    end subroutine

end module