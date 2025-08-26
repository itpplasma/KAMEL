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
        use IO_collection_m, only: write_complex_profile_abs

        implicit none

        type(plasma_t), intent(in) :: plasma_in
        type(EBdat_t), intent(inout) :: EBdat

        integer :: sp, j
        complex(dp), allocatable :: H(:), F2(:)
        real(dp), allocatable :: rhoL(:, :), B0_intp(:), n(:, :), Er(:)
        real(dp), allocatable :: om_E(:), kp(:), dpdr(:, :), lambda(:, :)
        real(dp), allocatable :: inv_lambda_tot_squared(:)
        complex(dp), allocatable :: offdiag(:)
        real(dp), allocatable :: A1(:), A2(:), I11(:), I13(:), nu(:), vT(:)

        integer :: nlagr = 4
        integer :: nder = 0
        integer :: ibeg, iend, ir
        real(dp), dimension(:,:), allocatable :: coef

        if (.not. allocated(coef)) allocate(coef(0:nder, nlagr))

        allocate(H(size(EBdat%Br)))
        allocate(F2(size(EBdat%Br)))
        allocate(offdiag(size(EBdat%Br)))

        allocate(rhoL(0:plasma_in%n_species-1, size(EBdat%Br)))
        allocate(dpdr(0:plasma_in%n_species-1, size(EBdat%Br)))
        allocate(lambda(0:plasma_in%n_species-1, size(EBdat%Br)))
        allocate(n(0:plasma_in%n_species-1, size(EBdat%Br)))
        allocate(A1(size(EBdat%Br)))
        allocate(A2(size(EBdat%Br)))
        allocate(I11(size(EBdat%Br)))
        allocate(I13(size(EBdat%Br)))
        allocate(nu(size(EBdat%Br)))
        allocate(vT(size(EBdat%Br)))

        allocate(B0_intp(size(EBdat%Br)))
        allocate(Er(size(EBdat%Br)))
        allocate(om_E(size(EBdat%Br)))
        allocate(kp(size(EBdat%Br)))
        allocate(inv_lambda_tot_squared(size(EBdat%Br)))

        allocate(EBdat%Phi_MA_asymptotic(size(EBdat%Br)))

        H = 0.0d0
        inv_lambda_tot_squared = 0.0d0
        
        do j = 1, size(EBdat%r_grid) ! r_grid of EBdat is identical to xl_grid%xc
            
            call binsrc(plasma_in%r_grid, 1, size(plasma_in%r_grid), EBdat%r_grid(j), ir) 
            ibeg = max(1, ir - nlagr/2)
            iend = ibeg + nlagr - 1
            if (iend .gt. size(plasma_in%r_grid)) then
                iend = size(plasma_in%r_grid)
                ibeg = iend -nlagr + 1
            end if

            call plag_coeff(nlagr, nder, EBdat%r_grid(j), plasma_in%r_grid(ibeg:iend), coef)

            do sp = 0, plasma_in%n_species-1
                lambda(sp, j) = sum(coef(0,:) * plasma_in%spec(sp)%lambda_D(ibeg:iend))
                rhoL(sp, j) = sum(coef(0,:) * plasma_in%spec(sp)%rho_L(ibeg:iend))
                dpdr(sp, j) = sum(coef(0,:) * ev * (plasma_in%spec(sp)%dndr(ibeg:iend) * plasma_in%spec(sp)%T(ibeg:iend) &
                    + plasma_in%spec(sp)%n(ibeg:iend) * plasma_in%spec(sp)%dTdr(ibeg:iend)))                
                n(sp, j) = sum(coef(0,:) * plasma_in%spec(sp)%n(ibeg:iend))
                inv_lambda_tot_squared(j) = inv_lambda_tot_squared(j) + 1.0d0 / lambda(sp, j)**2.0d0

                A1(j) = sum(coef(0,:) * plasma_in%spec(sp)%A1(ibeg:iend))
                A2(j) = sum(coef(0,:) * plasma_in%spec(sp)%A2(ibeg:iend))
                I11(j) = sum(coef(0,:) * plasma_in%spec(sp)%I11(ibeg:iend))
                I13(j) = sum(coef(0,:) * plasma_in%spec(sp)%I13(ibeg:iend))
                nu(j) = sum(coef(0,:) * plasma_in%spec(sp)%nu(ibeg:iend))
                vT(j) = sum(coef(0,:) * plasma_in%spec(sp)%vT(ibeg:iend))

                F2(j) = F2(j) + plasma_in%spec(sp)%Zspec* e_charge * n(sp, j) * vT(j)**2.0d0 &
                    * rhoL(sp, j)**2.0d0 / (2.0d0 * nu(j)) &
                    * ((A1(j) + 2.0d0 * A2(j)) * I11(j) + 0.5d0 * A2(j) * I13(j))
            end do

            B0_intp(j) = sum(coef(0,:) * B0(ibeg:iend))
            Er(j) = sum(coef(0,:) * plasma_in%Er(ibeg:iend))
            om_E(j) = sum(coef(0,:) * plasma_in%om_E(ibeg:iend))
            kp(j) = sum(coef(0,:) * plasma_in%kp(ibeg:iend))

        end do

        ! this uses expressions from FLR2 for the asymptotics far away from the resonant surface
        do j = 1, size(EBdat%Br)
            do sp = 0, plasma_in%n_species-1
                H(j) = H(j) + rhoL(sp, j)**2.0d0 / lambda(sp, j)**2.0d0 &
                    * (1.0d0 - (&
                        dpdr(sp, j) / (plasma_in%spec(sp)%Zspec * e_charge * n(sp, j) * Er(j)) &
                        ))
            end do
            H(j) = H(j) / (8.0d0 * pi)
        end do

        EBdat%Phi_MA_asymptotic = 0.0d0

        do j = 1, xl_grid%npts_c
            EBdat%Phi_MA_asymptotic(j) = - 4.0d0 * pi * com_unit / inv_lambda_tot_squared(j) &
                * (&
                    H(j) * Er(j) / kp(j) &
                    * sum(EBdat%Br(xl_grid%ipbeg(j):xl_grid%ipend(j)) / B0_intp(xl_grid%ipbeg(j):xl_grid%ipend(j)) &
                        * xl_grid%deriv2_coef(:,j)) &
                    + sum(EBdat%Br(xl_grid%ipbeg(j):xl_grid%ipend(j)) * Er(xl_grid%ipbeg(j):xl_grid%ipend(j)) &
                        * H(xl_grid%ipbeg(j):xl_grid%ipend(j))  &
                        / (kp(xl_grid%ipbeg(j):xl_grid%ipend(j)) * B0_intp(xl_grid%ipbeg(j):xl_grid%ipend(j)))&
                        * xl_grid%deriv2_coef(:,j)) &
                )

            offdiag(j) = H(j) - com_unit * kp(j)**2.0d0/(om_E(j)*Er(j)) * F2(j)

        end do

        call write_complex_profile_abs(xl_grid%xb, offdiag, xl_grid%npts_b, trim("output/fields/offdiag.dat"))

    end subroutine

end module
