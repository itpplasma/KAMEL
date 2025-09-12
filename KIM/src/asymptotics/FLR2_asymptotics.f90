module flr2_asymptotics_m

    ! module to calculate asymptotic quantities from the finite Larmor radius expansion
    ! model to second order (FLR2) based on the formulation of KiLCA

    implicit none

    contains

    subroutine calc_flr2_asymptotic_Phi_MA(plasma_in, EBdat)
        ! calculate the asymptotic form of the misalignment potential Phi_MA
        ! in the ideal MHD bulk region

        use KIM_kinds_m, only: dp
        use species_m, only: plasma_t
        use fields_m, only: EBdat_t
        use constants_m, only: pi, com_unit, sol, e_charge, ev
        use grid_m, only: xl_grid
        use equilibrium_m, only: B0
        use IO_collection_m, only: write_complex_profile_abs
        use config_m, only: output_path

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

                if (sp == 1) cycle
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
                if (sp == 1) cycle
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

        call write_complex_profile_abs(xl_grid%xb, offdiag, xl_grid%npts_b, trim(output_path)//"/fields/offdiag.dat")

    end subroutine


    subroutine calc_hatK_Phi_in_Fourier(plasma_in)

        use KIM_kinds_m, only: dp
        use species_m, only: plasma_t
        use constants_m, only: pi, com_unit, sol, e_charge, ev
        use grid_m, only: rg_grid
        use equilibrium_m, only: B0
        use IO_collection_m, only: write_complex_profile_abs
        use config_m, only: output_path
        use gsl_mod, only: gsl_sf_bessel_In
        use config_m, only: turn_off_ions, artificial_debye_case, turn_off_electrons

        implicit none

        type(plasma_t), intent(in) :: plasma_in
        integer :: j, sp, i
        complex(dp), allocatable :: kernel_phi(:)
        complex(dp), allocatable :: kernel_B(:)
        complex(dp) :: kernel_phi_temp
        real(dp) :: b
        real(dp) :: ks
        real(dp) :: kr
        real(dp) :: kr_arr(3)
        character(256) :: filename

        complex(dp) :: besselI ! complex bessel function from bessel.f90
        allocate(kernel_phi(rg_grid%npts_b))
        allocate(kernel_B(rg_grid%npts_b))

        kr_arr = [1.0d0, 5.0d0, 10.0d0]

        do i = 1, size(kr_arr)
            kr = kr_arr(i)
            print *, "Calculating hatK_Phi for kr = ", kr
            kernel_phi = 0.0d0
            kernel_B = 0.0d0
            
            do j = 1, size(rg_grid%xb)
                kernel_phi_temp = 0.0d0
                do sp = 0, plasma_in%n_species-1
                    if (turn_off_ions .and. sp >= 1) cycle
                    if (turn_off_electrons .and. sp == 0) cycle

                    ! Include full perpendicular wavenumber in FLR parameter: b = (k_r^2 + k_s^2) * rho_T^2
                    ks = plasma_in%ks(j)
                    ! b = (kr**2.0d0 + ks**2.0d0) * plasma_in%spec(sp)%rho_L(j)**2.0d0
                    b = kr**2.0d0 * plasma_in%spec(sp)%rho_L(j)**2.0d0

                    if (artificial_debye_case <= 1) then
                        kernel_phi_temp = - 1.0d0 / plasma_in%spec(sp)%lambda_D(j)**2.0d0
                    end if

                    if (artificial_debye_case == 0 .or. artificial_debye_case == 2) then
                        kernel_phi_temp = kernel_phi_temp + 1.0d0 / plasma_in%spec(sp)%lambda_D(j)**2.0d0 * com_unit * plasma_in%spec(sp)%vT(j)**2.0d0 * plasma_in%ks(j) &
                            / (plasma_in%spec(sp)%omega_c(j) * plasma_in%spec(sp)%nu(j)) * exp(-b) * &
                            (&
                                plasma_in%spec(sp)%I00(j) * (&
                                    gsl_sf_bessel_In(0, b) * (plasma_in%spec(sp)%A1(j) + plasma_in%spec(sp)%A2(j) * (1-b)) &
                                    + 0.5d0 * plasma_in%spec(sp)%A2(j) * b * gsl_sf_bessel_In(-1, b) &
                                )&
                                + 0.5d0 * plasma_in%spec(sp)%I20(j) * plasma_in%spec(sp)%A2(j) * gsl_sf_bessel_In(0, b) &
                            )
                        kernel_B(j) = kernel_B(j) - 1.0d0 / plasma_in%spec(sp)%lambda_D(j)**2.0d0 * plasma_in%spec(sp)%vT(j)**3.0d0 &
                            / (plasma_in%spec(sp)%omega_c(j) * plasma_in%spec(sp)%nu(j) * sol) * exp(-b) * &
                            (&
                                plasma_in%spec(sp)%I01(j) * (&
                                    gsl_sf_bessel_In(0, b) * (plasma_in%spec(sp)%A1(j) + plasma_in%spec(sp)%A2(j) * (1-b)) &
                                    + 0.5d0 * plasma_in%spec(sp)%A2(j) * b * gsl_sf_bessel_In(-1, b) &
                                )&
                                + 0.5d0 * plasma_in%spec(sp)%I21(j) * plasma_in%spec(sp)%A2(j) * gsl_sf_bessel_In(0, b) &
                            )
                    end if

                    kernel_phi(j) = kernel_phi(j) + kernel_phi_temp
                end do
            end do

            kernel_phi = kernel_phi / (4.0d0 * pi)
            kernel_B = kernel_B / (4.0d0 * pi)

            write(filename, '(A,I0,A)') trim(output_path)//"/fields/hatK_Phi_kr", int(kr), ".dat"
            call write_complex_profile_abs(rg_grid%xb, kernel_phi, rg_grid%npts_b, filename)
            write(filename, '(A,I0,A)') trim(output_path)//"/fields/hatK_B_kr", int(kr), ".dat"
            call write_complex_profile_abs(rg_grid%xb, kernel_B, rg_grid%npts_b, filename)

        end do

    end subroutine

end module
