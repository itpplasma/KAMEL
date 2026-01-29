MODULE Function_Input_Module
    ! File copied from the ZEAL package. Adapted for KIM.
    !
    ! Module for dispersion function input.
    ! Keep this name, since ZEAL uses it.

    USE Precision_Module, only: DP

    IMPLICIT NONE

    abstract interface
        subroutine dispersion_func_interface(kr, f)
            use KIM_kinds_m, only: dp
            import
            complex(kind=dp), intent(in)  :: kr
            complex(kind=dp), intent(out) :: f
        end subroutine dispersion_func_interface
    end interface

    procedure(dispersion_func_interface), pointer :: f_ptr => null()

    ! Function pointer for dispersion implementation (KIM or FLRE)
    ! Set once via init_dispersion_mode, avoids string comparison per call
    procedure(dispersion_func_interface), pointer :: dispersion_impl => null()

    !-----------------------------------------------------------------------
    ! Documentation from ZEAL:
    !
    ! > If ICON = 3 or 4 (as specified in Zeal_Input_Module), then specify
    ! > the values of NEWTONZ and NEWTONF.
    ! > These variables are used as follows. The modified Newton's method,
    ! > which takes into account the multiplicity of a zero and converges
    ! > quadratically, is used to refine the calculated approximations for
    ! > the zeros. The iteration stops if
    ! >   - the relative distance between two successive approximations is
    ! >     at most NEWTONZ
    ! > or
    ! >   - the absolute value of the function at the last approximation is
    ! >     at most NEWTONF
    ! > or if a maximum number of iterations is exceeded.
    !
     real(dp) :: NEWTONZ = 5.0d-08
     real(dp) :: NEWTONF = 1.0d-14

     integer :: rg_index

     CONTAINS

        subroutine init_dispersion_mode()
            ! Initialize the dispersion function pointer based on WKB_dispersion_mode.
            ! Call once before running the dispersion solver.
            use config_m, only: WKB_dispersion_mode

            select case (trim(WKB_dispersion_mode))
                case ('KIM')
                    dispersion_impl => dispersion_KIM
                case ('FLRE')
                    dispersion_impl => dispersion_FLRE
                case default
                    print *, 'Error: Unknown WKB_dispersion_mode: ', trim(WKB_dispersion_mode)
                    stop
            end select
        end subroutine

        subroutine FDF(kr, D_of_kr, dD_dkr)
            ! function that returns the dispersion function and its derivative with complex k_r as the input
            ! Uses numerical derivative (finite difference) since analytical formula had errors

            use kim_kinds_m, only: dp

            implicit none
            complex(dp), intent(in) :: kr
            complex(dp), intent(out) :: D_of_kr, dD_dkr
            complex(dp) :: D_plus, D_minus
            real(dp), parameter :: eps = 1.0d-7

            ! Compute function value at kr
            call dispersion_function(kr, D_of_kr)

            ! Compute numerical derivative using central difference
            call dispersion_function(kr + cmplx(eps, 0.0d0, dp), D_plus)
            call dispersion_function(kr - cmplx(eps, 0.0d0, dp), D_minus)
            dD_dkr = (D_plus - D_minus) / (2.0d0 * eps)

        end subroutine


        subroutine dispersion_function(kr, D_of_kr)
            ! Calls the dispersion implementation via function pointer.
            ! init_dispersion_mode() must be called first.

            use kim_kinds_m, only: dp

            implicit none
            complex(dp), intent(in) :: kr
            complex(dp), intent(out) :: D_of_kr

            if (.not. associated(dispersion_impl)) then
                print *, 'Error: dispersion_impl not initialized. Call init_dispersion_mode() first.'
                stop
            end if

            call dispersion_impl(kr, D_of_kr)

        end subroutine


        subroutine dispersion_KIM(kr, D_of_kr)
            ! KIM dispersion: full modified Bessel function formulation

            use species_m, only: plasma
            use constants_m, only: com_unit
            use kim_kinds_m, only: dp

            implicit none
            complex(dp), intent(in) :: kr
            complex(dp), intent(out) :: D_of_kr
            integer :: sp
            integer :: ifail_bess, nz
            complex(dp) :: bess0, bessm1
            real(dp) :: bess0_re(2), bess0_im(2)
            complex(dp) :: kr_rho_squared
            integer :: j

            if (rg_index .lt. 1 .or. rg_index .gt. plasma%grid_size) then
                print *, 'Error: rg_index out of bounds in dispersion function calculation'
                print *, '  rg_index =', rg_index, ' grid_size =', plasma%grid_size
                stop
            end if

            j = rg_index

            bess0_re = 0.0d0
            bess0_im = 0.0d0
            D_of_kr = (0.0d0, 0.0d0)

            do sp = 0, plasma%n_species-1

                kr_rho_squared = kr**2.0d0 * plasma%spec(sp)%rho_L(j)**2.0d0

                ! calculate modified Bessel functions I0 and I1 with complex argument
                call zbesi(real(kr_rho_squared, kind=dp), &
                           dimag(kr_rho_squared), &
                           0.0d0, &! initial order
                           2, &! KODE, 1= no exponential scaling, 2= with scaling
                           2, &! n  number of terms (I0 and I1)
                           bess0_re, &
                           bess0_im, &
                           nz, & ! number of underflows set to zero
                           ifail_bess)

                if (ifail_bess == 1) then
                    print *, 'Warning: Bessel function calculation, I0 and I1, bad input - no computation'
                else if (ifail_bess == 2) then
                    print *, 'Warning: Bessel function calculation, I0 and I1, overflow occurred - no computation'
                else if (ifail_bess == 3) then
                    print *, 'Warning: Bessel function calculation, I0 and I1, precision warning - computation completed'
                else if (ifail_bess == 4) then
                    print *, 'Error: Bessel function calculation, I0 and I1, precision error - no computation'
                else if (ifail_bess == 5) then
                    print *, 'Error: Bessel function calculation, I0 and I1, algorithmic error - no computation'
                end if

                bess0 = cmplx(bess0_re(1), bess0_im(1), kind=dp) ! zero order Bessel function
                bessm1 = cmplx(bess0_re(2), bess0_im(2), kind=dp) ! first order Bessel function, symmetric for integer order

                D_of_kr = D_of_kr + 1.0d0 / plasma%spec(sp)%lambda_D(j)**2.0d0 * ( -1.0d0 &
                    + com_unit * plasma%spec(sp)%vT(j)**2.0d0 * plasma%ks(j) &
                    / (plasma%spec(sp)%omega_c(j) * plasma%spec(sp)%nu(j)) * &
                    exp(abs(real(kr_rho_squared, kind=dp)) - kr_rho_squared) * & ! correct the scaling of the bessel functions
                    (&
                        plasma%spec(sp)%I00(j, 0) * (&
                            bess0 * (plasma%spec(sp)%A1(j) + plasma%spec(sp)%A2(j) * (1.0d0 - kr_rho_squared)) &
                            + plasma%spec(sp)%A2(j) * kr_rho_squared * bessm1 &
                        )&
                         + 0.5d0 * plasma%spec(sp)%I20(j, 0) * plasma%spec(sp)%A2(j) * bess0 &
                    )&
                )

            end do

            D_of_kr = kr**2.0d0 + plasma%kp(j)**2.0d0 - D_of_kr

        end subroutine


        subroutine dispersion_FLRE(kr, D_of_kr)
            ! FLRE dispersion: second-order finite Larmor radius expansion
            ! Replaces Bessel functions with Taylor expansion: I0(x)*exp(-x) ≈ 1 - x + ...

            use species_m, only: plasma
            use constants_m, only: com_unit
            use kim_kinds_m, only: dp

            implicit none
            complex(dp), intent(in) :: kr
            complex(dp), intent(out) :: D_of_kr
            integer :: sp
            complex(dp) :: kr_rho_squared
            integer :: j

            if (rg_index .lt. 1 .or. rg_index .gt. plasma%grid_size) then
                print *, 'Error: rg_index out of bounds in dispersion function calculation'
                print *, '  rg_index =', rg_index, ' grid_size =', plasma%grid_size
                stop
            end if

            j = rg_index

            D_of_kr = (0.0d0, 0.0d0)

            do sp = 0, plasma%n_species-1

                kr_rho_squared = kr**2.0d0 * plasma%spec(sp)%rho_L(j)**2.0d0

                ! FLRE: replace exp(-x)*I0(x) ≈ 1 - x, exp(-x)*I1(x) ≈ x/2
                ! This simplifies the Bessel terms to polynomial expressions
                D_of_kr = D_of_kr + 1.0d0 / plasma%spec(sp)%lambda_D(j)**2.0d0 * ( -1.0d0 &
                    + com_unit * plasma%spec(sp)%vT(j)**2.0d0 * plasma%ks(j) &
                    / (plasma%spec(sp)%omega_c(j) * plasma%spec(sp)%nu(j)) * &
                    (&
                        plasma%spec(sp)%I00(j, 0) * ((1.0d0 - kr_rho_squared) * &
                            (plasma%spec(sp)%A1(j) + plasma%spec(sp)%A2(j) * (1.0d0 - kr_rho_squared)) &
                        )&
                        + 0.5d0 * plasma%spec(sp)%I20(j, 0) * plasma%spec(sp)%A2(j) * (1.0d0 - kr_rho_squared) &
                    )&
                )
            end do

            D_of_kr = kr**2.0d0 + plasma%kp(j)**2.0d0 - D_of_kr

        end subroutine


        subroutine test_FDF_derivative(kr_test)
            !-----------------------------------------------------------------------
            ! Test that FDF derivative matches independent numerical calculation.
            ! Call this after rg_index is set to a valid grid point.
            !-----------------------------------------------------------------------
            use kim_kinds_m, only: dp

            implicit none

            complex(dp), intent(in) :: kr_test
            complex(dp) :: D_center, dD_from_FDF
            complex(dp) :: D_plus_re, D_minus_re
            complex(dp) :: D_plus_im, D_minus_im
            complex(dp) :: dD_numerical_re, dD_numerical_im
            real(dp) :: eps
            real(dp) :: rel_error_re, rel_error_im

            eps = 1.0d-6

            ! Get derivative from FDF
            call FDF(kr_test, D_center, dD_from_FDF)

            ! Independent numerical derivative using dispersion_function directly
            call dispersion_function(kr_test + cmplx(eps, 0.0d0, dp), D_plus_re)
            call dispersion_function(kr_test - cmplx(eps, 0.0d0, dp), D_minus_re)
            dD_numerical_re = (D_plus_re - D_minus_re) / (2.0d0 * eps)

            ! Numerical derivative in imaginary direction: dD/d(Im(kr))
            call dispersion_function(kr_test + cmplx(0.0d0, eps, dp), D_plus_im)
            call dispersion_function(kr_test - cmplx(0.0d0, eps, dp), D_minus_im)
            dD_numerical_im = (D_plus_im - D_minus_im) / (2.0d0 * eps)

            ! For a holomorphic function: dD/dkr = dD/d(Re(kr)) = -i * dD/d(Im(kr))

            print *
            print *, '=== FDF Derivative Test ==='
            print *, 'Test point kr = (', real(kr_test), ',', aimag(kr_test), ')'
            print *, 'D(kr) = (', real(D_center), ',', aimag(D_center), ')'
            print *
            print *, 'dD/dkr from FDF:'
            print *, '  Real part:      ', real(dD_from_FDF)
            print *, '  Imaginary part: ', aimag(dD_from_FDF)
            print *
            print *, 'Numerical dD/dkr (independent check):'
            print *, '  Real part:      ', real(dD_numerical_re)
            print *, '  Imaginary part: ', aimag(dD_numerical_re)
            print *

            ! Relative errors
            rel_error_re = 0.0d0
            if (abs(dD_from_FDF) > 1.0d-10) then
                rel_error_re = abs(dD_from_FDF - dD_numerical_re) / abs(dD_from_FDF)
                print *, 'Relative error: ', rel_error_re
            else
                print *, 'Absolute error: ', abs(dD_from_FDF - dD_numerical_re)
            end if

            print *
            print *, 'Cauchy-Riemann check (imag perturbation should give i * dD/dkr):'
            print *, '  dD/d(Im(kr)):   (', real(dD_numerical_im), ',', aimag(dD_numerical_im), ')'
            print *, '  i * dD/dkr:     (', -aimag(dD_from_FDF), ',', real(dD_from_FDF), ')'

            if (abs(dD_from_FDF) > 1.0d-10) then
                rel_error_im = abs(dD_numerical_im - cmplx(0.0d0, 1.0d0, dp) * dD_from_FDF) / abs(dD_from_FDF)
                print *, 'Cauchy-Riemann error: ', rel_error_im
            end if

            print *
            if (rel_error_re > 0.01d0) then
                print *, '*** WARNING: Derivative error > 1% ***'
            else if (rel_error_re > 0.001d0) then
                print *, '*** Note: Derivative error > 0.1% ***'
            else
                print *, 'Derivative OK (error < 0.1%)'
            end if
            print *, '=== End Derivative Test ==='
            print *

        end subroutine test_FDF_derivative


        FUNCTION VALREG(LV,H)
            !-----------------------------------------------------------------------
            !**PURPOSE
            !  Given a rectangular region specified by LV and H (cf. the module
            !  Zeal_Input_Module), decide whether the function is analytic inside
            !  this region or not.
            !-----------------------------------------------------------------------
            !
            implicit none

            logical :: VALREG
            REAL(KIND=dp), INTENT(IN) :: LV(2), H(2)

            VALREG = .TRUE.

            !  The following statement can be used for functions that have a
            !  branch cut along the non-positive real axis.
            !
            !    VALREG = .NOT. ( LV(2)*(LV(2)+H(2)) <= ZERO .AND. LV(1) <= ZERO )

        end function VALREG

END MODULE Function_Input_Module
