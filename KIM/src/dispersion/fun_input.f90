MODULE Function_Input_Module
    ! Dispersion function input for the WKB solver.
    !
    ! Holds the KIM and FLRE dispersion relations, the active-mode pointer used
    ! by the Muller solver, and the fortnum complex_root_fn_t adapter used by the
    ! region-root solver.

    use KIM_kinds_m, only: dp

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
            ! When WKB_solve_for_kr_squared=.false., kr is the unknown (solve D(kr)=0)
            ! When WKB_solve_for_kr_squared=.true., kr represents kr^2 (solve D(kr^2)=0)

            use species_m, only: plasma
            use constants_m, only: com_unit
            use kim_kinds_m, only: dp
            use config_m, only: WKB_solve_for_kr_squared

            implicit none
            complex(dp), intent(in) :: kr
            complex(dp), intent(out) :: D_of_kr
            integer :: sp
            integer :: ifail_bess, nz
            complex(dp) :: bess0, bessm1
            real(dp) :: bess0_re(2), bess0_im(2)
            complex(dp) :: kr_rho_squared
            complex(dp) :: kr_squared  ! Either kr^2 or kr depending on mode
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

            ! Determine kr^2 based on solver mode
            if (WKB_solve_for_kr_squared) then
                kr_squared = kr  ! Input is already kr^2
            else
                kr_squared = kr**2.0d0  ! Input is kr, compute kr^2
            end if

            do sp = 0, plasma%n_species-1

                kr_rho_squared = kr_squared * plasma%spec(sp)%rho_L(j)**2.0d0

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

            D_of_kr = kr_squared + plasma%kp(j)**2.0d0 - D_of_kr

        end subroutine


        subroutine dispersion_FLRE(kr, D_of_kr)
            ! FLRE dispersion: second-order finite Larmor radius expansion
            ! Replaces Bessel functions with Taylor expansion: I0(x)*exp(-x) ≈ 1 - x + ...
            ! When WKB_solve_for_kr_squared=.false., kr is the unknown (solve D(kr)=0)
            ! When WKB_solve_for_kr_squared=.true., kr represents kr^2 (solve D(kr^2)=0)

            use species_m, only: plasma
            use constants_m, only: com_unit
            use kim_kinds_m, only: dp
            use config_m, only: WKB_solve_for_kr_squared

            implicit none
            complex(dp), intent(in) :: kr
            complex(dp), intent(out) :: D_of_kr
            integer :: sp
            complex(dp) :: kr_rho_squared
            complex(dp) :: kr_squared  ! Either kr^2 or kr depending on mode
            integer :: j

            if (rg_index .lt. 1 .or. rg_index .gt. plasma%grid_size) then
                print *, 'Error: rg_index out of bounds in dispersion function calculation'
                print *, '  rg_index =', rg_index, ' grid_size =', plasma%grid_size
                stop
            end if

            j = rg_index

            D_of_kr = (0.0d0, 0.0d0)

            ! Determine kr^2 based on solver mode
            if (WKB_solve_for_kr_squared) then
                kr_squared = kr  ! Input is already kr^2
            else
                kr_squared = kr**2.0d0  ! Input is kr, compute kr^2
            end if

            do sp = 0, plasma%n_species-1

                kr_rho_squared = kr_squared * plasma%spec(sp)%rho_L(j)**2.0d0

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

            D_of_kr = kr_squared + plasma%kp(j)**2.0d0 - D_of_kr

        end subroutine


        subroutine dispersion_region_fn(kr, fk, ctx)
            ! fortnum complex_root_fn_t adapter: evaluates the active dispersion
            ! relation at the grid point selected by rg_index. The region-root
            ! finder ignores ctx; the grid index rides on module state, set by
            ! the caller before each search (as the Muller path does via rg_index).

            use kim_kinds_m, only: dp

            implicit none
            complex(dp), intent(in)            :: kr
            complex(dp), intent(out)           :: fk
            class(*),    intent(in), optional  :: ctx

            call dispersion_function(kr, fk)

        end subroutine dispersion_region_fn

END MODULE Function_Input_Module
