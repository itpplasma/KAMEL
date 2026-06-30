!> Incompressible, flowless ideal-MHD basis calculation for imhd_zone version 0
!> (formerly incompressible.{h,cpp}). The "_gsl" suffix in the oracle's
!> function names is stale: the C++ already used fortnum's callback-based
!> rk8pd/deriv_central, not GSL. This translation upgrades those calls to
!> fortnum's native Fortran interfaces (kilca_background_data_m's own
!> calculate_equilibrium already set this precedent for the same rk8pd
!> stepper).
!>
!> Ffunc/Gfunc are defined here (rather than in kilca_imhd_zone_m, where the
!> C++ oracle keeps them) because kilca_compressible_flow_m also needs them
!> and kilca_imhd_zone_m needs kilca_compressible_flow_m's
!> calculate_basis_flow: putting Ffunc/Gfunc in kilca_imhd_zone_m would make
!> that module and this one depend on each other in both directions, which
!> Fortran's `use` graph cannot express. All routines here take a bare
!> class(zone_t) since neither Ffunc/Gfunc nor the basis calculation touch
!> anything imhd_zone_t adds over the base type (imhd_zone has no data
!> members of its own in the C++ oracle either).
module kilca_incompressible_m
    use, intrinsic :: iso_c_binding, only: c_double, c_null_ptr
    use constants, only: dp, pi
    use kilca_zone_m, only: zone_t, BOUNDARY_CENTER, BOUNDARY_IDEALWALL
    use kilca_background_data_m, only: eval_Bt_Bz, eval_Bt, eval_mass_density, &
        eval_Bt_dBt_Bz_dBz
    use fortnum_status, only: fortnum_status_t, FORTNUM_OK
    use fortnum_ode_rk8pd, only: rk8pd_state_t, rk8pd_evolve_init, rk8pd_evolve_apply
    use fortnum_multiroot, only: deriv_central
    implicit none
    private

    public :: zone_ctx_t, zone_part_ctx_t
    public :: Ffunc, Gfunc, signum
    public :: calculate_basis_incompressible, state_to_EB_incompressible
    public :: rhs_incompressible

    complex(dp), parameter :: cmplx_zero = (0.0_dp, 0.0_dp)
    complex(dp), parameter :: cmplx_i = (0.0_dp, 1.0_dp)
    real(dp), parameter :: deriv_h = 1.0e-3_dp

    !> deriv_central/rk8pd callback context: the zone being integrated.
    type :: zone_ctx_t
        class(zone_t), pointer :: zone => null()
    end type zone_ctx_t

    !> deriv_central context for the real/imag parts of a complex coefficient
    !> function (matches the oracle's diff_params).
    type :: zone_part_ctx_t
        class(zone_t), pointer :: zone => null()
        integer :: part = 0
    end type zone_part_ctx_t

    interface
        function get_background_rtor() bind(C, name="get_background_rtor_") result(rtor)
            import :: c_double
            real(c_double) :: rtor
        end function get_background_rtor
    end interface

    external :: normalize_imhd_basis

contains

    real(dp) function Ffunc(zone, r) result(F)
        class(zone_t), pointer, intent(in) :: zone
        real(dp), intent(in) :: r
        real(dp) :: kt, kz, BtBz(2)
        kt = real(zone%wd%m, dp)/r
        kz = real(zone%wd%n, dp)/get_background_rtor()
        call eval_Bt_Bz(r, c_null_ptr, BtBz)
        F = kt*BtBz(1) + kz*BtBz(2)
    end function Ffunc

    real(dp) function Gfunc(zone, r) result(G)
        class(zone_t), pointer, intent(in) :: zone
        real(dp), intent(in) :: r
        real(dp) :: kt, kz, BtBz(2)
        kt = real(zone%wd%m, dp)/r
        kz = real(zone%wd%n, dp)/get_background_rtor()
        call eval_Bt_Bz(r, c_null_ptr, BtBz)
        G = kt*BtBz(2) - kz*BtBz(1)
    end function Gfunc

    integer function signum(x) result(s)
        real(dp), intent(in) :: x
        if (x < 0.0_dp) then
            s = -1
        else if (x == 0.0_dp) then
            s = 0
        else
            s = 1
        end if
    end function signum

    complex(dp) function Afunc_hi_inc(zone, r) result(A)
        class(zone_t), pointer, intent(in) :: zone
        real(dp), intent(in) :: r
        real(dp) :: kz, mdens(1), F, rk2
        complex(dp) :: omega2, omega2_a

        kz = real(zone%wd%n, dp)/get_background_rtor()
        omega2 = zone%wd%olab*zone%wd%olab
        call eval_mass_density(r, c_null_ptr, mdens)
        F = Ffunc(zone, r)
        omega2_a = cmplx(F*F/(4.0_dp*pi*mdens(1)), 0.0_dp, dp)
        rk2 = r*(kz*kz + real(zone%wd%m, dp)**2/(r*r))
        A = mdens(1)*(omega2 - omega2_a)/rk2
    end function Afunc_hi_inc

    complex(dp) function Bfunc_hi_inc(zone, r) result(B)
        class(zone_t), pointer, intent(in) :: zone
        real(dp), intent(in) :: r
        real(dp) :: kz, mdens(1), F, Bt(1), r3k2, t3, abserr
        complex(dp) :: omega2, omega2_a, t1, t2
        type(zone_ctx_t) :: ctx
        type(fortnum_status_t) :: status

        kz = real(zone%wd%n, dp)/get_background_rtor()
        omega2 = zone%wd%olab*zone%wd%olab
        call eval_mass_density(r, c_null_ptr, mdens)
        F = Ffunc(zone, r)
        call eval_Bt(r, c_null_ptr, Bt)

        omega2_a = cmplx(F*F/(4.0_dp*pi*mdens(1)), 0.0_dp, dp)
        r3k2 = r*r*r*(kz*kz + real(zone%wd%m, dp)**2/(r*r))

        t1 = -mdens(1)*(omega2 - omega2_a)/r
        t2 = 4.0_dp*kz*kz*Bt(1)*Bt(1)/(4.0_dp*pi*r3k2)*(omega2_a/(omega2 - omega2_a))

        ctx%zone => zone
        call deriv_central(func_der, r, deriv_h, t3, abserr, status, ctx)

        B = t1 + t2 + t3
    end function Bfunc_hi_inc

    real(dp) function Afunc_hi_inc_part(r, ctx) result(res)
        real(dp), intent(in) :: r
        class(*), intent(in), optional :: ctx
        complex(dp) :: A
        select type (ctx)
        type is (zone_part_ctx_t)
            A = Afunc_hi_inc(ctx%zone, r)
            if (ctx%part == 0) then
                res = real(A, dp)
            else
                res = aimag(A)
            end if
        end select
    end function Afunc_hi_inc_part

    real(dp) function Bfunc_hi_inc_part(r, ctx) result(res)
        real(dp), intent(in) :: r
        class(*), intent(in), optional :: ctx
        complex(dp) :: B
        select type (ctx)
        type is (zone_part_ctx_t)
            B = Bfunc_hi_inc(ctx%zone, r)
            if (ctx%part == 0) then
                res = real(B, dp)
            else
                res = aimag(B)
            end if
        end select
    end function Bfunc_hi_inc_part

    real(dp) function dAfunc_hi_inc_part(r, ctx) result(der)
        real(dp), intent(in) :: r
        class(*), intent(in), optional :: ctx
        real(dp) :: abserr
        type(fortnum_status_t) :: status
        call deriv_central(Afunc_hi_inc_part, r, deriv_h, der, abserr, status, ctx)
    end function dAfunc_hi_inc_part

    real(dp) function ddAfunc_hi_inc_part(r, ctx) result(der)
        real(dp), intent(in) :: r
        class(*), intent(in), optional :: ctx
        real(dp) :: abserr
        type(fortnum_status_t) :: status
        call deriv_central(dAfunc_hi_inc_part, r, deriv_h, der, abserr, status, ctx)
    end function ddAfunc_hi_inc_part

    real(dp) function dBfunc_hi_inc_part(r, ctx) result(der)
        real(dp), intent(in) :: r
        class(*), intent(in), optional :: ctx
        real(dp) :: abserr
        type(fortnum_status_t) :: status
        call deriv_central(Bfunc_hi_inc_part, r, deriv_h, der, abserr, status, ctx)
    end function dBfunc_hi_inc_part

    real(dp) function func_der(r, ctx) result(res)
        real(dp), intent(in) :: r
        class(*), intent(in), optional :: ctx
        real(dp) :: kz, G, Bt(1), r2k2
        select type (ctx)
        type is (zone_ctx_t)
            kz = real(ctx%zone%wd%n, dp)/get_background_rtor()
            G = Gfunc(ctx%zone, r)
            call eval_Bt(r, c_null_ptr, Bt)
            r2k2 = r*r*(kz*kz + real(ctx%zone%wd%m, dp)**2/(r*r))
            res = Bt(1)*Bt(1)/(4.0_dp*pi*r*r) + 2.0_dp*kz*Bt(1)*G/(4.0_dp*pi*r2k2)
        end select
    end function func_der

    !> Folds the oracle's rhs_incompressible (GSL-style int return) and its
    !> rhs_incompressible_fn fortnum_ode_rhs adapter into one ode_rhs_t-shaped
    !> routine; the GSL_SUCCESS return code carried no information once the
    !> adapter dropped it.
    subroutine rhs_incompressible(r, y, dydt, ctx)
        real(dp), intent(in) :: r
        real(dp), intent(in) :: y(:)
        real(dp), intent(out) :: dydt(:)
        class(*), intent(in), optional :: ctx

        class(zone_t), pointer :: zone
        complex(dp) :: A, B, dA, f, g, df, dg
        real(dp) :: der_re, der_im, abserr
        type(fortnum_status_t) :: status
        type(zone_part_ctx_t) :: dctx

        zone => null()
        select type (ctx)
        type is (zone_ctx_t)
            zone => ctx%zone
        end select

        A = Afunc_hi_inc(zone, r)
        B = Bfunc_hi_inc(zone, r)

        dctx%zone => zone

        dctx%part = 0
        call deriv_central(Afunc_hi_inc_part, r, deriv_h, der_re, abserr, status, dctx)
        dctx%part = 1
        call deriv_central(Afunc_hi_inc_part, r, deriv_h, der_im, abserr, status, dctx)
        dA = cmplx(der_re, der_im, dp)

        f = cmplx(y(1), y(2), dp) ! r*zeta
        g = cmplx(y(3), y(4), dp) ! (r*zeta)'

        df = g
        dg = -(dA*g + B*f)/A

        dydt(1) = real(df, dp)
        dydt(2) = aimag(df)
        dydt(3) = real(dg, dp)
        dydt(4) = aimag(dg)
    end subroutine rhs_incompressible

    !> Mirrors imhd_zone::calculate_basis_incompressible_gsl exactly,
    !> including the magic numbers (h0 = 1e-6 in the integration direction,
    !> max_steps = 100000) and the continuous rk8pd evolve recording every
    !> accepted step (single_step=.true., matching the GSL-faithful stepper
    !> driving fortnum_rk8pd_step_to in the oracle).
    subroutine calculate_basis_incompressible(self)
        class(zone_t), intent(inout), target :: self

        class(zone_t), pointer :: zone_ptr
        integer, parameter :: Neq = 4
        real(dp) :: t, tfinal
        real(dp) :: y(Neq)
        integer :: ind
        real(dp), allocatable :: grid(:)
        complex(dp), allocatable :: syst(:, :, :)
        real(dp) :: h0
        type(rk8pd_state_t) :: ode_state
        type(fortnum_status_t) :: status
        type(zone_ctx_t) :: ctx
        integer :: nfev, iter, j, k
        complex(dp) :: EB(8)
        logical :: valid_bc

        valid_bc = (self%bc1 == BOUNDARY_CENTER .or. self%bc2 == BOUNDARY_IDEALWALL)
        if (.not. valid_bc) then
            write (*, '(a,i0,a)') &
                'imhd::calculate_basis_incompressible_gsl: zone=', self%index, &
                ': error!'
            stop 1
        end if

        zone_ptr => self
        ctx%zone => zone_ptr

        if (self%bc1 == BOUNDARY_CENTER) then
            t = self%r1
            tfinal = self%r2

            y(1) = t**abs(self%wd%m)
            y(2) = 0.0_dp
            y(3) = real(abs(self%wd%m), dp)*t**(abs(self%wd%m) - 1)
            y(4) = 0.0_dp

            ind = 0
        else if (self%bc2 == BOUNDARY_IDEALWALL) then
            t = self%r2
            tfinal = self%r1

            y(1) = 0.0_dp
            y(2) = 0.0_dp
            y(3) = 1.0_dp
            y(4) = 0.0_dp

            ind = 0
        else
            write (*, '(a)') 'imhd::calculate_basis_incompressible_gsl: error!'
            stop 1
        end if

        allocate (grid(self%max_dim))
        allocate (syst(self%Ncomps, self%Nwaves, self%max_dim))
        syst = cmplx_zero

        h0 = 1.0e-6_dp*real(signum(tfinal - t), dp)

        call rk8pd_evolve_init(ode_state, Neq, h0, status)
        if (status%code /= FORTNUM_OK) then
            write (*, '(a)') &
                'calculate_basis_incompressible_gsl: failed to create ODE evolver'
            stop 1
        end if

        iter = 1
        grid(iter) = t
        call state_to_EB_incompressible(zone_ptr, grid(iter), y, EB)
        syst(:, ind + 1, iter) = EB

        nfev = 0
        do while (t /= tfinal)
            call rk8pd_evolve_apply(rhs_incompressible, ode_state, t, tfinal, y, &
                self%eps_abs, self%eps_rel, 100000, nfev, status, ctx, &
                single_step=.true.)

            if (status%code /= FORTNUM_OK) then
                write (*, '(a,es12.4e3)') &
                    'calculate_basis_incompressible_gsl: ODE solver failed at r = ', t
                stop 1
            end if

            iter = iter + 1
            grid(iter) = t
            call state_to_EB_incompressible(zone_ptr, grid(iter), y, EB)
            syst(:, ind + 1, iter) = EB

            if (iter == self%max_dim) then
                write (*, '(a,i0,a,es12.4e3)') &
                    'Maximum allowed iteration number is reached: iter=', iter, ' r=', t
                stop 1
            end if
        end do

        self%dim = iter

        if (allocated(self%r)) deallocate (self%r)
        allocate (self%r(self%dim))
        if (allocated(self%basis)) deallocate (self%basis)
        allocate (self%basis(self%Ncomps, self%Nwaves, self%dim))
        self%basis = cmplx_zero

        if (self%bc1 == BOUNDARY_CENTER) then
            do j = 1, self%dim
                self%r(j) = grid(j)
                self%basis(:, ind + 1, j) = syst(:, ind + 1, j)
            end do
            call normalize_imhd_basis(self%Ncomps, self%Nwaves, self%dim, self%basis, &
                ind, self%dim - 1)
        else if (self%bc2 == BOUNDARY_IDEALWALL) then
            k = 1
            do j = self%dim, 1, -1
                self%r(k) = grid(j)
                self%basis(:, ind + 1, k) = syst(:, ind + 1, j)
                k = k + 1
            end do
            call normalize_imhd_basis(self%Ncomps, self%Nwaves, self%dim, self%basis, &
                ind, 0)
        else
            write (*, '(a)') 'imhd::calculate_basis_incompressible_gsl: error!'
            stop 1
        end if

        deallocate (grid, syst)
    end subroutine calculate_basis_incompressible

    !> Mirrors imhd_zone::state_to_EB_incompressible exactly. state holds
    !> (r*zeta), (r*zeta)' as (re, im) pairs; EB returns Er, Et, Ez, Br, Bt,
    !> Bz, (r*zeta), (r*zeta)' as genuine complex values (the oracle's flat
    !> 16-double EB[] is just the same 8 complex values written as re/im
    !> pairs - this routine has no external/bind(C) caller, so the native
    !> complex representation is used directly instead).
    subroutine state_to_EB_incompressible(zone, r, state, EB)
        class(zone_t), pointer, intent(in) :: zone
        real(dp), intent(in) :: r
        real(dp), intent(in) :: state(4)
        complex(dp), intent(out) :: EB(8)

        complex(dp) :: rzeta, drzeta, zeta_r, dzeta_r
        real(dp) :: kt, kz, k0k0
        real(dp) :: Rb(4), B0t, dB0t, B0z, dB0z, F, rho(1)
        complex(dp) :: tB, omega2, omegaa2, comfac1, comfac2, zeta_t, zeta_z, alpha
        complex(dp) :: Br, Bt, Bz, Er, Et, Ez
        real(dp) :: divzeta, tst

        rzeta = cmplx(state(1), state(2), dp)
        drzeta = cmplx(state(3), state(4), dp)

        zeta_r = rzeta/r
        dzeta_r = (drzeta - zeta_r)/r

        kt = real(zone%wd%m, dp)/r
        kz = real(zone%wd%n, dp)/get_background_rtor()
        k0k0 = kt*kt + kz*kz

        call eval_Bt_dBt_Bz_dBz(r, c_null_ptr, Rb)
        B0t = Rb(1); dB0t = Rb(2); B0z = Rb(3); dB0z = Rb(4)

        F = kt*B0t + kz*B0z

        call eval_mass_density(r, c_null_ptr, rho)

        tB = 2.0_dp*cmplx_i*kz*B0t/(4.0_dp*pi*r)

        omega2 = zone%wd%olab*zone%wd%olab
        omegaa2 = cmplx(F*F/(4.0_dp*pi*rho(1)), 0.0_dp, dp)

        comfac1 = zeta_r*F*tB/(k0k0*rho(1)*(omega2 - omegaa2))
        comfac2 = drzeta*cmplx_i/(r*k0k0)

        zeta_t = -comfac1*kz + comfac2*kt
        zeta_z = comfac1*kt + comfac2*kz

        alpha = zeta_t*B0z - zeta_z*B0t

        Br = cmplx_i*F*zeta_r
        Bt = -(dzeta_r*B0t + zeta_r*dB0t) + cmplx_i*kz*alpha
        Bz = -(drzeta*B0z + rzeta*dB0z)/r - cmplx_i*kt*alpha

        Er = cmplx_zero
        Et = cmplx_zero
        Ez = cmplx_zero

        EB(1) = Er; EB(2) = Et; EB(3) = Ez
        EB(4) = Br; EB(5) = Bt; EB(6) = Bz
        EB(7) = rzeta; EB(8) = drzeta

        divzeta = abs(drzeta/r + cmplx_i*(kt*zeta_t + kz*zeta_z))
        tst = abs(drzeta/r) + abs(kt*zeta_t) + abs(kz*zeta_z)

        if (divzeta/tst > 1.0e-15_dp) then
            write (*, '(a,es12.4e3,a,es12.4e3,a,es12.4e3)') &
                'imhd_zone::state_to_EB_incompressible: warning: r = ', r, &
                char(9)//'|divzeta| = ', divzeta, char(9)//'tst = ', tst
        end if
    end subroutine state_to_EB_incompressible

end module kilca_incompressible_m
