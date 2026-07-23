!> Compressible ideal-MHD basis calculation with background flows, for
!> imhd_zone version 1 (formerly compressible_flow.{h,cpp}, the RWM
!> formulation of Bondenson et al. 1996). Same fortnum upgrade (native
!> rk8pd/deriv_central instead of the C-ABI callback shims) as
!> kilca_incompressible_m, whose Ffunc/zone_ctx_t/signum are reused here.
module kilca_compressible_flow_m
    use, intrinsic :: iso_c_binding, only: c_double, c_null_ptr
    use constants, only: dp, pi
    use kilca_zone_m, only: zone_t, BOUNDARY_CENTER, BOUNDARY_IDEALWALL
    use kilca_background_data_m, only: eval_Bt, eval_Bz, eval_Vt, eval_Vz, &
        eval_mass_density, eval_p_dp, eval_Bt_dBt_Bz_dBz
    use kilca_incompressible_m, only: zone_ctx_t, Ffunc, signum
    use fortnum_status, only: fortnum_status_t, FORTNUM_OK
    use fortnum_ode_rk8pd, only: rk8pd_state_t, rk8pd_evolve_init, rk8pd_evolve_apply
    use fortnum_multiroot, only: deriv_central
    implicit none
    private

    public :: calculate_basis_flow, state_to_EB_compressible_flow, rhs_flow

    complex(dp), parameter :: cmplx_zero = (0.0_dp, 0.0_dp)
    complex(dp), parameter :: cmplx_i = (0.0_dp, 1.0_dp)
    real(dp), parameter :: deriv_h = 1.0e-3_dp
    real(dp), parameter :: adiabat = 5.0_dp/3.0_dp

    interface
        function get_background_rtor() bind(C, name="get_background_rtor_") result(rtor)
            import :: c_double
            real(c_double) :: rtor
        end function get_background_rtor
    end interface

    external :: normalize_imhd_basis

contains

    !> Doppler-shifted frequency.
    complex(dp) function omega_D(zone, r) result(w)
        class(zone_t), pointer, intent(in) :: zone
        real(dp), intent(in) :: r
        real(dp) :: kz, Vt(1), Vz(1)
        kz = real(zone%wd%n, dp)/get_background_rtor()
        call eval_Vt(r, c_null_ptr, Vt)
        call eval_Vz(r, c_null_ptr, Vz)
        w = zone%wd%olab - (real(zone%wd%m, dp)/r)*Vt(1) - kz*Vz(1)
    end function omega_D

    complex(dp) function T_flow(zone, r) result(T)
        class(zone_t), pointer, intent(in) :: zone
        real(dp), intent(in) :: r
        real(dp) :: mdens(1), Bt(1), F, Vt(1)
        complex(dp) :: w
        call eval_mass_density(r, c_null_ptr, mdens)
        call eval_Bt(r, c_null_ptr, Bt)
        F = Ffunc(zone, r)
        w = omega_D(zone, r)
        call eval_Vt(r, c_null_ptr, Vt)
        T = F*Bt(1)/(4.0_dp*pi) + mdens(1)*w*Vt(1)
    end function T_flow

    complex(dp) function A_flow(zone, r) result(A)
        class(zone_t), pointer, intent(in) :: zone
        real(dp), intent(in) :: r
        real(dp) :: mdens(1), F
        complex(dp) :: w
        call eval_mass_density(r, c_null_ptr, mdens)
        F = Ffunc(zone, r)
        w = omega_D(zone, r)
        A = mdens(1)*w*w - F*F/(4.0_dp*pi)
    end function A_flow

    !> Q_flow's `A` local (the oracle's A_flow(r,this)) is computed and
    !> discarded in the C++ source; A_flow has no side effects, so dropping
    !> that dead call here changes nothing observable.
    complex(dp) function Q_flow(zone, r) result(Q)
        class(zone_t), pointer, intent(in) :: zone
        real(dp), intent(in) :: r
        real(dp) :: mdens(1), F, Vt(1), Bt(1)
        complex(dp) :: T, w
        call eval_mass_density(r, c_null_ptr, mdens)
        F = Ffunc(zone, r)
        call eval_Vt(r, c_null_ptr, Vt)
        call eval_Bt(r, c_null_ptr, Bt)
        T = T_flow(zone, r)
        w = omega_D(zone, r)
        Q = mdens(1)*(w*w*(Bt(1)*Bt(1)/(4.0_dp*pi) - mdens(1)*Vt(1)*Vt(1)) + &
            1.0_dp/(4.0_dp*pi)*(Bt(1)*w + F*Vt(1))*(Bt(1)*w + F*Vt(1)))
    end function Q_flow

    complex(dp) function S_flow(zone, r) result(S)
        class(zone_t), pointer, intent(in) :: zone
        real(dp), intent(in) :: r
        real(dp) :: mdens(1), F, Bt(1), Bz(1), press(2)
        complex(dp) :: w
        call eval_mass_density(r, c_null_ptr, mdens)
        F = Ffunc(zone, r)
        call eval_Bt(r, c_null_ptr, Bt)
        call eval_Bz(r, c_null_ptr, Bz)
        w = omega_D(zone, r)
        call eval_p_dp(r, c_null_ptr, press)
        S = ((Bt(1)*Bt(1) + Bz(1)*Bz(1))/(4.0_dp*pi) + adiabat*press(1))*mdens(1)*w*w &
            - adiabat*press(1)*F*F/(4.0_dp*pi)
    end function S_flow

    complex(dp) function C11_flow(zone, r) result(C11)
        class(zone_t), pointer, intent(in) :: zone
        real(dp), intent(in) :: r
        real(dp) :: mdens(1)
        complex(dp) :: w, Q, S, T
        call eval_mass_density(r, c_null_ptr, mdens)
        w = omega_D(zone, r)
        Q = Q_flow(zone, r)
        S = S_flow(zone, r)
        T = T_flow(zone, r)
        C11 = mdens(1)*w*w*Q/(r*r) - 2.0_dp*real(zone%wd%m, dp)*S*T/(r*r*r)
    end function C11_flow

    complex(dp) function C12_flow(zone, r) result(C12)
        class(zone_t), pointer, intent(in) :: zone
        real(dp), intent(in) :: r
        real(dp) :: mdens(1), kz
        complex(dp) :: S, w
        call eval_mass_density(r, c_null_ptr, mdens)
        kz = real(zone%wd%n, dp)/get_background_rtor()
        S = S_flow(zone, r)
        w = omega_D(zone, r)
        C12 = mdens(1)*mdens(1)*w*w*w*w - (kz*kz + real(zone%wd%m, dp)**2/(r*r))*S
    end function C12_flow

    complex(dp) function C21_flow(zone, r) result(C21)
        class(zone_t), pointer, intent(in) :: zone
        real(dp), intent(in) :: r
        complex(dp) :: S, A, T, Q, C4
        S = S_flow(zone, r)
        A = A_flow(zone, r)
        T = T_flow(zone, r)
        Q = Q_flow(zone, r)
        C4 = C4_flow(zone, r)
        C21 = A*S*C4/r - 4.0_dp*S*T*T/(r*r*r) + Q*Q/(r*r*r)
    end function C21_flow

    complex(dp) function C22_flow(zone, r) result(C22)
        class(zone_t), pointer, intent(in) :: zone
        real(dp), intent(in) :: r
        C22 = r*C11_flow(zone, r)
    end function C22_flow

    real(dp) function func_deriv_in_C4(r, ctx) result(res)
        real(dp), intent(in) :: r
        class(*), intent(in), optional :: ctx
        real(dp) :: mdens(1), Vt(1), Bt(1)
        select type (ctx)
        type is (zone_ctx_t)
            call eval_mass_density(r, c_null_ptr, mdens)
            call eval_Vt(r, c_null_ptr, Vt)
            call eval_Bt(r, c_null_ptr, Bt)
            res = (Bt(1)*Bt(1)/(4.0_dp*pi) - mdens(1)*Vt(1)*Vt(1))/(r*r)
        end select
    end function func_deriv_in_C4

    complex(dp) function C4_flow(zone, r) result(C4)
        class(zone_t), pointer, intent(in) :: zone
        real(dp), intent(in) :: r
        real(dp) :: der, abserr
        type(fortnum_status_t) :: status
        type(zone_ctx_t) :: ctx
        complex(dp) :: A

        ctx%zone => zone
        call deriv_central(func_deriv_in_C4, r, deriv_h, der, abserr, status, ctx)

        A = A_flow(zone, r)
        C4 = A + r*der
    end function C4_flow

    !> Folds the oracle's rhs_flow (GSL-style int return) and its rhs_flow_fn
    !> fortnum_ode_rhs adapter into one ode_rhs_t-shaped routine.
    subroutine rhs_flow(r, y, dydt, ctx)
        real(dp), intent(in) :: r
        real(dp), intent(in) :: y(:)
        real(dp), intent(out) :: dydt(:)
        class(*), intent(in), optional :: ctx

        class(zone_t), pointer :: zone
        complex(dp) :: C11, C12, C21, C22, A, S, f, g, df, dg

        zone => null()
        select type (ctx)
        type is (zone_ctx_t)
            zone => ctx%zone
        end select

        C11 = C11_flow(zone, r)
        C12 = C12_flow(zone, r)
        C21 = C21_flow(zone, r)
        C22 = C22_flow(zone, r)
        A = A_flow(zone, r)
        S = S_flow(zone, r)

        f = cmplx(y(1), y(2), dp) ! r*zeta
        g = cmplx(y(3), y(4), dp) ! pressure

        df = (C11*f - C12*g)/(A*S)*r
        dg = (C21*f - C22*g)/(A*S)

        dydt(1) = real(df, dp)
        dydt(2) = aimag(df)
        dydt(3) = real(dg, dp)
        dydt(4) = aimag(dg)
    end subroutine rhs_flow

    !> Mirrors imhd_zone::calculate_basis_flow_gsl exactly, including the
    !> oracle's own copy-paste artifact: the bc1/bc2 sanity-check error
    !> message names "calculate_basis_incompressible_gsl" even though this is
    !> the flow path.
    subroutine calculate_basis_flow(self)
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
        complex(dp) :: EB(8), rzeta, drzeta, pstar
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

            rzeta = cmplx(t**abs(self%wd%m), 0.0_dp, dp)
            drzeta = cmplx(real(abs(self%wd%m), dp)*t**(abs(self%wd%m) - 1), 0.0_dp, dp)

            pstar = (C11_flow(zone_ptr, t)*rzeta - &
                S_flow(zone_ptr, t)*A_flow(zone_ptr, t)*drzeta/t)/C12_flow(zone_ptr, t)
            y(3) = real(pstar, dp)
            y(4) = aimag(pstar)

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
            write (*, '(a)') 'imhd::calculate_basis_compressible_gsl: error!'
            stop 1
        end if

        allocate (grid(self%max_dim))
        allocate (syst(self%Ncomps, self%Nwaves, self%max_dim))
        syst = cmplx_zero

        h0 = 1.0e-6_dp*real(signum(tfinal - t), dp)

        call rk8pd_evolve_init(ode_state, Neq, h0, status)
        if (status%code /= FORTNUM_OK) then
            write (*, '(a)') 'calculate_basis_flow_gsl: failed to create ODE evolver'
            stop 1
        end if

        iter = 1
        grid(iter) = t
        call state_to_EB_compressible_flow(zone_ptr, grid(iter), y, EB)
        syst(:, ind + 1, iter) = EB

        nfev = 0
        do while (t /= tfinal)
            call rk8pd_evolve_apply(rhs_flow, ode_state, t, tfinal, y, &
                self%eps_abs, self%eps_rel, 100000, nfev, status, ctx, &
                single_step=.true.)

            if (status%code /= FORTNUM_OK) then
                write (*, '(a,es12.4e3)') &
                    'calculate_basis_flow_gsl: ODE solver failed at r = ', t
                stop 1
            end if

            iter = iter + 1
            grid(iter) = t
            call state_to_EB_compressible_flow(zone_ptr, grid(iter), y, EB)
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
            write (*, '(a)') 'imhd::calculate_basis_compressible_gsl: error!'
            stop 1
        end if

        deallocate (grid, syst)
    end subroutine calculate_basis_flow

    !> Mirrors imhd_zone::state_to_EB_compressible_flow exactly. Note the
    !> oracle's own subtlety, preserved here: EB's last complex slot is
    !> assigned the locally recomputed `drzeta` (from the C11/C12/S/A
    !> formula below), NOT the raw input `press` state component that feeds
    !> into that formula - the header comment calling it "(pressure)" is
    !> stale.
    subroutine state_to_EB_compressible_flow(zone, r, state, EB)
        class(zone_t), pointer, intent(in) :: zone
        real(dp), intent(in) :: r
        real(dp), intent(in) :: state(4)
        complex(dp), intent(out) :: EB(8)

        complex(dp) :: rzeta, press, zeta_r, drzeta, dzeta_r
        real(dp) :: kt, kz, k0k0
        real(dp) :: Rb(4), B0t, dB0t, B0z, dB0z, F
        real(dp) :: mdens(1), Vt(1), Vz(1), Kp(2), press_0, dpress_0
        complex(dp) :: omega, omega_Dop
        complex(dp) :: Acoef, M, G, N, H, zeta_t, zeta_z, alpha
        complex(dp) :: Br, Bt, Bz, Er, Et, Ez

        rzeta = cmplx(state(1), state(2), dp)
        press = cmplx(state(3), state(4), dp)

        zeta_r = rzeta/r

        kt = real(zone%wd%m, dp)/r
        kz = real(zone%wd%n, dp)/get_background_rtor()
        k0k0 = kt*kt + kz*kz

        call eval_Bt_dBt_Bz_dBz(r, c_null_ptr, Rb)
        B0t = Rb(1); dB0t = Rb(2); B0z = Rb(3); dB0z = Rb(4)

        F = kt*B0t + kz*B0z

        call eval_mass_density(r, c_null_ptr, mdens)
        call eval_Vt(r, c_null_ptr, Vt)
        call eval_Vz(r, c_null_ptr, Vz)

        omega = zone%wd%olab
        omega_Dop = omega - kt*Vt(1) - kz*Vz(1)

        call eval_p_dp(r, c_null_ptr, Kp)
        press_0 = Kp(1); dpress_0 = Kp(2)

        drzeta = r*(C11_flow(zone, r)*rzeta - C12_flow(zone, r)*press)/ &
            (S_flow(zone, r)*A_flow(zone, r))

        Acoef = -kz*adiabat*press_0*kt + B0t*B0z*k0k0/4.0_dp/pi

        M = zeta_r*cmplx_i*(F/4.0_dp/pi*dB0z + kz*dpress_0 - kt/4.0_dp/pi*dB0z*B0t - &
            kz/4.0_dp/pi*B0t*(B0t/r - dB0t)) + &
            drzeta*cmplx_i*(kz*B0t*B0t/4.0_dp/pi/r - kt/4.0_dp/pi*B0t*B0z/r + &
            kz*adiabat*press_0/r)

        G = zeta_r*cmplx_i*(F/4.0_dp/pi*(B0t/r + dB0t) + &
            2.0_dp*omega_Dop*mdens(1)*Vt(1)/r + dpress_0*kt + &
            kt/4.0_dp/pi*B0z*dB0z + kz/4.0_dp/pi*B0z*(B0t/r - dB0t)) + &
            drzeta*cmplx_i*(kt/4.0_dp/pi*B0z*B0z/r - kz/4.0_dp/pi/r*B0z*B0t + &
            kt*adiabat*press_0/r)

        N = -mdens(1)*omega_Dop*omega_Dop + kz*kz*adiabat*press_0 &
            + B0t*B0t/4.0_dp/pi*k0k0
        H = -mdens(1)*omega_Dop*omega_Dop + kt*kt*adiabat*press_0 &
            + B0z*B0z/4.0_dp/pi*k0k0

        dzeta_r = (drzeta - zeta_r)/r

        zeta_t = (M*Acoef + G*N)/(N*H - Acoef*Acoef)
        zeta_z = (G*Acoef + M*H)/(N*H - Acoef*Acoef)

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
    end subroutine state_to_EB_compressible_flow

end module kilca_compressible_flow_m
