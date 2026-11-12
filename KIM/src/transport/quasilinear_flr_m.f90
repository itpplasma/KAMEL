module quasilinear_flr_m
    !! One-harmonic, two-wave building block for KIM's finite-Larmor-radius ion tensor.
    !!
    !! The implementation follows the restricted field-insertion operators
    !! derived in mathematica/verify_quasilinear_flr_transport.wl.  The
    !! Bparallel derivatives act only on the Gaussian--Bessel gyro-factor;
    !! kinetic detuning and equilibrium quantities are deliberately absent
    !! from this interface.
    use KIM_kinds_m, only: dp
    implicit none
    private

    integer, parameter, public :: QL_PHI = 1
    integer, parameter, public :: QL_BR = 2
    integer, parameter, public :: QL_BPAR = 3
    integer, parameter, public :: QL_FLR_ALGEBRA_VERSION = 1
    public :: calc_ion_flr_harmonic
    public :: calc_ion_flr_harmonic_debug
    public :: build_transverse_moments
    public :: scaled_besseli_large_x

    type :: jet_t
        complex(dp) :: v = (0.0_dp, 0.0_dp)
        complex(dp) :: ds = (0.0_dp, 0.0_dp)
        complex(dp) :: dobs = (0.0_dp, 0.0_dp)
        complex(dp) :: dso = (0.0_dp, 0.0_dp)
    end type jet_t

contains

    subroutine calc_ion_flr_harmonic(ell, ks_s, kr_s, ks_o, kr_o, &
            vT, omega_c, omega_signed, c_light, B0, nu, fields_s, fields_o, &
            ifunc, tensor)
        integer, intent(in) :: ell
        real(dp), intent(in) :: ks_s, kr_s, ks_o, kr_o
        real(dp), intent(in) :: vT, omega_c, omega_signed, c_light, B0, nu
        complex(dp), intent(in) :: fields_s(3), fields_o(3)
        complex(dp), intent(in) :: ifunc(0:3,0:3)
        real(dp), intent(out) :: tensor(2,2)
        call calc_ion_flr_harmonic_impl(ell, ks_s, kr_s, ks_o, kr_o, &
            vT, omega_c, omega_signed, c_light, B0, nu, fields_s, fields_o, &
            ifunc, tensor)
    end subroutine calc_ion_flr_harmonic

    subroutine calc_ion_flr_harmonic_debug(ell, ks_s, kr_s, ks_o, kr_o, &
            vT, omega_c, omega_signed, c_light, B0, nu, fields_s, fields_o, &
            ifunc, tensor, channels)
        integer, intent(in) :: ell
        real(dp), intent(in) :: ks_s, kr_s, ks_o, kr_o
        real(dp), intent(in) :: vT, omega_c, omega_signed, c_light, B0, nu
        complex(dp), intent(in) :: fields_s(3), fields_o(3)
        complex(dp), intent(in) :: ifunc(0:3,0:3)
        real(dp), intent(out) :: tensor(2,2)
        complex(dp), intent(out) :: channels(2,2,3,3)
        call calc_ion_flr_harmonic_impl(ell, ks_s, kr_s, ks_o, kr_o, &
            vT, omega_c, omega_signed, c_light, B0, nu, fields_s, fields_o, &
            ifunc, tensor, channels)
    end subroutine calc_ion_flr_harmonic_debug

    subroutine calc_ion_flr_harmonic_impl(ell, ks_s, kr_s, ks_o, kr_o, &
            vT, omega_c, omega_signed, c_light, B0, nu, fields_s, fields_o, &
            ifunc, tensor, channels)
        integer, intent(in) :: ell
        real(dp), intent(in) :: ks_s, kr_s, ks_o, kr_o
        real(dp), intent(in) :: vT, omega_c, omega_signed, c_light, B0, nu
        complex(dp), intent(in) :: fields_s(3), fields_o(3)
        complex(dp), intent(in) :: ifunc(0:3,0:3)
        real(dp), intent(out) :: tensor(2,2)
        complex(dp), intent(out), optional :: channels(2,2,3,3)

        complex(dp) :: transverse(0:2,3,3), transverse_swapped(0:2,3,3)
        complex(dp) :: raw(2,2,3,3), raw_swapped(2,2,3,3)
        complex(dp) :: channel_work(2,2,3,3)
        real(dp) :: bessel_table(0:abs(ell)+3)
        logical :: table_ready, table_scaled
        integer :: i, j, observation, source

        if (vT <= 0.0_dp .or. omega_c <= 0.0_dp .or. &
                abs(B0) <= tiny(1.0_dp) .or. nu <= 0.0_dp) &
            error stop 'calc_ion_flr_harmonic: invalid physical scale'

        ! The Bessel argument bx is symmetric in source and observation, so
        ! both orientations share one table of I_n(bx).
        table_ready = .false.
        call build_transverse_moments(ell, ks_s, kr_s, ks_o, kr_o, vT, &
            omega_c, omega_signed, c_light, B0, transverse, &
            bessel_table, table_ready, table_scaled)
        call build_transverse_moments(ell, ks_o, kr_o, ks_s, kr_s, vT, &
            omega_c, omega_signed, c_light, B0, transverse_swapped, &
            bessel_table, table_ready, table_scaled)

        call build_raw_channels(transverse, fields_s, fields_o, ifunc, nu, raw)
        call build_raw_channels(transverse_swapped, fields_o, fields_s, ifunc, &
            nu, raw_swapped)
        do i = 1, 2
            do j = 1, 2
                do observation = 1, 3
                    do source = 1, 3
                        channel_work(i,j,observation,source) = &
                            cmplx(0.5_dp,0.0_dp,dp) * &
                            (raw(i,j,observation,source) &
                            + conjg(raw_swapped(j,i,source,observation)))
                    end do
                end do
                tensor(i,j) = real(sum(channel_work(i,j,:,:)), dp)
            end do
        end do
        if (present(channels)) channels = channel_work
    end subroutine calc_ion_flr_harmonic_impl

    subroutine build_raw_channels(transverse, fields_s, fields_o, ifunc, nu, raw)
        complex(dp), intent(in) :: transverse(0:2,3,3)
        complex(dp), intent(in) :: fields_s(3), fields_o(3)
        complex(dp), intent(in) :: ifunc(0:3,0:3)
        real(dp), intent(in) :: nu
        complex(dp), intent(out) :: raw(2,2,3,3)
        complex(dp) :: contraction
        integer :: i, j, observation, source

        do i = 1, 2
            do j = 1, 2
                do observation = 1, 3
                    do source = 1, 3
                        contraction = contract_energy(i, j, observation, source, &
                            transverse(:,observation,source), ifunc)
                        raw(i,j,observation,source) = &
                            conjg(fields_o(observation))*fields_s(source) * &
                            contraction/cmplx(2.0_dp*nu,0.0_dp,dp)
                    end do
                end do
            end do
        end do
    end subroutine build_raw_channels

    subroutine build_transverse_moments(ell, ks_s, kr_s, ks_o, kr_o, vT, &
            omega_c, omega_signed, c_light, B0, transverse, &
            bessel_table, table_ready, table_scaled)
        integer, intent(in) :: ell
        real(dp), intent(in) :: ks_s, kr_s, ks_o, kr_o
        real(dp), intent(in) :: vT, omega_c, omega_signed, c_light, B0
        complex(dp), intent(out) :: transverse(0:2,3,3)
        real(dp), intent(inout), optional :: bessel_table(0:)
        logical, intent(inout), optional :: table_ready, table_scaled

        type(jet_t) :: w(0:2)
        real(dp) :: inv_b02
        complex(dp), parameter :: imaginary_unit = (0.0_dp,1.0_dp)
        integer :: q

        if (omega_c <= 0.0_dp .or. abs(B0) <= tiny(1.0_dp)) &
            error stop 'build_transverse_moments: invalid magnetic scale'

        call gaussian_bessel_jets(ell, ks_s, kr_s, ks_o, kr_o, vT, omega_c, w, &
            bessel_table, table_ready, table_scaled)
        inv_b02 = 1.0_dp/B0**2

        do q = 0, 2
            transverse(q,QL_PHI,QL_PHI) = &
                cmplx(c_light**2*ks_o*ks_s*inv_b02,0.0_dp,dp)*w(q)%v
            transverse(q,QL_PHI,QL_BR) = &
                imaginary_unit*cmplx(c_light*ks_o*vT*inv_b02,0.0_dp,dp)*w(q)%v
            transverse(q,QL_BR,QL_PHI) = &
                -imaginary_unit*cmplx(c_light*ks_s*vT*inv_b02,0.0_dp,dp)*w(q)%v
            transverse(q,QL_BR,QL_BR) = cmplx(vT**2*inv_b02,0.0_dp,dp)*w(q)%v

            transverse(q,QL_PHI,QL_BPAR) = &
                cmplx(-c_light*ks_o*omega_signed*inv_b02,0.0_dp,dp)*w(q)%ds
            transverse(q,QL_BPAR,QL_PHI) = &
                cmplx(-c_light*ks_s*omega_signed*inv_b02,0.0_dp,dp)*w(q)%dobs
            transverse(q,QL_BR,QL_BPAR) = &
                imaginary_unit*cmplx(vT*omega_signed*inv_b02,0.0_dp,dp)*w(q)%ds
            transverse(q,QL_BPAR,QL_BR) = &
                -imaginary_unit*cmplx(vT*omega_signed*inv_b02,0.0_dp,dp)*w(q)%dobs
            transverse(q,QL_BPAR,QL_BPAR) = &
                cmplx(omega_signed**2*inv_b02,0.0_dp,dp)*w(q)%dso
        end do
    end subroutine build_transverse_moments

    function contract_energy(i, j, observation, source, transverse, ifunc) result(value)
        integer, intent(in) :: i, j, observation, source
        complex(dp), intent(in) :: transverse(0:2)
        complex(dp), intent(in) :: ifunc(0:3,0:3)
        complex(dp) :: value
        integer :: po, ps

        po = merge(1, 0, observation == QL_BR)
        ps = merge(1, 0, source == QL_BR)

        select case (10*i+j)
        case (11)
            value = transverse(0)*ifunc(po,ps)
        case (12)
            value = cmplx(0.5_dp,0.0_dp,dp)*transverse(0)*ifunc(po,ps+2) &
                + transverse(1)*ifunc(po,ps)
        case (21)
            value = cmplx(0.5_dp,0.0_dp,dp)*transverse(0)*ifunc(po+2,ps) &
                + transverse(1)*ifunc(po,ps)
        case (22)
            value = cmplx(0.25_dp,0.0_dp,dp)*transverse(0)*ifunc(po+2,ps+2) &
                + cmplx(0.5_dp,0.0_dp,dp)*transverse(1) &
                    *(ifunc(po+2,ps)+ifunc(po,ps+2)) &
                + transverse(2)*ifunc(po,ps)
        case default
            error stop 'contract_energy: invalid tensor index'
        end select
    end function contract_energy

    subroutine gaussian_bessel_jets(ell, ks_s, kr_s, ks_o, kr_o, vT, omega_c, w, &
            bessel_table, table_ready, table_scaled)
        integer, intent(in) :: ell
        real(dp), intent(in) :: ks_s, kr_s, ks_o, kr_o, vT, omega_c
        type(jet_t), intent(out) :: w(0:2)
        real(dp), intent(inout), optional :: bessel_table(0:)
        logical, intent(inout), optional :: table_ready, table_scaled
        real(dp) :: local_table(0:abs(ell)+3)
        logical :: scaled

        type(jet_t) :: bp, bx, phase, damping, il, ilm1
        type(jet_t) :: coeff, term1, term2
        real(dp) :: kp_s, kp_o, a, alpha_s, alpha_o
        complex(dp) :: ls, lo

        kp_s = hypot(ks_s, kr_s)
        kp_o = hypot(ks_o, kr_o)
        if (kp_s <= tiny(1.0_dp) .or. kp_o <= tiny(1.0_dp)) &
            then
            call gaussian_bessel_jets_zero_k(ell, ks_s, kr_s, ks_o, kr_o, &
                vT, omega_c, w)
            return
        end if

        a = (vT/omega_c)**2
        bp = jet_constant(0.5_dp*a*(kp_s**2+kp_o**2))
        bp%ds = cmplx(a*ks_s,0.0_dp,dp)
        bp%dobs = cmplx(a*ks_o,0.0_dp,dp)

        bx = jet_constant(a*kp_s*kp_o)
        bx%ds = cmplx(a*(ks_s/kp_s)*kp_o,0.0_dp,dp)
        bx%dobs = cmplx(a*kp_s*(ks_o/kp_o),0.0_dp,dp)
        bx%dso = cmplx(a*(ks_s/kp_s)*(ks_o/kp_o),0.0_dp,dp)

        alpha_s = atan2(kr_s, ks_s)
        alpha_o = atan2(kr_o, ks_o)
        phase = jet_constant(0.0_dp)
        phase%v = exp(cmplx(0.0_dp,real(ell,dp)*(alpha_o-alpha_s),dp))
        if (ell /= 0) then
            ls = cmplx(0.0_dp,real(ell,dp)*kr_s/kp_s**2,dp)
            lo = cmplx(0.0_dp,-real(ell,dp)*kr_o/kp_o**2,dp)
            phase%ds = phase%v*ls
            phase%dobs = phase%v*lo
            phase%dso = phase%v*ls*lo
        end if

        ! Orders ell-3 .. ell+2 enter il, ilm1 and their two derivatives.
        if (present(bessel_table)) then
            if (.not. table_ready) then
                table_scaled = real(bx%v,dp) > 500.0_dp
                call fill_bessel_table(real(bx%v,dp), table_scaled, bessel_table)
                table_ready = .true.
            end if
            scaled = table_scaled
            local_table = bessel_table(0:abs(ell)+3)
        else
            scaled = real(bx%v,dp) > 500.0_dp
            call fill_bessel_table(real(bx%v,dp), scaled, local_table)
        end if

        if (scaled) then
            ! Keep exp(-bPlus) fused with I_n(bCross); evaluating the two
            ! factors separately would produce 0*Inf for valid large FLR.
            damping = jet_constant(1.0_dp)
            il = jet_scaled_besseli(ell,bp,bx,local_table)
            ilm1 = jet_scaled_besseli(ell-1,bp,bx,local_table)
        else
            damping = jet_exp(jet_scale(bp,-1.0_dp))
            il = jet_besseli(ell,bx,local_table)
            ilm1 = jet_besseli(ell-1,bx,local_table)
        end if
        w(0) = jet_mul(phase,jet_mul(damping,il))

        coeff = jet_add(jet_scale(bp,-1.0_dp),jet_constant(1.0_dp-real(ell,dp)))
        term1 = jet_mul(coeff,il)
        term2 = jet_mul(bx,ilm1)
        w(1) = jet_mul(phase,jet_mul(damping,jet_add(term1,term2)))

        term1 = jet_mul(jet_add(jet_constant(3.0_dp),jet_scale(bp,-2.0_dp)), &
            jet_mul(bx,ilm1))
        coeff = jet_add(jet_constant(2.0_dp-3.0_dp*real(ell,dp)+real(ell,dp)**2), &
            jet_add(jet_mul(bp,bp),jet_mul(bx,bx)))
        coeff = jet_add(coeff,jet_scale(bp,2.0_dp*(-2.0_dp+real(ell,dp))))
        term2 = jet_mul(coeff,il)
        w(2) = jet_mul(phase,jet_mul(damping,jet_add(term1,term2)))
    end subroutine gaussian_bessel_jets

    subroutine gaussian_bessel_jets_zero_k(ell, ks_s, kr_s, ks_o, kr_o, &
            vT, omega_c, w)
        !! Cartesian regular branch for k_perp=0.  It uses
        !!   exp(i ell(alpha_o-alpha_s)) I_ell(bx)
        !! = a^|ell| z_o^|ell| z_s^|ell| H_|ell|(bx^2),
        !! where H_n(y)=sum y^k/[2^(n+2k) k! (n+k)!].  No polar angle or
        !! ratio of separately vanishing quantities is evaluated.
        integer, intent(in) :: ell
        real(dp), intent(in) :: ks_s, kr_s, ks_o, kr_o, vT, omega_c
        type(jet_t), intent(out) :: w(0:2)
        type(jet_t) :: bp, y, source_cart, observation_cart, amplitude
        type(jet_t) :: h0, h1, h2, term, power_y, exponential, u, bracket
        real(dp) :: a, coefficient, harmonic_sign
        integer :: n, k

        n = abs(ell)
        a = (vT/omega_c)**2
        harmonic_sign = merge(-1.0_dp,1.0_dp,ell < 0)
        bp = jet_constant(0.5_dp*a*((ks_s**2+kr_s**2)+(ks_o**2+kr_o**2)))
        bp%ds = cmplx(a*ks_s,0.0_dp,dp)
        bp%dobs = cmplx(a*ks_o,0.0_dp,dp)

        y = jet_constant(a**2*(ks_s**2+kr_s**2)*(ks_o**2+kr_o**2))
        y%ds = cmplx(2.0_dp*a**2*ks_s*(ks_o**2+kr_o**2),0.0_dp,dp)
        y%dobs = cmplx(2.0_dp*a**2*ks_o*(ks_s**2+kr_s**2),0.0_dp,dp)
        y%dso = cmplx(4.0_dp*a**2*ks_s*ks_o,0.0_dp,dp)

        source_cart = jet_complex_constant(cmplx(ks_s, &
            -harmonic_sign*kr_s,dp))
        observation_cart = jet_complex_constant(cmplx(ks_o, &
            harmonic_sign*kr_o,dp))
        source_cart%ds = (1.0_dp,0.0_dp)
        observation_cart%dobs = (1.0_dp,0.0_dp)
        amplitude = jet_scale(jet_mul(jet_power(source_cart,n), &
            jet_power(observation_cart,n)),a**n)

        coefficient = 1.0_dp
        do k = 1, n
            coefficient = coefficient/(2.0_dp*real(k,dp))
        end do
        h0 = jet_constant(0.0_dp)
        h1 = jet_constant(0.0_dp)
        h2 = jet_constant(0.0_dp)
        power_y = jet_constant(1.0_dp)
        do k = 0, 4
            if (k > 0) then
                coefficient = coefficient/(4.0_dp*real(k*(n+k),dp))
                power_y = jet_mul(power_y,y)
            end if
            h0 = jet_add(h0,jet_scale(power_y,coefficient))
            if (k >= 1) then
                term = jet_power(y,k-1)
                h1 = jet_add(h1,jet_scale(term,real(k,dp)*coefficient))
            end if
            if (k >= 2) then
                term = jet_power(y,k-2)
                h2 = jet_add(h2,jet_scale(term, &
                    real(k*(k-1),dp)*coefficient))
            end if
        end do

        exponential = jet_exp(jet_scale(bp,-1.0_dp))
        u = jet_add(jet_constant(real(n+1,dp)),jet_scale(bp,-1.0_dp))
        w(0) = jet_mul(amplitude,jet_mul(exponential,h0))
        bracket = jet_add(jet_mul(u,h0), &
            jet_scale(jet_mul(y,h1),2.0_dp))
        w(1) = jet_mul(amplitude,jet_mul(exponential,bracket))
        bracket = jet_mul(jet_add(jet_add(jet_mul(u,u),jet_scale(u,2.0_dp)), &
            jet_constant(-real(n+1,dp))),h0)
        bracket = jet_add(bracket,jet_mul(jet_scale(y,6.0_dp),h1))
        bracket = jet_add(bracket,jet_mul(jet_scale(jet_mul(u,y),4.0_dp),h1))
        bracket = jet_add(bracket,jet_scale(jet_mul(jet_mul(y,y),h2),4.0_dp))
        w(2) = jet_mul(amplitude,jet_mul(exponential,bracket))
    end subroutine gaussian_bessel_jets_zero_k

    pure function jet_constant(value) result(out)
        real(dp), intent(in) :: value
        type(jet_t) :: out
        out%v = cmplx(value,0.0_dp,dp)
    end function jet_constant

    pure function jet_complex_constant(value) result(out)
        complex(dp), intent(in) :: value
        type(jet_t) :: out
        out%v = value
    end function jet_complex_constant

    pure function jet_power(a,n) result(out)
        type(jet_t), intent(in) :: a
        integer, intent(in) :: n
        type(jet_t) :: out
        integer :: k
        out = jet_constant(1.0_dp)
        do k = 1, n
            out = jet_mul(out,a)
        end do
    end function jet_power

    pure function jet_add(a,b) result(out)
        type(jet_t), intent(in) :: a,b
        type(jet_t) :: out
        out%v=a%v+b%v; out%ds=a%ds+b%ds; out%dobs=a%dobs+b%dobs; out%dso=a%dso+b%dso
    end function jet_add

    pure function jet_scale(a,s) result(out)
        type(jet_t), intent(in) :: a
        real(dp), intent(in) :: s
        type(jet_t) :: out
        out%v=cmplx(s,0.0_dp,dp)*a%v
        out%ds=cmplx(s,0.0_dp,dp)*a%ds
        out%dobs=cmplx(s,0.0_dp,dp)*a%dobs
        out%dso=cmplx(s,0.0_dp,dp)*a%dso
    end function jet_scale

    pure function jet_mul(a,b) result(out)
        type(jet_t), intent(in) :: a,b
        type(jet_t) :: out
        out%v=a%v*b%v
        out%ds=a%ds*b%v+a%v*b%ds
        out%dobs=a%dobs*b%v+a%v*b%dobs
        out%dso=a%dso*b%v+a%ds*b%dobs+a%dobs*b%ds+a%v*b%dso
    end function jet_mul

    pure function jet_exp(a) result(out)
        type(jet_t), intent(in) :: a
        type(jet_t) :: out
        out%v=exp(a%v)
        out%ds=out%v*a%ds
        out%dobs=out%v*a%dobs
        out%dso=out%v*(a%dso+a%ds*a%dobs)
    end function jet_exp

    subroutine fill_bessel_table(x, scaled, table)
        !! table(n) = I_n(x), or exp(-x) I_n(x) when scaled, for n = 0..ubound.
        !! Each order is evaluated once instead of once per use. The scaled
        !! branch is used only for x > 500, where the asymptotic expansion is
        !! exact to rounding; the unscaled one uses bessel_in per order, which
        !! itself switches to its asymptotic expansion at large x.
        use fortnum_special_bessel, only: bessel_in
        real(dp), intent(in) :: x
        logical, intent(in) :: scaled
        real(dp), intent(out) :: table(0:)
        integer :: n

        if (scaled) then
            do n = 0, ubound(table,1)
                table(n) = scaled_besseli_large_x(n,x)
            end do
        else
            do n = 0, ubound(table,1)
                table(n) = bessel_in(n,x)
            end do
        end if
    end subroutine fill_bessel_table

    pure function scaled_besseli_large_x(n,x) result(value)
        !! exp(-x) I_n(x) from the large-argument expansion (DLMF 10.40.1),
        !!   exp(-x) I_n(x) ~ (2 pi x)^(-1/2) sum_k (-1)^k a_k(n) / x^k,
        !!   a_k(n) = prod_{i=1..k} (4 n^2 - (2i-1)^2) / (k! 8^k).
        !! For x > 500 and small n the terms fall below double-precision
        !! resolution within about ten terms (error < 1e-15 against SciPy's
        !! ive for n <= 8), long before the series starts to diverge.
        use constants_m, only: pi
        integer, intent(in) :: n
        real(dp), intent(in) :: x
        real(dp) :: value
        real(dp) :: mu, term, total
        integer :: k

        mu = 4.0_dp*real(n,dp)**2
        term = 1.0_dp
        total = 1.0_dp
        do k = 1, 40
            term = -term*(mu-real(2*k-1,dp)**2)/(real(k,dp)*8.0_dp*x)
            total = total+term
            if (abs(term) <= 1.0e-17_dp*abs(total)) exit
        end do
        value = total/sqrt(2.0_dp*pi*x)
    end function scaled_besseli_large_x

    function jet_scaled_besseli(order,bplus,x,table) result(out)
        integer, intent(in) :: order
        type(jet_t), intent(in) :: bplus, x
        real(dp), intent(in) :: table(0:)
        type(jet_t) :: out, scaled, prefactor
        real(dp) :: f, fp, fpp

        ! table holds S_n(x)=exp(-x) I_n(x) (I_{-n}=I_n).  Reconstruct
        ! exp(-bplus) I_n(x) with the bounded factor exp(x-bplus), where
        ! bplus >= x for the two-wave geometry.
        f = table(abs(order))
        fp = 0.5_dp*(table(abs(order-1))+table(abs(order+1)))-f
        fpp = 0.25_dp*(table(abs(order-2))+2.0_dp*f+table(abs(order+2))) &
            -(table(abs(order-1))+table(abs(order+1)))+f
        scaled%v = cmplx(f,0.0_dp,dp)
        scaled%ds = cmplx(fp,0.0_dp,dp)*x%ds
        scaled%dobs = cmplx(fp,0.0_dp,dp)*x%dobs
        scaled%dso = cmplx(fpp,0.0_dp,dp)*x%ds*x%dobs &
            + cmplx(fp,0.0_dp,dp)*x%dso
        prefactor = jet_exp(jet_add(x,jet_scale(bplus,-1.0_dp)))
        out = jet_mul(prefactor,scaled)
    end function jet_scaled_besseli

    function jet_besseli(order,x,table) result(out)
        integer, intent(in) :: order
        type(jet_t), intent(in) :: x
        real(dp), intent(in) :: table(0:)
        type(jet_t) :: out
        real(dp) :: f, fp, fpp

        ! table holds I_n(x) (I_{-n}=I_n).
        f = table(abs(order))
        fp = 0.5_dp*(table(abs(order-1))+table(abs(order+1)))
        fpp = 0.25_dp*(table(abs(order-2))+2.0_dp*f+table(abs(order+2)))
        out%v = cmplx(f,0.0_dp,dp)
        out%ds = cmplx(fp,0.0_dp,dp)*x%ds
        out%dobs = cmplx(fp,0.0_dp,dp)*x%dobs
        out%dso = cmplx(fpp,0.0_dp,dp)*x%ds*x%dobs &
            + cmplx(fp,0.0_dp,dp)*x%dso
    end function jet_besseli

end module quasilinear_flr_m
