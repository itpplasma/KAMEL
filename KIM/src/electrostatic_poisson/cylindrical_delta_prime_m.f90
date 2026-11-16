module cylindrical_delta_prime_m
    !! Cylindrical ideal outer-region tearing parameter Delta' of one (m, n)
    !! mode, in KIM's field-aligned conventions (Br = i ks A_par, Ampere
    !! nabla_perp^2 A_par = -(4 pi / c) j_par).
    !!
    !! Outer equation (low beta, cylinder, ideal: E_par = 0, xi_r = Br/(i kp B0)):
    !!     (1/r)(r A')' - kperp2 A - C A = 0,
    !!     C = ks J'/(kp B0),  J = (4 pi/c) j_par0 = hz (r B0th)'/r - hth B0z',
    !!     ks = hz m/r - hth n/R0,  kp = hth m/r + hz n/R0,
    !!     kperp2 = m^2/r^2 + n^2/R0^2 (= ks^2 + kp^2; exact for the vacuum A_z).
    !! The solution is regular on the axis (A ~ I_m(k r), k = n/R0, started at
    !! the first grid point) and, beyond the last grid point, either decays in
    !! vacuum (A ~ K_m(k r), no wall) or vanishes on an ideal wall at r_wall.
    !! Near the resonance kp(r_res) = 0, C ~ c0/x with x = r - r_res, and each
    !! side is A = a (1 + c0 x ln|x|) + b x; Delta' = b_R/a_R - b_L/a_L [1/cm].
    !! The one-sided values delta_prime_right = b_R/a_R and delta_prime_left =
    !! b_L/a_L are returned as well (x in cm: their mean depends on the length
    !! unit through c0 ln|x|; only the difference is unit-free).
    !!
    !! The equilibrium enters only through natural cubic splines of r*B0th and
    !! B0z on the given grid. Module state holds the splines during one call;
    !! the routine is not thread-safe.
    use KIM_kinds_m, only: dp
    use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
    implicit none
    private

    integer, parameter, public :: DELTA_PRIME_OK = 0
    integer, parameter, public :: DELTA_PRIME_INVALID_INPUT = 1
    integer, parameter, public :: DELTA_PRIME_NO_UNIQUE_RESONANCE = 2
    integer, parameter, public :: DELTA_PRIME_ODE_FAILED = 3

    type, public :: delta_prime_result_t
        real(dp) :: delta_prime = 0.0_dp   !! [1/cm]
        real(dp) :: delta_prime_right = 0.0_dp  !! b_R/a_R [1/cm], x in cm
        real(dp) :: delta_prime_left = 0.0_dp   !! b_L/a_L [1/cm], x in cm
        real(dp) :: r_res = 0.0_dp         !! kp(r_res) = 0 [cm]
        real(dp) :: ks_res = 0.0_dp        !! ks(r_res) [1/cm]
        real(dp) :: c0 = 0.0_dp            !! residue of C at r_res [1/cm]
        integer :: status = DELTA_PRIME_INVALID_INPUT
    end type delta_prime_result_t

    public :: cylindrical_delta_prime, delta_prime_current

    ! spline state of the current call
    real(dp), allocatable :: sr(:), sf(:), sf2(:), sz(:), sz2(:)
    real(dp) :: s_R0 = 0.0_dp, s_k = 0.0_dp, s_r_edge = 0.0_dp
    integer :: s_m = 0, s_n = 0
    logical :: s_current = .true.

contains

    subroutine cylindrical_delta_prime(r, B0z, B0th, R0, m, n, result, r_wall, eps, &
            include_current)
        use fortnum_ode_dop853, only: ode_solve_dop
        use fortnum_status, only: fortnum_status_t, FORTNUM_OK
        use fortnum_special_bessel, only: bessel_in, bessel_kn
        real(dp), intent(in) :: r(:), B0z(:), B0th(:), R0
        integer, intent(in) :: m, n
        type(delta_prime_result_t), intent(out) :: result
        real(dp), intent(in), optional :: r_wall    !! <= 0 or absent: no wall
        real(dp), intent(in), optional :: eps       !! Frobenius distance [cm], default 1e-5
        logical, intent(in), optional :: include_current  !! false drops C (tests)
        real(dp) :: e, x, rb, y0(2), yl(2), yr(2), d(5), wall
        real(dp), allocatable :: t_out(:), y_out(:,:)
        type(fortnum_status_t) :: ode_status
        integer :: npts

        result = delta_prime_result_t()
        npts = size(r)
        e = 1.0e-5_dp
        if (present(eps)) e = eps
        wall = -1.0_dp
        if (present(r_wall)) wall = r_wall
        s_current = .true.
        if (present(include_current)) s_current = include_current
        if (npts < 4 .or. size(B0z) /= npts .or. size(B0th) /= npts) return
        if (.not. all(ieee_is_finite(r)) .or. .not. all(ieee_is_finite(B0z)) .or. &
            .not. all(ieee_is_finite(B0th)) .or. .not. ieee_is_finite(R0)) return
        if (any(r(2:) <= r(:npts - 1)) .or. r(1) <= 0.0_dp .or. R0 <= 0.0_dp) return
        if (m == 0 .or. n == 0 .or. e <= 0.0_dp) return
        if (wall > 0.0_dp .and. wall <= r(npts)) return

        sr = r
        sf = r*B0th
        sz = B0z
        call natural_spline(sr, sf, sf2)
        call natural_spline(sr, sz, sz2)
        s_R0 = R0
        s_m = m
        s_n = n
        s_k = real(abs(n), dp)/R0
        s_r_edge = r(npts)

        call find_resonance(result%r_res, result%status)
        if (result%status /= DELTA_PRIME_OK) return
        if (result%r_res - e <= r(1) .or. result%r_res + e >= r(npts)) then
            result%status = DELTA_PRIME_NO_UNIQUE_RESONANCE
            return
        end if
        call local_geometry(result%r_res, d)
        result%ks_res = d(1)
        result%c0 = 0.0_dp
        if (s_current) result%c0 = d(1)*d(4)/(d(3)*d(5))   ! ks J' / (B0 kp')

        ! inner: regular solution, normalised to A = 1
        x = s_k*r(1)
        y0 = [1.0_dp, s_k*(bessel_in(abs(m) - 1, x)/bessel_in(abs(m), x) - abs(m)/x)]
        call ode_solve_dop(outer_rhs, r(1), result%r_res - e, y0, t_out, y_out, &
            ode_status, rtol=1.0e-11_dp, atol=1.0e-30_dp)
        if (ode_status%code /= FORTNUM_OK .or. size(t_out) < 2) then
            result%status = DELTA_PRIME_ODE_FAILED
            return
        end if
        yl = y_out(:, size(t_out))

        ! outer: vacuum decay or ideal wall, integrated inward
        if (wall > 0.0_dp) then
            rb = wall
            y0 = [0.0_dp, -1.0_dp]
        else
            rb = r(npts)
            x = s_k*rb
            y0 = [1.0_dp, s_k*(-bessel_kn(abs(m) - 1, x)/bessel_kn(abs(m), x) - abs(m)/x)]
        end if
        call ode_solve_dop(outer_rhs, rb, result%r_res + e, y0, t_out, y_out, &
            ode_status, rtol=1.0e-11_dp, atol=1.0e-30_dp)
        if (ode_status%code /= FORTNUM_OK .or. size(t_out) < 2) then
            result%status = DELTA_PRIME_ODE_FAILED
            return
        end if
        yr = y_out(:, size(t_out))

        result%delta_prime_right = small_over_large(yr, e, result%c0)
        result%delta_prime_left = small_over_large(yl, -e, result%c0)
        result%delta_prime = result%delta_prime_right - result%delta_prime_left
        if (.not. ieee_is_finite(result%delta_prime_right) .or. &
                .not. ieee_is_finite(result%delta_prime_left)) then
            result%status = DELTA_PRIME_ODE_FAILED
            return
        end if
        result%status = DELTA_PRIME_OK
    end subroutine cylindrical_delta_prime

    pure complex(dp) function delta_prime_current(delta_prime, r_res, ks_res, c_light) &
            result(u_delta)
        !! Outer-region current per unit Br at r_res in the periodic
        !! normalization's units (int j_par r dr, unit drive Br = 1 G):
        !! scale = I c / (2 pi (u_unit - u_delta)). Derivation: Ampere jump
        !! [A'] = -(4 pi/c) int j_par dr with int j_par dr = u/r_res,
        !! A(r_res) = Br/(i ks), and the ideal limit int j_par r dr = I c/(2 pi).
        real(dp), intent(in) :: delta_prime, r_res, ks_res, c_light
        u_delta = cmplx(0.0_dp, c_light*delta_prime*r_res/(4.0_dp*acos(-1.0_dp)*ks_res), dp)
    end function delta_prime_current

    real(dp) function small_over_large(y, x, c0) result(ratio)
        !! Solve [1 + c0 x ln|x|, x; c0 (ln|x| + 1), 1] (a, b) = (A, A') for b/a.
        real(dp), intent(in) :: y(2), x, c0
        real(dp) :: l, m11, m12, m21, m22, det, a, b
        l = log(abs(x))
        m11 = 1.0_dp + c0*x*l
        m12 = x
        m21 = c0*(l + 1.0_dp)
        m22 = 1.0_dp
        det = m11*m22 - m12*m21
        a = (m22*y(1) - m12*y(2))/det
        b = (-m21*y(1) + m11*y(2))/det
        ratio = b/a
    end function small_over_large

    subroutine outer_rhs(t, y, dydt, ctx)
        real(dp), intent(in) :: t
        real(dp), intent(in) :: y(:)
        real(dp), intent(out) :: dydt(:)
        class(*), intent(in), optional :: ctx
        real(dp) :: coef, d(5)
        coef = (real(s_m, dp)/t)**2 + s_k**2
        if (s_current .and. t <= s_r_edge) then
            call local_geometry(t, d)
            coef = coef + d(1)*d(4)/(d(2)*d(3))   ! ks J' / (kp B0)
        end if
        dydt(1) = y(2)
        dydt(2) = -y(2)/t + coef*y(1)
    end subroutine outer_rhs

    subroutine local_geometry(rr, d)
        !! d = [ks, kp, B0, J', kp'] at rr from the splines (J = (4 pi/c) j_par0).
        real(dp), intent(in) :: rr
        real(dp), intent(out) :: d(5)
        real(dp) :: f(0:3), z(0:3), bt, bt1, b, b1, hz, hth, hz1, hth1, mm, nn
        call spline_eval(sr, sf, sf2, rr, f)
        call spline_eval(sr, sz, sz2, rr, z)
        mm = real(s_m, dp)
        nn = real(s_n, dp)
        bt = f(0)/rr
        bt1 = f(1)/rr - f(0)/rr**2
        b = hypot(z(0), bt)
        b1 = (z(0)*z(1) + bt*bt1)/b
        hz = z(0)/b
        hth = bt/b
        hz1 = (z(1)*b - z(0)*b1)/b**2
        hth1 = (bt1*b - bt*b1)/b**2
        d(1) = hz*mm/rr - hth*nn/s_R0
        d(2) = hth*mm/rr + hz*nn/s_R0
        d(3) = b
        d(4) = hz1*f(1)/rr + hz*(f(2)/rr - f(1)/rr**2) - hth1*z(1) - hth*z(2)
        d(5) = hth1*mm/rr - hth*mm/rr**2 + hz1*nn/s_R0
    end subroutine local_geometry

    subroutine find_resonance(r_res, status)
        real(dp), intent(out) :: r_res
        integer, intent(out) :: status
        real(dp) :: d(5), a, b, c, fa, fc
        real(dp), allocatable :: kp(:)
        integer :: i, found, it
        allocate(kp(size(sr)))
        do i = 1, size(sr)
            call local_geometry(sr(i), d)
            kp(i) = d(2)
        end do
        found = 0
        status = DELTA_PRIME_NO_UNIQUE_RESONANCE
        r_res = 0.0_dp
        do i = 1, size(sr) - 1
            if (kp(i) == 0.0_dp .or. kp(i)*kp(i + 1) < 0.0_dp) then
                found = found + 1
                a = sr(i)
                b = sr(i + 1)
            end if
        end do
        if (found /= 1) return
        call local_geometry(a, d)
        fa = d(2)
        do it = 1, 200
            c = 0.5_dp*(a + b)
            call local_geometry(c, d)
            fc = d(2)
            if (fa*fc <= 0.0_dp) then
                b = c
            else
                a = c
                fa = fc
            end if
            if (b - a < 1.0e-13_dp*c) exit
        end do
        r_res = 0.5_dp*(a + b)
        status = DELTA_PRIME_OK
    end subroutine find_resonance

    subroutine natural_spline(x, y, y2)
        real(dp), intent(in) :: x(:), y(:)
        real(dp), allocatable, intent(out) :: y2(:)
        real(dp), allocatable :: u(:)
        real(dp) :: sig, p
        integer :: i, nn
        nn = size(x)
        allocate(y2(nn), u(nn))
        y2(1) = 0.0_dp
        u(1) = 0.0_dp
        do i = 2, nn - 1
            sig = (x(i) - x(i - 1))/(x(i + 1) - x(i - 1))
            p = sig*y2(i - 1) + 2.0_dp
            y2(i) = (sig - 1.0_dp)/p
            u(i) = (y(i + 1) - y(i))/(x(i + 1) - x(i)) - (y(i) - y(i - 1))/(x(i) - x(i - 1))
            u(i) = (6.0_dp*u(i)/(x(i + 1) - x(i - 1)) - sig*u(i - 1))/p
        end do
        y2(nn) = 0.0_dp
        do i = nn - 1, 1, -1
            y2(i) = y2(i)*y2(i + 1) + u(i)
        end do
    end subroutine natural_spline

    subroutine spline_eval(x, y, y2, xx, v)
        !! v = [s, s', s'', s'''] of the natural cubic spline at xx (clamped to the grid).
        real(dp), intent(in) :: x(:), y(:), y2(:), xx
        real(dp), intent(out) :: v(0:3)
        integer :: lo, hi, mid
        real(dp) :: h, a, b, xc
        xc = min(max(xx, x(1)), x(size(x)))
        lo = 1
        hi = size(x)
        do while (hi - lo > 1)
            mid = (lo + hi)/2
            if (x(mid) > xc) then
                hi = mid
            else
                lo = mid
            end if
        end do
        h = x(hi) - x(lo)
        a = (x(hi) - xc)/h
        b = (xc - x(lo))/h
        v(0) = a*y(lo) + b*y(hi) + ((a**3 - a)*y2(lo) + (b**3 - b)*y2(hi))*h**2/6.0_dp
        v(1) = (y(hi) - y(lo))/h - (3.0_dp*a**2 - 1.0_dp)/6.0_dp*h*y2(lo) &
            + (3.0_dp*b**2 - 1.0_dp)/6.0_dp*h*y2(hi)
        v(2) = a*y2(lo) + b*y2(hi)
        v(3) = (y2(hi) - y2(lo))/h
    end subroutine spline_eval

end module cylindrical_delta_prime_m
