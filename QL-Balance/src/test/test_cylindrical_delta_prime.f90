program test_cylindrical_delta_prime
    !! Cylindrical outer-region Delta' (cylindrical_delta_prime_m).
    !! Synthetic equilibrium: R0 = 165 cm, B0z = -1.8e4 G, q = -(1 + 3 (r/60)^2),
    !! (m, n) = (6, 2), resonance q = -3 at r = 60 sqrt(2/3) cm.
    !! References: analytic Bessel solutions without the current term; with it,
    !! the independent Python implementation (scipy DOP853, natural splines,
    !! rmp_integral_ion_tensor_kamel/scripts/delta_prime_cyl.py).
    use QLBalance_kinds, only: dp
    use cylindrical_delta_prime_m
    use fortnum_special_bessel, only: bessel_in, bessel_kn
    implicit none
    integer, parameter :: npts = 4000, m = 6, n = 2
    real(dp), parameter :: R0 = 165.0_dp, bz = -1.8e4_dp, c_light = 2.99792458e10_dp
    real(dp) :: r(npts), q(npts), B0z(npts), B0th(npts), k, x, analytic, rw, f, df
    type(delta_prime_result_t) :: res, res2
    complex(dp) :: u_delta
    integer :: i, failures

    failures = 0
    do i = 1, npts
        r(i) = 3.0_dp + 67.0_dp*real(i - 1, dp)/real(npts - 1, dp)
    end do
    q = -(1.0_dp + 3.0_dp*(r/60.0_dp)**2)
    B0z = bz
    B0th = r*B0z/(q*R0)
    k = real(n, dp)/R0

    ! resonance and local geometry
    call cylindrical_delta_prime(r, B0z, B0th, R0, m, n, res)
    call check('status', real(res%status, dp), 0.0_dp, 0.0_dp)
    call check('r_res', res%r_res, 60.0_dp*sqrt(2.0_dp/3.0_dp), 1.0e-11_dp)
    call check('ks_res', res%ks_res, -1.230728393403168e-01_dp, 1.0e-9_dp)
    call check('c0', res%c0, 2.718355275252230e-02_dp, 1.0e-6_dp)

    ! no current term, no wall: k (K_m'/K_m - I_m'/I_m) at r_res
    x = k*res%r_res
    analytic = k*((-bessel_kn(m - 1, x)/bessel_kn(m, x) - m/x) &
        - (bessel_in(m - 1, x)/bessel_in(m, x) - m/x))
    call cylindrical_delta_prime(r, B0z, B0th, R0, m, n, res2, include_current=.false.)
    call check('vacuum, no wall', res2%delta_prime, analytic, 5.0e-6_dp)

    ! no current term, ideal wall at 80 cm
    rw = 80.0_dp
    f = bessel_in(m, k*rw)*bessel_kn(m, x) - bessel_kn(m, k*rw)*bessel_in(m, x)
    df = k*(bessel_in(m, k*rw)*(-bessel_kn(m - 1, x) - m/x*bessel_kn(m, x)) &
        - bessel_kn(m, k*rw)*(bessel_in(m - 1, x) - m/x*bessel_in(m, x)))
    analytic = df/f - k*(bessel_in(m - 1, x)/bessel_in(m, x) - m/x)
    call cylindrical_delta_prime(r, B0z, B0th, R0, m, n, res2, r_wall=rw, include_current=.false.)
    call check('vacuum, wall 80', res2%delta_prime, analytic, 5.0e-6_dp)

    ! with the equilibrium current (Python reference, same eps)
    call check('current, no wall', res%delta_prime, -2.282741958666485e-01_dp, 1.0e-7_dp)
    call cylindrical_delta_prime(r, B0z, B0th, R0, m, n, res2, eps=1.0e-4_dp)
    call check('current, eps 1e-4', res2%delta_prime, -2.282722764238130e-01_dp, 1.0e-7_dp)
    call check('eps convergence', res2%delta_prime, res%delta_prime, 2.0e-5_dp)

    ! matching current: u_delta = i c Delta' r_res / (4 pi ks_res)
    u_delta = delta_prime_current(res%delta_prime, res%r_res, res%ks_res, c_light)
    call check('u_delta real part', real(u_delta, dp), 0.0_dp, 0.0_dp)
    call check('u_delta', aimag(u_delta), &
        c_light*res%delta_prime*res%r_res/(4.0_dp*acos(-1.0_dp)*res%ks_res), 1.0e-14_dp)

    ! invalid input and missing resonance
    call cylindrical_delta_prime(r, B0z, B0th, R0, m, n, res2, r_wall=60.0_dp)
    call check('wall inside grid', real(res2%status, dp), real(DELTA_PRIME_INVALID_INPUT, dp), 0.0_dp)
    call cylindrical_delta_prime(r, B0z, r*B0z/(-1.5_dp*R0), R0, m, n, res2)
    call check('no resonance', real(res2%status, dp), &
        real(DELTA_PRIME_NO_UNIQUE_RESONANCE, dp), 0.0_dp)
    call cylindrical_delta_prime(r(npts:1:-1), B0z, B0th, R0, m, n, res2)
    call check('descending grid', real(res2%status, dp), real(DELTA_PRIME_INVALID_INPUT, dp), 0.0_dp)

    if (failures > 0) then
        write(*, *) failures, ' failures'
        error stop 1
    end if
    write(*, *) 'test_cylindrical_delta_prime: all checks passed'

contains

    subroutine check(name, got, want, rtol)
        character(*), intent(in) :: name
        real(dp), intent(in) :: got, want, rtol
        real(dp) :: err
        err = abs(got - want)
        if (want /= 0.0_dp) err = err/abs(want)
        if (err > rtol .or. err /= err) then
            failures = failures + 1
            write(*, '(A,A,3ES24.15)') 'FAIL ', name, got, want, err
        else
            write(*, '(A,A,3ES24.15)') 'ok   ', name, got, want, err
        end if
    end subroutine check

end program test_cylindrical_delta_prime
