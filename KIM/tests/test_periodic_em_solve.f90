program test_periodic_em_solve
    !! periodic_solve_m::solve_periodic_em (coupled Poisson-Ampere on the period).
    !!
    !! (a) No plasma current (Kjphi = KjB = 0): Br_m = Br_ref delta_{m,0} and Phi_m
    !!     equal the electrostatic solve_periodic result, for both psi definitions.
    !! (b) Thin local conducting layer j = sigma(x) Br, sigma Gaussian of width
    !!     delta at r_m, weak coupling eps: to first order in eps (ks L -> 0)
    !!       A(x) = psi - (2 pi/c) i ks psi int |x - x'| sigma(x') dx'
    !!     (kink of the local current; the compensating sheet at the period
    !!     boundary removes the far field), with the 'extrapolated' psi. Checked
    !!     on the core, error O(eps^2): halving eps must reduce it ~4x.
    !! (c) 'resonance' psi: A(r_m) = psi exactly.
    use KIM_kinds_m, only: dp
    use constants_m, only: pi, sol, com_unit
    use periodic_solve_m, only: solve_periodic, solve_periodic_em, reconstruct_delta_phi, &
        core_integral_weights, PERIODIC_PSI_EXTRAPOLATED, PERIODIC_PSI_RESONANCE
    implicit none

    integer :: failures
    failures = 0

    call test_no_current()
    call test_thin_layer()
    call test_core_weights()

    if (failures > 0) then
        print *, failures, ' failures'
        error stop 1
    end if
    print *, 'test_periodic_em_solve: all checks passed'

contains

    subroutine check(name, err, tol)
        character(*), intent(in) :: name
        real(dp), intent(in) :: err, tol
        if (.not. (err <= tol)) then
            failures = failures + 1
            write(*, '(A,A,2ES12.4)') 'FAIL ', name, err, tol
        else
            write(*, '(A,A,2ES12.4)') 'ok   ', name, err, tol
        end if
    end subroutine check

    subroutine test_no_current()
        integer, parameter :: M = 8, dim = 2*M + 1
        complex(dp) :: Kphi(dim, dim), KB(dim, dim), Z(dim, dim)
        complex(dp), allocatable :: Phi_es(:), Phi_em(:), Br_m(:)
        complex(dp) :: Br_ref, expect(dim)
        real(dp) :: re(dim, dim), im(dim, dim), L, rm
        integer :: info, i, def

        call random_seed()
        call random_number(re)
        call random_number(im)
        Kphi = cmplx(re - 0.5_dp, im - 0.5_dp, dp)
        do i = 1, dim
            Kphi(i, i) = Kphi(i, i) + 5.0_dp
        end do
        call random_number(re)
        call random_number(im)
        KB = cmplx(re, im, dp)
        Z = (0.0_dp, 0.0_dp)
        L = 3.0_dp
        rm = 55.0_dp
        Br_ref = (0.7_dp, -0.3_dp)
        call solve_periodic(Kphi, KB, L, M, Br_ref, Phi_es, info)
        call check('ES solve info', real(abs(info), dp), 0.0_dp)
        expect = (0.0_dp, 0.0_dp)
        expect(M + 1) = Br_ref
        do def = PERIODIC_PSI_EXTRAPOLATED, PERIODIC_PSI_RESONANCE
            call solve_periodic_em(Kphi, KB, Z, Z, L, M, rm, 0.5_dp, -0.11_dp, Br_ref, &
                def, Phi_em, Br_m, info)
            call check('EM solve info', real(abs(info), dp), 0.0_dp)
            call check('no current: Br_m = Br_ref delta', maxval(abs(Br_m - expect)), &
                1.0e-14_dp)
            call check('no current: Phi_m = electrostatic', &
                maxval(abs(Phi_em - Phi_es))/maxval(abs(Phi_es)), 1.0e-12_dp)
        end do
    end subroutine test_no_current

    !> Error of the first-order kink solution on the core for coupling eps.
    real(dp) function thin_layer_error(eps, def, at_rm) result(err)
        real(dp), intent(in) :: eps
        integer, intent(in) :: def
        complex(dp), intent(out) :: at_rm
        integer, parameter :: M = 64, dim = 2*M + 1, nchk = 81
        real(dp), parameter :: L = 8.0_dp, rm = 55.0_dp, x_a = 2.0_dp, delta = 0.25_dp
        real(dp), parameter :: ks = 1.0e-3_dp
        complex(dp) :: Kphi(dim, dim), KB(dim, dim), Kjphi(dim, dim), KjB(dim, dim)
        complex(dp), allocatable :: Phi_m(:), Br_m(:)
        complex(dp) :: Br_ref, br(nchk), shat, analytic
        real(dp) :: sigma0, x(nchk), f, kq
        integer :: i, j, info

        ! eps = (2 pi/c) ks sigma0 sqrt(2 pi) delta x_a
        sigma0 = eps*sol/(2.0_dp*pi*ks*sqrt(2.0_dp*pi)*delta*x_a)
        Kphi = (0.0_dp, 0.0_dp)
        do i = 1, dim
            Kphi(i, i) = (1.0_dp, 0.0_dp)
        end do
        KB = (0.0_dp, 0.0_dp)
        Kjphi = (0.0_dp, 0.0_dp)
        ! KjB(m, m') = sigma_hat(m - m'), sigma_hat(q) = (1/L) int sigma exp(-i k_q r) dr
        do i = -M, M
            do j = -M, M
                kq = 2.0_dp*pi*real(i - j, dp)/L
                shat = sigma0*sqrt(2.0_dp*pi)*delta/L*exp(-0.5_dp*(kq*delta)**2) &
                    *exp(-com_unit*kq*rm)
                KjB(i + M + 1, j + M + 1) = shat
            end do
        end do
        Br_ref = (1.0_dp, 0.0_dp)
        call solve_periodic_em(Kphi, KB, Kjphi, KjB, L, M, rm, x_a, ks, Br_ref, def, &
            Phi_m, Br_m, info)
        if (info /= 0) error stop 'thin layer: solve failed'
        do i = 1, nchk
            x(i) = -x_a + 2.0_dp*x_a*real(i - 1, dp)/real(nchk - 1, dp)
        end do
        br = reconstruct_delta_phi(Br_m, L, M, rm + x)
        err = 0.0_dp
        do i = 1, nchk
            f = x(i)*erf(x(i)/(sqrt(2.0_dp)*delta)) &
                + delta*sqrt(2.0_dp/pi)*exp(-0.5_dp*(x(i)/delta)**2)
            ! Br/Br_ref = A/psi = 1 - (2 pi/c) i ks sigma0 sqrt(2 pi) delta f(x)
            analytic = 1.0_dp - com_unit*eps*f/x_a
            err = max(err, abs(br(i)/Br_ref - analytic))
        end do
        at_rm = sum(reconstruct_delta_phi(Br_m, L, M, [rm]))/Br_ref
    end function thin_layer_error

    subroutine test_thin_layer()
        real(dp) :: e1, e2
        complex(dp) :: at_rm, dummy

        e1 = thin_layer_error(1.0e-2_dp, PERIODIC_PSI_EXTRAPOLATED, dummy)
        e2 = thin_layer_error(0.5e-2_dp, PERIODIC_PSI_EXTRAPOLATED, dummy)
        ! first-order term is eps (|f|/x_a <= 1); the remainder must be O(eps^2)
        call check('thin layer eps=1e-2 vs first order', e1, 3.0e-4_dp)
        call check('thin layer error ratio (eps halved) ~ 4', abs(e1/e2 - 4.0_dp), 0.3_dp)
        e1 = thin_layer_error(1.0e-2_dp, PERIODIC_PSI_RESONANCE, at_rm)
        call check('resonance psi: Br(r_m) = Br_ref', abs(at_rm - 1.0_dp), 1.0e-12_dp)
    end subroutine test_thin_layer

    subroutine test_core_weights()
        integer, parameter :: M = 5
        real(dp), parameter :: L = 4.0_dp, rm = 10.0_dp, x_a = 0.7_dp
        complex(dp) :: f(2*M + 1), w(2*M + 1), exact
        f = (0.0_dp, 0.0_dp)
        f(M + 1 + 3) = (2.0_dp, 1.0_dp)
        f(M + 1) = (0.5_dp, 0.0_dp)
        w = core_integral_weights(L, M, rm, x_a)
        ! int exp(i k r) over [rm - x_a, rm + x_a] = (exp(i k (rm+x_a)) - exp(i k (rm-x_a)))/(i k)
        exact = 0.5_dp*2.0_dp*x_a + (2.0_dp, 1.0_dp)* &
            (exp(com_unit*2.0_dp*pi*3.0_dp/L*(rm + x_a)) &
            - exp(com_unit*2.0_dp*pi*3.0_dp/L*(rm - x_a)))/(com_unit*2.0_dp*pi*3.0_dp/L)
        call check('core integral weights', abs(sum(w*f) - exact), 1.0e-13_dp)
    end subroutine test_core_weights

end program test_periodic_em_solve
