program test_ampere_coupling
    ! Normalization of the Ampere rows against the current kernels.
    !
    ! A local current j = -(c/4 pi) q^2 A_par enters the Ampere rows as the
    ! kernel K_jB = -(c/4 pi) q^2 M / alpha, with M the P1 mass matrix of the
    ! same dr hat moments (Br = alpha A_par). The weak equation then is
    !   A'' + A'/r - q^2 A = 0,  with ks = 0,
    ! solved exactly by c1 I0(q r) + c2 K0(q r). The discrete solution must
    ! converge to it at second order; a mass-inverted kernel, a missing 1/c
    ! or 4 pi, or a dropped alpha leaves an O(1) error.
    use KIM_kinds_m, only: dp
    use constants_m, only: pi, sol
    use ampere_operator_m, only: assemble_ampere_weak, add_ampere_current_coupling
    use fortnum_special_bessel, only: bessel_in, bessel_kn
    use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
    implicit none

    real(dp), parameter :: r_inner = 1.0_dp, r_outer = 2.5_dp, q = 1.7_dp
    complex(dp), parameter :: c1 = (1.0_dp, 0.0_dp), c2 = (0.5_dp, -0.3_dp)
    real(dp), allocatable :: r(:), ks(:), laplace(:, :), mass(:, :)
    complex(dp), allocatable :: block(:, :), K_jB(:, :), alpha(:), rhs(:), exact(:)
    integer, allocatable :: pivots(:)
    real(dp) :: h, error, previous
    integer :: level, n, i, j, info
    external :: zgesv

    previous = huge(1.0_dp)
    do level = 4, 8
        n = 2**level + 1
        allocate(r(n), ks(n), laplace(n, n), mass(n, n), block(n, n), K_jB(n, n), &
                 alpha(n), rhs(n), exact(n), pivots(n))
        h = (r_outer - r_inner) / real(n - 1, dp)
        do i = 1, n
            r(i) = r_inner + h * real(i - 1, dp)
        end do
        ks = 0.0_dp
        alpha = cmplx(1.0_dp + 0.3_dp * r, 0.2_dp * r, dp)

        ! P1 mass matrix of the dr hat moments, element by element.
        mass = 0.0_dp
        do i = 1, n - 1
            mass(i, i) = mass(i, i) + h / 3.0_dp
            mass(i + 1, i + 1) = mass(i + 1, i + 1) + h / 3.0_dp
            mass(i, i + 1) = mass(i, i + 1) + h / 6.0_dp
            mass(i + 1, i) = mass(i + 1, i) + h / 6.0_dp
        end do
        do j = 1, n
            K_jB(:, j) = -sol / (4.0_dp * pi) * q**2 * mass(:, j) / alpha(j)
        end do

        call assemble_ampere_weak(r, ks, laplace)
        block = cmplx(laplace, 0.0_dp, dp)
        call add_ampere_current_coupling(block, K_jB, alpha)

        do i = 1, n
            exact(i) = c1 * bessel_in(0, q * r(i)) + c2 * bessel_kn(0, q * r(i))
        end do
        rhs = (0.0_dp, 0.0_dp)
        block(1, :) = (0.0_dp, 0.0_dp)
        block(1, 1) = (1.0_dp, 0.0_dp)
        rhs(1) = exact(1)
        block(n, :) = (0.0_dp, 0.0_dp)
        block(n, n) = (1.0_dp, 0.0_dp)
        rhs(n) = exact(n)

        call zgesv(n, 1, block, n, pivots, rhs, n, info)
        if (info /= 0) error stop 'Coupled Ampere solve failed'
        if (.not. all(ieee_is_finite(real(rhs, dp))) .or. .not. all(ieee_is_finite(aimag(rhs)))) &
            error stop 'Nonfinite coupled Ampere solution'
        error = maxval(abs(rhs - exact)) / maxval(abs(exact))
        print '(a,i5,a,es10.3)', 'nodes', n, '  max rel error', error
        if (level > 4 .and. error > 0.3_dp * previous) &
            error stop 'Coupled Ampere rows do not converge at second order'
        previous = error
        deallocate(r, ks, laplace, mass, block, K_jB, alpha, rhs, exact, pivots)
    end do
    if (previous > 1.0e-5_dp) error stop 'Coupled Ampere rows miss the modified Bessel solution'
    print *, 'PASS: Ampere rows and current kernel share one normalization'
end program test_ampere_coupling
