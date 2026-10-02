program test_ampere_operator
    use KIM_kinds_m, only: dp
    use ampere_operator_m, only: assemble_ampere_weak
    use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
    implicit none
    real(dp), allocatable :: r(:), ks(:), mat(:, :), rhs(:), exact(:)
    integer, allocatable :: pivots(:)
    real(dp) :: scale, t, h, a, b, lr, previous(2), error, k
    integer :: level, n, i, mode, scaling, info
    external :: dgesv

    do scaling = 1, 2
        scale = 100.0_dp**(scaling-1)
        previous = huge(1.0_dp)
        do level = 5, 8
            n = 2**level + 1
            allocate(r(n), ks(n), mat(n, n), rhs(n), exact(n), pivots(n))
            do i = 1, n
                t = real(i-1, dp)/real(n-1, dp)
                r(i) = scale*(1.0_dp+t**1.4_dp)
            end do
            do mode = 1, 2
                if (mode == 1) then
                    ks = 0.0_dp
                    exact = log(r/scale)/log(2.0_dp)
                else
                    ks = 2.0_dp/r
                    exact = ((r/scale)**2-(scale/r)**2)/3.75_dp
                end if
                call assemble_ampere_weak(r, ks, mat)
                rhs = 0.0_dp
                call boundaries(0.0_dp, 1.0_dp)
                call dgesv(n, 1, mat, n, pivots, rhs, n, info)
                if (info /= 0) error stop 'Ampere vacuum solve failed'
                if (.not. all(ieee_is_finite(rhs))) error stop 'Nonfinite Ampere vacuum field'
                error = maxval(abs(rhs-exact))
                print *, 'vacuum scale, nodes, mode, error: ', scale, n, mode, error
                if (level > 5) then
                    if (error > 0.3_dp*previous(mode)) &
                        error stop 'Ampere vacuum does not converge at second order'
                end if
                previous(mode) = error
                if (level == 8) then
                    if (error > 2.0e-5_dp) error stop 'Ampere vacuum accuracy'
                end if
            end do

            ! Analytic load for A=r, L A=1/r-k**2*r. It represents a
            ! physical prescribed current, independently of the assembler.
            k = 0.3_dp/scale
            ks = k
            call assemble_ampere_weak(r, ks, mat)
            rhs = 0.0_dp
            do i = 1, n-1
                a = r(i)
                b = r(i+1)
                h = b-a
                lr = log(b/a)
                rhs(i) = rhs(i)+(b*lr-h)/h-k*k*h*(2*a+b)/6.0_dp
                rhs(i+1) = rhs(i+1)+(h-a*lr)/h-k*k*h*(a+2*b)/6.0_dp
            end do
            call boundaries(r(1), r(n))
            call dgesv(n, 1, mat, n, pivots, rhs, n, info)
            if (info /= 0) error stop 'Ampere manufactured-current solve failed'
            if (.not. all(ieee_is_finite(rhs))) error stop 'Nonfinite Ampere current field'
            error = maxval(abs(rhs-r))/scale
            if (error > 2.0e-8_dp) error stop 'Ampere current normalization failed'
            deallocate(r, ks, mat, rhs, exact, pivots)
        end do
    end do
    print *, 'Ampere analytic vacuum and current tests passed'

contains

    subroutine boundaries(left, right)
        real(dp), intent(in) :: left, right
        mat(1, :) = 0.0_dp
        mat(1, 1) = 1.0_dp
        rhs(1) = left
        mat(n, :) = 0.0_dp
        mat(n, n) = 1.0_dp
        rhs(n) = right
    end subroutine boundaries
end program test_ampere_operator
