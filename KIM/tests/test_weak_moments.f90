program test_weak_moments
    use KIM_kinds_m, only: dp
    use weak_moments_m, only: weak_moment_projector_t
    use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
    implicit none
    real(dp), parameter :: nodes(3) = [-sqrt(3.0_dp/5.0_dp), 0.0_dp, &
        sqrt(3.0_dp/5.0_dp)]
    real(dp), parameter :: weights(3) = [5.0_dp/9.0_dp, 8.0_dp/9.0_dp, &
        5.0_dp/9.0_dp]
    real(dp) :: r(7), shape(2), h, x
    complex(dp) :: exact(7), load(7), at_x
    complex(dp), allocatable :: recovered(:)
    type(weak_moment_projector_t) :: projector
    integer :: i, k, scaling

    do scaling = 1, 2
        r = [1.0_dp, 1.1_dp, 1.35_dp, 1.8_dp, 2.0_dp, 2.9_dp, 4.0_dp] &
            *100.0_dp**(scaling-1)
        exact = cmplx(sin(r/100.0_dp**(scaling-1)), &
            cos(r/100.0_dp**(scaling-1)), dp)
        ! Independent quadrature of piecewise-linear physical moments.
        ! Both endpoints are nonzero: a Dirichlet mass mutant must fail.
        load = cmplx(0.0_dp, 0.0_dp, dp)
        do i = 1, 6
            h = r(i+1)-r(i)
            do k = 1, 3
                x = (nodes(k)+1.0_dp)/2.0_dp
                shape(1) = 1.0_dp-x
                shape(2) = x
                at_x = shape(1)*exact(i)+shape(2)*exact(i+1)
                load(i:i+1) = load(i:i+1)+shape*at_x*weights(k)*h/2.0_dp
            end do
        end do
        call projector%init(r)
        call projector%project(load, recovered)
        if (.not. all(ieee_is_finite(real(recovered, dp)))) &
            error stop 'Nonfinite real physical moment'
        if (.not. all(ieee_is_finite(aimag(recovered)))) &
            error stop 'Nonfinite imaginary physical moment'
        if (maxval(abs(recovered-exact)) > 2.0e-14_dp) &
            error stop 'Physical moment projection failed quadrature oracle'
    end do
    print *, 'Raw weak moment projection passed independent quadrature'
end program test_weak_moments
