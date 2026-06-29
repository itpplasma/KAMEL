program test_plag_coeff
    !> Unit tests for the Lagrange-coefficient and binary-search utilities in
    !> KIM/src/util/plag_coeff.f90 (binsrc, plag_coeff). Both are external
    !> (non-module) subroutines, called here through an implicit interface.
    !> These routines underpin grid generation and finite-difference operators.
    use KIM_kinds_m, only: dp
    implicit none

    logical :: all_passed
    all_passed = .true.

    call test_binsrc_interior(all_passed)
    call test_binsrc_upper(all_passed)
    call test_plag_partition_of_unity(all_passed)
    call test_plag_interpolation_cubic(all_passed)
    call test_plag_first_derivative(all_passed)
    call test_plag_second_derivative(all_passed)

    if (all_passed) then
        print *, 'All plag_coeff tests PASSED'
        stop 0
    else
        print *, 'Some plag_coeff tests FAILED'
        stop 1
    end if

contains

    !> Test polynomial: f(x) = 2x^3 - 3x^2 + x + 5
    !> f'(x) = 6x^2 - 6x + 1 ; f''(x) = 12x - 6
    pure function f(x) result(y)
        real(dp), intent(in) :: x
        real(dp) :: y
        y = 2.0_dp*x**3 - 3.0_dp*x**2 + x + 5.0_dp
    end function f

    subroutine test_binsrc_interior(passed)
        logical, intent(inout) :: passed
        real(dp) :: p(5)
        integer :: idx
        p = [10.0_dp, 20.0_dp, 30.0_dp, 40.0_dp, 50.0_dp]
        call binsrc(p, 1, 5, 25.0_dp, idx)
        call report('binsrc: 25 lies in (p2,p3] -> idx 3', idx == 3, passed)
    end subroutine test_binsrc_interior

    subroutine test_binsrc_upper(passed)
        logical, intent(inout) :: passed
        real(dp) :: p(5)
        integer :: idx
        p = [10.0_dp, 20.0_dp, 30.0_dp, 40.0_dp, 50.0_dp]
        call binsrc(p, 1, 5, 45.0_dp, idx)
        call report('binsrc: 45 lies in (p4,p5] -> idx 5', idx == 5, passed)
    end subroutine test_binsrc_upper

    subroutine test_plag_partition_of_unity(passed)
        logical, intent(inout) :: passed
        real(dp) :: xp(4), coef(0:2, 4)
        xp = [0.0_dp, 1.0_dp, 2.0_dp, 3.0_dp]
        call plag_coeff(4, 2, 1.5_dp, xp, coef)
        call report('plag: value weights sum to 1', &
                    abs(sum(coef(0, :)) - 1.0_dp) < 1.0e-12_dp, passed)
    end subroutine test_plag_partition_of_unity

    subroutine test_plag_interpolation_cubic(passed)
        logical, intent(inout) :: passed
        real(dp) :: xp(4), coef(0:2, 4), interp
        integer :: i
        xp = [0.0_dp, 1.0_dp, 2.0_dp, 3.0_dp]
        call plag_coeff(4, 2, 1.5_dp, xp, coef)
        interp = 0.0_dp
        do i = 1, 4
            interp = interp + coef(0, i) * f(xp(i))
        end do
        call report('plag: exact value for cubic at x=1.5', &
                    abs(interp - f(1.5_dp)) < 1.0e-9_dp, passed)
    end subroutine test_plag_interpolation_cubic

    subroutine test_plag_first_derivative(passed)
        logical, intent(inout) :: passed
        real(dp) :: xp(4), coef(0:2, 4), d1
        integer :: i
        xp = [0.0_dp, 1.0_dp, 2.0_dp, 3.0_dp]
        call plag_coeff(4, 2, 1.5_dp, xp, coef)
        d1 = 0.0_dp
        do i = 1, 4
            d1 = d1 + coef(1, i) * f(xp(i))
        end do
        ! f'(1.5) = 6*2.25 - 6*1.5 + 1 = 5.5
        call report('plag: exact 1st derivative at x=1.5', &
                    abs(d1 - 5.5_dp) < 1.0e-9_dp, passed)
    end subroutine test_plag_first_derivative

    subroutine test_plag_second_derivative(passed)
        logical, intent(inout) :: passed
        real(dp) :: xp(4), coef(0:2, 4), d2
        integer :: i
        xp = [0.0_dp, 1.0_dp, 2.0_dp, 3.0_dp]
        call plag_coeff(4, 2, 1.5_dp, xp, coef)
        d2 = 0.0_dp
        do i = 1, 4
            d2 = d2 + coef(2, i) * f(xp(i))
        end do
        ! f''(1.5) = 12*1.5 - 6 = 12
        call report('plag: exact 2nd derivative at x=1.5', &
                    abs(d2 - 12.0_dp) < 1.0e-9_dp, passed)
    end subroutine test_plag_second_derivative

    subroutine report(name, ok, passed)
        character(*), intent(in) :: name
        logical, intent(in) :: ok
        logical, intent(inout) :: passed
        if (ok) then
            print *, 'PASS: ', name
        else
            print *, 'FAIL: ', name
            passed = .false.
        end if
    end subroutine report

end program test_plag_coeff
