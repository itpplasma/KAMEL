program test_grid_equidistant
    !> Unit test for grid_m equidistant grid generation: uniform spacing,
    !> first boundary at min_val, and cell centres at boundary midpoints.
    !> Exercises grid_init_equidistant + grid_generate_equidistant, which in
    !> turn drive binsrc and plag_coeff for the difference-operator stencils.
    use KIM_kinds_m, only: dp
    use grid_m, only: grid_type
    use kim_resonances_m, only: r_res
    implicit none

    type(grid_type) :: g
    integer :: npts, i
    real(dp) :: h
    logical :: all_passed, uniform
    real(dp), parameter :: tol = 1.0e-12_dp

    all_passed = .true.
    r_res = 0.0_dp      ! deterministic: avoids reading uninitialised module state in binsrc

    npts = 11
    call g%grid_init_equidistant(npts, 1.0_dp, 6.0_dp, 'test')
    call g%grid_generate_equidistant()

    h = (6.0_dp - 1.0_dp) / 11.0_dp

    call report('npts_b set to requested count', g%npts_b == 11, all_passed)
    call report('npts_c = npts_b - 1', g%npts_c == 10, all_passed)
    call report('first boundary == min_val', abs(g%xb(1) - 1.0_dp) < tol, all_passed)

    uniform = .true.
    do i = 2, g%npts_b
        if (abs((g%xb(i) - g%xb(i-1)) - h) > tol) uniform = .false.
    end do
    call report('uniform boundary spacing h', uniform, all_passed)

    call report('cell centre is boundary midpoint', &
                abs(g%xc(1) - 0.5_dp * (g%xb(1) + g%xb(2))) < tol, all_passed)

    if (all_passed) then
        print *, 'All grid tests PASSED'
        stop 0
    else
        print *, 'Some grid tests FAILED'
        stop 1
    end if

contains

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

end program test_grid_equidistant
