program test_grid_equidistant
    !> Unit test for grid_m equidistant grid generation: uniform spacing,
    !> first boundary at min_val, and cell centres at boundary midpoints.
    !> Exercises grid_init_equidistant + grid_generate_equidistant, which in
    !> turn drive binsrc and plag_coeff for the difference-operator stencils.
    use KIM_kinds_m, only: dp
    use grid_m, only: grid_type, rg_grid, xl_grid, rg_space_dim, l_space_dim, &
        r_min, r_plas, grid_spacing_rg, grid_spacing_xl
    use config_m, only: output_path, hdf5_output
    use kim_resonances_m, only: r_res
    implicit none

    type(grid_type) :: g
    integer :: npts, i, j, mode, fixture
    real(dp) :: h, lower, upper, integral
    complex(dp) :: fourier_sum
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

    ! Independent periodic oracle: distinct samples integrate nonzero Fourier
    ! harmonics to zero. A duplicated physical endpoint violates this identity.
    do mode = 1, 5
        fourier_sum = sum(exp(cmplx(0.0_dp, 1.0_dp, dp)* &
            2.0_dp*acos(-1.0_dp)*real(mode, dp)*(g%xb-1.0_dp)/5.0_dp))
        call report('periodic Fourier orthogonality', &
            abs(fourier_sum) < 1.0e-10_dp, all_passed)
    end do
    do fixture = 1, 3
        select case (fixture)
        case (1)
            lower = -2.0_dp
            upper = 5.0_dp
            npts = 11
        case (2)
            lower = 3.0_dp
            upper = 67.0_dp
            npts = 17
        case (3)
            ! Arithmetic reconstruction of this endpoint overshoots by one ULP.
            lower = 0.1_dp
            upper = 0.3_dp
            npts = 13
        end select
        call g%grid_init_equidistant(npts, lower, upper, 'closed')
        call g%grid_generate_equidistant(endpoint_inclusive=.true.)
        call report('closed interval reaches both bounds', &
            g%xb(1) == lower .and. g%xb(npts) == upper, all_passed)
        integral = sum(g%xb(2:)-g%xb(:npts-1))
        call report('integrated constant over declared interval', &
            abs(integral-(upper-lower)) < tol, all_passed)
        integral = sum((g%xb(2:)-g%xb(:npts-1))* &
            (2.0_dp+3.0_dp*g%xc))
        call report('integrated affine function over declared interval', &
            abs(integral-(2.0_dp*(upper-lower)+ &
            1.5_dp*(upper**2-lower**2))) < 1.0e-9_dp, all_passed)
        do j = 1, g%npts_c
            integral = dot_product(g%deriv_coef(:, j), &
                g%xc(g%ipbeg(j):g%ipend(j))**2)
            call report('quadratic derivative at physical boundary node', &
                abs(integral-2.0_dp*g%xb(j)) < 1.0e-9_dp, all_passed)
        end do
    end do

    ! Exercise production routing with unequal field/background resolutions.
    ! Endpoint correctness is a physical-domain oracle, not source inspection.
    rg_space_dim = 13
    l_space_dim = 17
    r_min = 0.1_dp
    r_plas = 0.3_dp
    grid_spacing_rg = 'equidistant'
    grid_spacing_xl = 'equidistant'
    output_path = './grid-contract-output/'
    hdf5_output = .false.
    call execute_command_line('mkdir -p grid-contract-output/grid')
    call generate_grids()
    call report('production background spans declared domain', &
        rg_grid%xb(1) == r_min .and. &
        rg_grid%xb(rg_grid%npts_b) == r_plas, all_passed)
    call report('production field spans declared domain', &
        xl_grid%xb(1) == r_min .and. &
        xl_grid%xb(xl_grid%npts_b) == r_plas, all_passed)

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
