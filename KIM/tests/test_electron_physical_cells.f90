program test_electron_physical_cells
    use KIM_kinds_m, only: dp
    use constants_m, only: pi, e_charge, ev
    use grid_m, only: rg_grid, xl_grid
    use setup_m, only: spline_base
    use species_m, only: plasma, plasma_t
    use FP_kernel_plasma_prefacs_m, only: FP_G0_rho_phi
    use integrals_gauss_m, only: gauss_config_t
    use kernel_m, only: FP_calc_kernel_element_electrons, &
        pref_rho_phi_g0, pref_rho_phi_g1, pref_rho_B_g1, pref_j_phi_g1, pref_j_B_g1
    implicit none
    integer, parameter :: n = 6, nc = 5
    real(dp), parameter :: density = 1.0e12_dp, temperature = 800.0_dp*ev
    real(dp), parameter :: positions(n) = [0.0_dp, 0.1_dp, 0.45_dp, 0.55_dp, 0.8_dp, 1.0_dp]
    real(dp), parameter :: centers(nc+1) = [0.0_dp, 0.2_dp, 0.35_dp, 0.62_dp, 0.7_dp, 1.0_dp]
    real(dp), parameter :: tolerance = 512.0_dp*epsilon(1.0_dp)
    complex(dp) :: coefficients(nc, 4), matrix(n, n, 4), expected(n, n, 4)
    complex(dp) :: source_coefficients(nc, 4), source(n, n, 4), total, i00, i10, i01, i11, a
    real(dp) :: debye, scale, affine(n), length, x1, x2, vt, offset, width
    complex(dp) :: action, expected_action
    type(gauss_config_t) :: gauss
    integer :: scenario, j, l, lp, channel, checks, degree
    character(32) :: argument

    checks = 0
    spline_base = 1
    call get_command_argument(1, argument)
    plasma = plasma_t(n_species=1, grid_size=nc+1)
    allocate(plasma%spec(0:0))
    allocate(plasma%spec(0)%lambda_D_cc(nc), plasma%spec(0)%rho_L_cc(nc))
    plasma%spec(0)%rho_L_cc = 0.002_dp
    plasma%spec(0)%lambda_D_cc = sqrt(temperature/(4.0_dp*pi*density*e_charge**2))
    debye = density*e_charge**2/temperature
    allocate(pref_rho_phi_g0(1, nc), pref_rho_phi_g1(1, nc, 0:0), &
        pref_rho_B_g1(1, nc, 0:0), pref_j_phi_g1(1, nc, 0:0), pref_j_B_g1(1, nc, 0:0))
    do j = 1, nc
        pref_rho_phi_g0(1, j) = FP_G0_rho_phi(j, plasma%spec(0))
        call expect(abs(pref_rho_phi_g0(1, j)/(4.0_dp*pi)+debye) &
            < tolerance*debye, 'Independent physical -n e^2/T Debye normalization')
    end do
    xl_grid%npts_b = n
    xl_grid%npts_c = n-1
    rg_grid%npts_b = nc+1
    rg_grid%npts_c = nc

    do scenario = 1, 3
        offset = 1.0_dp
        width = 1.0_dp
        if (scenario == 2) then
            offset = 1.0e6_dp
            width = 1.0e-3_dp
        end if
        xl_grid%xb = offset+width*positions
        rg_grid%xb = offset+width*centers
        if (scenario == 3) then
            ! Preserve the inherited virtual endpoint hats inside a larger rg domain.
            rg_grid%xb(1) = 0.95_dp
            rg_grid%xb(nc+1) = 2.04_dp
        end if
        rg_grid%xc = rg_grid%xb(:nc)+0.5_dp*(rg_grid%xb(2:)-rg_grid%xb(:nc))
        length = xl_grid%xb(n)-xl_grid%xb(1)
        affine = (xl_grid%xb-xl_grid%xb(1))/length
        coefficients(:, 1) = cmplx(-debye, 0.0_dp, dp)
        coefficients(:, 2) = (0.4_dp, -0.3_dp)
        coefficients(:, 3) = (1.2_dp, 0.5_dp)
        coefficients(:, 4) = (-0.2_dp, -0.7_dp)
        call set_prefactors()
        if (trim(argument) == 'unsupported_basis') then
            spline_base = 3
            call active_matrix()
            print *, 'Unsupported electron basis was accepted'
            stop 0
        end if
        call active_matrix()
        call independent_gl2(coefficients, expected)
        call matrix_oracles('Homogeneous full physical Debye mass')
        if (scenario /= 3) then
            do channel = 1, 4
                do degree = 0, 1
                    total = (0.0_dp, 0.0_dp)
                    do lp = 1, n
                        scale = 1.0_dp
                        if (degree == 1) scale = affine(lp)
                        total = total+sum(matrix(:, lp, channel))*scale
                    end do
                    scale = length
                    if (degree == 1) scale = 0.5_dp*length
                    call expect(abs(total/coefficients(1, channel)-scale) &
                        < tolerance*length, 'Absolute constant/affine physical integral')
                end do
            end do
        end if

        ! Nonconstant complex cell coefficients cross unequal background/field knots.
        do j = 1, nc
            coefficients(j, 1) = -debye*cmplx(1.0_dp+0.13_dp*j, 0.07_dp*j, dp)
            coefficients(j, 2) = cmplx(0.2_dp*j, -0.3_dp+0.04_dp*j, dp)
            coefficients(j, 3) = cmplx(-0.4_dp+0.21_dp*j, 0.5_dp*j, dp)
            coefficients(j, 4) = cmplx(0.1_dp*j*j, -0.2_dp*j, dp)
        end do
        call set_prefactors()
        call active_matrix()
        call independent_gl2(coefficients, expected)
        call matrix_oracles('Cellwise complex prefactors and kink crossings')
        do channel = 1, 4
            do l = 1, n
                action = (0.0_dp, 0.0_dp)
                expected_action = (0.0_dp, 0.0_dp)
                do lp = 1, n
                    action = action+matrix(l, lp, channel)*affine(lp)
                    expected_action = expected_action+expected(l, lp, channel)*affine(lp)
                end do
                call expect(abs(action-expected_action) &
                    < tolerance*maxval(abs(expected(:, :, channel))), &
                    'Affine source action with cellwise coefficients')
            end do
        end do

        ! Manufactured conserved moment columns: x1 I10-x2 I00=-i and
        ! x1 I11-x2 I01=0. This tests common geometry, not an FP susceptibility fit.
        x1 = 0.7_dp
        x2 = 0.3_dp
        vt = 2.5_dp
        i00 = (0.4_dp, 0.2_dp)
        i01 = (-0.1_dp, 0.6_dp)
        i10 = (x2*i00-(0.0_dp, 1.0_dp))/x1
        i11 = x2*i01/x1
        source_coefficients = (0.0_dp, 0.0_dp)
        do j = 1, nc
            a = cmplx(0.0_dp, 0.13_dp*j, dp)
            coefficients(j, 1) = debye*(-1.0_dp+a*i00)
            coefficients(j, 2) = 0.2_dp*j*i01
            coefficients(j, 3) = debye*vt*a*i10
            coefficients(j, 4) = 0.2_dp*j*vt*i11
            source_coefficients(j, 1) = debye*(x2-(0.0_dp, 1.0_dp)*a)
        end do
        call set_prefactors()
        call active_matrix()
        call independent_gl2(source_coefficients, source)
        scale = maxval(abs(source(:, :, 1)))
        call expect(maxval(abs(x1/vt*matrix(:, :, 3)-x2*matrix(:, :, 1) &
            -source(:, :, 1))) < tolerance*scale, 'Geometric particle Ward source identity')
        scale = maxval(abs(matrix(:, :, 2)))
        call expect(maxval(abs(x1/vt*matrix(:, :, 4)-x2*matrix(:, :, 2))) &
            < tolerance*scale, 'Geometric magnetic-column particle Ward identity')
    end do
    print *, 'Electron full physical-cell independent oracles passed:', checks

contains

    subroutine expect(condition, label)
        logical, intent(in) :: condition
        character(*), intent(in) :: label
        if (.not. condition) then
            print *, 'FAIL: ', label
            error stop 'Electron physical-cell behavioral oracle failed'
        end if
        checks = checks+1
    end subroutine expect

    subroutine set_prefactors()
        pref_rho_phi_g1(1, :, 0) = 4.0_dp*pi*coefficients(:, 1)-pref_rho_phi_g0(1, :)
        pref_rho_B_g1(1, :, 0) = 4.0_dp*pi*coefficients(:, 2)
        pref_j_phi_g1(1, :, 0) = 4.0_dp*pi*coefficients(:, 3)
        pref_j_B_g1(1, :, 0) = 4.0_dp*pi*coefficients(:, 4)
    end subroutine set_prefactors

    subroutine active_matrix()
        do lp = 1, n
            do l = 1, n
                call FP_calc_kernel_element_electrons(l, lp, matrix(l, lp, 1), &
                    matrix(l, lp, 2), matrix(l, lp, 3), matrix(l, lp, 4), gauss)
            end do
        end do
    end subroutine active_matrix

    subroutine matrix_oracles(label)
        character(*), intent(in) :: label
        integer :: c, row, column
        real(dp) :: norm
        do c = 1, 4
            norm = maxval(abs(expected(:, :, c)))
            do column = 1, n
                do row = 1, n
                    call expect(abs(matrix(row, column, c)-expected(row, column, c)) &
                        < tolerance*norm, label)
                    if (abs(row-column) > 1) then
                        call expect(matrix(row, column, c) == (0.0_dp, 0.0_dp), &
                            'Disjoint P1 supports contribute exactly zero')
                    end if
                end do
            end do
            call expect(maxval(abs(matrix(:, :, c)-transpose(matrix(:, :, c)))) &
                < tolerance*norm, 'Symmetric bilinear complex kernel without conjugation')
        end do
    end subroutine matrix_oracles

    subroutine independent_gl2(cell_coefficients, result)
        complex(dp), intent(in) :: cell_coefficients(nc, 4)
        complex(dp), intent(out) :: result(n, n, 4)
        real(dp), parameter :: points(2) = &
            [0.5_dp-0.5_dp/sqrt(3.0_dp), 0.5_dp+0.5_dp/sqrt(3.0_dp)]
        real(dp) :: edges(0:n+1), left, right, h, coordinate, shapes(2)
        integer :: cell, element, q, row, column, r, c
        result = (0.0_dp, 0.0_dp)
        edges(1:n) = xl_grid%xb
        edges(0) = 2.0_dp*edges(1)-edges(2)
        edges(n+1) = 2.0_dp*edges(n)-edges(n-1)
        ! Assemble on background/field interval intersections, independently of
        ! production hat evaluation. Relative coordinates avoid large-offset loss.
        do cell = 1, nc
            do element = 0, n
                left = max(rg_grid%xb(cell), edges(element))
                right = min(rg_grid%xb(cell+1), edges(element+1))
                if (right <= left) cycle
                h = right-left
                do q = 1, 2
                    coordinate = ((left-edges(element))+points(q)*h) &
                        /(edges(element+1)-edges(element))
                    shapes(1) = 1.0_dp-coordinate
                    shapes(2) = coordinate
                    do c = 1, 2
                        column = element+c-1
                        if (column < 1) cycle
                        if (column > n) cycle
                        do r = 1, 2
                            row = element+r-1
                            if (row < 1) cycle
                            if (row > n) cycle
                            result(row, column, :) = result(row, column, :) &
                                +0.5_dp*h*shapes(r)*shapes(c)*cell_coefficients(cell, :)
                        end do
                    end do
                end do
            end do
        end do
    end subroutine independent_gl2
end program test_electron_physical_cells
