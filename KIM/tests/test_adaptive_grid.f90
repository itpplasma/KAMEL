program test_adaptive_grid
    use KIM_kinds_m, only: dp
    use grid_m, only: grid_type, width_res, ampl_res, hrmax_scaling, &
        xl_grid, rg_grid, calc_mass_matrix, rg_space_dim, l_space_dim, &
        grid_spacing_rg, grid_spacing_xl, r_min, r_plas
    use config_m, only: output_path, hdf5_output
    use kim_resonances_m, only: r_res, prop
    use, intrinsic :: ieee_arithmetic, only: ieee_is_finite, ieee_value, ieee_quiet_nan
    implicit none
    type(grid_type) :: grid
    integer :: nominal, checks, level
    real(dp) :: left, right, gap
    character(32) :: argument

    checks = 0
    prop = .false.
    r_res = 1.53_dp
    width_res = 0.06_dp
    ampl_res = 20.0_dp
    hrmax_scaling = 1.0_dp
    nominal = 32
    left = 1.0_dp
    right = 2.0_dp
    call get_command_argument(1, argument)
    if (len_trim(argument) > 0) then
        select case (trim(argument))
        case ('nan_min')
            left = ieee_value(0.0_dp, ieee_quiet_nan)
        case ('nan_max')
            right = ieee_value(0.0_dp, ieee_quiet_nan)
        case ('unordered')
            right = left
        case ('one_node')
            nominal = 1
        case ('nan_width')
            width_res = ieee_value(0.0_dp, ieee_quiet_nan)
        case ('zero_width')
            width_res = 0.0_dp
        case ('nan_amplitude')
            ampl_res = ieee_value(0.0_dp, ieee_quiet_nan)
        case ('nonpositive_factor')
            ampl_res = -1.0_dp
        case ('nan_scale')
            hrmax_scaling = ieee_value(0.0_dp, ieee_quiet_nan)
        case ('zero_scale')
            hrmax_scaling = 0.0_dp
        case ('stalled')
            right = nearest(left, 1.0_dp)
        case default
            error stop 'Unknown adaptive-grid rejection case'
        end select
        call grid%grid_init(nominal, left, right, 'invalid')
        call grid%grid_generate()
        print *, 'Invalid adaptive-grid controls were accepted'
        stop 0
    end if

    call grid%grid_init(nominal, left, right, 'manufactured')
    call grid%grid_generate()
    call geometry(left, right, nominal)
    call weak_oracles()
    call gaussian_density(left, right)

    ! Independent field/background budgets share the physical domain and monitor.
    r_res = 59.22635072676093_dp
    width_res = 0.15_dp
    ampl_res = 32.0_dp
    do level = 1, 2
        nominal = 256*level
        call grid%grid_init(nominal, 3.0_dp, 61.45547565601116_dp, 'physical')
        call grid%grid_generate()
        call geometry(3.0_dp, 61.45547565601116_dp, nominal)
        gap = minval(grid%xb(2:)-grid%xb(:grid%npts_c))
        call expect(gap < 0.007_dp, 'Targeted physical layer refinement')
        call gaussian_density(3.0_dp, 61.45547565601116_dp)
        print *, 'nominal, actual nodes, minimum h:', nominal, grid%npts_b, gap
    end do

    ! Uniform spacing nearly divides the interval: avoid a duplicate terminal cell.
    nominal = 16
    ampl_res = 0.0_dp
    hrmax_scaling = nearest(1.0_dp, -1.0_dp)
    call grid%grid_init(nominal, 1.0_dp, 2.0_dp, 'near-terminal')
    call grid%grid_generate()
    call geometry(1.0_dp, 2.0_dp, nominal)
    call expect(grid%npts_b == nominal+1, 'No numerically duplicate terminal node')
    call expect(minval(grid%xb(2:)-grid%xb(:grid%npts_c)) > 0.06_dp, &
        'Terminal cell keeps a physical width')

    call flat_monitor_oracles()

    ! A valid coarsening amplitude remains supported if its monitor is positive.
    ampl_res = -0.5_dp
    r_res = 1.5_dp
    width_res = 0.2_dp
    call grid%grid_init(nominal, 1.0_dp, 2.0_dp, 'coarsening')
    call grid%grid_generate()
    call geometry(1.0_dp, 2.0_dp, nominal)
    call production_grid_oracles()
    print *, 'Adaptive grid independent geometry/weak oracles passed:', checks

contains

    subroutine production_grid_oracles()
        integer :: fixture

        ! Exercise the public caller with unequal field/background budgets and
        ! both adaptive routing names, rather than constructing grid types alone.
        rg_space_dim = 128
        l_space_dim = 256
        r_min = 3.0_dp
        r_plas = 61.45547565601116_dp
        grid_spacing_rg = 'adaptive'
        grid_spacing_xl = 'non-equidistant'
        output_path = './adaptive-grid-contract-output/'
        hdf5_output = .false.
        call execute_command_line('mkdir -p adaptive-grid-contract-output/grid')
        r_res = 59.22635072676093_dp
        width_res = 0.15_dp
        hrmax_scaling = 1.0_dp
        do fixture = 1, 2
            ampl_res = 0.0_dp
            if (fixture == 2) ampl_res = 32.0_dp
            call generate_grids()
            grid = rg_grid
            call geometry(r_min, r_plas, rg_space_dim)
            if (fixture == 1) then
                call expect(grid%npts_b == rg_space_dim+1, &
                    'Production background keeps its flat-monitor cell budget')
            else
                call gaussian_density(r_min, r_plas)
            end if
            grid = xl_grid
            call geometry(r_min, r_plas, l_space_dim)
            if (fixture == 1) then
                call expect(grid%npts_b == l_space_dim+1, &
                    'Production field keeps its distinct flat-monitor cell budget')
            else
                call gaussian_density(r_min, r_plas)
            end if
        end do
    end subroutine production_grid_oracles

    subroutine expect(condition, label)
        logical, intent(in) :: condition
        character(*), intent(in) :: label
        if (.not. condition) then
            print *, 'FAIL: ', label
            error stop 'Adaptive grid behavioral oracle failed'
        end if
        checks = checks+1
    end subroutine expect

    subroutine geometry(a, b, expected_nominal)
        real(dp), intent(in) :: a, b
        integer, intent(in) :: expected_nominal
        call expect(grid%xb(1) == a, 'Exact requested lower physical endpoint')
        call expect(grid%xb(grid%npts_b) == b, 'Exact requested upper physical endpoint')
        call expect(grid%npts == expected_nominal, 'Nominal budget metadata')
        call expect(grid%npts_c == grid%npts_b-1, 'Staggered cell count')
        call expect(size(grid%xb) == grid%npts_b, 'Physical node count')
        call expect(size(grid%xc) == grid%npts_c, 'Physical cell-center count')
        call expect(all(ieee_is_finite(grid%xb)), 'Finite physical nodes')
        call expect(all(grid%xb(2:) > grid%xb(:grid%npts_c)), 'Positive local cell widths')
        call expect(all(grid%xc > grid%xb(:grid%npts_c)), 'Centers above left cell boundary')
        call expect(all(grid%xc < grid%xb(2:)), 'Centers below right cell boundary')
        call expect(maxval(abs(grid%xc-grid%xb(:grid%npts_c) &
            -0.5_dp*(grid%xb(2:)-grid%xb(:grid%npts_c)))) &
            < 8.0_dp*epsilon(b)*max(1.0_dp, abs(a), abs(b)), 'Accepted staggered midpoints')
    end subroutine geometry

    subroutine flat_monitor_oracles()
        integer, parameter :: budgets(4) = [16, 128, 512, 1024]
        real(dp), parameter :: bounds(2, 3) = reshape([ &
            3.0_dp, 61.45547565601116_dp, &
            -61.45547565601116_dp, 0.0_dp, &
            -61.45547565601116_dp, -3.0_dp], [2, 3])
        integer :: b, j
        real(dp) :: a, z, expected_h, expected_remainder

        ampl_res = 0.0_dp
        hrmax_scaling = 1.0_dp
        do b = 1, size(bounds, 2)
            a = bounds(1, b)
            z = bounds(2, b)
            do j = 1, size(budgets)
                nominal = budgets(j)
                call grid%grid_init(nominal, a, z, 'nondyadic-flat')
                call grid%grid_generate()
                call geometry(a, z, nominal)
                call expect(grid%npts_b == nominal + 1, &
                    'Flat monitor closes exactly the nominal number of cells')
                expected_h = (z - a)/real(nominal, dp)
                call expect(minval(grid%xb(2:) - grid%xb(:grid%npts_c)) &
                    > 0.99_dp*expected_h, 'No roundoff-sized flat terminal cell')
                call weak_oracles()
            end do
        end do

        ! A small but resolved remainder is physical and must survive endpoint closure.
        nominal = 16
        a = bounds(1, 1)
        z = bounds(2, 1)
        hrmax_scaling = real(nominal, dp)/(real(nominal, dp) + 1.0e-8_dp)
        call grid%grid_init(nominal, a, z, 'resolved-small-terminal')
        call grid%grid_generate()
        call geometry(a, z, nominal)
        call expect(grid%npts_b == nominal + 2, 'Retain a resolved short final cell')
        expected_remainder = (z - a)*1.0e-8_dp/(real(nominal, dp) + 1.0e-8_dp)
        call expect(abs(grid%xb(grid%npts_b) - grid%xb(grid%npts_c) &
            - expected_remainder) < 1.0e-5_dp*expected_remainder, &
            'Resolved short final cell has its analytic constant-monitor width')
        hrmax_scaling = 1.0_dp
    end subroutine flat_monitor_oracles

    subroutine gaussian_density(a, b)
        real(dp), intent(in) :: a, b
        real(dp) :: pi, expected_cells, smallest, center, h
        integer :: i
        pi = acos(-1.0_dp)
        expected_cells = ((b-a)+ampl_res*width_res*sqrt(pi)*0.5_dp &
            *(erf((b-r_res)/width_res)-erf((a-r_res)/width_res)))/grid%hrmax
        call expect(abs(real(grid%npts_c, dp)-expected_cells) < 3.0_dp, &
            'Analytic integrated Gaussian mesh density')
        smallest = minval(grid%xb(2:)-grid%xb(:grid%npts_c))
        call expect(smallest < 0.06_dp*grid%hrmax, 'Gaussian density concentrates nodes')
        smallest = huge(1.0_dp)
        center = 0.0_dp
        do i = 1, grid%npts_c
            if (abs(grid%xc(i)-r_res) > width_res) cycle
            h = grid%xb(i+1)-grid%xb(i)
            if (h >= smallest) cycle
            smallest = h
            center = grid%xc(i)
        end do
        call expect(abs(center-r_res) < 2.0_dp*smallest, &
            'Refinement remains centered on the physical rational surface')
    end subroutine gaussian_density

    subroutine weak_oracles()
        real(dp), parameter :: points(2) = &
            [0.5_dp-0.5_dp/sqrt(3.0_dp), 0.5_dp+0.5_dp/sqrt(3.0_dp)]
        real(dp), allocatable :: mass(:, :), qmass(:, :), radius(:)
        complex(dp), allocatable :: exact(:), load(:), active_load(:)
        complex(dp) :: value, integral
        real(dp) :: h, shape(2), at, scale
        integer :: i, j, k, q, n, degree

        radius = grid%xb
        n = size(radius)
        allocate(mass(n, n), qmass(n, n), load(n), active_load(n))
        xl_grid%xb = radius
        xl_grid%npts_b = n
        call calc_mass_matrix(mass)
        qmass = 0.0_dp
        do i = 1, n-1
            h = radius(i+1)-radius(i)
            do q = 1, 2
                shape(1) = 1.0_dp-points(q)
                shape(2) = points(q)
                do k = 1, 2
                    do j = 1, 2
                        qmass(i+j-1, i+k-1) = qmass(i+j-1, i+k-1) &
                            +0.5_dp*h*shape(j)*shape(k)
                    end do
                end do
            end do
        end do
        call expect(maxval(abs(mass(:n-1, :n-1)-qmass(:n-1, :n-1))) &
            < 100.0_dp*epsilon(h), 'Native P1 mass on generated nonuniform mesh')
        call expect(all(mass(n, :n-1) == 0.0_dp), 'Right Dirichlet mass row')
        call expect(all(mass(:n-1, n) == 0.0_dp), 'Right Dirichlet mass column')
        call expect(mass(n, n) == 1.0_dp, 'Right Dirichlet mass unit diagonal')
        do degree = 0, 1
            if (degree == 0) then
                exact = cmplx(2.0_dp+0.0_dp*radius, 1.0_dp+0.0_dp*radius, dp)
                integral = cmplx(2.0_dp, 1.0_dp, dp)*(radius(n)-radius(1))
            else
                exact = cmplx(2.0_dp+3.0_dp*radius, 1.0_dp-2.0_dp*radius, dp)
                integral = cmplx(2.0_dp, 1.0_dp, dp)*(radius(n)-radius(1)) &
                    +cmplx(3.0_dp, -2.0_dp, dp)*(radius(n)**2-radius(1)**2)/2.0_dp
            end if
            load = (0.0_dp, 0.0_dp)
            do i = 1, n-1
                h = radius(i+1)-radius(i)
                do j = 1, 2
                    shape(1) = 1.0_dp-points(j)
                    shape(2) = points(j)
                    at = radius(i)+h*points(j)
                    value = cmplx(2.0_dp, 1.0_dp, dp)
                    if (degree == 1) value = value+cmplx(3.0_dp, -2.0_dp, dp)*at
                    load(i:i+1) = load(i:i+1)+0.5_dp*h*shape*value
                end do
            end do
            scale = maxval(abs(exact))
            call expect(abs(sum(load)-integral) &
                < 400.0_dp*epsilon(scale)*max(scale, abs(integral)), &
                'Analytic constant/linear physical-domain integral')
            active_load = load-qmass(:, n)*exact(n)
            exact(n) = (0.0_dp, 0.0_dp)
            do i = 1, n-1
                value = (0.0_dp, 0.0_dp)
                do j = 1, n
                    value = value+mass(i, j)*exact(j)
                end do
                call expect(abs(value-active_load(i)) &
                    < 400.0_dp*epsilon(scale)*scale, &
                    'Native constant/linear P1 loads with explicit Dirichlet data')
            end do
        end do
    end subroutine weak_oracles
end program test_adaptive_grid
