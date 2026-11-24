program test_host_geometry
    ! Host flux-surface geometry for the forced-periodicity solver:
    !   - the metric-weighted radial Poisson operator d/dr (M d/dr),
    !   - loading and evaluating a host geometry table,
    !   - an end-to-end periodic solve whose host table reproduces KIM's own
    !     screw-pinch geometry, and one with a non-unit metric.
    use KIM_kinds_m, only: dp

    implicit none

    logical :: all_passed

    all_passed = .true.
    call test_uniform_metric_operator(all_passed)
    call test_cosine_metric_operator(all_passed)
    call test_uniform_metric_solve(all_passed)
    call test_flr_arguments(all_passed)
    call test_table_evaluation(all_passed)
    call test_periodic_solve_with_host_geometry(all_passed)

    if (all_passed) then
        print *, 'All host geometry tests PASSED'
        stop 0
    else
        print *, 'Host geometry tests FAILED'
        stop 1
    end if

contains

    subroutine check(name, ok, passed)
        character(*), intent(in) :: name
        logical, intent(in) :: ok
        logical, intent(inout) :: passed

        if (ok) then
            print *, 'PASS: ', name
        else
            print *, 'FAIL: ', name
            passed = .false.
        end if
    end subroutine check

    subroutine window_nodes(n_nodes, r_lo, L, r_nodes)
        integer, intent(in) :: n_nodes
        real(dp), intent(in) :: r_lo, L
        real(dp), intent(out) :: r_nodes(n_nodes)
        integer :: j

        do j = 1, n_nodes
            r_nodes(j) = r_lo + L * real(j - 1, dp) / real(n_nodes, dp)
        end do
    end subroutine window_nodes

    subroutine test_uniform_metric_operator(passed)
        use constants_m, only: pi
        use periodic_solve_m, only: metric_radial_operator
        logical, intent(inout) :: passed

        integer, parameter :: M = 5, n_nodes = 32
        real(dp), parameter :: L = 2.0_dp, metric = 1.7_dp
        real(dp) :: r_nodes(n_nodes), metric_nodes(n_nodes), k_m, error
        complex(dp) :: D(2*M + 1, 2*M + 1)
        integer :: m_row, m_col

        call window_nodes(n_nodes, 54.0_dp, L, r_nodes)
        metric_nodes = metric
        D = metric_radial_operator(metric_nodes, r_nodes, L, M)
        error = 0.0_dp
        do m_col = -M, M
            do m_row = -M, M
                k_m = 2.0_dp * pi * real(m_row, dp) / L
                if (m_row == m_col) then
                    error = max(error, abs(D(m_row + M + 1, m_col + M + 1) + metric * k_m**2))
                else
                    error = max(error, abs(D(m_row + M + 1, m_col + M + 1)))
                end if
            end do
        end do
        call check('uniform metric gives the diagonal -M k^2', &
            error < 1.0e-12_dp * metric * (2.0_dp * pi * M / L)**2, passed)
    end subroutine test_uniform_metric_operator

    subroutine test_cosine_metric_operator(passed)
        use constants_m, only: pi, com_unit
        use periodic_solve_m, only: metric_radial_operator
        logical, intent(inout) :: passed

        ! M(r) = 1 + a cos(2 pi (r - r_lo)/L) has Mhat_0 = 1 and
        ! Mhat_{+-1} = (a/2) exp(-+ i 2 pi r_lo / L).
        integer, parameter :: M = 4, n_nodes = 40
        real(dp), parameter :: L = 3.0_dp, a = 0.3_dp, r_lo = 50.0_dp
        real(dp) :: r_nodes(n_nodes), metric_nodes(n_nodes), k_m, k_mp, error
        complex(dp) :: D(2*M + 1, 2*M + 1), expected, metric_hat
        integer :: m_row, m_col, shift

        call window_nodes(n_nodes, r_lo, L, r_nodes)
        metric_nodes = 1.0_dp + a * cos(2.0_dp * pi * (r_nodes - r_lo) / L)
        D = metric_radial_operator(metric_nodes, r_nodes, L, M)
        error = 0.0_dp
        do m_col = -M, M
            k_mp = 2.0_dp * pi * real(m_col, dp) / L
            do m_row = -M, M
                k_m = 2.0_dp * pi * real(m_row, dp) / L
                shift = m_row - m_col
                select case (shift)
                case (0)
                    metric_hat = (1.0_dp, 0.0_dp)
                case (1, -1)
                    metric_hat = 0.5_dp * a * exp(-com_unit * 2.0_dp * pi * shift * r_lo / L)
                case default
                    metric_hat = (0.0_dp, 0.0_dp)
                end select
                expected = -k_m * k_mp * metric_hat
                error = max(error, abs(D(m_row + M + 1, m_col + M + 1) - expected))
            end do
        end do
        call check('cosine metric couples neighbouring Fourier modes', &
            error < 1.0e-11_dp * (2.0_dp * pi * M / L)**2, passed)
    end subroutine test_cosine_metric_operator

    subroutine test_uniform_metric_solve(passed)
        use periodic_solve_m, only: solve_periodic
        logical, intent(inout) :: passed

        ! (-c k^2 + 4 pi K) Phi = -4 pi KB Br  equals
        ! (-k^2 + 4 pi K/c) Phi = -4 pi (KB/c) Br.
        integer, parameter :: M = 3, n_nodes = 24, dim = 2*M + 1
        real(dp), parameter :: L = 1.5_dp, metric = 2.3_dp
        complex(dp), parameter :: Br = (0.7_dp, -0.2_dp)
        real(dp) :: r_nodes(n_nodes), metric_nodes(n_nodes)
        complex(dp) :: Kphi(dim, dim), KB(dim, dim)
        complex(dp), allocatable :: Phi_metric(:), Phi_scaled(:)
        integer :: i, k, info_metric, info_scaled

        do k = 1, dim
            do i = 1, dim
                Kphi(i, k) = cmplx(0.01_dp * real(i - k, dp), 0.02_dp * real(i + k, dp), dp)
                KB(i, k) = cmplx(0.3_dp / real(i + k, dp), -0.1_dp * real(i * k, dp) / dim, dp)
            end do
            Kphi(k, k) = Kphi(k, k) - (1.0_dp, 0.0_dp)
        end do
        call window_nodes(n_nodes, 40.0_dp, L, r_nodes)
        metric_nodes = metric

        call solve_periodic(Kphi, KB, L, M, Br, Phi_metric, info_metric, &
            radial_metric=metric_nodes, r_nodes=r_nodes)
        call solve_periodic(Kphi / metric, KB / metric, L, M, Br, Phi_scaled, info_scaled)
        call check('uniform-metric solve succeeds', info_metric == 0 .and. info_scaled == 0, passed)
        call check('uniform metric equals rescaled kernels', &
            maxval(abs(Phi_metric - Phi_scaled)) <= 1.0e-12_dp * maxval(abs(Phi_scaled)), passed)
    end subroutine test_uniform_metric_solve

    subroutine test_flr_arguments(passed)
        use species_m, only: plasma_t, species_t
        use flr2_fourier_kernel_m, only: flr_arg_pair_sp, boltzmann_gyroaverage, &
            kern_include_ks2
        logical, intent(inout) :: passed

        ! A metric M at node j acts like the radial wavenumbers sqrt(M) k_r.
        type(plasma_t) :: with_metric, without_metric
        real(dp), parameter :: metric = 1.44_dp, kr = 3.0_dp, krp = -1.0_dp
        real(dp) :: bplus_m, bcross_m, bplus_0, bcross_0
        logical :: saved_ks2

        saved_ks2 = kern_include_ks2
        allocate(without_metric%spec(0:0))
        allocate(without_metric%spec(0)%rho_L(1), without_metric%ks(1))
        without_metric%spec(0)%rho_L = 0.2_dp
        without_metric%ks = 0.5_dp
        with_metric = without_metric
        allocate(with_metric%radial_metric(1))
        with_metric%radial_metric = metric

        kern_include_ks2 = .true.
        call flr_arg_pair_sp(with_metric, 0, kr, krp, 1, bplus_m, bcross_m)
        call flr_arg_pair_sp(without_metric, 0, sqrt(metric) * kr, sqrt(metric) * krp, 1, &
            bplus_0, bcross_0)
        call check('metric scales k_r in the Larmor arguments (with k_s)', &
            abs(bplus_m - bplus_0) < 1.0e-14_dp .and. abs(bcross_m - bcross_0) < 1.0e-14_dp, &
            passed)

        kern_include_ks2 = .false.
        call flr_arg_pair_sp(with_metric, 0, kr, krp, 1, bplus_m, bcross_m)
        call flr_arg_pair_sp(without_metric, 0, sqrt(metric) * kr, sqrt(metric) * krp, 1, &
            bplus_0, bcross_0)
        call check('metric scales k_r in the Larmor arguments (without k_s)', &
            abs(bplus_m - bplus_0) < 1.0e-14_dp .and. abs(bcross_m - bcross_0) < 1.0e-14_dp, &
            passed)
        kern_include_ks2 = saved_ks2

        call check('metric scales k_r in the Boltzmann gyroaverage', &
            abs(boltzmann_gyroaverage(with_metric, 0, kr, krp, 1) - &
                boltzmann_gyroaverage(without_metric, 0, sqrt(metric) * kr, &
                    sqrt(metric) * krp, 1)) < 1.0e-14_dp, passed)
    end subroutine test_flr_arguments

    subroutine test_table_evaluation(passed)
        use host_geometry_m, only: load_host_geometry, evaluate_host_geometry, &
            reset_host_geometry, host_geometry_active, host_geometry_has_metric
        logical, intent(inout) :: passed

        ! Cubic Lagrange interpolation reproduces linear tables exactly.
        character(*), parameter :: path = './host_geometry_linear_test.dat'
        real(dp) :: B0, kp, ks, metric, r
        integer :: unit, i

        open(newunit=unit, file=path, status='replace', action='write')
        write(unit, '(a)') '# r B0 kp_m kp_n ks_m ks_n radial_metric'
        do i = 0, 10
            r = 10.0_dp + 5.0_dp * real(i, dp)
            write(unit, '(7es24.16)') r, 2.0e4_dp - 10.0_dp * r, 1.0e-3_dp * r, 2.0e-3_dp, &
                1.0_dp / 30.0_dp, -1.0e-4_dp * r, 1.0_dp + 0.01_dp * r
        end do
        close(unit)

        call load_host_geometry(path)
        call check('host geometry is active after loading', host_geometry_active(), passed)
        call check('non-unit metric is detected', host_geometry_has_metric(), passed)
        r = 37.3_dp
        call evaluate_host_geometry(r, -6, 2, B0, kp, ks, metric)
        call check('B0 interpolates linearly', abs(B0 - (2.0e4_dp - 10.0_dp * r)) < 1.0e-9_dp, passed)
        call check('kp = m kp_m + n kp_n', &
            abs(kp - (-6.0_dp * 1.0e-3_dp * r + 2.0_dp * 2.0e-3_dp)) < 1.0e-13_dp, passed)
        call check('ks = m ks_m + n ks_n', &
            abs(ks - (-6.0_dp / 30.0_dp - 2.0_dp * 1.0e-4_dp * r)) < 1.0e-13_dp, passed)
        call check('metric interpolates linearly', abs(metric - (1.0_dp + 0.01_dp * r)) < 1.0e-13_dp, &
            passed)
        call evaluate_host_geometry(70.0_dp, -6, 2, B0, kp, ks, metric)
        call check('geometry is held constant beyond the table', &
            abs(B0 - (2.0e4_dp - 10.0_dp * 60.0_dp)) < 1.0e-9_dp .and. &
            abs(metric - 1.6_dp) < 1.0e-13_dp, passed)
        call reset_host_geometry()
        call check('reset deactivates the host geometry', .not. host_geometry_active(), passed)
    end subroutine test_table_evaluation

    subroutine write_config(path, geometry_path)
        character(*), intent(in) :: path, geometry_path
        character(len=512) :: line
        integer :: unit_in, unit_out, ios

        open(newunit=unit_in, file='KIM_config_em_small.nml', status='old', action='read')
        open(newunit=unit_out, file=path, status='replace', action='write')
        do
            read(unit_in, '(a)', iostat=ios) line
            if (ios /= 0) exit
            write(unit_out, '(a)') trim(line)
            if (index(line, "q_file = 'q.dat'") > 0 .and. len_trim(geometry_path) > 0) then
                write(unit_out, '(a)') "    geometry_file = '" // trim(geometry_path) // "'"
            end if
        end do
        close(unit_in)
        write(unit_out, '(a)') '&KIM_PERIODIC'
        write(unit_out, '(a)') '    periodic_calculate_radial_current = .false.'
        write(unit_out, '(a)') '    periodic_calculate_ion_tensor = .false.'
        write(unit_out, '(a)') '/'
        close(unit_out)
    end subroutine write_config

    subroutine write_geometry(path, background, m, metric)
        use kim_solver_m, only: kim_results_t
        character(*), intent(in) :: path
        type(kim_results_t), intent(in) :: background
        integer, intent(in) :: m
        real(dp), intent(in) :: metric
        integer :: unit, i

        ! Mode-specific table: kp = m kp_m and ks = m ks_m reproduce KIM's own
        ! screw-pinch values for this (m, n).
        open(newunit=unit, file=path, status='replace', action='write')
        do i = 1, size(background%r_plasma)
            write(unit, '(7es26.17e3)') background%r_plasma(i), background%B0(i), &
                background%kp(i) / real(m, dp), 0.0_dp, &
                background%ks(i) / real(m, dp), 0.0_dp, metric
        end do
        close(unit)
    end subroutine write_geometry

    subroutine test_periodic_solve_with_host_geometry(passed)
        use kim_solver_m, only: kim_solver_t, kim_results_t, kim_profiles_t, KIM_OK
        use host_geometry_m, only: reset_host_geometry
        logical, intent(inout) :: passed

        ! The profile grid equals KIM's global grid (64 equidistant points on
        ! [r_min, r_plas] = [3, 67]), so the host table built from the global
        ! background is sampled exactly at its nodes.
        integer, parameter :: npts = 64, m_mode = -6, n_mode = 2
        type(kim_solver_t) :: reference, same_geometry, stretched
        type(kim_results_t) :: res_reference, res_same, res_stretched, background
        type(kim_profiles_t) :: prof
        character(len=1024) :: cwd
        character(len=:), allocatable :: geometry_same, geometry_stretched
        real(dp) :: frac
        integer :: i, ierr

        allocate(prof%r(npts), prof%n(npts), prof%Te(npts), prof%Ti(npts), &
                 prof%q(npts), prof%Er(npts))
        do i = 1, npts
            prof%r(i) = 3.0_dp + 64.0_dp * real(i - 1, dp) / real(npts - 1, dp)
            frac = (prof%r(i) - 3.0_dp) / 64.0_dp
            prof%n(i) = 2.0e13_dp * (1.1_dp - prof%r(i) / 100.0_dp)
            prof%Te(i) = 1.0e3_dp * (1.2_dp - prof%r(i) / 100.0_dp)
            prof%Ti(i) = prof%Te(i)
            prof%q(i) = 1.5_dp + 2.5_dp * frac
            prof%Er(i) = -0.5_dp
        end do

        call getcwd(cwd)
        geometry_same = trim(cwd) // '/host_geometry_same_test.dat'
        geometry_stretched = trim(cwd) // '/host_geometry_stretched_test.dat'

        call write_config('KIM_config_host_geometry_none.nml', '')
        call reference%init('KIM_config_host_geometry_none.nml', &
            run_type='electrostatic_periodic', profiles=prof, stat=ierr)
        if (ierr == KIM_OK) call reference%solve(m=m_mode, n=n_mode, stat=ierr)
        call check('reference periodic solve succeeds', ierr == KIM_OK, passed)
        if (ierr /= KIM_OK) return
        res_reference = reference%results()
        background = reference%background()
        call reference%finalize()
        call check('reference background has the global geometry', &
            allocated(background%r_plasma) .and. allocated(background%B0) .and. &
            allocated(background%kp) .and. allocated(background%ks), passed)
        call check('global grid equals the profile grid', &
            size(background%r_plasma) == npts, passed)
        if (size(background%r_plasma) /= npts) return
        call check('global grid nodes equal the profile nodes', &
            maxval(abs(background%r_plasma - prof%r)) < 1.0e-12_dp, passed)

        call write_geometry(geometry_same, background, m_mode, 1.0_dp)
        call write_config('KIM_config_host_geometry_same.nml', geometry_same)
        call same_geometry%init('KIM_config_host_geometry_same.nml', &
            run_type='electrostatic_periodic', profiles=prof, stat=ierr)
        if (ierr == KIM_OK) call same_geometry%solve(m=m_mode, n=n_mode, stat=ierr)
        call check('periodic solve with host geometry succeeds', ierr == KIM_OK, passed)
        if (ierr /= KIM_OK) return
        res_same = same_geometry%results()
        call same_geometry%finalize()
        call check('host copy of the screw pinch reproduces Phi', &
            maxval(abs(res_same%Phi - res_reference%Phi)) <= &
            1.0e-9_dp * maxval(abs(res_reference%Phi)), passed)
        call check('host copy of the screw pinch reproduces jpar', &
            maxval(abs(res_same%jpar - res_reference%jpar)) <= &
            1.0e-9_dp * maxval(abs(res_reference%jpar)), passed)

        call write_geometry(geometry_stretched, background, m_mode, 1.5_dp)
        call write_config('KIM_config_host_geometry_stretched.nml', geometry_stretched)
        call stretched%init('KIM_config_host_geometry_stretched.nml', &
            run_type='electrostatic_periodic', profiles=prof, stat=ierr)
        if (ierr == KIM_OK) call stretched%solve(m=m_mode, n=n_mode, stat=ierr)
        call check('periodic solve with a radial metric succeeds', ierr == KIM_OK, passed)
        if (ierr /= KIM_OK) return
        res_stretched = stretched%results()
        call stretched%finalize()
        call reset_host_geometry()
        call check('radial metric changes the potential', &
            maxval(abs(res_stretched%Phi - res_reference%Phi)) > &
            1.0e-3_dp * maxval(abs(res_reference%Phi)), passed)
        call check('radial metric keeps the potential finite', &
            all(abs(res_stretched%Phi) < huge(1.0_dp)), passed)
    end subroutine test_periodic_solve_with_host_geometry

end program test_host_geometry
