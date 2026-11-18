program test_equilibrium_source_stability
    !! Stability of the unperturbed equilibrium against a temperature offset.
    !!
    !! The equilibrium sources of particles and energy keep the initial profiles
    !! stationary [Markl2023 (22), (24), (25)]. Particles supplied by the particle
    !! source arrive at the source temperature, so a plasma that is hotter than the
    !! equilibrium is cooled by dilution. A radially uniform temperature offset must
    !! therefore decay on the particle confinement time, without any RMP transport.
    !!
    !! If the sources are applied as a fixed rate of temperature change instead,
    !! the dilution term is lost and the offset is damped only through the outer
    !! boundary. With the temperature-dependent neoclassical ion heat conductivity
    !! and the one-sided convective difference it then grows.
    !!
    !! Asserted: after 4 s (about 1.2 radial diffusion times a^2/D^a of the fixture)
    !! the offset of T_e and T_i has decayed to less than 10 % of its initial value
    !! everywhere, and the density stays at its equilibrium profile.
    use QLBalance_kinds, only: dp
    use grid_mod, only: rb, rc, npoib, npoi, rmin, rmax, npoimin, iboutype, gg_factor, &
                        gg_width, set_boundary_condition, mwind, neqset, nbaleqs, &
                        dery_equisource, source_term, rb_cut_in, rb_cut_out, re_cut_out, &
                        dqle11, dqle12, dqle21, dqle22, dqli11, dqli12, dqli21, dqli22
    use plasma_parameters, only: params, init_background_profiles
    use baseparam_mod, only: am, rtor, Z_i, btor, rsepar, ev
    use paramscan_mod, only: viscosity_factor
    use resonances_mod, only: prop, numres
    use control_mod, only: ihdf5IO
    use wave_code_data, only: r, n, Te, Ti, q, Vth, Vz
    use recstep_mod, only: timstep_arr
    use h5mod, only: path2inp, h5_id, group_id_1, h5_create, h5_define_group, &
                     h5_close_group, h5_close, h5_deinit, h5_add

    implicit none

    real(dp), parameter :: time_step = 5.0e-2_dp
    integer, parameter :: num_steps = 80
    real(dp), parameter :: temperature_offset = 1.0_dp  ! eV
    real(dp), parameter :: remaining_offset_tol = 0.1_dp
    real(dp), parameter :: density_tol = 1.0e-6_dp
    real(dp), allocatable :: input_r(:), equilibrium(:, :)
    real(dp) :: remaining(3:4), density_change
    integer :: num_failed, istep, ieq
    external :: gengrid, allocate_wave_code_data, calc_geometric_parameter_profiles
    external :: det_balance_eqs_source_terms, evolvestep

    num_failed = 0
    call setup_plasma()
    allocate (timstep_arr(neqset))
    timstep_arr = time_step

    call det_balance_eqs_source_terms()
    equilibrium = params(:, 1:npoi)

    params(3:4, 1:npoi) = equilibrium(3:4, :) + temperature_offset * ev
    do istep = 1, num_steps
        call evolvestep(time_step, 1.0e-8_dp)
    end do

    do ieq = 3, 4
        remaining(ieq) = maxval(abs(params(ieq, 1:npoi) - equilibrium(ieq, :))) / &
                         (temperature_offset * ev)
    end do
    density_change = maxval(abs(params(1, 1:npoi) - equilibrium(1, :)) / equilibrium(1, :))

    print '(A,F5.2,A,2ES10.2)', '  remaining fraction of the offset after ', &
        num_steps * time_step, ' s (Te, Ti):', remaining
    print '(A,ES10.2)', '  relative change of the density:', density_change
    call require(remaining(3) < remaining_offset_tol, 'T_e offset does not decay')
    call require(remaining(4) < remaining_offset_tol, 'T_i offset does not decay')
    call require(density_change < density_tol, 'density left its equilibrium profile')

    if (num_failed > 0) then
        print *, 'TESTS FAILED: ', num_failed
        error stop 1
    end if
    print *, 'PASS: the unperturbed equilibrium is stable against a temperature offset'

contains

    subroutine setup_plasma()
        real(dp), allocatable :: diffusion(:), pedestal(:)

        rmin = 2.0_dp
        rmax = 60.0_dp
        npoimin = 1500
        iboutype = 1
        gg_factor = 0.0_dp
        gg_width = 58.0_dp
        prop = .false.
        numres = 0
        call gengrid()
        call set_boundary_condition()

        allocate (input_r(npoib))
        input_r = rb
        call allocate_wave_code_data(npoib, input_r)
        ! H-mode like profiles with a pedestal at r = 54 cm.
        pedestal = 0.5_dp * (1.0_dp - tanh((r - 54.0_dp) / 2.0_dp))
        n = 1.0e12_dp + 2.0e13_dp * pedestal * (1.2_dp - r / 100.0_dp)
        Te = 50.0_dp + 1.0e3_dp * pedestal * (1.2_dp - r / 100.0_dp)
        Ti = 80.0_dp + 0.8e3_dp * pedestal * (1.2_dp - r / 100.0_dp)
        q = 1.5_dp + 3.0_dp * (r - 2.0_dp) / 58.0_dp
        Vth = 3.0e5_dp * r / 60.0_dp
        Vz = -8.0e4_dp * 165.0_dp * (1.0_dp - r / 120.0_dp)
        rtor = 165.0_dp
        am = 2.0_dp
        Z_i = 1.0_dp
        btor = -18000.0_dp
        rsepar = 60.0_dp
        mwind = 1
        dery_equisource = 0.0_dp
        source_term = 0.0_dp
        call init_background_profiles()
        call calc_geometric_parameter_profiles()

        viscosity_factor = 1.0_dp

        ! No RMP: the quasilinear transport coefficients vanish.
        dqle11 = 0.0_dp
        dqle12 = 0.0_dp
        dqle21 = 0.0_dp
        dqle22 = 0.0_dp
        dqli11 = 0.0_dp
        dqli12 = 0.0_dp
        dqli21 = 0.0_dp
        dqli22 = 0.0_dp

        ihdf5IO = 1
        rb_cut_in = 0.0_dp
        rb_cut_out = maxval(rb) + 1.0_dp
        re_cut_out = rb_cut_out + 1.0_dp
        diffusion = spread(1.0e3_dp, 1, npoib)
        path2inp = 'equilibrium_source_stability_input.h5'
        call h5_create(trim(path2inp), h5_id)
        call h5_define_group(h5_id, 'da_estimation', group_id_1)
        call h5_close_group(group_id_1)
        call h5_add(h5_id, '/da_estimation/r', rb, [1], [npoib])
        call h5_add(h5_id, '/da_estimation/Da', diffusion, [1], [npoib])
        call h5_close(h5_id)
        call h5_deinit()
    end subroutine setup_plasma

    subroutine require(condition, message)
        logical, intent(in) :: condition
        character(*), intent(in) :: message
        if (.not. condition) then
            print *, 'FAIL: ', message
            num_failed = num_failed + 1
        end if
    end subroutine require

end program test_equilibrium_source_stability
