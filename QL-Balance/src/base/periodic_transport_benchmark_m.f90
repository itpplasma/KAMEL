module periodic_transport_benchmark_m
    !! Optional comparison of the actual QL drift path and embedded KIM tensor.
    !! Values are per mode, before antenna weighting, summation and smoothing.
    use QLBalance_kinds, only: dp
    use control_mod, only: kim_transport_benchmark, kim_ion_transport_model, &
        ion_transport_model_id, ION_TRANSPORT_FLR, ION_TRANSPORT_DRIFT_KINETIC
    use kim_wave_code_adapter_m, only: kim_D_ion_modes
    use kim_qldiff_m, only: calc_dqli_benchmark_residuals
    implicit none
    private
    public :: select_periodic_ion_transport, reset_transport_benchmark, write_transport_benchmark
    real(dp), allocatable :: drift(:, :, :, :), flr(:, :, :, :)
    logical, allocatable :: recorded(:)

contains

    subroutine reset_transport_benchmark()
        if (allocated(recorded)) deallocate(recorded, drift, flr)
    end subroutine reset_transport_benchmark

    subroutine select_periodic_ion_transport(mode, vt, nu, selected)
        use wave_code_data, only: dim_mn
        use QLBalance_diag, only: i_mn_loop
        integer, intent(in) :: mode
        real(dp), intent(in) :: vt(:), nu(:)
        real(dp), intent(out) :: selected(2, 2, size(vt))
        real(dp) :: reference(2, 2, size(vt))
        integer :: model

        if (mode < 1 .or. mode > dim_mn) error stop 'periodic transport: invalid mode index'
        model = ion_transport_model_id(kim_ion_transport_model)
        if (mode /= i_mn_loop) error stop 'periodic transport: inconsistent mode index'
        if (size(nu) /= size(vt)) error stop 'periodic transport: inconsistent radial grid'
        if (model == ION_TRANSPORT_DRIFT_KINETIC .or. kim_transport_benchmark) then
            call calc_transport_coeffs_ornuhl(size(vt), vt, nu, &
                reference(1, 1, :), reference(1, 2, :), reference(2, 1, :), reference(2, 2, :))
        end if
        if (model == ION_TRANSPORT_FLR .or. kim_transport_benchmark) then
            if (.not. allocated(kim_D_ion_modes)) &
                error stop 'finite-Larmor-radius ion transport tensor is unavailable'
            if (size(kim_D_ion_modes, 3) /= size(vt) .or. &
                    mode < 1 .or. mode > size(kim_D_ion_modes, 4)) &
                error stop 'finite-Larmor-radius ion transport tensor has an invalid shape'
        end if
        select case (model)
        case (ION_TRANSPORT_FLR)
            selected = kim_D_ion_modes(:, :, :, mode)
        case (ION_TRANSPORT_DRIFT_KINETIC)
            selected = reference
        case default
            error stop 'invalid kim_ion_transport_model'
        end select
        if (.not. kim_transport_benchmark) return
        if (allocated(recorded)) then
            if (size(recorded) /= dim_mn .or. size(drift, 3) /= size(vt)) &
                call reset_transport_benchmark()
        end if
        if (.not. allocated(recorded)) then
            allocate(recorded(dim_mn), drift(2, 2, size(vt), dim_mn), &
                flr(2, 2, size(vt), dim_mn))
            recorded = .false.
        end if
        drift(:, :, :, mode) = reference
        flr(:, :, :, mode) = kim_D_ion_modes(:, :, :, mode)
        recorded(mode) = .true.
    end subroutine select_periodic_ion_transport

    subroutine write_transport_benchmark(time_index)
        use h5mod, only: h5_init, h5_open_rw, h5_close, h5_deinit, h5_add, h5_add_string, &
            h5_id, path2out, h5_mode_groupname, create_group_if_not_existent, h5_delete, &
            h5_obj_exists
        use wave_code_data, only: r, m_vals, n_vals
        use grid_mod, only: gg_width
        use kim_wave_code_adapter_m, only: kim_transition_weights, kim_embedding_metadata
        use config_m, only: resolved_ion_ifunc_conservation_model
        use getIfunc_config_m, only: boole_energy_conservation
        integer, intent(in) :: time_index
        real(dp), allocatable :: absolute(:, :, :), relative(:, :, :), antisymmetric(:, :, :)
        character(1024) :: root, group
        integer :: mode, point
        logical :: exists

        if (.not. kim_transport_benchmark) return
        if (.not. allocated(recorded)) return
        if (.not. all(recorded)) error stop 'transport benchmark: incomplete mode collection'
        allocate(absolute(2, 2, size(r)), relative(2, 2, size(r)), antisymmetric(2, 2, size(r)))
        call h5_init()
        call h5_open_rw(path2out, h5_id)
        write(root, '(A,A,A,I0)') '/', trim(h5_mode_groupname), '/TransportBenchmark/', time_index
        call h5_obj_exists(h5_id, trim(root), exists)
        if (exists) call h5_delete(h5_id, trim(root))
        call create_group_if_not_existent(trim(root))
        call h5_add_string(h5_id, trim(root)//'/selected_ion_model', trim(kim_ion_transport_model))
        call h5_add_string(h5_id, trim(root)//'/convention', &
            'CGS; per mode before weighting/smoothing; Fortran tensor axes (2,2,r)')
        call h5_add_string(h5_id, trim(root)//'/relative_residual_definition', &
            'abs(FLR-drift)/max(abs(FLR),abs(drift)); zero/zero=0')
        call h5_add_string(h5_id, trim(root)//'/interpretation', &
            'Full production difference includes FLR, spectral and embedding errors; '// &
            'not an isolated asymptotic-limit error. FLR tensor is symmetric.')
        call h5_add(h5_id, trim(root)//'/kim_ion_conservation_model', &
            resolved_ion_ifunc_conservation_model)
        call h5_add(h5_id, trim(root)//'/drift_energy_conservation', &
            merge(1, 0, boole_energy_conservation))
        call h5_add(h5_id, trim(root)//'/drift_cutoff_halfwidth', 2.0_dp*gg_width)
        do mode = 1, size(recorded)
            write(group, '(A,A,I0,A)') trim(root), '/mode_', mode, '/'
            call create_group_if_not_existent(trim(group))
            do point = 1, size(r)
                call calc_dqli_benchmark_residuals(drift(:, :, point, mode), &
                    flr(:, :, point, mode), absolute(:, :, point), relative(:, :, point))
                antisymmetric(:, :, point) = 0.5_dp*(drift(:, :, point, mode) - &
                    transpose(drift(:, :, point, mode)))
            end do
            call h5_add(h5_id, trim(group)//'m', m_vals(mode))
            call h5_add(h5_id, trim(group)//'n', n_vals(mode))
            call h5_add(h5_id, trim(group)//'r', r, lbound(r), ubound(r))
            call h5_add(h5_id, trim(group)//'transition_weight', kim_transition_weights(:, mode), &
                [1], [size(r)])
            call h5_add(h5_id, trim(group)//'embedding_bounds_and_widths', &
                kim_embedding_metadata(:, mode), [1], [4])
            call h5_add(h5_id, trim(group)//'drift_kinetic', drift(:, :, :, mode), &
                [1, 1, 1], [2, 2, size(r)])
            call h5_add(h5_id, trim(group)//'finite_larmor_radius', flr(:, :, :, mode), &
                [1, 1, 1], [2, 2, size(r)])
            call h5_add(h5_id, trim(group)//'absolute_residual', absolute, &
                lbound(absolute), ubound(absolute))
            call h5_add(h5_id, trim(group)//'relative_residual', relative, &
                lbound(relative), ubound(relative))
            call h5_add(h5_id, trim(group)//'drift_antisymmetric', antisymmetric, &
                lbound(antisymmetric), ubound(antisymmetric))
        end do
        call h5_close(h5_id)
        call h5_deinit()
    end subroutine write_transport_benchmark
end module periodic_transport_benchmark_m
