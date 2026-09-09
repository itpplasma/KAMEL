module periodic_current_diagnostics_m
    !! Saved normalization uses the native trusted current layer, before the
    !! response is interpolated or tapered onto the QL-Balance grid.
    use kim_wave_code_adapter_m, only: kim_current_records
    implicit none
    private
    public :: write_periodic_current_diagnostics
contains
    subroutine write_periodic_current_diagnostics(time_index)
        use h5mod, only: h5_init, h5_open_rw, h5_close, h5_deinit, h5_add, h5_add_string, &
            h5_id, path2out, h5_mode_groupname, create_group_if_not_existent, h5_delete, &
            h5_obj_exists
        use control_mod, only: ihdf5IO
        use wave_code_data, only: m_vals, n_vals
        integer, intent(in) :: time_index
        integer :: mode, nr, entry
        character(24), parameter :: dataset_names(*) = [character(24) :: &
            'm', 'n', 'target_current', 'unit_current', 'achieved_current', 'residual', &
            'scale', 'relative_residual', 'relaxation', 'current_floor', 'max_scale_ratio', &
            'status', 'core_bounds', 'r', 'unit_jpar', 'normalized_jpar']
        real(8) :: relative, denominator
        real(8), allocatable :: profile(:, :)
        character(1024) :: root, group
        logical :: exists, mode_exists

        if (ihdf5IO /= 1) return
        if (.not. allocated(kim_current_records)) return
        call h5_init()
        call h5_open_rw(path2out, h5_id)
        write(root, '(A,A,A,I0)') '/', trim(h5_mode_groupname), &
            '/CurrentNormalization/', time_index
        call h5_obj_exists(h5_id, trim(root), exists)
        ! fortio's h5_delete removes datasets, not groups. Explicitly clear
        ! every old mode's data, including modes absent from the new solve.
        if (exists) then
            call h5_delete(h5_id, trim(root)//'/active')
            mode = 1
            do
                write(group, '(A,A,I0,A)') trim(root), '/mode_', mode, '/'
                call h5_obj_exists(h5_id, trim(group), mode_exists)
                if (.not. mode_exists) exit
                do entry = 1, size(dataset_names)
                    call h5_delete(h5_id, trim(group)//trim(dataset_names(entry)))
                end do
                mode = mode + 1
            end do
        end if
        if (.not. any(kim_current_records%active)) then
            ! No new group for a manual run; mark a reused index explicitly
            ! inactive after clearing its old target and profile datasets.
            if (exists) call h5_add(h5_id, trim(root)//'/active', 0)
            call h5_close(h5_id)
            call h5_deinit()
            return
        end if
        call create_group_if_not_existent(trim(root))
        call h5_add(h5_id, trim(root)//'/active', 1)
        call h5_add_string(h5_id, trim(root)//'/convention', &
            'Complex data: [real,imag]. Unit current=int(r*Jpar dr) in full CGS; '// &
            'target, achieved and residual use c=1 CGS, including 2*pi.')
        call h5_add_string(h5_id, trim(root)//'/layer', &
            'Native KIM grid, resonance +/- dx_asis, before embedding; linear r*J interpolation')
        call h5_add_string(h5_id, trim(root)//'/status_codes', &
            '0=success; 1=invalid configuration; 2=current floor; 3=non-finite/excessive response')
        call h5_add_string(h5_id, trim(root)//'/relaxation_convention', &
            'One stationary update from unit amplitude: s=(1-alpha)+alpha*s_target')
        do mode = 1, size(kim_current_records)
            if (.not. kim_current_records(mode)%active) cycle
            associate(rec => kim_current_records(mode))
                write(group, '(A,A,I0,A)') trim(root), '/mode_', mode, '/'
                call create_group_if_not_existent(trim(group))
                call h5_add(h5_id, trim(group)//'m', m_vals(mode))
                call h5_add(h5_id, trim(group)//'n', n_vals(mode))
                call h5_add(h5_id, trim(group)//'target_current', rec%target)
                call write_complex('unit_current', rec%unit)
                call write_complex('achieved_current', rec%achieved)
                call write_complex('residual', rec%residual)
                call write_complex('scale', rec%scale)
                denominator = max(abs(rec%target), abs(rec%achieved))
                relative = 0.0d0
                if (denominator > 0.0d0) relative = abs(rec%residual)/denominator
                call h5_add(h5_id, trim(group)//'relative_residual', relative)
                call h5_add(h5_id, trim(group)//'relaxation', rec%relaxation)
                call h5_add(h5_id, trim(group)//'current_floor', rec%current_floor)
                call h5_add(h5_id, trim(group)//'max_scale_ratio', rec%max_scale)
                call h5_add(h5_id, trim(group)//'status', rec%status)
                call h5_add(h5_id, trim(group)//'core_bounds', rec%core, [1], [2])
                nr = size(rec%r)
                call h5_add(h5_id, trim(group)//'r', rec%r, [1], [nr])
                allocate(profile(nr, 2))
                profile(:, 1) = real(rec%unit_jpar)
                profile(:, 2) = aimag(rec%unit_jpar)
                call h5_add(h5_id, trim(group)//'unit_jpar', profile, [1, 1], [nr, 2])
                profile(:, 1) = real(rec%normalized_jpar)
                profile(:, 2) = aimag(rec%normalized_jpar)
                call h5_add(h5_id, trim(group)//'normalized_jpar', profile, [1, 1], [nr, 2])
                deallocate(profile)
            end associate
        end do
        call h5_close(h5_id)
        call h5_deinit()
    contains
        subroutine write_complex(name, value)
            character(*), intent(in) :: name
            complex(8), intent(in) :: value
            call h5_add(h5_id, trim(group)//name, [real(value), aimag(value)], [1], [2])
        end subroutine write_complex
    end subroutine write_periodic_current_diagnostics
end module periodic_current_diagnostics_m
