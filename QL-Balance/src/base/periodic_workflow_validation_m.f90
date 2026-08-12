module periodic_workflow_validation_m
    use QLBalance_kinds, only: dp
    implicit none
    private
    public :: validate_periodic_workflow

contains

    subroutine validate_periodic_workflow(wave_code, kim_run_type, type_of_run, profiles_from_balance, &
            n_modes, m_modes, n_mode_values, target_current, jpar_method, electron_model, ion_model, &
            bparallel_source, benchmark_mode)
        character(*), intent(in) :: wave_code, kim_run_type, type_of_run, jpar_method
        character(*), intent(in) :: electron_model, ion_model, bparallel_source, benchmark_mode
        logical, intent(in) :: profiles_from_balance
        integer, intent(in) :: n_modes, m_modes(:), n_mode_values(:)
        real(dp), intent(in) :: target_current
        integer :: i

        if (trim(wave_code) /= 'KIM' .and. trim(wave_code) /= 'KiLCA') then
            error stop 'wave_code must be KIM or KiLCA'
        end if
        if (trim(wave_code) /= 'KIM') return

        select case (trim(kim_run_type))
        case ('electrostatic_periodic', 'electrostatic', 'electromagnetic', 'flr2', &
              'flr2_benchmark', 'WKB_dispersion')
            continue
        case default
            error stop 'unsupported KIM run type; select an explicit supported solver'
        end select

        if (trim(kim_run_type) /= 'electrostatic_periodic') return
        if (.not. profiles_from_balance) then
            error stop 'periodic KIM production workflow requires profiles_from_balance=.true.'
        end if
        if (trim(type_of_run) /= 'SingleStep' .and. trim(type_of_run) /= 'TimeEvolution' &
            .and. trim(type_of_run) /= 'ParameterScan') then
            error stop 'unsupported QL-Balance run type for periodic KIM workflow'
        end if
        if (n_modes < 1 .or. n_modes > size(m_modes)) error stop 'kim_n_modes is outside its declared bounds'
        if (target_current < 0.0_dp) error stop 'I_par_toroidal must be non-negative'
        do i = 1, n_modes
            if (m_modes(i) == 0 .or. n_mode_values(i) == 0) then
                error stop 'periodic KIM modes require nonzero m and n'
            end if
        end do
        if (trim(jpar_method) /= 'conductivity' .and. trim(jpar_method) /= 'curlB') then
            error stop 'jpar_method must be conductivity or curlB'
        end if
        if (trim(electron_model) /= 'drift_kinetic') then
            error stop 'periodic workflow electrons must use drift_kinetic transport'
        end if
        if (trim(ion_model) /= 'integral') then
            error stop 'periodic workflow ions must use integral transport'
        end if
        if (trim(bparallel_source) /= 'periodic') then
            error stop 'periodic workflow requires the periodic Bparallel source'
        end if
        if (trim(benchmark_mode) /= 'none' .and. trim(benchmark_mode) /= 'drift_kinetic_limit') then
            error stop 'kim_benchmark_mode must be none or drift_kinetic_limit'
        end if
    end subroutine validate_periodic_workflow

end module periodic_workflow_validation_m
