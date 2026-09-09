module periodic_checkpoint_m
    !! Accepted-step continuation, independent of solver and wave-code globals.
    use QLBalance_kinds, only: dp
    use periodic_amplitude_state_m, only: periodic_amplitude_state_t, &
                                          periodic_normalization_version, periodic_phase_policy
    use KAMEL_hdf5_tools
    use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
    implicit none
    private
    public :: periodic_restart_t, pending_periodic_restart
    public :: write_periodic_checkpoint, read_periodic_checkpoint, reset_periodic_restart
    public :: periodic_mode_group

    type :: periodic_restart_t
        logical :: active = .false.
        integer :: accepted_step = 0
        integer :: nstack = 0
        integer :: ramp_mode = 0
        integer :: ramp_up_down = 0
        integer :: iexit = 0
        integer :: hyst_mod_stage = 0
        real(dp) :: timscal_dql = 0.0_dp
        real(dp) :: timscal_dqli = 0.0_dp
        real(dp) :: rate_dql = 0.0_dp
        real(dp) :: time = 0.0_dp
        real(dp) :: next_dt = 0.0_dp
        real(dp) :: tolerance = 0.0_dp
        real(dp) :: br_beta = 0.0_dp
        real(dp) :: br_predicted = 0.0_dp
        real(dp) :: t_hysteresis_turn = 0.0_dp
        real(dp) :: antenna_factor = 1.0_dp
        real(dp) :: antenna_factor_max = 0.0_dp
        real(dp) :: antenna_max_stopping = 0.0_dp
        real(dp) :: t_max_ramp_up = 0.0_dp
        real(dp) :: target_current = 0.0_dp
        real(dp) :: accepted_floor = 1.0e-30_dp
        real(dp) :: accepted_max_scale = 1.0e12_dp
        real(dp) :: current_floor = 1.0e-30_dp
        real(dp) :: max_scale = 1.0e12_dp
        real(dp) :: relaxation = 1.0_dp
        real(dp) :: constant_dt = 0.0_dp
        real(dp) :: tmax = 0.0_dp
        real(dp) :: timescale = 0.0_dp
        real(dp) :: tmax_factor = 0.0_dp
        real(dp) :: stop_time_step = 0.0_dp
        real(dp) :: timstep_min = 0.0_dp
        real(dp) :: t_flattop_begin = 0.0_dp
        real(dp) :: t_flattop_end = 0.0_dp
        real(dp) :: ant_fac_flattop = 0.0_dp
        real(dp) :: delta_t_flattop = 0.0_dp
        real(dp) :: hyst_mod_amp_fac = 0.0_dp
        real(dp) :: hyst_mod_freq = 0.0_dp
        real(dp) :: hyst_mod_phase = 0.0_dp
        real(dp) :: bifcrit_flattop = 0.0_dp
        logical :: constant_dt_enabled = .false.
        logical :: br_stopping = .false.
        logical :: discr_reached = .false.
        logical :: scratch = .false.
        complex(dp) :: boundary_br = (1.0_dp, 0.0_dp)
        real(dp), allocatable :: rc(:)
        real(dp), allocatable :: rb(:)
        real(dp), allocatable :: q(:)
        real(dp), allocatable :: er(:)
        real(dp), allocatable :: vth(:)
        real(dp), allocatable :: dery_equisource(:)
        real(dp), allocatable :: timstep_arr(:)
        real(dp), allocatable :: br_abs(:)
        real(dp), allocatable :: br_abs_time(:)
        real(dp), allocatable :: br_abs_antenna_factor(:)
        real(dp), allocatable :: dqle22_res_time(:)
        real(dp), allocatable :: dae22_res_time(:)
        real(dp), allocatable :: bif_criterion(:)
        complex(dp), allocatable :: Ipar_time(:)
        complex(dp), allocatable :: br_formfactor(:)
        complex(dp), allocatable :: br_vac_res(:)
        real(dp), allocatable :: params(:, :), tim_stack(:), y_stack(:, :)
    end type periodic_restart_t
    type(periodic_restart_t), save :: pending_periodic_restart

contains
    subroutine reset_periodic_restart()
        pending_periodic_restart = periodic_restart_t()
    end subroutine reset_periodic_restart

    function checkpoint_prefix(group) result(prefix)
        character(*), intent(in) :: group
        character(:), allocatable :: prefix
        prefix = trim(group)
        if (len(prefix) == 0) prefix = '/'
        if (prefix(len(prefix):) /= '/') prefix = prefix//'/'
    end function checkpoint_prefix

    function periodic_mode_group(m, n) result(name)
        integer, intent(in) :: m(:), n(:)
        character(128) :: name
        if (size(m) /= size(n) .or. size(m) < 1) error stop 'invalid checkpoint mode list'
        if (size(m) == 1) then
            write (name, '(A,I0,A,I0)') 'f_', m(1), '_', n(1)
        else
            name = 'multi_mode'
        end if
    end function periodic_mode_group

    subroutine write_periodic_checkpoint(file_id, group, m, n, state, payload)
        integer(HID_T), intent(in) :: file_id
        character(*), intent(in) :: group
        integer, intent(in) :: m(:), n(:)
        type(periodic_amplitude_state_t), intent(in) :: state
        type(periodic_restart_t), intent(in) :: payload
        character(:), allocatable :: root
        logical :: exists
        integer :: stack_present, leaf
        character(32), parameter :: leaves(*) = [character(32) :: &
                                                 'version', 'phase', 'residual_convention', &
                                                 'm', 'n', 'amplitude', &
                                                 'unit_current', 'residual', 'status', &
                                              'accepted_target', 'accepted_relaxation', 'scalars', &
                                                 'integers', 'boundary_br', 'params', &
                                                 'rc', 'rb', 'q', &
                                                 'er', 'vth', 'timstep_arr', &
                                                 'br_abs', 'br_abs_time', 'br_abs_antenna_factor', &
                                             'dqle22_res_time', 'dae22_res_time', 'bif_criterion', &
                                                 'Ipar_time', 'br_formfactor', 'br_vac_res', &
                                         'stack_present', 'tim_stack', 'y_stack', 'dery_equisource']
        real(dp) :: scalars(33), boundary(2)
        integer :: integers(10)
        call validate_modes(m, n)
        call validate_payload(payload)
        if (.not. state%initialized) error stop 'checkpoint amplitude state is uninitialized'
        if (.not. allocated(state%accepted)) error stop 'checkpoint accepted amplitude is absent'
        if (size(state%accepted) /= size(m)) error stop 'checkpoint amplitude mode count mismatch'
        if (state%trial_ready) error stop 'checkpoint contains an uncommitted amplitude trial'
        if (.not. allocated(state%accepted_current_unit) .or. &
            .not. allocated(state%accepted_residual) .or. &
            .not. allocated(state%accepted_status) .or. .not. allocated(state%trial) .or. &
            .not. allocated(state%trial_current_unit) .or. &
            .not. allocated(state%trial_residual) .or. .not. allocated(state%trial_status)) &
            error stop 'checkpoint incomplete amplitude state'
        if (size(state%accepted_current_unit) /= size(m) .or. &
            size(state%accepted_residual) /= size(m) .or. size(state%accepted_status) /= size(m) &
            .or. size(state%trial) /= size(m) .or. size(state%trial_current_unit) /= size(m) &
            .or. size(state%trial_residual) /= size(m) .or. size(state%trial_status) /= size(m)) &
            error stop 'checkpoint amplitude record dimensions mismatch'
        if (.not. all(ieee_is_finite([state%accepted_target_current, &
                   state%accepted_relaxation]))) error stop 'checkpoint nonfinite accepted controls'
        if (state%accepted_relaxation <= 0.0_dp .or. state%accepted_relaxation > 1.0_dp) &
            error stop 'checkpoint invalid accepted relaxation'
        if (state%normalization_version /= periodic_normalization_version .or. &
            trim(state%phase_policy) /= periodic_phase_policy) &
            error stop 'checkpoint incompatible amplitude contract'
        if (any(state%trial /= state%accepted) .or. &
            any(state%trial_current_unit /= state%accepted_current_unit) .or. &
            any(state%trial_residual /= state%accepted_residual) .or. &
            any(state%trial_status /= state%accepted_status) .or. &
            state%trial_target_current /= state%accepted_target_current .or. &
            state%trial_relaxation /= state%accepted_relaxation) &
            error stop 'checkpoint requires a committed accepted state'
        if (any(state%accepted_status /= -1 .and. state%accepted_status /= 0)) &
            error stop 'checkpoint contains a failed normalization'
        call check_complex(state%accepted)
        call check_complex(state%accepted_current_unit)
        call check_complex(state%accepted_residual)
        root = checkpoint_prefix(group)//'PeriodicCheckpoint/'
        do leaf = 1, size(leaves)
            call h5_obj_exists(file_id, root//trim(leaves(leaf)), exists)
            if (exists) call h5_delete(file_id, root//trim(leaves(leaf)))
        end do
        call h5_add(file_id, root//'version', periodic_normalization_version)
        call h5_add_string(file_id, root//'phase', periodic_phase_policy)
        call h5_add_string(file_id, root//'residual_convention', &
                           'achieved-target; achieved=(2*pi/c)*scale*integral(r*Jparallel dr)')
        call h5_add(file_id, root//'m', m, [1], [size(m)])
        call h5_add(file_id, root//'n', n, [1], [size(n)])
        call put_complex(file_id, root//'amplitude', state%accepted)
        call put_complex(file_id, root//'unit_current', state%accepted_current_unit)
        call put_complex(file_id, root//'residual', state%accepted_residual)
        call h5_add(file_id, root//'status', state%accepted_status, [1], [size(m)])
        call h5_add(file_id, root//'accepted_target', state%accepted_target_current)
        call h5_add(file_id, root//'accepted_relaxation', state%accepted_relaxation)
        scalars(1) = payload%time
        scalars(2) = payload%next_dt
        scalars(3) = payload%tolerance
        scalars(4) = payload%br_beta
        scalars(5) = payload%br_predicted
        scalars(6) = payload%t_hysteresis_turn
        scalars(7) = payload%antenna_factor
        scalars(8) = payload%antenna_factor_max
        scalars(9) = payload%antenna_max_stopping
        scalars(10) = payload%t_max_ramp_up
        scalars(11) = payload%target_current
        scalars(12) = payload%current_floor
        scalars(13) = payload%max_scale
        scalars(14) = payload%relaxation
        scalars(15) = payload%constant_dt
        scalars(16) = payload%tmax
        scalars(17) = payload%timescale
        scalars(18) = payload%tmax_factor
        scalars(19) = payload%stop_time_step
        scalars(20) = payload%timstep_min
        scalars(21) = payload%t_flattop_begin
        scalars(22) = payload%t_flattop_end
        scalars(23) = payload%ant_fac_flattop
        scalars(24) = payload%delta_t_flattop
        scalars(25) = payload%hyst_mod_amp_fac
        scalars(26) = payload%hyst_mod_freq
        scalars(27) = payload%hyst_mod_phase
        scalars(28) = payload%bifcrit_flattop
        scalars(29) = payload%timscal_dql
        scalars(30) = payload%timscal_dqli
        scalars(31) = payload%rate_dql
        scalars(32) = payload%accepted_floor
        scalars(33) = payload%accepted_max_scale
        integers(1) = payload%accepted_step
        integers(2) = payload%nstack
        integers(3) = payload%ramp_mode
        integers(4) = payload%ramp_up_down
        integers(5) = payload%iexit
        integers(6) = payload%hyst_mod_stage
        integers(7) = merge(1, 0, payload%constant_dt_enabled)
        integers(8) = merge(1, 0, payload%br_stopping)
        integers(9) = merge(1, 0, payload%discr_reached)
        integers(10) = merge(1, 0, payload%scratch)
        boundary = [real(payload%boundary_br, dp), aimag(payload%boundary_br)]
        call h5_add(file_id, root//'scalars', scalars, [1], [size(scalars)])
        call h5_add(file_id, root//'integers', integers, [1], [size(integers)])
        call h5_add(file_id, root//'boundary_br', boundary, [1], [2])
        call h5_add(file_id, root//'dery_equisource', payload%dery_equisource, &
                    [1], [size(payload%dery_equisource)])
        call h5_add(file_id, root//'params', payload%params, [1, 1], shape(payload%params))
        call h5_add(file_id, root//'rc', payload%rc, [1], [size(payload%rc)])
        call h5_add(file_id, root//'rb', payload%rb, [1], [size(payload%rb)])
        call h5_add(file_id, root//'q', payload%q, [1], [size(payload%q)])
        call h5_add(file_id, root//'er', payload%er, [1], [size(payload%er)])
        call h5_add(file_id, root//'vth', payload%vth, [1], [size(payload%vth)])
        call h5_add(file_id, root//'timstep_arr', payload%timstep_arr, [1], &
                    [size(payload%timstep_arr)])
        call h5_add(file_id, root//'br_abs', payload%br_abs, [1], [size(payload%br_abs)])
        call h5_add(file_id, root//'br_abs_time', payload%br_abs_time, [1], &
                    [size(payload%br_abs_time)])
        call h5_add(file_id, root//'br_abs_antenna_factor', payload%br_abs_antenna_factor, &
                    [1], [size(payload%br_abs_antenna_factor)])
        call h5_add(file_id, root//'dqle22_res_time', payload%dqle22_res_time, [1], &
                    [size(payload%dqle22_res_time)])
        call h5_add(file_id, root//'dae22_res_time', payload%dae22_res_time, [1], &
                    [size(payload%dae22_res_time)])
        call h5_add(file_id, root//'bif_criterion', payload%bif_criterion, [1], &
                    [size(payload%bif_criterion)])
        call put_complex(file_id, root//'Ipar_time', payload%Ipar_time)
        call put_complex(file_id, root//'br_formfactor', payload%br_formfactor)
        call put_complex(file_id, root//'br_vac_res', payload%br_vac_res)
        stack_present = merge(1, 0, allocated(payload%tim_stack)) &
                        + merge(2, 0, allocated(payload%y_stack))
        call h5_add(file_id, root//'stack_present', stack_present)
        if (btest(stack_present, 0)) then
            call h5_add(file_id, root//'tim_stack', payload%tim_stack, [1], &
                        [size(payload%tim_stack)])
        end if
        if (btest(stack_present, 1)) then
            call h5_add(file_id, root//'y_stack', payload%y_stack, [1, 1], shape(payload%y_stack))
        end if
    end subroutine write_periodic_checkpoint

    subroutine read_periodic_checkpoint(file_id, group, m, n, state, payload, found)
        integer(HID_T), intent(in) :: file_id
        character(*), intent(in) :: group
        integer, intent(in) :: m(:), n(:)
        type(periodic_amplitude_state_t), intent(inout) :: state
        type(periodic_restart_t), intent(out) :: payload
        logical, intent(out) :: found
        integer, allocatable :: saved_m(:), saved_n(:), statuses(:), mapping(:), integers(:)
        real(dp), allocatable :: scalars(:), boundary(:)
        complex(dp), allocatable :: amplitudes(:), currents(:), residuals(:)
        real(dp) :: target, relaxation
        integer :: version, i, j, stack_present
        character(64) :: phase
        character(:), allocatable :: root
        logical :: old_state

        call state%reset()
        payload = periodic_restart_t()
        found = .false.
        call validate_modes(m, n)
        root = checkpoint_prefix(group)//'PeriodicCheckpoint/'
        call h5_obj_exists(file_id, root, found)
        if (.not. found) then
            call h5_obj_exists(file_id, &
                            checkpoint_prefix(group)//'periodic_amplitude_accepted_real', old_state)
            if (old_state) error stop 'checkpoint v1 lacks a compatible continuation contract'
            call h5_obj_exists(file_id, &
                              checkpoint_prefix(group)//'periodic_normalization_version', old_state)
            if (old_state) error stop 'checkpoint v1 lacks a compatible continuation contract'
            call h5_obj_exists(file_id, &
                               checkpoint_prefix(group)//'periodic_amplitude_trial_real', old_state)
            if (old_state) error stop 'checkpoint partial legacy amplitude state'
            call h5_obj_exists(file_id, &
                               checkpoint_prefix(group)//'periodic_phase_policy', old_state)
            if (old_state) error stop 'checkpoint partial legacy amplitude state'
            return
        end if
        call require_dataset(file_id, root//'version')
        call h5_get(file_id, root//'version', version)
        if (version /= periodic_normalization_version) error stop 'checkpoint version mismatch'
        call require_dataset(file_id, root//'phase')
        call h5_get(file_id, root//'phase', phase)
        if (trim(phase) /= periodic_phase_policy) error stop 'checkpoint phase policy mismatch'
        call get_integer_vector(file_id, root//'m', saved_m)
        call get_integer_vector(file_id, root//'n', saved_n)
        call validate_modes(saved_m, saved_n)
        if (size(saved_m) /= size(m)) error stop 'checkpoint mode count mismatch'
        allocate (mapping(size(m)))
        do i = 1, size(m)
            mapping(i) = 0
            do j = 1, size(saved_m)
                if (m(i) == saved_m(j) .and. n(i) == saved_n(j)) mapping(i) = j
            end do
            if (mapping(i) == 0) error stop 'checkpoint signed mode identity mismatch'
        end do
        call get_complex(file_id, root//'amplitude', amplitudes)
        call get_complex(file_id, root//'unit_current', currents)
        call get_complex(file_id, root//'residual', residuals)
        call get_integer_vector(file_id, root//'status', statuses)
        if (size(amplitudes) /= size(m) .or. size(currents) /= size(m) .or. &
            size(residuals) /= size(m) .or. size(statuses) /= size(m)) &
            error stop 'checkpoint amplitude dimensions mismatch'
        if (any(statuses /= -1 .and. statuses /= 0)) &
            error stop 'checkpoint contains a failed normalization'
        call require_dataset(file_id, root//'accepted_target')
        call require_dataset(file_id, root//'accepted_relaxation')
        call h5_get(file_id, root//'accepted_target', target)
        call h5_get(file_id, root//'accepted_relaxation', relaxation)
        if (.not. all(ieee_is_finite([target, relaxation]))) &
            error stop 'checkpoint amplitude controls are nonfinite'
        if (relaxation <= 0.0_dp .or. relaxation > 1.0_dp) &
            error stop 'checkpoint amplitude relaxation is invalid'
        call get_vector(file_id, root//'scalars', scalars)
        call get_integer_vector(file_id, root//'integers', integers)
        call get_vector(file_id, root//'boundary_br', boundary)
        if (size(scalars) /= 33 .or. size(integers) /= 10 &
            .or. size(boundary) /= 2) error stop 'checkpoint continuation scalar shape mismatch'
        payload%time = scalars(1)
        payload%next_dt = scalars(2)
        payload%tolerance = scalars(3)
        payload%br_beta = scalars(4)
        payload%br_predicted = scalars(5)
        payload%t_hysteresis_turn = scalars(6)
        payload%antenna_factor = scalars(7)
        payload%antenna_factor_max = scalars(8)
        payload%antenna_max_stopping = scalars(9)
        payload%t_max_ramp_up = scalars(10)
        payload%target_current = scalars(11)
        payload%current_floor = scalars(12)
        payload%max_scale = scalars(13)
        payload%relaxation = scalars(14)
        payload%constant_dt = scalars(15)
        payload%tmax = scalars(16)
        payload%timescale = scalars(17)
        payload%tmax_factor = scalars(18)
        payload%stop_time_step = scalars(19)
        payload%timstep_min = scalars(20)
        payload%t_flattop_begin = scalars(21)
        payload%t_flattop_end = scalars(22)
        payload%ant_fac_flattop = scalars(23)
        payload%delta_t_flattop = scalars(24)
        payload%hyst_mod_amp_fac = scalars(25)
        payload%hyst_mod_freq = scalars(26)
        payload%hyst_mod_phase = scalars(27)
        payload%bifcrit_flattop = scalars(28)
        payload%timscal_dql = scalars(29)
        payload%timscal_dqli = scalars(30)
        payload%rate_dql = scalars(31)
        payload%accepted_floor = scalars(32)
        payload%accepted_max_scale = scalars(33)
        payload%accepted_step = integers(1)
        payload%nstack = integers(2)
        payload%ramp_mode = integers(3)
        payload%ramp_up_down = integers(4)
        payload%iexit = integers(5)
        payload%hyst_mod_stage = integers(6)
        if (integers(7) /= 0 .and. integers(7) /= 1) error stop 'invalid checkpoint logical'
        payload%constant_dt_enabled = integers(7) == 1
        if (integers(8) /= 0 .and. integers(8) /= 1) error stop 'invalid checkpoint logical'
        payload%br_stopping = integers(8) == 1
        if (integers(9) /= 0 .and. integers(9) /= 1) error stop 'invalid checkpoint logical'
        payload%discr_reached = integers(9) == 1
        if (integers(10) /= 0 .and. integers(10) /= 1) error stop 'invalid checkpoint logical'
        payload%scratch = integers(10) == 1
        payload%boundary_br = cmplx(boundary(1), boundary(2), dp)
        call get_vector(file_id, root//'dery_equisource', payload%dery_equisource)
        call get_matrix(file_id, root//'params', payload%params)
        call get_vector(file_id, root//'rc', payload%rc)
        call get_vector(file_id, root//'rb', payload%rb)
        call get_vector(file_id, root//'q', payload%q)
        call get_vector(file_id, root//'er', payload%er)
        call get_vector(file_id, root//'vth', payload%vth)
        call get_vector(file_id, root//'timstep_arr', payload%timstep_arr)
        call get_vector(file_id, root//'br_abs', payload%br_abs)
        call get_vector(file_id, root//'br_abs_time', payload%br_abs_time)
        call get_vector(file_id, root//'br_abs_antenna_factor', payload%br_abs_antenna_factor)
        call get_vector(file_id, root//'dqle22_res_time', payload%dqle22_res_time)
        call get_vector(file_id, root//'dae22_res_time', payload%dae22_res_time)
        call get_vector(file_id, root//'bif_criterion', payload%bif_criterion)
        call get_complex(file_id, root//'Ipar_time', payload%Ipar_time)
        call get_complex(file_id, root//'br_formfactor', payload%br_formfactor)
        call get_complex(file_id, root//'br_vac_res', payload%br_vac_res)
        call require_dataset(file_id, root//'stack_present')
        call h5_get(file_id, root//'stack_present', stack_present)
        if (stack_present < 0 .or. stack_present > 3) &
            error stop 'checkpoint adaptive-stack marker is invalid'
        if (btest(stack_present, 0)) call get_vector(file_id, root//'tim_stack', payload%tim_stack)
        if (btest(stack_present, 1)) call get_matrix(file_id, root//'y_stack', payload%y_stack)
        call validate_payload(payload)
        call state%initialize(amplitudes(mapping), currents(mapping), residuals(mapping), &
                              statuses(mapping), target, relaxation)
        payload%active = .true.
    end subroutine read_periodic_checkpoint

    subroutine validate_modes(m, n)
        integer, intent(in) :: m(:), n(:)
        integer :: i, j
        if (size(m) < 1 .or. size(m) /= size(n)) error stop 'checkpoint mode dimensions mismatch'
        do i = 1, size(m)
            do j = 1, i - 1
                if (m(i) == m(j) .and. n(i) == n(j)) error stop 'checkpoint duplicate signed modes'
            end do
        end do
    end subroutine validate_modes

    subroutine validate_payload(payload)
        type(periodic_restart_t), intent(in) :: payload
        if (.not. allocated(payload%rc)) error stop 'checkpoint missing rc'
        if (.not. all(ieee_is_finite(payload%rc))) &
            error stop 'checkpoint nonfinite rc'
        if (.not. allocated(payload%rb)) error stop 'checkpoint missing rb'
        if (.not. all(ieee_is_finite(payload%rb))) &
            error stop 'checkpoint nonfinite rb'
        if (.not. allocated(payload%q)) error stop 'checkpoint missing q'
        if (.not. all(ieee_is_finite(payload%q))) &
            error stop 'checkpoint nonfinite q'
        if (.not. allocated(payload%er)) error stop 'checkpoint missing er'
        if (.not. all(ieee_is_finite(payload%er))) &
            error stop 'checkpoint nonfinite er'
        if (.not. allocated(payload%vth)) error stop 'checkpoint missing vth'
        if (.not. all(ieee_is_finite(payload%vth))) &
            error stop 'checkpoint nonfinite vth'
        if (.not. allocated(payload%timstep_arr)) error stop 'checkpoint missing timstep_arr'
        if (.not. all(ieee_is_finite(payload%timstep_arr))) &
            error stop 'checkpoint nonfinite timstep_arr'
        if (.not. allocated(payload%br_abs)) error stop 'checkpoint missing br_abs'
        if (.not. all(ieee_is_finite(payload%br_abs))) &
            error stop 'checkpoint nonfinite br_abs'
        if (.not. allocated(payload%br_abs_time)) error stop 'checkpoint missing br_abs_time'
        if (.not. all(ieee_is_finite(payload%br_abs_time))) &
            error stop 'checkpoint nonfinite br_abs_time'
        if (.not. allocated(payload%br_abs_antenna_factor)) &
            error stop 'checkpoint missing br_abs_antenna_factor'
        if (.not. all(ieee_is_finite(payload%br_abs_antenna_factor))) &
            error stop 'checkpoint nonfinite br_abs_antenna_factor'
        if (.not. allocated(payload%dqle22_res_time)) &
            error stop 'checkpoint missing dqle22_res_time'
        if (.not. all(ieee_is_finite(payload%dqle22_res_time))) &
            error stop 'checkpoint nonfinite dqle22_res_time'
        if (.not. allocated(payload%dae22_res_time)) error stop 'checkpoint missing dae22_res_time'
        if (.not. all(ieee_is_finite(payload%dae22_res_time))) &
            error stop 'checkpoint nonfinite dae22_res_time'
        if (.not. allocated(payload%bif_criterion)) error stop 'checkpoint missing bif_criterion'
        if (.not. all(ieee_is_finite(payload%bif_criterion))) &
            error stop 'checkpoint nonfinite bif_criterion'
        if (.not. allocated(payload%Ipar_time)) error stop 'checkpoint missing Ipar_time'
        call check_complex(payload%Ipar_time)
        if (.not. allocated(payload%br_formfactor)) error stop 'checkpoint missing br_formfactor'
        call check_complex(payload%br_formfactor)
        if (.not. allocated(payload%br_vac_res)) error stop 'checkpoint missing br_vac_res'
        call check_complex(payload%br_vac_res)
        if (.not. allocated(payload%params)) error stop 'checkpoint missing params'
        if (.not. all(ieee_is_finite(payload%params))) &
            error stop 'checkpoint nonfinite params'
        if (.not. all(ieee_is_finite([ &
                                     payload%time, payload%next_dt, payload%tolerance]))) &
            error stop 'checkpoint nonfinite continuation scalar'
        if (.not. all(ieee_is_finite([ &
                              payload%br_beta, payload%br_predicted, payload%t_hysteresis_turn]))) &
            error stop 'checkpoint nonfinite continuation scalar'
        if (.not. all(ieee_is_finite([ &
              payload%antenna_factor, payload%antenna_factor_max, payload%antenna_max_stopping]))) &
            error stop 'checkpoint nonfinite continuation scalar'
        if (.not. all(ieee_is_finite([ &
                          payload%t_max_ramp_up, payload%target_current, payload%current_floor]))) &
            error stop 'checkpoint nonfinite continuation scalar'
        if (.not. all(ieee_is_finite([ &
                                    payload%max_scale, payload%relaxation, payload%constant_dt]))) &
            error stop 'checkpoint nonfinite continuation scalar'
        if (.not. all(ieee_is_finite([ &
                                     payload%tmax, payload%timescale, payload%tmax_factor]))) &
            error stop 'checkpoint nonfinite continuation scalar'
        if (.not. all(ieee_is_finite([ &
                          payload%stop_time_step, payload%timstep_min, payload%t_flattop_begin]))) &
            error stop 'checkpoint nonfinite continuation scalar'
        if (.not. all(ieee_is_finite([ &
                       payload%t_flattop_end, payload%ant_fac_flattop, payload%delta_t_flattop]))) &
            error stop 'checkpoint nonfinite continuation scalar'
        if (.not. all(ieee_is_finite([ &
                       payload%hyst_mod_amp_fac, payload%hyst_mod_freq, payload%hyst_mod_phase]))) &
            error stop 'checkpoint nonfinite continuation scalar'
        if (.not. all(ieee_is_finite([ &
                                     payload%bifcrit_flattop]))) &
            error stop 'checkpoint nonfinite continuation scalar'
        if (.not. all(ieee_is_finite([payload%timscal_dql, payload%timscal_dqli, &
                               payload%rate_dql]))) error stop 'checkpoint nonfinite QL time scales'
        if (.not. all(ieee_is_finite([payload%accepted_floor, &
                    payload%accepted_max_scale]))) error stop 'checkpoint nonfinite accepted guards'
        if (payload%accepted_floor <= 0.0_dp .or. payload%accepted_max_scale <= 0.0_dp) &
            error stop 'checkpoint invalid accepted guards'
        if (.not. allocated(payload%dery_equisource)) error stop 'checkpoint missing fixed sources'
        if (.not. all(ieee_is_finite(payload%dery_equisource))) &
            error stop 'checkpoint nonfinite fixed sources'
        if (size(payload%dery_equisource) /= size(payload%timstep_arr)) &
            error stop 'checkpoint fixed source dimensions mismatch'
        call check_complex([payload%boundary_br])
        if (payload%accepted_step < 0 .or. payload%time < 0.0_dp .or. &
            payload%next_dt <= 0.0_dp .or. payload%tolerance <= 0.0_dp) &
            error stop 'checkpoint invalid continuation clock'
        if (payload%current_floor <= 0.0_dp .or. payload%max_scale <= 0.0_dp .or. &
            payload%relaxation <= 0.0_dp .or. payload%relaxation > 1.0_dp) &
            error stop 'checkpoint invalid normalization controls'
        if (size(payload%params, 1) /= 4 .or. size(payload%params, 2) /= size(payload%rc)) &
            error stop 'checkpoint profile grid mismatch'
        if (size(payload%rc) < 1 .or. size(payload%rb) < 2) &
            error stop 'checkpoint empty physical grid'
        if (any(payload%rc(2:) <= payload%rc(:size(payload%rc) - 1)) .or. &
            any(payload%rb(2:) <= payload%rb(:size(payload%rb) - 1))) &
            error stop 'checkpoint grids must increase'
        if (size(payload%q) /= size(payload%rb) .or. size(payload%er) /= size(payload%rb) &
            .or. size(payload%vth) /= size(payload%rb)) &
            error stop 'checkpoint background shape mismatch'
        if (payload%nstack < 0) error stop 'checkpoint negative stack size'
        if (allocated(payload%tim_stack)) then
            if (.not. all(ieee_is_finite(payload%tim_stack))) &
                error stop 'nonfinite adaptive time stack'
            if (size(payload%tim_stack) /= size(payload%timstep_arr)) &
                error stop 'checkpoint adaptive time stack shape'
        end if
        if (allocated(payload%y_stack)) then
            if (.not. all(ieee_is_finite(payload%y_stack))) &
                error stop 'nonfinite adaptive profile stack'
            if (payload%nstack > size(payload%y_stack, 2)) &
                error stop 'checkpoint adaptive profile stack shape'
        elseif (payload%nstack /= 0) then
            error stop 'checkpoint missing active adaptive profile stack'
        end if
        if (size(payload%br_abs) /= payload%accepted_step) &
            error stop 'checkpoint history length mismatch: br_abs'
        if (size(payload%br_abs_time) /= payload%accepted_step) &
            error stop 'checkpoint history length mismatch: br_abs_time'
        if (size(payload%br_abs_antenna_factor) /= payload%accepted_step) &
            error stop 'checkpoint history length mismatch: br_abs_antenna_factor'
        if (size(payload%dqle22_res_time) /= payload%accepted_step) &
            error stop 'checkpoint history length mismatch: dqle22_res_time'
        if (size(payload%dae22_res_time) /= payload%accepted_step) &
            error stop 'checkpoint history length mismatch: dae22_res_time'
        if (size(payload%bif_criterion) /= payload%accepted_step) &
            error stop 'checkpoint history length mismatch: bif_criterion'
        if (size(payload%Ipar_time) /= payload%accepted_step) &
            error stop 'checkpoint history length mismatch: Ipar_time'
        if (size(payload%br_formfactor) /= payload%accepted_step) &
            error stop 'checkpoint history length mismatch: br_formfactor'
        if (size(payload%br_vac_res) /= payload%accepted_step) &
            error stop 'checkpoint history length mismatch: br_vac_res'
    end subroutine validate_payload

    subroutine require_dataset(file_id, path)
        integer(HID_T), intent(in) :: file_id
        character(*), intent(in) :: path
        logical :: exists
        call h5_obj_exists(file_id, path, exists)
        if (.not. exists) error stop 'partial periodic checkpoint: missing '//path
    end subroutine require_dataset

    subroutine get_vector(file_id, path, values)
        integer(HID_T), intent(in) :: file_id
        character(*), intent(in) :: path
        real(dp), allocatable, intent(out) :: values(:)
        integer :: lo, hi
        call require_dataset(file_id, path)
        call h5_get_bounds_1(file_id, path, lo, hi)
        if (lo /= 1 .or. hi < 0) error stop 'checkpoint vector bounds mismatch'
        allocate (values(hi))
        call h5_get(file_id, path, values)
        if (.not. all(ieee_is_finite(values))) error stop 'checkpoint nonfinite vector'
    end subroutine get_vector

    subroutine get_integer_vector(file_id, path, values)
        integer(HID_T), intent(in) :: file_id
        character(*), intent(in) :: path
        integer, allocatable, intent(out) :: values(:)
        integer :: lo, hi
        call require_dataset(file_id, path)
        call h5_get_bounds_1(file_id, path, lo, hi)
        if (lo /= 1 .or. hi < 1) error stop 'checkpoint integer bounds mismatch'
        allocate (values(hi))
        call h5_get(file_id, path, values)
    end subroutine get_integer_vector

    subroutine get_matrix(file_id, path, values)
        integer(HID_T), intent(in) :: file_id
        character(*), intent(in) :: path
        real(dp), allocatable, intent(out) :: values(:, :)
        integer :: lo1, lo2, hi1, hi2
        call require_dataset(file_id, path)
        call h5_get_bounds(file_id, path, lo1, lo2, hi1, hi2)
        if (lo1 /= 1 .or. lo2 /= 1 .or. hi1 < 1 .or. hi2 < 0) &
            error stop 'checkpoint matrix bounds mismatch'
        allocate (values(hi1, hi2))
        call h5_get(file_id, path, values)
        if (.not. all(ieee_is_finite(values))) error stop 'checkpoint nonfinite matrix'
    end subroutine get_matrix

    subroutine put_complex(file_id, path, values)
        integer(HID_T), intent(in) :: file_id
        character(*), intent(in) :: path
        complex(dp), intent(in) :: values(:)
        real(dp) :: parts(2, size(values))
        parts(1, :) = real(values, dp)
        parts(2, :) = aimag(values)
        call h5_add(file_id, path, parts, [1, 1], shape(parts))
    end subroutine put_complex

    subroutine get_complex(file_id, path, values)
        integer(HID_T), intent(in) :: file_id
        character(*), intent(in) :: path
        complex(dp), allocatable, intent(out) :: values(:)
        real(dp), allocatable :: parts(:, :)
        call get_matrix(file_id, path, parts)
        if (size(parts, 1) /= 2) error stop 'checkpoint complex shape mismatch'
        values = cmplx(parts(1, :), parts(2, :), dp)
    end subroutine get_complex

    subroutine check_complex(values)
        complex(dp), intent(in) :: values(:)
        if (.not. all(ieee_is_finite(real(values, dp))) .or. &
            .not. all(ieee_is_finite(aimag(values)))) error stop 'checkpoint nonfinite complex data'
    end subroutine check_complex
end module periodic_checkpoint_m
