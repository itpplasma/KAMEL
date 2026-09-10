module periodic_amplitude_state_m
    use QLBalance_kinds, only: dp
    use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
    implicit none
    private
    public :: periodic_amplitude_state_t, periodic_amplitudes
    integer, parameter, public :: periodic_normalization_version = 2
    character(len=*), parameter, public :: periodic_phase_policy = 'complex-current'

    type :: periodic_amplitude_state_t
        complex(dp), allocatable :: accepted(:), trial(:)
        complex(dp), allocatable :: accepted_current_unit(:), trial_current_unit(:)
        complex(dp), allocatable :: accepted_residual(:), trial_residual(:)
        integer, allocatable :: accepted_status(:), trial_status(:)
        real(dp) :: target_current = 0.0_dp
        real(dp) :: relaxation = 1.0_dp
        real(dp) :: accepted_target_current = 0.0_dp, trial_target_current = 0.0_dp
        real(dp) :: accepted_relaxation = 1.0_dp, trial_relaxation = 1.0_dp
        integer :: normalization_version = periodic_normalization_version
        character(len=32) :: phase_policy = periodic_phase_policy
        logical :: initialized = .false.
        logical :: trial_ready = .false.
    contains
        procedure :: reset => amplitude_reset
        procedure :: initialize => amplitude_initialize
        procedure :: begin_trial => amplitude_begin_trial
        procedure :: accept => amplitude_accept
        procedure :: reject => amplitude_reject
    end type periodic_amplitude_state_t

    type(periodic_amplitude_state_t), save :: periodic_amplitudes

contains
    subroutine amplitude_reset(self)
        class(periodic_amplitude_state_t), intent(inout) :: self
        if (allocated(self%accepted)) deallocate (self%accepted)
        if (allocated(self%trial)) deallocate (self%trial)
        if (allocated(self%accepted_current_unit)) deallocate (self%accepted_current_unit)
        if (allocated(self%trial_current_unit)) deallocate (self%trial_current_unit)
        if (allocated(self%accepted_residual)) deallocate (self%accepted_residual)
        if (allocated(self%trial_residual)) deallocate (self%trial_residual)
        if (allocated(self%accepted_status)) deallocate (self%accepted_status)
        if (allocated(self%trial_status)) deallocate (self%trial_status)
        self%target_current = 0.0_dp
        self%relaxation = 1.0_dp
        self%accepted_target_current = 0.0_dp
        self%trial_target_current = 0.0_dp
        self%accepted_relaxation = 1.0_dp
        self%trial_relaxation = 1.0_dp
        self%normalization_version = periodic_normalization_version
        self%phase_policy = periodic_phase_policy
        self%initialized = .false.
        self%trial_ready = .false.
    end subroutine amplitude_reset

    subroutine amplitude_initialize(self, values, current_unit, residual, status, &
                                    target_current, relaxation)
        class(periodic_amplitude_state_t), intent(inout) :: self
        complex(dp), intent(in) :: values(:)
        complex(dp), intent(in), optional :: current_unit(:), residual(:)
        integer, intent(in), optional :: status(:)
        real(dp), intent(in), optional :: target_current, relaxation
        type(periodic_amplitude_state_t) :: candidate
        integer :: stat

        ! Validate a complete candidate before replacing an existing accepted state.
      call candidate%begin_trial(values, current_unit, residual, status, target_current, relaxation)
        call candidate%accept(stat)
        if (stat /= 0) error stop 'cannot initialize invalid amplitude state'
        call self%reset()
        call self%begin_trial(candidate%accepted, candidate%accepted_current_unit, &
                              candidate%accepted_residual, candidate%accepted_status, &
                              candidate%accepted_target_current, candidate%accepted_relaxation)
        call self%accept()
    end subroutine amplitude_initialize

    subroutine amplitude_begin_trial(self, values, current_unit, residual, status, &
                                     target_current, relaxation)
        class(periodic_amplitude_state_t), intent(inout) :: self
        complex(dp), intent(in) :: values(:)
        complex(dp), intent(in), optional :: current_unit(:), residual(:)
        integer, intent(in), optional :: status(:)
        real(dp), intent(in), optional :: target_current, relaxation
        integer :: n

        n = size(values)
        if (n == 0) error stop 'amplitude state requires at least one mode'
        if (self%initialized) then
            if (.not. allocated(self%accepted)) error stop 'accepted amplitude storage missing'
            if (size(self%accepted) /= n) error stop 'amplitude trial mode count changed'
        end if
        if (present(current_unit)) then
            if (size(current_unit) /= n) error stop 'amplitude current shape mismatch'
        end if
        if (present(residual)) then
            if (size(residual) /= n) error stop 'amplitude residual shape mismatch'
        end if
        if (present(status)) then
            if (size(status) /= n) error stop 'amplitude status shape mismatch'
        end if

        self%trial = values
        ! Omitted metadata belongs to the accepted state, never an earlier trial.
        if (self%initialized) then
            self%trial_current_unit = self%accepted_current_unit
            self%trial_residual = self%accepted_residual
            self%trial_status = self%accepted_status
            self%trial_target_current = self%accepted_target_current
            self%trial_relaxation = self%accepted_relaxation
        else
            self%trial_current_unit = spread(cmplx(0.0_dp, 0.0_dp, dp), 1, n)
            self%trial_residual = self%trial_current_unit
            self%trial_status = spread(0, 1, n)
            self%trial_target_current = 0.0_dp
            self%trial_relaxation = 1.0_dp
        end if
        if (present(current_unit)) self%trial_current_unit = current_unit
        if (present(residual)) self%trial_residual = residual
        if (present(status)) self%trial_status = status
        if (present(target_current)) self%trial_target_current = target_current
        if (present(relaxation)) self%trial_relaxation = relaxation
        self%target_current = self%trial_target_current
        self%relaxation = self%trial_relaxation
        self%trial_ready = .true.
    end subroutine amplitude_begin_trial

    subroutine amplitude_accept(self, stat)
        class(periodic_amplitude_state_t), intent(inout) :: self
        integer, intent(out), optional :: stat
        integer :: validation_status

        ! A failed candidate remains inspectable; accepted storage is never changed.
        validation_status = trial_validation_status(self)
        if (present(stat)) stat = validation_status
        if (validation_status /= 0) then
            if (present(stat)) return
            error stop 'cannot accept invalid amplitude state'
        end if
        self%accepted = self%trial
        self%accepted_current_unit = self%trial_current_unit
        self%accepted_residual = self%trial_residual
        self%accepted_status = self%trial_status
        self%accepted_target_current = self%trial_target_current
        self%accepted_relaxation = self%trial_relaxation
        self%target_current = self%accepted_target_current
        self%relaxation = self%accepted_relaxation
        self%initialized = .true.
        self%trial_ready = .false.
    end subroutine amplitude_accept

    integer function trial_validation_status(self) result(stat)
        class(periodic_amplitude_state_t), intent(in) :: self
        integer :: n

        stat = 1
        if (.not. self%trial_ready) return
        if (.not. allocated(self%trial)) return
        if (.not. allocated(self%trial_current_unit)) return
        if (.not. allocated(self%trial_residual)) return
        if (.not. allocated(self%trial_status)) return
        n = size(self%trial)
        if (n == 0) return
        if (size(self%trial_current_unit) /= n) return
        if (size(self%trial_residual) /= n) return
        if (size(self%trial_status) /= n) return
        if (self%initialized) then
            if (.not. allocated(self%accepted)) return
            if (size(self%accepted) /= n) return
        end if
        stat = 2
        if (.not. all(ieee_is_finite(real(self%trial, dp)))) return
        if (.not. all(ieee_is_finite(aimag(self%trial)))) return
        if (.not. all(ieee_is_finite(real(self%trial_current_unit, dp)))) return
        if (.not. all(ieee_is_finite(aimag(self%trial_current_unit)))) return
        if (.not. all(ieee_is_finite(real(self%trial_residual, dp)))) return
        if (.not. all(ieee_is_finite(aimag(self%trial_residual)))) return
        if (.not. ieee_is_finite(self%trial_target_current)) return
        if (.not. ieee_is_finite(self%trial_relaxation)) return
        if (self%trial_relaxation <= 0.0_dp .or. self%trial_relaxation > 1.0_dp) return
        stat = 3
        ! Zero is successful normalization; -1 denotes a manual-amplitude mode.
        if (any(self%trial_status /= 0 .and. self%trial_status /= -1)) return
        stat = 0
    end function trial_validation_status

    subroutine amplitude_reject(self)
        class(periodic_amplitude_state_t), intent(inout) :: self
        if (.not. self%initialized) error stop 'cannot reject uninitialized amplitude state'
        self%trial = self%accepted
        self%trial_current_unit = self%accepted_current_unit
        self%trial_residual = self%accepted_residual
        self%trial_status = self%accepted_status
        self%trial_target_current = self%accepted_target_current
        self%trial_relaxation = self%accepted_relaxation
        self%target_current = self%accepted_target_current
        self%relaxation = self%accepted_relaxation
        self%trial_ready = .false.
    end subroutine amplitude_reject
end module periodic_amplitude_state_m
