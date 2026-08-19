module periodic_amplitude_state_m
    use QLBalance_kinds, only: dp
    implicit none
    private
    public :: periodic_amplitude_state_t, periodic_amplitudes
    integer, parameter, public :: periodic_normalization_version = 1
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
    contains
        procedure :: initialize => amplitude_initialize
        procedure :: begin_trial => amplitude_begin_trial
        procedure :: accept => amplitude_accept
        procedure :: reject => amplitude_reject
    end type periodic_amplitude_state_t

    type(periodic_amplitude_state_t), save :: periodic_amplitudes

contains
    subroutine amplitude_initialize(self, values, current_unit, residual, status, target_current, relaxation)
        class(periodic_amplitude_state_t), intent(inout) :: self
        complex(dp), intent(in) :: values(:)
        complex(dp), intent(in), optional :: current_unit(:), residual(:)
        integer, intent(in), optional :: status(:)
        real(dp), intent(in), optional :: target_current, relaxation
        if (allocated(self%accepted_current_unit)) deallocate(self%accepted_current_unit)
        if (allocated(self%trial_current_unit)) deallocate(self%trial_current_unit)
        if (allocated(self%accepted_residual)) deallocate(self%accepted_residual)
        if (allocated(self%trial_residual)) deallocate(self%trial_residual)
        if (allocated(self%accepted_status)) deallocate(self%accepted_status)
        if (allocated(self%trial_status)) deallocate(self%trial_status)
        self%target_current = 0.0_dp
        self%relaxation = 1.0_dp
        self%accepted = values
        self%trial = values
        allocate(self%accepted_current_unit(size(values)), self%trial_current_unit(size(values)))
        allocate(self%accepted_residual(size(values)), self%trial_residual(size(values)))
        allocate(self%accepted_status(size(values)), self%trial_status(size(values)))
        self%accepted_current_unit = (0.0_dp, 0.0_dp)
        self%accepted_residual = (0.0_dp, 0.0_dp)
        self%accepted_status = 0
        if (present(current_unit)) self%accepted_current_unit = current_unit
        if (present(residual)) self%accepted_residual = residual
        if (present(status)) self%accepted_status = status
        self%trial_current_unit = self%accepted_current_unit
        self%trial_residual = self%accepted_residual
        self%trial_status = self%accepted_status
        if (present(target_current)) self%target_current = target_current
        if (present(relaxation)) self%relaxation = relaxation
        self%accepted_target_current = self%target_current
        self%trial_target_current = self%target_current
        self%accepted_relaxation = self%relaxation
        self%trial_relaxation = self%relaxation
        self%initialized = .true.
    end subroutine amplitude_initialize

    subroutine amplitude_begin_trial(self, values, current_unit, residual, status, target_current, relaxation)
        class(periodic_amplitude_state_t), intent(inout) :: self
        complex(dp), intent(in) :: values(:)
        complex(dp), intent(in), optional :: current_unit(:), residual(:)
        integer, intent(in), optional :: status(:)
        real(dp), intent(in), optional :: target_current, relaxation
        if (.not. self%initialized .or. size(self%accepted) /= size(values)) then
            call self%initialize(values, current_unit, residual, status, target_current, relaxation)
        else
            self%trial = values
            if (present(current_unit)) self%trial_current_unit = current_unit
            if (present(residual)) self%trial_residual = residual
            if (present(status)) self%trial_status = status
            if (present(target_current)) self%trial_target_current = target_current
            if (present(relaxation)) self%trial_relaxation = relaxation
            self%target_current = self%trial_target_current
            self%relaxation = self%trial_relaxation
        end if
    end subroutine amplitude_begin_trial

    subroutine amplitude_accept(self)
        class(periodic_amplitude_state_t), intent(inout) :: self
        if (.not. self%initialized) error stop 'cannot accept uninitialized amplitude state'
        self%accepted = self%trial
        self%accepted_current_unit = self%trial_current_unit
        self%accepted_residual = self%trial_residual
        self%accepted_status = self%trial_status
        self%accepted_target_current = self%trial_target_current
        self%accepted_relaxation = self%trial_relaxation
        self%target_current = self%accepted_target_current
        self%relaxation = self%accepted_relaxation
    end subroutine amplitude_accept

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
    end subroutine amplitude_reject
end module periodic_amplitude_state_m
