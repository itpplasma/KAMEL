module kim_solver_m
    !> In-process API seam for KIM.
    !>
    !> kim_solver_t is a state-owning handle over a KIM run: init -> [set_profiles]
    !> -> solve(m,n) -> results -> finalize. It hides the global module state
    !> (config_m, setup_m, grid_m, species_m, equilibrium_m, fields_m) and the
    !> call-ordering / state-reset that callers (e.g. the QL-Balance adapter)
    !> currently hand-code.
    !>
    !> Phase 1 (interface-first): the handle delegates to the existing free
    !> subroutine kim_init, the kim_mod_m factory, and the run-type init/run,
    !> copying results out of the globals. Behaviour is unchanged; only the
    !> contract becomes explicit. Later phases absorb the globals into the type
    !> behind this same public surface. See docs/plans/2026-06-15-kim-solver-api-design.md.

    use KIM_kinds_m, only: dp
    use kim_base_m, only: kim_t

    implicit none
    private

    public :: kim_solver_t, kim_results_t, kim_profiles_t
    public :: KIM_OK, KIM_NOT_SETUP, KIM_BAD_CONFIG, KIM_GRID_COVERAGE, KIM_SOLVE_FAILED

    ! Status codes reported at the seam (never halts the host process).
    integer, parameter :: KIM_OK            = 0
    integer, parameter :: KIM_NOT_SETUP     = 1
    integer, parameter :: KIM_BAD_CONFIG    = 2
    integer, parameter :: KIM_GRID_COVERAGE = 3
    integer, parameter :: KIM_SOLVE_FAILED  = 4

    !> In-memory profile inputs for one solve (optional alternative to file input).
    type :: kim_profiles_t
        real(dp), allocatable :: r(:), n(:), Te(:), Ti(:), q(:), Er(:)
    end type kim_profiles_t

    !> Plain-data results of one solve. Deep copy, valid after the next solve.
    !> Field solution lives on r_field; derived background on r_plasma.
    type :: kim_results_t
        integer :: m = 0, n = 0
        integer :: stat = KIM_OK

        ! field solution (field grid)
        real(dp),    allocatable :: r_field(:)
        complex(dp), allocatable :: Es(:), Ep(:), Er(:), Etheta(:), Ez(:), Br(:)
        complex(dp), allocatable :: jpar(:), jpar_e(:), jpar_i(:)

        ! derived background (plasma grid)
        real(dp), allocatable :: r_plasma(:)
        real(dp), allocatable :: kp(:), ks(:), om_E(:)
        real(dp), allocatable :: nu_e(:), nu_i(:)
        real(dp), allocatable :: B0(:), B0z(:), B0th(:)
    end type kim_results_t

    !> Stateful KIM solver handle.
    type :: kim_solver_t
        private
        class(kim_t), allocatable :: run_type
        logical :: is_setup = .false.
        logical :: has_solved = .false.
        integer :: status = KIM_OK
        type(kim_results_t) :: last
    contains
        procedure :: init        => solver_init
        procedure :: set_profiles => solver_set_profiles
        procedure :: solve       => solver_solve
        procedure :: results     => solver_results
        procedure :: finalize    => solver_finalize
        procedure :: status_code => solver_status_code
    end type kim_solver_t

contains

    !> Read config, optionally inject in-memory profiles, select the run-type
    !> and build grids + first equilibrium. run_type overrides the namelist
    !> type_of_run. profiles, when present, switch KIM to in-memory mode.
    subroutine solver_init(self, config_path, run_type, profiles, stat)
        use config_m, only: nml_config_path, type_of_run, profiles_in_memory
        use kim_mod_m, only: from_kim_factory_get_kim

        class(kim_solver_t), intent(inout) :: self
        character(*), intent(in) :: config_path
        character(*), intent(in), optional :: run_type
        type(kim_profiles_t), intent(in), optional :: profiles
        integer, intent(out), optional :: stat

        logical :: config_exists

        self%status = KIM_OK

        ! Handle-level guard: report a missing config rather than letting the
        ! delegated namelist open() abort the process.
        inquire(file=trim(config_path), exist=config_exists)
        if (.not. config_exists) then
            self%status = KIM_BAD_CONFIG
            if (present(stat)) stat = self%status
            return
        end if

        nml_config_path = trim(config_path)
        if (present(profiles)) profiles_in_memory = .true.

        ! Config read + plasma allocation/init (file profiles skipped when in-memory).
        call kim_init

        if (present(profiles)) call inject_profiles(profiles)

        if (present(run_type)) type_of_run = trim(run_type)

        call from_kim_factory_get_kim(trim(type_of_run), self%run_type)
        call self%run_type%init()

        self%is_setup = .true.
        if (present(stat)) stat = self%status
    end subroutine solver_init

    !> Update the in-memory profiles between solves (e.g. time evolution).
    subroutine solver_set_profiles(self, profiles, stat)
        class(kim_solver_t), intent(inout) :: self
        type(kim_profiles_t), intent(in) :: profiles
        integer, intent(out), optional :: stat

        if (.not. self%is_setup) then
            self%status = KIM_NOT_SETUP
            if (present(stat)) stat = self%status
            return
        end if

        call inject_profiles(profiles)
        self%status = KIM_OK
        if (present(stat)) stat = self%status
    end subroutine solver_set_profiles

    !> Solve for one (m, n) mode: set the mode, ensure the equilibrium is
    !> consistent for it, run, and store results. The per-mode equilibrium
    !> recompute and field reset are the orchestration the QL-Balance adapter
    !> currently hand-codes; here they are owned by the handle.
    subroutine solver_solve(self, m, n, stat)
        use setup_m, only: m_mode, n_mode

        class(kim_solver_t), intent(inout) :: self
        integer, intent(in) :: m, n
        integer, intent(out), optional :: stat

        if (.not. self%is_setup) then
            self%status = KIM_NOT_SETUP
            if (present(stat)) stat = self%status
            return
        end if

        m_mode = m
        n_mode = n

        ! The first solve reuses the equilibrium that init() built for the
        ! configured mode; later solves recompute it for the new (m, n).
        if (self%has_solved) call recompute_equilibrium_for_mode()

        call reset_fields()
        call self%run_type%run()

        call copy_results_from_globals(self%last, m, n)
        self%has_solved = .true.
        self%status = KIM_OK
        if (present(stat)) stat = self%status
    end subroutine solver_solve

    !> Return the last solve's results (deep copy).
    function solver_results(self) result(res)
        class(kim_solver_t), intent(in) :: self
        type(kim_results_t) :: res
        res = self%last
    end function solver_results

    !> Release run state: reset the field buffers and equilibrium background
    !> (so a later init starts clean), drop the run-type, clear the flags.
    subroutine solver_finalize(self)
        use equilibrium_m, only: B0, B0z, B0th

        class(kim_solver_t), intent(inout) :: self

        call reset_fields()
        if (allocated(B0))   deallocate(B0)
        if (allocated(B0z))  deallocate(B0z)
        if (allocated(B0th)) deallocate(B0th)
        if (allocated(self%run_type)) deallocate(self%run_type)
        self%is_setup = .false.
        self%has_solved = .false.
        self%status = KIM_OK
    end subroutine solver_finalize

    !> Last recorded status code.
    integer function solver_status_code(self) result(code)
        class(kim_solver_t), intent(in) :: self
        code = self%status
    end function solver_status_code

    ! ----------------------------------------------------------------------
    ! Internal helpers
    ! ----------------------------------------------------------------------

    subroutine inject_profiles(profiles)
        use species_m, only: set_profiles_from_arrays
        type(kim_profiles_t), intent(in) :: profiles
        call set_profiles_from_arrays(profiles%r, profiles%n, profiles%Te, &
                                      profiles%Ti, profiles%q, profiles%Er, &
                                      size(profiles%r))
    end subroutine inject_profiles

    !> Recompute the background equilibrium for the current (m_mode, n_mode).
    subroutine recompute_equilibrium_for_mode()
        use equilibrium_m, only: calculate_equil, interpolate_equil
        use species_m, only: plasma, set_plasma_quantities
        use grid_m, only: rg_grid

        call calculate_equil(.false.)
        call set_plasma_quantities(plasma)
        call interpolate_equil(rg_grid%xb)
    end subroutine recompute_equilibrium_for_mode

    !> Deallocate the global field buffers so the next run() re-allocates cleanly.
    subroutine reset_fields()
        use fields_m, only: EBdat

        if (allocated(EBdat%r_grid))            deallocate(EBdat%r_grid)
        if (allocated(EBdat%Br))                deallocate(EBdat%Br)
        if (allocated(EBdat%Apar))              deallocate(EBdat%Apar)
        if (allocated(EBdat%E_perp_psi))        deallocate(EBdat%E_perp_psi)
        if (allocated(EBdat%E_perp))            deallocate(EBdat%E_perp)
        if (allocated(EBdat%E_perp_MA))         deallocate(EBdat%E_perp_MA)
        if (allocated(EBdat%Er))                deallocate(EBdat%Er)
        if (allocated(EBdat%Etheta))            deallocate(EBdat%Etheta)
        if (allocated(EBdat%Ez))                deallocate(EBdat%Ez)
        if (allocated(EBdat%Es))                deallocate(EBdat%Es)
        if (allocated(EBdat%Ep))                deallocate(EBdat%Ep)
        if (allocated(EBdat%jpar))              deallocate(EBdat%jpar)
        if (allocated(EBdat%jpar_e))            deallocate(EBdat%jpar_e)
        if (allocated(EBdat%jpar_i))            deallocate(EBdat%jpar_i)
        if (allocated(EBdat%Phi))               deallocate(EBdat%Phi)
        if (allocated(EBdat%Phi_e))             deallocate(EBdat%Phi_e)
        if (allocated(EBdat%Phi_i))             deallocate(EBdat%Phi_i)
        if (allocated(EBdat%Phi_aligned))       deallocate(EBdat%Phi_aligned)
        if (allocated(EBdat%Phi_MA))            deallocate(EBdat%Phi_MA)
        if (allocated(EBdat%Phi_MA_ideal))      deallocate(EBdat%Phi_MA_ideal)
        if (allocated(EBdat%Phi_MA_asymptotic)) deallocate(EBdat%Phi_MA_asymptotic)
    end subroutine reset_fields

    !> Deep-copy the solve outputs from the global state into a results record.
    !> Every field is guarded: a run-type only fills the buffers it produces.
    subroutine copy_results_from_globals(res, m, n)
        use fields_m, only: EBdat
        use species_m, only: plasma
        use equilibrium_m, only: B0, B0z, B0th

        type(kim_results_t), intent(out) :: res
        integer, intent(in) :: m, n

        res%m = m
        res%n = n
        res%stat = KIM_OK

        ! field solution (field grid)
        if (allocated(EBdat%r_grid)) res%r_field = EBdat%r_grid
        if (allocated(EBdat%Es))     res%Es      = EBdat%Es
        if (allocated(EBdat%Ep))     res%Ep      = EBdat%Ep
        if (allocated(EBdat%Er))     res%Er      = EBdat%Er
        if (allocated(EBdat%Etheta)) res%Etheta  = EBdat%Etheta
        if (allocated(EBdat%Ez))     res%Ez      = EBdat%Ez
        if (allocated(EBdat%Br))     res%Br      = EBdat%Br
        if (allocated(EBdat%jpar))   res%jpar    = EBdat%jpar
        if (allocated(EBdat%jpar_e)) res%jpar_e  = EBdat%jpar_e
        if (allocated(EBdat%jpar_i)) res%jpar_i  = EBdat%jpar_i

        ! derived background (plasma grid)
        if (allocated(plasma%r_grid)) res%r_plasma = plasma%r_grid
        if (allocated(plasma%kp))     res%kp       = plasma%kp
        if (allocated(plasma%ks))     res%ks       = plasma%ks
        if (allocated(plasma%om_E))   res%om_E     = plasma%om_E
        if (allocated(plasma%spec)) then
            if (lbound(plasma%spec, 1) <= 0 .and. 0 <= ubound(plasma%spec, 1)) then
                if (allocated(plasma%spec(0)%nu)) res%nu_e = plasma%spec(0)%nu
            end if
            if (lbound(plasma%spec, 1) <= 1 .and. 1 <= ubound(plasma%spec, 1)) then
                if (allocated(plasma%spec(1)%nu)) res%nu_i = plasma%spec(1)%nu
            end if
        end if

        ! equilibrium background
        if (allocated(B0))   res%B0   = B0
        if (allocated(B0z))  res%B0z  = B0z
        if (allocated(B0th)) res%B0th = B0th
    end subroutine copy_results_from_globals

end module kim_solver_m
