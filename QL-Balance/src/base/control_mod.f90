module control_mod

    use QLBalance_kinds, only: dp

    implicit none

    character(100) :: type_of_run
    character(len=32) :: wave_code = 'KiLCA'  ! 'KiLCA' or 'KIM'
    character(len=1024) :: kim_config_path = './KIM_config.nml'
    character(len=64) :: kim_config_sha256 = 'unavailable'
    ! Preserve existing KIM decks; periodic coupling requires explicit selection.
    character(len=32) :: kim_run_type = 'electromagnetic'
    character(len=32) :: kim_electron_transport_model = 'drift_kinetic'
    character(len=32) :: kim_ion_transport_model = 'finite_larmor_radius'
    character(len=32) :: kim_bparallel_source = 'disabled'
    character(len=32) :: kim_benchmark_mode = 'none'
    real(dp) :: kim_current_floor = 1.0e-30_dp
    real(dp) :: kim_current_max_scale = 1.0e12_dp
    real(dp) :: kim_current_relaxation = 1.0_dp
    ! Periodic amplitude: 'target_current' fixes int(j_par) to I_par_toroidal;
    ! 'delta_prime' matches the constant-psi layer to the cylindrical ideal
    ! outer region (I_par_toroidal = fully shielded limit).
    character(len=32) :: kim_periodic_matching = 'target_current'
    ! First periodic TimeEvolution response: 'target' normalizes it like every
    ! later response; 'unit' keeps the constant-psi unit drive (Br_ref = 1 G)
    ! for the first transport step (behaviour before this switch).
    character(len=32) :: kim_periodic_initial_drive = 'target'
    ! Gaussian width (std, in radial grid cells) applied to the E x B
    ! frequency used by the electron drift-kinetic coefficients only;
    ! 0 = pointwise (default, previous behaviour). Diagnostic switch.
    real(dp) :: electron_er_smooth_cells = 0.0_dp
    ! Periodic KIM TimeEvolution: ramp the target current with the coil
    ! current, I = I_par_toroidal*sqrt(antenna_factor/antenna_factor_max),
    ! antenna_factor from ramp_up_mode as for KiLCA. .false.: fixed target
    ! (previous behaviour).
    logical :: kim_periodic_coil_ramp = .false.
    ! Ideal wall radius for the outer region [cm]; <= 0: no wall (vacuum decay).
    real(dp) :: kim_periodic_wall_radius = -1.0_dp
    logical :: kim_transport_benchmark = .false.
    logical :: kim_profiles_from_balance = .true.
    integer, parameter :: kim_max_modes = 100
    integer :: kim_n_modes = 0
    integer :: kim_m_list(100) = 0
    integer :: kim_n_list(100) = 0
    integer :: ihdf5IO ! added: Markus Markl
    logical :: paramscan ! added: Markus Markl, 03.03.2021
    logical :: timing_mode         ! added by Markus Markl 06.04.2021
    integer :: log_level = 3         ! maps to LVL_INFO
    logical :: suppression_mode    ! added by Markus Markl 13.04.2021
    !> Treatment of negative equilibrium source densities (local sinks):
    !> 'explicit' keeps them fixed, 'implicit' scales them with the local
    !> density and pressure relative to the initial state.
    character(len=16) :: equilibrium_sink_mode = 'explicit'
    logical :: misalign_diffusion ! trigger the calculation and addition of the diffusion due to misaligned equipotentials and flux surfaces
    character(len=32) :: jpar_method = 'conductivity'  ! 'conductivity' or 'curlB'
    integer :: data_verbosity = 1    ! standard output
    logical :: write_gyro_current
    integer :: irf
    integer :: readfromtimestep ! added by Markus Markl 02.06.2021. Reads the background profiles from hdf5 file in which profiles of a ql time evolution are stored.
    integer :: gyro_current_study
    character(len=1024) :: equil_path ! path to equil file containing q, psi, phi,...
    real(dp) :: eps
    real(dp) :: temperature_limit ! limits ion and electron temperatures from below, in eV

    integer, parameter :: ION_TRANSPORT_INVALID = 0
    integer, parameter :: ION_TRANSPORT_FLR = 1
    integer, parameter :: ION_TRANSPORT_DRIFT_KINETIC = 2

contains

    pure integer function ion_transport_model_id(model) result(model_id)
        character(*), intent(in) :: model

        select case (trim(model))
        case ('finite_larmor_radius')
            model_id = ION_TRANSPORT_FLR
        case ('drift_kinetic')
            model_id = ION_TRANSPORT_DRIFT_KINETIC
        case default
            model_id = ION_TRANSPORT_INVALID
        end select
    end function ion_transport_model_id

    pure logical function is_time_evolution_run(run_type)
        ! TimeEvolutionNTV extends TimeEvolution and shares its time stepping,
        ! histories and periodic-amplitude handling.
        character(*), intent(in) :: run_type

        select case (trim(run_type))
        case ('TimeEvolution', 'TimeEvolutionNTV')
            is_time_evolution_run = .true.
        case default
            is_time_evolution_run = .false.
        end select
    end function is_time_evolution_run

end module control_mod
