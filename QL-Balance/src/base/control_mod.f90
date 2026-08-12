module control_mod

    use QLBalance_kinds, only: dp

    implicit none

    character(100) :: type_of_run
    character(len=32) :: wave_code = 'KiLCA'  ! 'KiLCA' or 'KIM'
    character(len=1024) :: kim_config_path = './KIM_config.nml'
    character(len=32) :: kim_run_type = 'electrostatic_periodic'
    character(len=32) :: kim_electron_transport_model = 'drift_kinetic'
    character(len=32) :: kim_ion_transport_model = 'finite_larmor_radius'
    character(len=32) :: kim_bparallel_source = 'periodic'
    character(len=32) :: kim_benchmark_mode = 'none'
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

end module control_mod
