module control_mod

    use QLBalance_kinds, only: dp

    implicit none

    character(100) :: type_of_run
    character(len=32) :: wave_code = 'KiLCA'  ! 'KiLCA' or 'KIM'
    character(len=1024) :: kim_config_path = './KIM_config.nml'
    logical :: kim_profiles_from_balance = .true.
    integer, parameter :: kim_max_modes = 100
    integer :: kim_n_modes = 0
    integer :: kim_m_list(100) = 0
    integer :: kim_n_list(100) = 0
    integer :: ihdf5IO ! added: Markus Markl
    logical :: paramscan ! added: Markus Markl, 03.03.2021
    logical :: timing_mode         ! added by Markus Markl 06.04.2021
    logical :: debug_mode
    logical :: suppression_mode    ! added by Markus Markl 13.04.2021
    logical :: misalign_diffusion ! trigger the calculation and addition of the diffusion due to misaligned equipotentials and flux surfaces
    logical :: diagnostics_output
    logical :: write_gyro_current
    integer :: irf
    integer :: readfromtimestep ! added by Markus Markl 02.06.2021. Reads the background profiles from hdf5 file in which profiles of a ql time evolution are stored.
    integer :: gyro_current_study
    character(len=1024) :: equil_path ! path to equil file containing q, psi, phi,...
    real(dp) :: eps
    real(dp) :: temperature_limit ! limits ion and electron temperatures from below, in eV

end module control_mod
