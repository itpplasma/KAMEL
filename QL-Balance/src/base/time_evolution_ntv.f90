module time_evolution_ntv
    use iso_fortran_env, only: dp => real64
    use omp_lib
    use time_evolution, only: TimeEvolution_t
    use neort_lib

    implicit none

    private

    type, public, extends(TimeEvolution_t) :: TimeEvolutionNTV_t
    contains
        procedure :: init_balance => initTimeEvolutionNTV
        procedure :: run_balance => runTimeEvolutionNTV
    end type

    ! for NEO-RT (shared data prepared once)
    real(dp), dimension(:), allocatable :: r
    real(dp), dimension(:), allocatable :: s_tor
    real(dp), dimension(:), allocatable :: Omega_tE
    real(dp), dimension(:, :), allocatable :: plasma_data
    real(dp), dimension(:, :), allocatable :: profile_data

    ! from NEO-RT
    type(transport_data_t), dimension(:), allocatable :: transport_data

    ! NEO-RT cached parameters for parallel use
    real(dp) :: am1, am2, Z1, Z2

contains

    subroutine initTimeEvolutionNTV(this)
        use baseparam_mod, only: am, Z_i
        use grid_mod, only: rmin, rmax
        use logger, only: set_log_level
        use neort_interface, only: read_equil_file, calculate_s_tor, calculate_coarse_s_tor, &
                                   calculate_Omega_tE_splined, prepare_plasma_data_for_neort, &
                                   prepare_profile_data_for_neort
        use spline, only: spline_coeff, spline_val

        class(TimeEvolutionNTV_t), intent(inout) :: this

        real(dp), dimension(:), allocatable :: r_eff, psi_tor, s_tor_equil
        real(dp), dimension(:, :), allocatable :: r_of_s_coeffs, s_of_r_coeffs
        real(dp), dimension(:, :), allocatable :: s_splined, r_splined
        real(dp) :: s_min, s_max
        class(config_t), allocatable :: config

        integer, parameter :: S_SIZE = 100

        call this%TimeEvolution_t%init_balance
        this%runType = "TimeEvolutionNTV"

        ! NEO-RT
        ! note: profiles live on rc, derivatives on rb
        allocate (r_splined(S_SIZE, 3))
        allocate (r(S_SIZE))
        allocate (s_tor(S_SIZE))
        allocate (Omega_tE(S_SIZE))
        allocate (plasma_data(S_SIZE, 6))
        allocate (profile_data(S_SIZE, 2))
        allocate (transport_data(S_SIZE))

        ! prepare grid splines for NEO-RT
        call read_equil_file(r_eff=r_eff, psi_tor=psi_tor)
        allocate (s_tor_equil(size(psi_tor)))
        call calculate_s_tor(s_tor_equil, psi_tor)

        ! Find s values corresponding to rmin and rmax by splining s(r)
        allocate (s_of_r_coeffs(size(r_eff) - 1, 5))
        allocate (s_splined(2, 3))
        s_of_r_coeffs = spline_coeff(r_eff, s_tor_equil)
        s_splined = spline_val(s_of_r_coeffs, [rmin, rmax])
        s_min = s_splined(1, 1)
        s_max = s_splined(2, 1)

        call calculate_coarse_s_tor(s_tor, s_min, s_max, S_SIZE)
        allocate (r_of_s_coeffs(size(s_tor_equil) - 1, 5))
        r_of_s_coeffs = spline_coeff(s_tor_equil, r_eff)
        r_splined = spline_val(r_of_s_coeffs, s_tor)
        r = r_splined(:, 1)
        call calculate_Omega_tE_splined(Omega_tE, r)

        ! Cache species parameters for parallel use (same mass/charge for both species)
        am1 = am
        am2 = am
        Z1 = Z_i
        Z2 = Z_i

        ! prepare config struct for NEO-RT
        allocate(config)
        config%epsmn = 0.001
        config%m0 = 6 ! TODO: fix
        config%mph = 2 ! TODO: fix
        config%magdrift = .true.
        config%nopassing = .true.
        config%noshear = .false.
        config%pertfile = .false.
        config%nonlin = .false.
        config%comptorque = .true.
        config%inp_swi = 9  ! AUG
        config%vsteps = 512
        config%log_level = -1  ! silent

        ! NEO-RT initialization
        call neort_init(config, "neo-rt/g33353_2900_EQH_MARKL.bc", "neo-rt/in_file_pert")
        ! TODO: overwrite modes from driftorbit input file with modes from ql balance input

        call prepare_plasma_data_for_neort(plasma_data, r, s_tor)
        call prepare_profile_data_for_neort(profile_data, r, s_tor, Omega_tE)
        call neort_prepare_splines(S_SIZE, am1, am2, Z1, Z2, plasma_data, profile_data)

        deallocate (r_eff)
        deallocate (psi_tor)
        deallocate (s_tor_equil)
        deallocate (r_of_s_coeffs)
        deallocate (s_of_r_coeffs)
        deallocate (s_splined)
        deallocate (r_splined)
    end subroutine initTimeEvolutionNTV

    subroutine runTimeEvolutionNTV(this)
        use time_evolution, only: time_ind, Nstorage

        class(TimeEvolutionNTV_t), intent(inout) :: this

        do time_ind = 1, Nstorage
            call doStep(this)
        end do
    end subroutine runTimeEvolutionNTV

    subroutine doStep(this)
        use neort_interface, only: prepare_plasma_data_for_neort, prepare_profile_data_for_neort
        use time_evolution, only: doStepBase => doStep, time_ind

        class(TimeEvolutionNTV_t), intent(inout) :: this

        integer :: s_size, s_idx

        call doStepBase(this%TimeEvolution_t)

        ! NEO-RT
        s_size = size(s_tor)

        call prepare_plasma_data_for_neort(plasma_data, r, s_tor)
        call prepare_profile_data_for_neort(profile_data, r, s_tor, Omega_tE)
        call neort_prepare_splines(s_size, am1, am2, Z1, Z2, plasma_data, profile_data)

        !$omp parallel do schedule(dynamic)
        do s_idx = 1, s_size
            call neort_compute_at_s(s_tor(s_idx), transport_data(s_idx))
            ! TODO: Apply NEO-RT transport coefficients back to KAMEL
            ! This would involve updating the transport coefficient arrays in grid_mod
            ! For example:
            ! call apply_ntv_transport(D11_ntv, D12_ntv, torque_ntv)
        end do
    end subroutine doStep

end module
