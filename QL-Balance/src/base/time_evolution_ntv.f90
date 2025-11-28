module time_evolution_ntv
    use iso_fortran_env, only: dp => real64
    use omp_lib
    use time_evolution, only: TimeEvolution_t
    use neort_datatypes, only: magfie_data_t, transport_data_t

    implicit none

    type, public, extends(TimeEvolution_t) :: TimeEvolutionNTV_t
    contains
        procedure :: init_balance => initTimeEvolutionNTV
        procedure :: run_balance => runTimeEvolutionNTV
    end type

    private :: initTimeEvolutionNTV
    private :: runTimeEvolutionNTV
    private :: doStep
    private :: neo_rt

    ! for NEO-RT (shared data prepared once)
    real(dp), dimension(:), allocatable :: r
    real(dp), dimension(:), allocatable :: s_tor
    real(dp), dimension(:), allocatable :: Omega_tE
    real(dp), dimension(:, :), allocatable :: plasma_data
    real(dp), dimension(:, :), allocatable :: profile_data

    ! from NEO-RT
    type(magfie_data_t) :: magfie_data
    type(transport_data_t), dimension(:), allocatable :: transport_data

    ! NEO-RT cached parameters for parallel use
    real(dp) :: am1, am2, Z1, Z2

contains

    subroutine initTimeEvolutionNTV(this)
        use baseparam_mod, only: am, Z_i
        use do_magfie_mod, only: read_boozer_file
        use do_magfie_pert_mod, only: read_boozer_pert_file
        use driftorbit, only: pertfile
        use grid_mod, only: rmin, rmax
        use logger, only: set_log_level
        use neort, only: read_and_set_control
        use neort_interface, only: read_equil_file, calculate_s_tor, calculate_coarse_s_tor, &
                                   calculate_Omega_tE_splined, prepare_plasma_data_for_neort, &
                                   prepare_profile_data_for_neort
        use neort_profiles, only: prepare_plasma_splines, prepare_profile_splines
        use spline, only: spline_coeff, spline_val

        class(TimeEvolutionNTV_t), intent(inout) :: this

        real(dp), dimension(:), allocatable :: r_eff, psi_tor, s_tor_equil
        real(dp), dimension(:, :), allocatable :: r_of_s_coeffs, s_of_r_coeffs
        real(dp), dimension(:, :), allocatable :: s_splined, r_splined
        real(dp) :: s_min, s_max

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

        ! NEO-RT initialization
        call read_and_set_control("neo-rt/driftorbit")
        call set_log_level(-1)

        call read_boozer_file("neo-rt/g33353_2900_EQH_MARKL.bc")
        if (pertfile) call read_boozer_pert_file("neo-rt/in_file_pert")

        call prepare_plasma_data_for_neort(plasma_data, r, s_tor)
        call prepare_profile_data_for_neort(profile_data, r, s_tor, Omega_tE)
        call prepare_plasma_splines(S_SIZE, am1, am2, Z1, Z2, plasma_data)
        call prepare_profile_splines(profile_data)
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
        use neort_profiles, only: prepare_plasma_splines, prepare_profile_splines
        use time_evolution, only: doStepBase => doStep, time_ind

        class(TimeEvolutionNTV_t), intent(inout) :: this

        integer :: s_size, s_idx

        call doStepBase(this%TimeEvolution_t)

        ! NEO-RT
        s_size = size(s_tor)

        call prepare_plasma_data_for_neort(plasma_data, r, s_tor)
        call prepare_profile_data_for_neort(profile_data, r, s_tor, Omega_tE)

        call prepare_plasma_splines(s_size, am1, am2, Z1, Z2, plasma_data)
        call prepare_profile_splines(profile_data)

        !$omp parallel do schedule(auto)
        do s_idx = 1, s_size
            call neo_rt(transport_data(s_idx), s_tor(s_idx))
            ! TODO: Apply NEO-RT transport coefficients back to KAMEL
            ! This would involve updating the transport coefficient arrays in grid_mod
            ! For example:
            ! call apply_ntv_transport(D11_ntv, D12_ntv, torque_ntv)
        end do
        !$omp end parallel do
    end subroutine doStep

    subroutine neo_rt(transport_data_, s_val)
        use do_magfie_mod, only: init_magfie_at_s, R0, bfac, psi_pr, q
        use do_magfie_pert_mod, only: init_magfie_pert_at_s, init_mph_from_shared
        use driftorbit, only: pertfile, efac, nopassing, sign_vpar, etamin, etamax
        use neort, only: set_s, compute_transport, set_to_trapped_region
        use neort_freq, only: init_canon_freq_trapped_spline, init_canon_freq_passing_spline
        use neort_magfie, only: init_flux_surface_average
        use neort_profiles, only: init_plasma_at_s, init_profile_at_s, init_thermodynamic_forces

        type(transport_data_t), intent(out) :: transport_data_
        real(8), intent(in) :: s_val

        call set_s(s_val)

        call init_magfie_at_s()
        if (pertfile) call init_magfie_pert_at_s()
        call init_mph_from_shared()

        call init_plasma_at_s()
        call init_profile_at_s(R0, efac, bfac)

        ! subroutine init
        call init_flux_surface_average(s_val)
        call init_canon_freq_trapped_spline()  ! sets etamin and etamax
        if (.not. nopassing) call init_canon_freq_passing_spline()
        sign_vpar = 1
        call set_to_trapped_region(etamin, etamax)  ! sets etamin and etamax again
        ! psi_pr is the torodial flux at plasma boundary, fixed for all s
        call init_thermodynamic_forces(psi_pr, q)  ! psi_pr=const., q is set in do_magfie

        call compute_transport(transport_data_)
    end subroutine neo_rt

end module
