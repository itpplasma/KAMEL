module time_evolution_ntv
    use iso_fortran_env, only: dp => real64
    use time_evolution, only: TimeEvolution_t
    use neort_lib, only: transport_data_t, equilibrium_t, flux_surface_t, neort_init, &
                         neort_get_equilibrium, neort_get_flux_surface, &
                         neort_prepare_splines, neort_compute_at_s

    implicit none

    private

    type, public, extends(TimeEvolution_t) :: TimeEvolutionNTV_t
    contains
        procedure :: init_balance => initTimeEvolutionNTV
        procedure :: run_balance => runTimeEvolutionNTV
    end type

    ! for NEO-RT
    real(dp), dimension(:), allocatable :: r
    real(dp), dimension(:), allocatable :: s_tor
    real(dp) :: am1
    real(dp) :: am2
    real(dp) :: Z1
    real(dp) :: Z2
    real(dp), dimension(:, :), allocatable :: plasma_data
    real(dp), dimension(:, :), allocatable :: profile_data
    real(dp) :: antenna_factor_exponent

    ! from NEO-RT
    type(equilibrium_t) :: equil
    type(transport_data_t), dimension(:), allocatable :: transport_data

contains

    subroutine initTimeEvolutionNTV(this)
        use baseparam_mod, only: am, Z_i, rsepar
        use grid_mod, only: rc, rmin
        use time_evolution, only: antenna_factor_max
        use neort_interface, only: meta_config_neort_t, read_neort_meta_config, &
                                   meta_config_neort_problem, read_equil_file, &
                                   calculate_s_tor, calculate_coarse_s_tor, dpsi_pol_dr, &
                                   outer_s_of_grid, poloidal_flux_problem, &
                                   boozer_s_range_problem, s_grid_problem
        use plasma_parameters, only: qsaf
        use spline, only: spline_coeff, spline_val

        class(TimeEvolutionNTV_t), intent(inout) :: this

        real(dp), dimension(:), allocatable :: r_eff, psi_tor, s_tor_equil
        real(dp), dimension(:, :), allocatable :: r_of_s_coeffs, s_of_r_coeffs
        real(dp), dimension(:, :), allocatable :: s_splined, r_splined, q_splined
        real(dp), dimension(:, :), allocatable :: q_of_r_coeffs
        real(dp), dimension(:), allocatable :: dpsi_pol_ds_neort, dpsi_pol_ds_kamel, s_evaluated
        real(dp) :: s_min, s_max
        type(meta_config_neort_t) :: meta_config
        type(flux_surface_t) :: surface
        integer :: s_size, i

        character(len=*), parameter :: balance_config_file = "balance_conf.nml"

        ! Reject a bad &NEORT namelist before the wave-code setup and the
        ! first output; the mode numbers are filled in after it.
        call read_neort_meta_config(balance_config_file, meta_config, with_modes=.false.)
        if (len_trim(meta_config_neort_problem(meta_config)) > 0) then
            error stop trim(meta_config_neort_problem(meta_config))
        end if

        call this%TimeEvolution_t%init_balance
        this%runType = "TimeEvolutionNTV"
        if (.not. (antenna_factor_max > 0.0_dp)) error stop &
            "TimeEvolutionNTV needs antenna_factor > 0: the reference amplitude of the NEO-RT "// &
            "perturbation"

        ! NEO-RT
        call read_neort_meta_config(balance_config_file, meta_config)
        antenna_factor_exponent = meta_config%antenna_factor_exponent
        s_size = meta_config%amount_of_s

        allocate (r_splined(s_size, 3))
        allocate (r(s_size))
        allocate (s_tor(s_size))
        allocate (plasma_data(s_size, 6))
        allocate (profile_data(s_size, 2))
        allocate (transport_data(s_size))

        ! NEO-RT initialization
        call neort_init(meta_config%config, meta_config%boozer_file, meta_config%boozer_pert_file)
        call neort_get_equilibrium(equil)

        ! Prepare grid splines for NEO-RT
        call read_equil_file(r_eff=r_eff, psi_tor=psi_tor)
        ! The s(r) spline must not extrapolate: NEO-RT needs s inside [0, 1].
        if (rmin < minval(r_eff) .or. rsepar > maxval(r_eff)) error stop &
            "TimeEvolutionNTV: the equilibrium file must cover rmin..rsepar"
        allocate (s_tor_equil(size(psi_tor)))
        call calculate_s_tor(s_tor_equil, psi_tor, equil%toroidal_flux_boundary)

        ! Find s values corresponding to rmin and rsepar by splining s(r)
        allocate (s_of_r_coeffs(size(r_eff) - 1, 5))
        allocate (s_splined(2, 3))
        s_of_r_coeffs = spline_coeff(r_eff, s_tor_equil)
        s_splined = spline_val(s_of_r_coeffs, [rmin, rsepar])
        s_min = s_splined(1, 1)
        ! The grid ends at min(s(rsepar), s_max), and only that range must lie in
        ! 0 < s <= 1. Beyond its outer end the NTV torque is zero (apply_ntv_transport).
        if (len_trim(s_grid_problem(s_min, s_splined(2, 1), meta_config%s_max)) > 0) then
            print *, "NEO-RT s(rmin), s(rsepar), NEORT s_max: ", s_min, s_splined(2, 1), &
                meta_config%s_max
            error stop trim(s_grid_problem(s_min, s_splined(2, 1), meta_config%s_max))
        end if
        s_max = outer_s_of_grid(s_splined(2, 1), meta_config%s_max)

        ! Create an equidistant grid in s from s(rmin) to the outer end
        call calculate_coarse_s_tor(s_tor, s_min, s_max, s_size)
        allocate (r_of_s_coeffs(size(s_tor_equil) - 1, 5))
        r_of_s_coeffs = spline_coeff(s_tor_equil, r_eff)
        r_splined = spline_val(r_of_s_coeffs, s_tor)
        r = r_splined(:, 1)

        ! The Boozer file and KAMEL's equilibrium must agree on the flux surfaces,
        ! including the orientation that fixes the sign of dpsi_pol/ds. NEO-RT
        ! defines Omega_tE through dPhi/ds = -sign_theta psi_pr / (q c) Omega_tE
        ! (neort_profiles), so poloidal_flux_derivative has to be
        ! sign_theta * iota * psi_pr: the denominator NEO-RT itself uses.
        allocate (q_of_r_coeffs(size(rc) - 1, 5), q_splined(s_size, 3))
        allocate (dpsi_pol_ds_neort(s_size), dpsi_pol_ds_kamel(s_size), s_evaluated(s_size))
        q_of_r_coeffs = spline_coeff(rc, qsaf)
        q_splined = spline_val(q_of_r_coeffs, r)
        dpsi_pol_ds_kamel = dpsi_pol_dr(r, q_splined(:, 1)) * r_splined(:, 2)
        do i = 1, s_size
            call neort_get_flux_surface(s_tor(i), surface)
            s_evaluated(i) = surface%s
            dpsi_pol_ds_neort(i) = surface%poloidal_flux_derivative
        end do
        if (len_trim(boozer_s_range_problem(s_tor, s_evaluated)) > 0) then
            error stop trim(boozer_s_range_problem(s_tor, s_evaluated))
        end if
        if (len_trim(poloidal_flux_problem(dpsi_pol_ds_neort, dpsi_pol_ds_kamel)) > 0) then
            error stop trim(poloidal_flux_problem(dpsi_pol_ds_neort, dpsi_pol_ds_kamel))
        end if

        ! Cache species parameters for parallel use (same mass/charge for both species)
        am1 = am
        am2 = am
        Z1 = Z_i
        Z2 = Z_i

        deallocate (r_eff)
        deallocate (psi_tor)
        deallocate (s_tor_equil)
        deallocate (r_of_s_coeffs)
        deallocate (s_of_r_coeffs)
        deallocate (s_splined)
        deallocate (r_splined)
        deallocate (q_of_r_coeffs, q_splined, dpsi_pol_ds_neort, dpsi_pol_ds_kamel, s_evaluated)
    end subroutine initTimeEvolutionNTV

    subroutine runTimeEvolutionNTV(this)
        use time_evolution, only: time_ind, time_step_range

        class(TimeEvolutionNTV_t), intent(inout) :: this

        integer :: first_step, last_step

        call time_step_range(first_step, last_step)
        do time_ind = first_step, last_step
            call doStep(this)
        end do
    end subroutine runTimeEvolutionNTV

    !> The NTV diagnostic written by neo_rt shares the output session of the base step.
    subroutine doStep(this)
        use time_evolution, only: doStepBase => doStep
        use control_mod, only: ihdf5IO
        use h5mod, only: path2out
        use KAMEL_hdf5_tools, only: h5_begin_output_session, h5_end_output_session

        class(TimeEvolutionNTV_t), intent(inout) :: this

        if (ihdf5IO == 1) call h5_begin_output_session(trim(path2out))
        call neo_rt
        call doStepBase(this%TimeEvolution_t)
        if (ihdf5IO == 1) call h5_end_output_session()
    end subroutine doStep

    subroutine neo_rt
        use loading_bar_m, only: updateLoadingBarWithETA
        use neort_interface, only: prepare_plasma_data_for_neort, &
                                   prepare_profile_data_for_neort, apply_ntv_transport

        integer :: s_size, work_idx, s_idx
        real(dp) :: s_val
        integer(kind=8) :: start_count, count_rate
        integer :: completed_steps

        s_size = size(s_tor)

        call prepare_plasma_data_for_neort(plasma_data, r, s_tor)
        call prepare_profile_data_for_neort(profile_data, r, s_tor, equil%major_radius)
        call neort_prepare_splines(s_size, am1, am2, Z1, Z2, plasma_data, profile_data)

        ! Initialize progress tracking
        call system_clock(start_count, count_rate)
        completed_steps = 0

        ! Higher s takes longer to compute, therefore schedule it in reverse
        !$omp parallel do private(s_idx, s_val) schedule(runtime)
        do work_idx = 1, s_size
            s_idx = s_size - work_idx + 1
            s_val = s_tor(s_idx)
            call neort_compute_at_s(s_val, transport_data(s_idx))

            !$omp critical
            completed_steps = completed_steps + 1
            call updateLoadingBarWithETA(completed_steps, s_size, start_count, count_rate, "NEO-RT")
            !$omp end critical
        end do
        !$omp end parallel do

        ! Add torque to global torque_density_ntv array to be used in rhs_balance
        call apply_ntv_transport(r, transport_data, antenna_factor_exponent)

    end subroutine neo_rt

end module
