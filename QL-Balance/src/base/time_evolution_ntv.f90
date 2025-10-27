module time_evolution_ntv
    use iso_fortran_env, only: dp => real64
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

    ! for NEO-RT
    real(dp), allocatable :: plasma_data(:, :)
    real(dp), allocatable :: profile_data(:, :)
    real(dp) :: am1, am2, Z1, Z2, s, efac, bfac

    ! from NEO-RT
    type(magfie_data_t) :: magfie_data
    type(transport_data_t), allocatable :: transport_data(:)

contains

    subroutine initTimeEvolutionNTV(this)
        use do_magfie_mod, only: do_magfie_init
        use do_magfie_pert_mod, only: do_magfie_pert_init
        use grid_mod, only: npoic
        use logger, only: set_log_level
        use neort, only: read_and_set_control, init, check_magfie
        use neort_interface, only: prepare_plasma_data_for_neort, prepare_profile_data_for_neort
        use neort_profiles, only: init_profiles, init_plasma_input, init_profile_input
        use baseparam_mod, only: R0 => rtor, am, Z_i

        class(TimeEvolutionNTV_t), intent(inout) :: this

        call this%TimeEvolution_t%init_balance
        this%runType = "TimeEvolutionNTV"

        ! NEO-RT
        allocate (plasma_data(npoic, 6))
        allocate (profile_data(npoic, 2))

        ! set parameters
        s = 0.5d0
        efac = 1d0
        bfac = 1d0
        ! pass same mass and charge for species 1 and 2
        am1 = am
        am2 = am
        Z1 = Z_i
        Z2 = Z_i

        call set_log_level(4)  ! for development purposes

        call read_and_set_control("neo-rt/driftorbit") ! NEO-RT config
        call do_magfie_init("neo-rt/in_file") ! Boozer field file
        ! call do_magfie_pert_init("neo-rt/in_file_pert") ! Boozer perturbed field file
        call read_and_init_profile_input("neo-rt/profile.in", s, R0, efac, bfac)
        call prepare_plasma_data_for_neort(this%plasma_data, this%am1, this%am2, this%Z1, this%Z2)
        call prepare_profile_data_for_neort(this%profile_data)
        call init_profiles(R0)

        call prepare_plasma_data_for_neort(plasma_data)
        call prepare_profile_data_for_neort(profile_data)

        call init_plasma_input(s, npoic, am1, am2, Z1, Z2, plasma_data)
        call init_profile_input(s, R0, efac, bfac, profile_data)

        call init
        call check_magfie(magfie_data)
    end subroutine initTimeEvolutionNTV

    subroutine runTimeEvolutionNTV(this)
        use baseparam_mod, only: R0 => rtor, btor
        use grid_mod, only: npoic, rc
        use neort, only: compute_transport
        use neort_interface, only: prepare_plasma_data_for_neort, prepare_profile_data_for_neort
        use neort_profiles, only: init_profiles, init_plasma_input, init_profile_input, &
                                  init_thermodynamic_forces
        use plasma_parameters, only: qsafb
        use time_evolution, only: time_ind, Nstorage, doStep

        class(TimeEvolutionNTV_t), intent(inout) :: this

        real(dp) :: q, r, psi_pr

        allocate (transport_data(Nstorage))

        q = sum(qsafb) / size(qsafb)
        r = sum(rc) / size(rc)
        psi_pr = r * btor / q

        do time_ind = 1, Nstorage
            call doStep(this%TimeEvolution_t)

            call prepare_plasma_data_for_neort(plasma_data)
            call prepare_profile_data_for_neort(profile_data)

            call init_profiles(R0)
            call init_plasma_input(s, npoic, am1, am2, Z1, Z2, plasma_data)
            call init_profile_input(s, R0, efac, bfac, profile_data)
            call init_thermodynamic_forces(psi_pr, q)

            call compute_transport(transport_data(time_ind))

            ! TODO: Apply NEO-RT transport coefficients back to KAMEL
            ! This would involve updating the transport coefficient arrays in grid_mod
            ! For example:
            ! call apply_ntv_transport(D11_ntv, D12_ntv, torque_ntv)
        end do
    end subroutine runTimeEvolutionNTV

end module
