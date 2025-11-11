!> @brief Module for plasma parameter storage and management
!> @details Stores plasma profiles (density, temperature, rotation) on radial grid.
!>          All units are in cgs unless explicitly noted otherwise.
module plasma_parameters
    use control_mod

    implicit none

    !> @brief Main plasma parameter array (nparams, npoic)
    !> @details params(1, :) = electron density [cm^-3]
    !>          params(2, :) = toroidal rotation frequency [rad/s]
    !>          params(3, :) = electron temperature [erg]
    !>          params(4, :) = ion temperature [erg]
    real(dp), dimension(:, :), allocatable :: params

    !> @brief Time derivative of plasma parameters [same units as params, per second]
    real(dp), dimension(:, :), allocatable :: dot_params

    !> @brief Radial derivatives of plasma parameters (linear and nonlinear contributions)
    !> @details Same units as params but per cm (e.g., [cm^-3/cm] = [cm^-4])
    real(dp), dimension(:, :), allocatable :: ddr_params, ddr_params_nl

    !> @brief Storage for parameters at previous timesteps (Runge-Kutta integration)
    real(dp), dimension(:, :), allocatable :: params_beg, params_begbeg

    !> @brief Storage for numerator and denominator in quasilinear transport calculations
    real(dp), dimension(:, :), allocatable :: params_num, params_denom

    !> @brief Backup storage for plasma parameters at current boundary
    real(dp), dimension(:, :), allocatable :: params_b

    !> @brief Initial plasma parameters (stored at t=0 for reference)
    real(dp), dimension(:, :), allocatable :: init_params

    !> @brief Linear response parameters and boundary values
    real(dp), dimension(:, :), allocatable :: params_lin, params_b_lin

    !> @brief Safety factor profiles
    !> @details qsafb(npoib) = q at boundary grid points (b-grid)
    !>          qsaf(npoic)  = q at cell centers (c-grid)
    !>          These are DIMENSIONLESS (cgs units but no conversion needed)
    real(dp), dimension(:), allocatable :: qsafb, qsaf

    !> @brief Storage arrays to hold initial background profiles
    !> @details hold_n     = initial density [cm^-3]
    !>          hold_Te    = initial electron temperature [erg]
    !>          hold_Ti    = initial ion temperature [erg]
    !>          hold_Vz    = initial toroidal rotation frequency [rad/s]
    !>          hold_dphi0 = initial electric potential gradient [statV/cm]
    real(dp), dimension(:), allocatable :: hold_n, hold_Te, hold_Ti, hold_Vz, hold_dphi0

contains

    subroutine alloc_hold_parameters
        use grid_mod, only: npoib
        use wave_code_data, only: idPhi0

        allocate (hold_n(npoib))
        allocate (hold_Vz(npoib))
        allocate (hold_Te(npoib))
        allocate (hold_Ti(npoib))
        allocate (hold_dphi0(npoib))

        hold_n = params(1, :)
        hold_Vz = params(2, :)
        hold_Te = params(3, :)
        hold_Ti = params(4, :)
        hold_dphi0 = idPhi0
    end subroutine

    subroutine limit_temps_from_below
        use baseparam_mod, only: ev, rsepar
        use grid_mod, only: npoic
        use wave_code_data, only: r

        integer :: ipoi

        do ipoi = 1, npoic
            if (.true.) then
                ! convert to erg
                params(3, ipoi) = max(params(3, ipoi), temperature_limit * ev)
                params(4, ipoi) = max(params(4, ipoi), temperature_limit * ev)
            else
                !> Quick fix of steady state solution. Keep boundary inside the separatrix.
                if (r(ipoi) > rsepar - 0.5d0) then
                    params(3, ipoi) = hold_Te(ipoi)
                    params(4, ipoi) = hold_Ti(ipoi)
                end if
            end if
        end do
    end subroutine limit_temps_from_below

    subroutine init_background_profiles
        use baseparam_mod, only: ev, rtor
        use grid_mod, only: npoic
        use wave_code_data, only: q, n, Vz, Te, Ti

        integer :: ipoi

        do ipoi = 1, npoic
            ! Safety factor: dimensionless, averaged to cell centers
            qsaf(ipoi) = 0.5 * (q(ipoi) + q(ipoi + 1))

            ! Electron density [cm^-3]
            params(1, ipoi) = 0.5 * (n(ipoi) + n(ipoi + 1))
            init_params(1, ipoi) = params(1, ipoi)

            ! Toroidal rotation frequency [rad/s]
            params(2, ipoi) = 0.5 * (Vz(ipoi) + Vz(ipoi + 1)) / rtor
            init_params(2, ipoi) = params(2, ipoi)

            ! Electron temperature [erg]
            ! Input Te is in eV, convert to erg
            params(3, ipoi) = 0.5 * (Te(ipoi) + Te(ipoi + 1)) * ev
            init_params(3, ipoi) = params(3, ipoi)

            ! Ion temperature [erg]
            ! Input Ti is in eV, convert to erg
            params(4, ipoi) = 0.5 * (Ti(ipoi) + Ti(ipoi + 1)) * ev
            init_params(4, ipoi) = params(4, ipoi)
        end do
    end subroutine

    !> @brief Write initial plasma profiles to HDF5 or binary file
    !> @details Writes profiles to output file for post-processing and analysis.
    !>          NOTE: Temperatures are converted back to [eV] for output by dividing by ev.
    !>          HDF5 output structure:
    !>          - /init_params/n:    density [cm^-3]
    !>          - /init_params/Vz:   toroidal rotation frequency [rad/s]
    !>          - /init_params/Te:   electron temperature [eV]
    !>          - /init_params/Ti:   ion temperature [eV]
    !>          - /init_params/qsaf: safety factor [1]
    !>          - /init_params/r:    radial coordinate [cm]
    !>          - /init_params/Er:   radial electric field [statV/cm]
    !>          - /init_params/Vth:  thermal velocity [cm/s]
    !> @author Markus Markl
    !> @date 05.10.2022
    subroutine write_initial_parameters
        use baseparam_mod, only: ev
        use control_mod, only: debug_mode, ihdf5IO
        use h5mod
        use wave_code_data, only: r, Vth, dPhi0

        if (debug_mode) write (*, *) "Debug: writing initial background profiles"
        if (ihdf5IO .eq. 1) then
            call h5_init()
            ! open hdf5 file
            call h5_open_rw(path2out, h5_id)
            call h5_obj_exists(h5_id, "/init_params/n", h5_exists_log)
            if (.not. h5_exists_log) then
                ! density [cm^-3]
                call h5_add_double_1(h5_id, "/init_params/n", params(1, :), &
                                     lbound(params(1, :)), ubound(params(1, :)))
                ! rotation frequency [rad/s]
                call h5_add_double_1(h5_id, "/init_params/Vz", params(2, :), &
                                     lbound(params(2, :)), ubound(params(2, :)))
                ! electron temperature [eV]
                call h5_add_double_1(h5_id, "/init_params/Te", params(3, :) / ev, &
                                     lbound(params(3, :)), ubound(params(3, :)))
                ! ion temperature [eV]
                call h5_add_double_1(h5_id, "/init_params/Ti", params(4, :) / ev, &
                                     lbound(params(4, :)), ubound(params(4, :)))
                ! safety factor [1]
                call h5_add_double_1(h5_id, "/init_params/qsaf", qsaf(:), &
                                     lbound(qsaf(:)), ubound(qsaf(:)))
                ! radial coordinate [cm]
                call h5_add_double_1(h5_id, "/init_params/r", r, lbound(r), ubound(r))
                ! radial electric field [statV/cm]
                ! Note: Er = -dPhi/dr, hence the minus sign
                call h5_add_double_1(h5_id, "/init_params/Er", -dPhi0, lbound(dPhi0), ubound(dPhi0))
                ! thermal velocity [cm/s]
                call h5_add_double_1(h5_id, "/init_params/Vth", Vth, lbound(Vth), ubound(Vth))
            else
                if (debug_mode) write (*, *) "Debug: they are already there -> skiping"
            end if

            call h5_close(h5_id)
            call h5_deinit()
            if (debug_mode) write (*, *) "Debug: finished writing initial background profiles"
            !stop ! for test purposes

        else
            ! Binary output: write params in internal units (temperatures in erg)
            open (123, form='unformatted', file='init_params.dat')
            write (123) params
            close (123)
        end if
    end subroutine

end module
