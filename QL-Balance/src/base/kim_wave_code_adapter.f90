module kim_wave_code_adapter_m
    !! Adapter module that wraps KIM library calls to populate
    !! wave_code_data module arrays, matching the contract expected
    !! by get_dql.f90 and diff_coeffs.f90.
    !!
    !! KIM modules are accessed via per-file INCLUDE_DIRECTORIES in
    !! CMake to avoid collision with the QL-Balance resonances_mod.
    !! All KIM symbols are renamed on import to avoid name clashes.

    use kim_solver_m, only: kim_solver_t, kim_results_t, kim_profiles_t, KIM_OK
    use species_m, only: kim_plasma => plasma     ! Path B reads KIM's file-loaded inputs
    use equilibrium_m, only: kim_B0 => B0         ! Path B background B0
    use config_m, only: kim_hdf5_output => hdf5_output
    use setup_m, only: kim_m_mode => m_mode, kim_n_mode => n_mode
    use grid_m, only: kim_xl_grid => xl_grid, &
                      kim_r_min => r_min, kim_r_plas => r_plas

    implicit none
    private

    public :: kim_initialize
    public :: kim_run_for_all_modes
    public :: kim_update_profiles
    public :: kim_get_wave_fields
    public :: kim_get_wave_vectors
    public :: kim_get_background_magnetic_fields
    public :: kim_get_collision_frequencies
    public :: kim_load_vacuum_fields
    public :: kim_check_domain_consistency
    public :: kim_get_current_densities
    public :: interp_complex_profile  ! exposed for testing

    !! Module-level KIM solver handle (reused across calls)
    type(kim_solver_t) :: kim_handle

    !! Actual KIM field grid boundary radius (xl_grid%xb(N)).
    !! This may differ from r_plas due to non-equidistant grid construction.
    !! Used for vacuum continuation cutoff and Br boundary condition.
    real(8) :: kim_r_boundary = 0.0d0

    !! Per-mode vacuum Br on QL-Balance grid (dim_r, dim_mn)
    complex(8), allocatable, public :: kim_vac_Br(:,:)

    !! Per-mode stored field results (nrad, dim_mn)
    !! Filled by kim_run_for_all_modes, read by kim_get_wave_fields
    complex(8), allocatable :: kim_Es_modes(:,:)
    complex(8), allocatable :: kim_Ep_modes(:,:)
    complex(8), allocatable :: kim_Er_modes(:,:)
    complex(8), allocatable :: kim_Et_modes(:,:)
    complex(8), allocatable :: kim_Ez_modes(:,:)
    complex(8), allocatable, public :: kim_Br_modes(:,:)

    !! Per-mode stored wave vectors (nrad, dim_mn)
    !! kp and ks depend on (m,n) via the equilibrium formulas.
    !! Filled by kim_run_for_all_modes, read by kim_get_wave_vectors
    real(8), allocatable :: kim_kp_modes(:,:)
    real(8), allocatable :: kim_ks_modes(:,:)

    !! Per-mode stored current densities (nrad, dim_mn)
    complex(8), allocatable :: kim_jpar_modes(:,:)
    complex(8), allocatable :: kim_jpar_e_modes(:,:)
    complex(8), allocatable :: kim_jpar_i_modes(:,:)

contains

    subroutine kim_initialize(nrad, r_grid)
        !! Initialize KIM backend: read config, profiles, grids.
        !! Populate wave_code_data module arrays with background
        !! quantities by running the electromagnetic init for the
        !! first mode to generate equilibrium and plasma backgrounds.
        !!
        !! Two paths:
        !!   kim_profiles_from_balance = .true.  (Path A, default):
        !!     Caller already loaded QL-Balance profiles into wave_code_data
        !!     and set mode arrays.  KIM reads config only; profiles injected
        !!     in memory via set_profiles_from_arrays.
        !!   kim_profiles_from_balance = .false. (Path B):
        !!     KIM reads its own files, adapter reads modes.in, extracts
        !!     everything — original behavior.
        use control_mod, only: kim_config_path, kim_profiles_from_balance
        use IO_collection_m, only: deinitialize_hdf5_output
        use wave_code_data, only: dim_mn, m_vals, n_vals, dim_r, &
            r => r, q => q, n => n, Te => Te, Ti => Ti, &
            Vth => Vth, Vz => Vz, dPhi0 => dPhi0, &
            kp => kp, ks => ks, om_E => om_E, &
            B0 => B0, nue => nue, nui => nui, &
            B0t => B0t, B0z => B0z, &
            flre_path

        implicit none

        integer, intent(in) :: nrad
        real(8), intent(in) :: r_grid(nrad)

        type(kim_profiles_t) :: prof
        integer :: ierr, kim_npts
        real(8), allocatable :: kim_r(:), work_old(:)

        ! The first mode sets the equilibrium that init() builds; the first
        ! solve in kim_run_for_all_modes reuses it (no recompute on solve 1).
        kim_m_mode = m_vals(1)
        kim_n_mode = n_vals(1)

        ! Re-init safe: clear any equilibrium/field state from a prior run.
        call kim_handle%finalize()

        if (kim_profiles_from_balance) then
            ! -----------------------------------------------------------
            ! Path A: inject QL-Balance profiles in memory (dPhi0 = -Er).
            ! The derived background is read from results() after the first
            ! solve (in kim_run_for_all_modes), not here.
            ! -----------------------------------------------------------
            if (r(1) > kim_r_min + 1.0d-6 .or. r(dim_r) < kim_r_plas - 1.0d-6) then
                write(*,*) 'WARNING: QL-Balance grid [', r(1), ',', r(dim_r), &
                           '] does not fully cover KIM range [', kim_r_min, ',', kim_r_plas, ']'
            end if

            allocate(prof%r(dim_r), prof%n(dim_r), prof%Te(dim_r), &
                     prof%Ti(dim_r), prof%q(dim_r), prof%Er(dim_r))
            prof%r = r; prof%n = n; prof%Te = Te
            prof%Ti = Ti; prof%q = q; prof%Er = -dPhi0

            call kim_handle%init(trim(kim_config_path), run_type='electromagnetic', &
                                 profiles=prof, stat=ierr)
        else
            ! -----------------------------------------------------------
            ! Path B (legacy): KIM reads its own profile files.
            ! -----------------------------------------------------------
            call read_antenna_modes(flre_path)
            call allocate_wave_code_data(nrad, r_grid)

            call kim_handle%init(trim(kim_config_path), run_type='electromagnetic', &
                                 stat=ierr)
        end if

        if (ierr /= KIM_OK) then
            write(*,*) 'ERROR: KIM init failed with status ', ierr
            stop 1
        end if

        ! Disable KIM HDF5 output for QL-Balance integration.
        ! We only need in-memory fields; HDF5 writes conflict on re-init.
        if (kim_hdf5_output) call deinitialize_hdf5_output()
        kim_hdf5_output = .false.

        ! Store actual KIM grid boundary (may differ from r_plas).
        kim_r_boundary = kim_xl_grid%xb(kim_xl_grid%npts_b)

        if (.not. kim_profiles_from_balance) then
            ! Path B: read all background quantities back from KIM's globals
            ! (populated by init()); per-mode fields come from the solve loop.
            kim_npts = kim_plasma%grid_size
            allocate(kim_r(kim_npts), work_old(kim_npts))
            kim_r = kim_plasma%r_grid

            work_old = kim_plasma%q
            call interp_profile(kim_npts, kim_r, work_old, nrad, r, q)

            work_old = kim_plasma%kp
            call interp_profile(kim_npts, kim_r, work_old, nrad, r, kp)

            work_old = kim_plasma%ks
            call interp_profile(kim_npts, kim_r, work_old, nrad, r, ks)

            work_old = kim_plasma%om_E
            call interp_profile(kim_npts, kim_r, work_old, nrad, r, om_E)

            ! dPhi0 = -Er
            work_old = -kim_plasma%Er
            call interp_profile(kim_npts, kim_r, work_old, nrad, r, dPhi0)

            work_old = kim_plasma%spec(0)%n
            call interp_profile(kim_npts, kim_r, work_old, nrad, r, n)

            work_old = kim_plasma%spec(0)%T
            call interp_profile(kim_npts, kim_r, work_old, nrad, r, Te)

            work_old = kim_plasma%spec(1)%T
            call interp_profile(kim_npts, kim_r, work_old, nrad, r, Ti)

            work_old = kim_plasma%spec(0)%nu
            call interp_profile(kim_npts, kim_r, work_old, nrad, r, nue)

            work_old = kim_plasma%spec(1)%nu
            call interp_profile(kim_npts, kim_r, work_old, nrad, r, nui)

            work_old = kim_B0
            call interp_profile(kim_npts, kim_r, work_old, nrad, r, B0)

            call interp_B0_components(kim_npts, kim_r, nrad, r)

            Vth = 0.0d0
            Vz = 0.0d0

            deallocate(kim_r, work_old)
        end if

        ! Vacuum Br placeholder (NaN until kim_load_vacuum_fields fills it)
        block
            use ieee_arithmetic, only: ieee_value, ieee_quiet_nan
            real(8) :: nan_val
            nan_val = ieee_value(1.0d0, ieee_quiet_nan)
            if (allocated(kim_vac_Br)) deallocate(kim_vac_Br)
            allocate(kim_vac_Br(nrad, dim_mn))
            kim_vac_Br = cmplx(nan_val, nan_val, kind=8)
        end block

        write(*, *) "KIM adapter: initialization complete"
        write(*, *) "  Balance grid points: ", nrad
        write(*, *) "  Number of modes:     ", dim_mn
        if (kim_profiles_from_balance) then
            write(*, *) "  Profile source:      QL-Balance (in-memory)"
        else
            write(*, *) "  Profile source:      KIM files"
        end if

    end subroutine kim_initialize

    ! ---------------------------------------------------------------
    ! Helper: interpolate B0 toroidal and poloidal components
    ! ---------------------------------------------------------------
    subroutine interp_B0_components(kim_npts, kim_r, nbal, bal_r)
        !! Interpolate B0z (toroidal) and B0th (poloidal) from
        !! equilibrium_m onto the QL-Balance grid.
        use equilibrium_m, only: B0z_equil => B0z, B0th_equil => B0th
        use wave_code_data, only: B0t, B0z

        implicit none

        integer, intent(in) :: kim_npts, nbal
        real(8), intent(in) :: kim_r(kim_npts), bal_r(nbal)
        real(8), allocatable :: work(:)

        allocate(work(kim_npts))

        ! B0z (toroidal component) -> wave_code_data::B0z
        work = B0z_equil
        call interp_profile(kim_npts, kim_r, work, nbal, bal_r, B0z)

        ! B0th (poloidal component) -> wave_code_data::B0t
        work = B0th_equil
        call interp_profile(kim_npts, kim_r, work, nbal, bal_r, B0t)

        deallocate(work)

    end subroutine interp_B0_components

    ! ---------------------------------------------------------------
    ! Helper: interpolate up to r_plas, clamp to the boundary beyond
    ! ---------------------------------------------------------------
    subroutine clamp_to_balance(n_old, r_old, f_old, nrad, nrad_inside, bal_r, f_new)
        !! Interpolate f_old onto bal_r up to nrad_inside points, then clamp
        !! the remaining points to the last in-domain value (the KIM solution
        !! is valid only inside the plasma domain).
        implicit none

        integer, intent(in) :: n_old, nrad, nrad_inside
        real(8), intent(in) :: r_old(n_old), f_old(n_old), bal_r(nrad)
        real(8), intent(out) :: f_new(nrad)

        call interp_profile(n_old, r_old, f_old, nrad_inside, &
            bal_r(1:nrad_inside), f_new(1:nrad_inside))
        f_new(nrad_inside+1:nrad) = f_new(nrad_inside)

    end subroutine clamp_to_balance

    subroutine kim_run_for_all_modes()
        !! Run KIM electrostatic solver for each (m,n) mode and
        !! store the resulting fields in per-mode arrays.
        !!
        !! After this call, kim_Es_modes(:,i_mn) etc. hold the
        !! field solutions interpolated onto the QL-Balance radial
        !! grid.  kim_get_wave_fields(i_mn) copies from these
        !! arrays into the wave_code_data module scalars.
        use wave_code_data, only: dim_mn, m_vals, n_vals, &
            dim_r, bal_r => r, &
            wcd_kp => kp, wcd_ks => ks, wcd_om_E => om_E, &
            wcd_B0 => B0, wcd_nue => nue, wcd_nui => nui, &
            wcd_B0t => B0t, wcd_B0z => B0z, wcd_Vth => Vth, wcd_Vz => Vz
        use control_mod, only: kim_profiles_from_balance

        implicit none

        type(kim_results_t) :: res
        integer :: i_mn, ierr, kim_npts, kim_plasma_npts, nrad_inside, i
        real(8), allocatable :: kim_r(:), kim_plasma_r(:)

        ! -------------------------------------------------------
        ! 1. (Re-)allocate per-mode storage
        ! -------------------------------------------------------
        if (allocated(kim_Es_modes)) deallocate(kim_Es_modes)
        if (allocated(kim_Ep_modes)) deallocate(kim_Ep_modes)
        if (allocated(kim_Er_modes)) deallocate(kim_Er_modes)
        if (allocated(kim_Et_modes)) deallocate(kim_Et_modes)
        if (allocated(kim_Ez_modes)) deallocate(kim_Ez_modes)
        if (allocated(kim_Br_modes)) deallocate(kim_Br_modes)
        if (allocated(kim_kp_modes)) deallocate(kim_kp_modes)
        if (allocated(kim_ks_modes)) deallocate(kim_ks_modes)
        if (allocated(kim_jpar_modes)) deallocate(kim_jpar_modes)
        if (allocated(kim_jpar_e_modes)) deallocate(kim_jpar_e_modes)
        if (allocated(kim_jpar_i_modes)) deallocate(kim_jpar_i_modes)

        allocate(kim_Es_modes(dim_r, dim_mn))
        allocate(kim_Ep_modes(dim_r, dim_mn))
        allocate(kim_Er_modes(dim_r, dim_mn))
        allocate(kim_Et_modes(dim_r, dim_mn))
        allocate(kim_Ez_modes(dim_r, dim_mn))
        allocate(kim_Br_modes(dim_r, dim_mn))
        allocate(kim_kp_modes(dim_r, dim_mn))
        allocate(kim_ks_modes(dim_r, dim_mn))
        allocate(kim_jpar_modes(dim_r, dim_mn))
        allocate(kim_jpar_e_modes(dim_r, dim_mn))
        allocate(kim_jpar_i_modes(dim_r, dim_mn))

        kim_Es_modes = (0.0d0, 0.0d0)
        kim_Ep_modes = (0.0d0, 0.0d0)
        kim_Er_modes = (0.0d0, 0.0d0)
        kim_Et_modes = (0.0d0, 0.0d0)
        kim_Ez_modes = (0.0d0, 0.0d0)
        kim_Br_modes = (0.0d0, 0.0d0)
        kim_kp_modes = 0.0d0
        kim_ks_modes = 0.0d0
        kim_jpar_modes = (0.0d0, 0.0d0)
        kim_jpar_e_modes = (0.0d0, 0.0d0)
        kim_jpar_i_modes = (0.0d0, 0.0d0)

        ! -------------------------------------------------------
        ! 2. Loop over modes: solve through the seam and store
        ! -------------------------------------------------------
        do i_mn = 1, dim_mn

            ! The seam owns mode setup, the per-mode equilibrium recompute
            ! (modes 2+), the field reset, and the run.
            call kim_handle%solve(m_vals(i_mn), n_vals(i_mn), stat=ierr)
            if (ierr /= KIM_OK) then
                write(*,*) 'ERROR: KIM solve failed for mode ', i_mn, &
                           ' status ', ierr
                stop 1
            end if
            res = kim_handle%results()

            ! Interpolate KIM fields (on res%r_field) onto the QL-Balance grid.
            kim_npts = size(res%r_field)
            allocate(kim_r(kim_npts))
            kim_r = res%r_field

            if (i_mn == 1) then
                write(*,*) '  KIM field grid: r_min=', kim_r(1), &
                           ' r_max=', kim_r(kim_npts), ' npts=', kim_npts
                write(*,*) '  KIM r_plas =', kim_r_plas
                write(*,*) '  KIM |Br| at last grid pt =', abs(res%Br(kim_npts))
            end if

            ! Es (perpendicular E field in rsp coordinates)
            call interp_complex_profile(kim_npts, kim_r, res%Es, &
                dim_r, bal_r, kim_Es_modes(:, i_mn))

            ! Ep (parallel E field in rsp coordinates)
            call interp_complex_profile(kim_npts, kim_r, res%Ep, &
                dim_r, bal_r, kim_Ep_modes(:, i_mn))

            ! Er (radial E field, cylindrical)
            call interp_complex_profile(kim_npts, kim_r, res%Er, &
                dim_r, bal_r, kim_Er_modes(:, i_mn))

            ! Etheta -> Et (poloidal E field, cylindrical)
            call interp_complex_profile(kim_npts, kim_r, res%Etheta, &
                dim_r, bal_r, kim_Et_modes(:, i_mn))

            ! Ez (axial E field, cylindrical)
            call interp_complex_profile(kim_npts, kim_r, res%Ez, &
                dim_r, bal_r, kim_Ez_modes(:, i_mn))

            ! Br (radial magnetic field perturbation)
            call interp_complex_profile(kim_npts, kim_r, res%Br, &
                dim_r, bal_r, kim_Br_modes(:, i_mn))

            ! jpar (parallel current density: total, electron, ion)
            if (allocated(res%jpar)) then
                call interp_complex_profile(kim_npts, kim_r, res%jpar, &
                    dim_r, bal_r, kim_jpar_modes(:, i_mn))
            end if
            if (allocated(res%jpar_e)) then
                call interp_complex_profile(kim_npts, kim_r, res%jpar_e, &
                    dim_r, bal_r, kim_jpar_e_modes(:, i_mn))
            end if
            if (allocated(res%jpar_i)) then
                call interp_complex_profile(kim_npts, kim_r, res%jpar_i, &
                    dim_r, bal_r, kim_jpar_i_modes(:, i_mn))
            end if

            deallocate(kim_r)

            ! Apply vacuum continuation beyond r_plas:
            ! KIM fields are only valid inside the plasma domain.
            ! Beyond r_plas, use the vacuum Br from KiLCA and zero E fields.
            call apply_vacuum_continuation(i_mn, dim_r, bal_r)

            ! kp and ks (mode-dependent wave vectors on plasma grid)
            kim_plasma_npts = size(res%r_plasma)
            allocate(kim_plasma_r(kim_plasma_npts))
            kim_plasma_r = res%r_plasma

            call interp_profile(kim_plasma_npts, kim_plasma_r, res%kp, &
                dim_r, bal_r, kim_kp_modes(:, i_mn))
            call interp_profile(kim_plasma_npts, kim_plasma_r, res%ks, &
                dim_r, bal_r, kim_ks_modes(:, i_mn))

            ! Path A: derived background from the first solve, interpolated up
            ! to r_plas and clamped to the boundary value beyond.  n/Te/Ti/q/
            ! dPhi0 are NOT touched -- they came from QL-Balance.
            if (i_mn == 1 .and. kim_profiles_from_balance) then
                nrad_inside = dim_r
                do i = 1, dim_r
                    if (bal_r(i) > kim_r_plas) then
                        nrad_inside = i - 1
                        exit
                    end if
                end do

                call clamp_to_balance(kim_plasma_npts, kim_plasma_r, res%B0, &
                    dim_r, nrad_inside, bal_r, wcd_B0)
                call clamp_to_balance(kim_plasma_npts, kim_plasma_r, res%B0z, &
                    dim_r, nrad_inside, bal_r, wcd_B0z)
                call clamp_to_balance(kim_plasma_npts, kim_plasma_r, res%B0th, &
                    dim_r, nrad_inside, bal_r, wcd_B0t)
                call clamp_to_balance(kim_plasma_npts, kim_plasma_r, res%kp, &
                    dim_r, nrad_inside, bal_r, wcd_kp)
                call clamp_to_balance(kim_plasma_npts, kim_plasma_r, res%ks, &
                    dim_r, nrad_inside, bal_r, wcd_ks)
                call clamp_to_balance(kim_plasma_npts, kim_plasma_r, res%om_E, &
                    dim_r, nrad_inside, bal_r, wcd_om_E)
                call clamp_to_balance(kim_plasma_npts, kim_plasma_r, res%nu_e, &
                    dim_r, nrad_inside, bal_r, wcd_nue)
                call clamp_to_balance(kim_plasma_npts, kim_plasma_r, res%nu_i, &
                    dim_r, nrad_inside, bal_r, wcd_nui)

                wcd_Vth = 0.0d0
                wcd_Vz = 0.0d0
            end if

            deallocate(kim_plasma_r)

            write(*, '(A,I3,A,I4,A,I4,A)') &
                "  KIM adapter: solved mode ", i_mn, &
                "  (m=", m_vals(i_mn), ", n=", n_vals(i_mn), ")"

        end do

        ! -------------------------------------------------------
        ! 3. Path A: refresh derived quantities from the final
        !    equilibrium state (non-clamped, full grid)
        ! -------------------------------------------------------
        if (kim_profiles_from_balance) then
            kim_plasma_npts = size(res%r_plasma)
            allocate(kim_plasma_r(kim_plasma_npts))
            kim_plasma_r = res%r_plasma

            call interp_profile(kim_plasma_npts, kim_plasma_r, res%om_E, &
                dim_r, bal_r, wcd_om_E)
            call interp_profile(kim_plasma_npts, kim_plasma_r, res%nu_e, &
                dim_r, bal_r, wcd_nue)
            call interp_profile(kim_plasma_npts, kim_plasma_r, res%nu_i, &
                dim_r, bal_r, wcd_nui)
            call interp_profile(kim_plasma_npts, kim_plasma_r, res%B0, &
                dim_r, bal_r, wcd_B0)

            deallocate(kim_plasma_r)
        end if

        write(*, *) "KIM adapter: all modes solved"

    end subroutine kim_run_for_all_modes

    subroutine kim_update_profiles()
        !! Transfer QL-Balance time-evolved profiles into the
        !! wave_code_data in-memory arrays, mirroring what
        !! update_background_files does for the KiLCA path.
        !!
        !! QL-Balance params_b layout (at boundary grid points):
        !!   params_b(1,:) = density      [1/cm^3]
        !!   params_b(2,:) = Vphi         [rad/s]
        !!   params_b(3,:) = Te           [erg]
        !!   params_b(4,:) = Ti           [erg]
        !!
        !! wave_code_data stores Te, Ti in [eV], Vz in [cm/s].
        !!
        !! Note: dPhi0 is only updated from Ercov during time evolution.
        !! For SingleStep runs, the input Er profile is preserved.
        use wave_code_data, only: dim_r, &
            wcd_n => n, wcd_Te => Te, wcd_Ti => Ti, &
            wcd_Vz => Vz, wcd_dPhi0 => dPhi0
        use plasma_parameters, only: params_b
        use baseparam_mod, only: ev, rtor
        use grid_mod, only: Ercov
        use control_mod, only: type_of_run

        implicit none

        integer :: k

        do k = 1, dim_r
            wcd_n(k)     = params_b(1, k)
            wcd_Te(k)    = params_b(3, k) / ev   ! erg -> eV
            wcd_Ti(k)    = params_b(4, k) / ev   ! erg -> eV
            wcd_Vz(k)    = params_b(2, k) * rtor  ! rad/s -> cm/s
        end do

        ! Only update dPhi0 from Ercov during time evolution.
        ! For SingleStep, preserve the input Er profile (from Er.dat).
        if (trim(type_of_run) /= 'SingleStep') then
            do k = 1, dim_r
                wcd_dPhi0(k) = -Ercov(k)
            end do
        end if

    end subroutine kim_update_profiles

    subroutine kim_get_wave_fields(i_mn)
        !! Copy per-mode stored fields from kim_*_modes arrays
        !! into the wave_code_data module scalars for mode i_mn.
        !! Called inside the mode loop in get_dql.
        use wave_code_data, only: Es, Br, Er, Ep, Et, Ez, &
            Bs, Bp, Bt, Bz

        implicit none

        integer, intent(in) :: i_mn

        Es = kim_Es_modes(:, i_mn)
        Br = kim_Br_modes(:, i_mn)
        Er = kim_Er_modes(:, i_mn)
        Ep = kim_Ep_modes(:, i_mn)
        Et = kim_Et_modes(:, i_mn)
        Ez = kim_Ez_modes(:, i_mn)

        ! KIM does not compute magnetic field perturbation components
        ! other than Br. Set remaining B components to zero.
        Bs = (0.0d0, 0.0d0)
        Bp = (0.0d0, 0.0d0)
        Bt = (0.0d0, 0.0d0)
        Bz = (0.0d0, 0.0d0)

    end subroutine kim_get_wave_fields

    subroutine kim_get_wave_vectors(i_mn)
        !! Copy per-mode stored wave vectors kp and ks from
        !! kim_kp_modes / kim_ks_modes into wave_code_data for
        !! mode i_mn.  kp and ks depend on (m,n) through
        !! the equilibrium geometry and are recomputed by KIM
        !! for each mode in kim_run_for_all_modes.
        use wave_code_data, only: kp, ks

        implicit none

        integer, intent(in) :: i_mn

        kp = kim_kp_modes(:, i_mn)
        ks = kim_ks_modes(:, i_mn)

    end subroutine kim_get_wave_vectors

    subroutine kim_get_background_magnetic_fields()
        !! No-op: B0, B0t, B0z were already populated during
        !! kim_initialize and stored directly in wave_code_data.
        !! This subroutine exists to satisfy the adapter interface
        !! contract expected by get_dql.

        implicit none

        ! Nothing to do -- B0, B0t, B0z are already set.

    end subroutine kim_get_background_magnetic_fields

    subroutine kim_get_collision_frequencies()
        !! No-op: nue and nui were already populated during
        !! kim_initialize and stored directly in wave_code_data.
        !! This subroutine exists to satisfy the adapter interface
        !! contract expected by get_dql.

        implicit none

        ! Nothing to do -- nue, nui are already set.

    end subroutine kim_get_collision_frequencies

    ! ---------------------------------------------------------------
    ! Helper: interpolate a complex profile onto a new grid
    ! ---------------------------------------------------------------
    subroutine interp_complex_profile(n_old, r_old, z_old, &
                                      n_new, r_new, z_new)
        !! Interpolate a complex(8) profile from one radial grid
        !! to another.  Real and imaginary parts are interpolated
        !! independently via the existing interp_profile routine.
        implicit none

        integer, intent(in)  :: n_old, n_new
        real(8), intent(in)  :: r_old(n_old), r_new(n_new)
        complex(8), intent(in)  :: z_old(n_old)
        complex(8), intent(out) :: z_new(n_new)

        real(8) :: re_old(n_old), im_old(n_old)
        real(8) :: re_new(n_new), im_new(n_new)

        re_old = real(z_old)
        im_old = aimag(z_old)

        call interp_profile(n_old, r_old, re_old, n_new, r_new, re_new)
        call interp_profile(n_old, r_old, im_old, n_new, r_new, im_new)

        z_new = cmplx(re_new, im_new, kind=8)

    end subroutine interp_complex_profile

    ! ---------------------------------------------------------------
    ! Vacuum field loading and domain consistency
    ! ---------------------------------------------------------------

    subroutine kim_load_vacuum_fields()
        !! Extract vacuum Br from KiLCA vacuum solution (vac_cd_ptr)
        !! onto the balance grid and store in kim_vac_Br for each mode.
        use wave_code_data, only: dim_r, r, dim_mn, m_vals, n_vals, &
            vac_cd_ptr, Br, Bz

        implicit none

        integer :: k

        ! Extract vacuum Br for each mode.
        ! get_wave_fields_from_wave_code fills Br; we use Bz as dummy
        ! for all other field components.
        do k = 1, dim_mn
            call get_wave_fields_from_wave_code(vac_cd_ptr(k), dim_r, r, &
                m_vals(k), n_vals(k), Bz, Bz, Bz, Bz, Bz, Br, Bz, Bz, Bz, Bz)
            kim_vac_Br(:, k) = Br
        end do

        write(*,*) 'KIM adapter: loaded vacuum Br from KiLCA for ', dim_mn, ' modes'
        do k = 1, dim_mn
            write(*,*) '  mode ', k, ': |Br_vac| at r_max = ', &
                abs(kim_vac_Br(dim_r, k)), ' at r=', r(dim_r)
        end do

        call kim_check_domain_consistency()

    end subroutine kim_load_vacuum_fields

    subroutine kim_check_domain_consistency()
        !! Set KIM Br boundary condition from KiLCA vacuum solution
        !! at the actual KIM grid boundary (kim_r_boundary = xl_grid%xb(N)),
        !! which may differ from r_plas due to non-equidistant grid construction.
        !! Uses linear interpolation on the balance grid to get vacuum Br
        !! at exactly kim_r_boundary.
        use wave_code_data, only: dim_r, r, vac_cd_ptr
        use setup_m, only: Br_boundary_re, Br_boundary_im

        implicit none

        integer :: i, i_lo, i_hi
        real(8) :: w, r_bc
        complex(8) :: Br_vac_at_bc

        ! Use actual KIM grid boundary for the boundary condition
        r_bc = kim_r_boundary

        ! 1. Report domain info
        write(*,*) 'KIM adapter: domain check'
        write(*,*) '  KIM r_plas       = ', kim_r_plas
        write(*,*) '  KIM r_boundary   = ', r_bc, ' (actual grid boundary)'
        write(*,*) '  Balance grid max = ', r(dim_r)

        if (r_bc < r(dim_r)) then
            write(*,*) '  Vacuum continuation beyond r_boundary to r=', r(dim_r), ' cm'
        end if

        ! 2. Check vacuum solution is available
        if (vac_cd_ptr(1) == 0) then
            write(*,*) 'ERROR: KIM mode requires vacuum solution but vac_cd_ptr is not loaded.'
            write(*,*) '       Set vac_path in balance_conf.nml and run KiLCA vacuum solver.'
            stop 1
        end if

        ! 3. Set KIM Br boundary from vacuum Br at exactly r_bc
        ! Find bracketing indices for linear interpolation
        i_lo = 1
        do i = 1, dim_r - 1
            if (r(i) <= r_bc .and. r(i+1) >= r_bc) then
                i_lo = i
                exit
            end if
        end do
        i_hi = i_lo + 1

        ! Linear interpolation weight
        if (abs(r(i_hi) - r(i_lo)) > 1.0d-30) then
            w = (r_bc - r(i_lo)) / (r(i_hi) - r(i_lo))
        else
            w = 0.0d0
        end if

        Br_vac_at_bc = (1.0d0 - w) * kim_vac_Br(i_lo, 1) &
                     + w * kim_vac_Br(i_hi, 1)

        write(*,*) '  Vacuum Br at r_boundary (linear interp, mode 1):'
        write(*,*) '    r_lo=', r(i_lo), ' r_hi=', r(i_hi), ' w=', w
        write(*,*) '    |Br_vac(r_lo)| = ', abs(kim_vac_Br(i_lo, 1))
        write(*,*) '    |Br_vac(r_hi)| = ', abs(kim_vac_Br(i_hi, 1))
        write(*,*) '    |Br_vac(r_bc)| = ', abs(Br_vac_at_bc)
        write(*,*) '    Re(Br)    = ', real(Br_vac_at_bc)
        write(*,*) '    Im(Br)    = ', aimag(Br_vac_at_bc)

        ! Override KIM boundary condition with vacuum value
        write(*,*) '  Old Br_boundary: Re=', Br_boundary_re, ' Im=', Br_boundary_im
        Br_boundary_re = real(Br_vac_at_bc)
        Br_boundary_im = aimag(Br_vac_at_bc)
        write(*,*) '  New Br_boundary: Re=', Br_boundary_re, ' Im=', Br_boundary_im

    end subroutine kim_check_domain_consistency

    subroutine kim_get_current_densities(i_mn)
        !! Copy per-mode stored parallel current densities from
        !! kim_jpar_e_modes / kim_jpar_i_modes into the wave_code_data
        !! J arrays for mode i_mn. KIM computes species-resolved
        !! currents via separate electron and ion kernels.
        use wave_code_data, only: Jpi, Jpe, Jri, Jre, Jsi, Jse

        implicit none

        integer, intent(in) :: i_mn

        Jpe = kim_jpar_e_modes(:, i_mn)
        Jpi = kim_jpar_i_modes(:, i_mn)
        ! KIM does not compute radial/perpendicular current components
        Jri = (0.0d0, 0.0d0)
        Jre = (0.0d0, 0.0d0)
        Jsi = (0.0d0, 0.0d0)
        Jse = (0.0d0, 0.0d0)

    end subroutine kim_get_current_densities

    subroutine apply_vacuum_continuation(i_mn, nrad, bal_r)
        !! For grid points beyond kim_r_boundary (actual KIM grid extent),
        !! replace KIM fields with vacuum values: E fields = 0,
        !! Br = kim_vac_Br, wave vectors = 0.
        implicit none

        integer, intent(in) :: i_mn, nrad
        real(8), intent(in) :: bal_r(nrad)
        integer :: i, i_last_kim

        ! Find last balance grid point inside KIM domain
        ! Use kim_r_boundary (actual grid extent), not kim_r_plas
        i_last_kim = nrad
        do i = 1, nrad
            if (bal_r(i) > kim_r_boundary) then
                i_last_kim = i - 1
                exit
            end if
        end do

        ! Diagnostic: print transition values at the stitching boundary
        if (i_mn == 1 .and. i_last_kim < nrad) then
            write(*,*) '  Vacuum continuation transition (at r_boundary=', &
                       kim_r_boundary, '):'
            write(*,*) '    Last KIM point:  r=', bal_r(i_last_kim), &
                       ' |Br_KIM|=', abs(kim_Br_modes(i_last_kim, i_mn))
            write(*,*) '    First vac point: r=', bal_r(i_last_kim+1), &
                       ' |Br_vac|=', abs(kim_vac_Br(i_last_kim+1, i_mn))
            write(*,*) '    Vac Br at last KIM r: |Br_vac|=', &
                       abs(kim_vac_Br(i_last_kim, i_mn))
            write(*,*) '    Value jump: |Br_KIM - Br_vac| = ', &
                       abs(kim_Br_modes(i_last_kim, i_mn) - kim_vac_Br(i_last_kim, i_mn))
        end if

        do i = i_last_kim + 1, nrad
            ! No plasma response in vacuum: E fields = 0
            kim_Es_modes(i, i_mn) = (0.0d0, 0.0d0)
            kim_Ep_modes(i, i_mn) = (0.0d0, 0.0d0)
            kim_Er_modes(i, i_mn) = (0.0d0, 0.0d0)
            kim_Et_modes(i, i_mn) = (0.0d0, 0.0d0)
            kim_Ez_modes(i, i_mn) = (0.0d0, 0.0d0)
            ! Total Br = vacuum Br beyond plasma
            kim_Br_modes(i, i_mn) = kim_vac_Br(i, i_mn)
            ! Wave vectors not physical in vacuum
            kim_kp_modes(i, i_mn) = 0.0d0
            kim_ks_modes(i, i_mn) = 0.0d0
        end do

    end subroutine apply_vacuum_continuation

end module kim_wave_code_adapter_m
