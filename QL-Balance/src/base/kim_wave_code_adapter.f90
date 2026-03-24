module kim_wave_code_adapter_m
    !! Adapter module that wraps KIM library calls to populate
    !! wave_code_data module arrays, matching the contract expected
    !! by get_dql.f90 and diff_coeffs.f90.
    !!
    !! KIM modules are accessed via per-file INCLUDE_DIRECTORIES in
    !! CMake to avoid collision with the QL-Balance resonances_mod.
    !! All KIM symbols are renamed on import to avoid name clashes.

    use species_m, only: kim_plasma => plasma, &
                        set_profiles_from_arrays, deallocate_plasma_derived, &
                        set_plasma_quantities
    use equilibrium_m, only: kim_B0 => B0, &
        eq_B0z => B0z, eq_B0th => B0th, &
        eq_hz => hz, eq_hth => hth, eq_equil_grid => equil_grid, &
        eq_u => u, eq_dpress_prof => dpress_prof, &
        calculate_equil, interpolate_equil
    use config_m, only: kim_type_of_run => type_of_run, &
                        kim_nml_config_path => nml_config_path, &
                        kim_profiles_in_memory => profiles_in_memory, &
                        kim_hdf5_output => hdf5_output
    use setup_m, only: kim_m_mode => m_mode, kim_n_mode => n_mode
    use grid_m, only: kim_xl_grid => xl_grid, kim_rg_grid => rg_grid, &
                      kim_r_min => r_min, kim_r_plas => r_plas
    use fields_m, only: EBdat
    use kim_base_m, only: kim_t
    use kim_mod_m, only: from_kim_factory_get_kim

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

    !! Module-level KIM solver instance (reused across calls)
    class(kim_t), allocatable :: kim_instance

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

        integer :: kim_npts, nrad_inside, i
        real(8), allocatable :: kim_r(:)
        real(8), allocatable :: work_old(:), work_new(:)

        ! -----------------------------------------------------------
        ! 1. Set KIM config path and initialize KIM backend
        ! -----------------------------------------------------------
        kim_nml_config_path = trim(kim_config_path)

        if (kim_profiles_from_balance) then
            ! Path A: profiles come from QL-Balance in memory
            kim_profiles_in_memory = .true.
        end if

        call kim_init()

        ! Disable KIM HDF5 output for QL-Balance integration.
        ! We only need in-memory fields; HDF5 writes conflict on re-init.
        if (kim_hdf5_output) then
            call deinitialize_hdf5_output()
        end if
        kim_hdf5_output = .false.

        if (kim_profiles_from_balance) then
            ! -----------------------------------------------------------
            ! Path A: inject QL-Balance profiles into KIM's plasma struct
            ! -----------------------------------------------------------
            ! Grid coverage check
            if (r(1) > kim_r_min + 1.0d-6 .or. r(dim_r) < kim_r_plas - 1.0d-6) then
                write(*,*) 'WARNING: QL-Balance grid [', r(1), ',', r(dim_r), &
                           '] does not fully cover KIM range [', kim_r_min, ',', kim_r_plas, ']'
            end if

            ! Inject profiles: note dPhi0 = -Er, so pass -dPhi0 as Er
            call deallocate_equilibrium_arrays()
            call set_profiles_from_arrays(r, n, Te, Ti, q, -dPhi0, dim_r)

            ! Initialize KIM grids and equilibrium for the first mode
            kim_type_of_run = 'electromagnetic'
            kim_m_mode = m_vals(1)
            kim_n_mode = n_vals(1)

            call from_kim_factory_get_kim(kim_type_of_run, kim_instance)
            call kim_instance%init()

            ! Store actual KIM grid boundary (may differ from r_plas)
            kim_r_boundary = kim_xl_grid%xb(kim_xl_grid%npts_b)
            write(*,*) '  KIM xl_grid boundary: r_plas=', kim_r_plas, &
                       ' xl_grid(N)=', kim_r_boundary

            ! Extract derived quantities from KIM back to wave_code_data
            ! Do NOT overwrite n, Te, Ti, q, dPhi0 (they came from QL-Balance)
            ! Interpolate only up to r_plas; clamp to boundary value beyond.
            kim_npts = kim_plasma%grid_size
            allocate(kim_r(kim_npts), work_old(kim_npts))

            kim_r = kim_plasma%r_grid

            ! Find last balance grid index within KIM domain
            nrad_inside = nrad
            do i = 1, nrad
                if (r(i) > kim_r_plas) then
                    nrad_inside = i - 1
                    exit
                end if
            end do

            ! Background magnetic field B0
            work_old = kim_B0
            call interp_profile(kim_npts, kim_r, work_old, nrad_inside, r(1:nrad_inside), B0(1:nrad_inside))
            B0(nrad_inside+1:nrad) = B0(nrad_inside)

            ! B0 toroidal and poloidal components
            call interp_B0_components_clamped(kim_npts, kim_r, nrad, r, nrad_inside)

            ! Wave vectors kp, ks
            work_old = kim_plasma%kp
            call interp_profile(kim_npts, kim_r, work_old, nrad_inside, r(1:nrad_inside), kp(1:nrad_inside))
            kp(nrad_inside+1:nrad) = kp(nrad_inside)

            work_old = kim_plasma%ks
            call interp_profile(kim_npts, kim_r, work_old, nrad_inside, r(1:nrad_inside), ks(1:nrad_inside))
            ks(nrad_inside+1:nrad) = ks(nrad_inside)

            ! ExB drift frequency
            work_old = kim_plasma%om_E
            call interp_profile(kim_npts, kim_r, work_old, nrad_inside, r(1:nrad_inside), om_E(1:nrad_inside))
            om_E(nrad_inside+1:nrad) = om_E(nrad_inside)

            ! Collision frequencies
            work_old = kim_plasma%spec(0)%nu
            call interp_profile(kim_npts, kim_r, work_old, nrad_inside, r(1:nrad_inside), nue(1:nrad_inside))
            nue(nrad_inside+1:nrad) = nue(nrad_inside)

            work_old = kim_plasma%spec(1)%nu
            call interp_profile(kim_npts, kim_r, work_old, nrad_inside, r(1:nrad_inside), nui(1:nrad_inside))
            nui(nrad_inside+1:nrad) = nui(nrad_inside)

            ! Vth and Vz: KIM does not store these
            Vth = 0.0d0
            Vz = 0.0d0

            deallocate(kim_r, work_old)

        else
            ! -----------------------------------------------------------
            ! Path B: KIM reads its own files (original behavior)
            ! -----------------------------------------------------------
            call read_antenna_modes(flre_path)
            call allocate_wave_code_data(nrad, r_grid)

            ! Initialize KIM grids and equilibrium for the first mode
            kim_type_of_run = 'electromagnetic'
            kim_m_mode = m_vals(1)
            kim_n_mode = n_vals(1)

            call from_kim_factory_get_kim(kim_type_of_run, kim_instance)
            call kim_instance%init()

            ! Store actual KIM grid boundary (may differ from r_plas)
            kim_r_boundary = kim_xl_grid%xb(kim_xl_grid%npts_b)

            ! Extract everything from KIM onto QL-Balance grid
            kim_npts = kim_plasma%grid_size
            allocate(kim_r(kim_npts), work_old(kim_npts), work_new(nrad))

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

            deallocate(kim_r, work_old, work_new)
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
        write(*, *) "  KIM grid points:     ", kim_npts
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
    ! Helper: interpolate B0 components, clamped at r_plas
    ! ---------------------------------------------------------------
    subroutine interp_B0_components_clamped(kim_npts, kim_r, nbal, bal_r, nbal_inside)
        !! Like interp_B0_components but only interpolates up to
        !! nbal_inside points; clamps the rest to the boundary value.
        use equilibrium_m, only: B0z_equil => B0z, B0th_equil => B0th
        use wave_code_data, only: B0t, B0z

        implicit none

        integer, intent(in) :: kim_npts, nbal, nbal_inside
        real(8), intent(in) :: kim_r(kim_npts), bal_r(nbal)
        real(8), allocatable :: work(:)

        allocate(work(kim_npts))

        work = B0z_equil
        call interp_profile(kim_npts, kim_r, work, nbal_inside, bal_r(1:nbal_inside), B0z(1:nbal_inside))
        B0z(nbal_inside+1:nbal) = B0z(nbal_inside)

        work = B0th_equil
        call interp_profile(kim_npts, kim_r, work, nbal_inside, bal_r(1:nbal_inside), B0t(1:nbal_inside))
        B0t(nbal_inside+1:nbal) = B0t(nbal_inside)

        deallocate(work)

    end subroutine interp_B0_components_clamped

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
            wcd_n => n, wcd_Te => Te, wcd_Ti => Ti, &
            wcd_q => q, wcd_dPhi0 => dPhi0, &
            wcd_kp => kp, wcd_ks => ks, wcd_om_E => om_E, &
            wcd_B0 => B0, wcd_nue => nue, wcd_nui => nui
        use control_mod, only: kim_profiles_from_balance

        implicit none

        integer :: i_mn, kim_npts, kim_plasma_npts
        real(8), allocatable :: kim_r(:)
        real(8), allocatable :: kim_plasma_r(:), work_real(:)

        ! -------------------------------------------------------
        ! 0. Re-inject time-evolved QL-Balance profiles into KIM
        ! -------------------------------------------------------
        ! NOTE: Profile re-injection is only needed in TimeEvolution mode
        ! where profiles change between calls. For SingleStep, skip it
        ! because set_profiles_from_arrays deallocates derived quantities
        ! (B0, kp, ks, vT, nu, x1, x2, etc.) that the solver needs.
        ! Re-injection for TimeEvolution would require also calling
        ! calculate_equil and set_plasma_quantities to recompute derived
        ! quantities on the rg_grid.

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
        ! 2. Loop over modes: solve and store
        ! -------------------------------------------------------
        do i_mn = 1, dim_mn

            ! Set mode numbers for this solve
            kim_m_mode = m_vals(i_mn)
            kim_n_mode = n_vals(i_mn)

            ! Clean up EBdat from previous solve
            call deallocate_EBdat()

            ! For modes 2+, recompute equilibrium quantities (kp, ks, om_E)
            ! which depend on the (m, n) mode numbers.
            if (i_mn > 1) then
                call calculate_equil(.false.)
                call set_plasma_quantities(kim_plasma)
                call interpolate_equil(kim_rg_grid%xb)
            end if

            call kim_instance%run()

            ! run() fills kernels, solves coupled Poisson-Ampere, and calls
            ! postprocess_electric_field which computes Es, Ep, Er, Etheta, Ez.
            ! EBdat%Br = alpha * A_par (self-consistent).

            ! Interpolate KIM fields (on xl_grid%xb) onto
            ! the QL-Balance radial grid and store per-mode.
            kim_npts = size(EBdat%r_grid)
            allocate(kim_r(kim_npts))
            kim_r = EBdat%r_grid

            if (i_mn == 1) then
                write(*,*) '  KIM field grid: r_min=', kim_r(1), &
                           ' r_max=', kim_r(kim_npts), ' npts=', kim_npts
                write(*,*) '  KIM r_plas =', kim_r_plas
                write(*,*) '  KIM |Br| at last grid pt =', abs(EBdat%Br(kim_npts))
                write(*,*) '  KIM |Br| at 2nd-last     =', abs(EBdat%Br(kim_npts-1))
            end if

            ! Es (perpendicular E field in rsp coordinates)
            call interp_complex_profile(kim_npts, kim_r, EBdat%Es, &
                dim_r, bal_r, kim_Es_modes(:, i_mn))

            ! Ep (parallel E field in rsp coordinates)
            call interp_complex_profile(kim_npts, kim_r, EBdat%Ep, &
                dim_r, bal_r, kim_Ep_modes(:, i_mn))

            ! Er (radial E field, cylindrical)
            call interp_complex_profile(kim_npts, kim_r, EBdat%Er, &
                dim_r, bal_r, kim_Er_modes(:, i_mn))

            ! Etheta -> Et (poloidal E field, cylindrical)
            call interp_complex_profile(kim_npts, kim_r, EBdat%Etheta, &
                dim_r, bal_r, kim_Et_modes(:, i_mn))

            ! Ez (axial E field, cylindrical)
            call interp_complex_profile(kim_npts, kim_r, EBdat%Ez, &
                dim_r, bal_r, kim_Ez_modes(:, i_mn))

            ! Br (radial magnetic field perturbation)
            call interp_complex_profile(kim_npts, kim_r, EBdat%Br, &
                dim_r, bal_r, kim_Br_modes(:, i_mn))

            ! jpar (parallel current density: total, electron, ion)
            if (allocated(EBdat%jpar)) then
                call interp_complex_profile(kim_npts, kim_r, EBdat%jpar, &
                    dim_r, bal_r, kim_jpar_modes(:, i_mn))
            end if
            if (allocated(EBdat%jpar_e)) then
                call interp_complex_profile(kim_npts, kim_r, EBdat%jpar_e, &
                    dim_r, bal_r, kim_jpar_e_modes(:, i_mn))
            end if
            if (allocated(EBdat%jpar_i)) then
                call interp_complex_profile(kim_npts, kim_r, EBdat%jpar_i, &
                    dim_r, bal_r, kim_jpar_i_modes(:, i_mn))
            end if

            deallocate(kim_r)

            ! Apply vacuum continuation beyond r_plas:
            ! KIM fields are only valid inside the plasma domain.
            ! Beyond r_plas, use the vacuum Br from KiLCA and zero E fields.
            call apply_vacuum_continuation(i_mn, dim_r, bal_r)

            ! kp and ks (mode-dependent wave vectors on plasma grid)
            kim_plasma_npts = kim_plasma%grid_size
            allocate(kim_plasma_r(kim_plasma_npts))
            allocate(work_real(kim_plasma_npts))
            kim_plasma_r = kim_plasma%r_grid

            work_real = kim_plasma%kp
            call interp_profile(kim_plasma_npts, kim_plasma_r, &
                work_real, dim_r, bal_r, kim_kp_modes(:, i_mn))

            work_real = kim_plasma%ks
            call interp_profile(kim_plasma_npts, kim_plasma_r, &
                work_real, dim_r, bal_r, kim_ks_modes(:, i_mn))

            deallocate(kim_plasma_r, work_real)

            write(*, '(A,I3,A,I4,A,I4,A)') &
                "  KIM adapter: solved mode ", i_mn, &
                "  (m=", m_vals(i_mn), ", n=", n_vals(i_mn), ")"

        end do

        ! -------------------------------------------------------
        ! 3. Re-extract updated derived quantities when profiles
        !    came from QL-Balance (they change with equilibrium)
        ! -------------------------------------------------------
        if (kim_profiles_from_balance) then
            kim_plasma_npts = kim_plasma%grid_size
            allocate(kim_plasma_r(kim_plasma_npts), work_real(kim_plasma_npts))
            kim_plasma_r = kim_plasma%r_grid

            work_real = kim_plasma%om_E
            call interp_profile(kim_plasma_npts, kim_plasma_r, &
                work_real, dim_r, bal_r, wcd_om_E)

            work_real = kim_plasma%spec(0)%nu
            call interp_profile(kim_plasma_npts, kim_plasma_r, &
                work_real, dim_r, bal_r, wcd_nue)

            work_real = kim_plasma%spec(1)%nu
            call interp_profile(kim_plasma_npts, kim_plasma_r, &
                work_real, dim_r, bal_r, wcd_nui)

            work_real = kim_B0
            call interp_profile(kim_plasma_npts, kim_plasma_r, &
                work_real, dim_r, bal_r, wcd_B0)

            deallocate(kim_plasma_r, work_real)
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
    ! Helper: deallocate equilibrium_m arrays for re-initialization
    ! ---------------------------------------------------------------
    subroutine deallocate_equilibrium_arrays()
        !! Deallocate arrays from equilibrium_m that calculate_equil
        !! allocates with bare allocate() (no guard).
        !! Separated from deallocate_plasma_derived to avoid circular
        !! module dependency (equilibrium_m uses species_m).
        implicit none

        if (allocated(kim_B0)) deallocate(kim_B0)
        if (allocated(eq_B0z)) deallocate(eq_B0z)
        if (allocated(eq_B0th)) deallocate(eq_B0th)
        if (allocated(eq_hz)) deallocate(eq_hz)
        if (allocated(eq_hth)) deallocate(eq_hth)
        if (allocated(eq_equil_grid)) deallocate(eq_equil_grid)
        if (allocated(eq_u)) deallocate(eq_u)
        if (allocated(eq_dpress_prof)) deallocate(eq_dpress_prof)

    end subroutine deallocate_equilibrium_arrays

    ! ---------------------------------------------------------------
    ! Helper: deallocate EBdat fields between mode solves
    ! ---------------------------------------------------------------
    subroutine deallocate_EBdat()
        !! Deallocate all allocated components of the global EBdat
        !! so that the next call to run() can re-allocate them.
        use fields_m, only: EBdat_t

        implicit none

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
        if (allocated(EBdat%Phi))               deallocate(EBdat%Phi)
        if (allocated(EBdat%Phi_e))             deallocate(EBdat%Phi_e)
        if (allocated(EBdat%Phi_i))             deallocate(EBdat%Phi_i)
        if (allocated(EBdat%Phi_aligned))       deallocate(EBdat%Phi_aligned)
        if (allocated(EBdat%Phi_MA))            deallocate(EBdat%Phi_MA)
        if (allocated(EBdat%Phi_MA_ideal))      deallocate(EBdat%Phi_MA_ideal)
        if (allocated(EBdat%Phi_MA_asymptotic)) deallocate(EBdat%Phi_MA_asymptotic)
        if (allocated(EBdat%jpar))              deallocate(EBdat%jpar)
        if (allocated(EBdat%jpar_e))            deallocate(EBdat%jpar_e)
        if (allocated(EBdat%jpar_i))            deallocate(EBdat%jpar_i)

    end subroutine deallocate_EBdat

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
