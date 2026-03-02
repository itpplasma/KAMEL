module kim_wave_code_adapter_m
    !! Adapter module that wraps KIM library calls to populate
    !! wave_code_data module arrays, matching the contract expected
    !! by get_dql.f90 and diff_coeffs.f90.
    !!
    !! KIM modules are accessed via per-file INCLUDE_DIRECTORIES in
    !! CMake to avoid collision with the QL-Balance resonances_mod.
    !! All KIM symbols are renamed on import to avoid name clashes.

    use species_m, only: kim_plasma => plasma
    use equilibrium_m, only: kim_B0 => B0
    use config_m, only: kim_type_of_run => type_of_run, &
                        nml_config_path
    use setup_m, only: kim_m_mode => m_mode, kim_n_mode => n_mode
    use grid_m, only: kim_xl_grid => xl_grid, kim_rg_grid => rg_grid
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

    !! Module-level KIM solver instance (reused across calls)
    class(kim_t), allocatable :: kim_instance

    !! Vacuum Br placeholder (filled properly in Task 7)
    complex(8), allocatable, public :: kim_vac_Br(:)

contains

    subroutine kim_initialize(nrad, r_grid)
        !! Initialize KIM backend: read config, profiles, grids.
        !! Populate wave_code_data module arrays with background
        !! quantities by running the electrostatic solver for the
        !! first mode to generate equilibrium and plasma backgrounds.
        use wave_code_data, only: dim_mn, m_vals, n_vals, &
            r => r, q => q, n => n, Te => Te, Ti => Ti, &
            Vth => Vth, Vz => Vz, dPhi0 => dPhi0, &
            kp => kp, ks => ks, om_E => om_E, &
            B0 => B0, nue => nue, nui => nui, &
            flre_path

        implicit none

        integer, intent(in) :: nrad
        real(8), intent(in) :: r_grid(nrad)

        integer :: kim_npts
        real(8), allocatable :: kim_r(:)
        real(8), allocatable :: work_old(:), work_new(:)

        ! -----------------------------------------------------------
        ! 1. Set KIM config path and initialize KIM backend
        ! -----------------------------------------------------------
        ! nml_config_path (from config_m) tells KIM where its
        ! namelist file lives.  Default is ./KIM_config.nml which
        ! the user can override via Task 8 (balance namelist).
        nml_config_path = "./KIM_config.nml"

        call kim_init()

        ! -----------------------------------------------------------
        ! 2. Read (m,n) modes from modes.in via existing utility
        ! -----------------------------------------------------------
        call read_antenna_modes(flre_path)

        ! -----------------------------------------------------------
        ! 3. Allocate wave_code_data arrays on QL-Balance grid
        ! -----------------------------------------------------------
        call allocate_wave_code_data(nrad, r_grid)

        ! -----------------------------------------------------------
        ! 4. Run KIM electrostatic solver for the first mode to
        !    generate equilibrium, grids, and plasma backgrounds
        ! -----------------------------------------------------------
        kim_type_of_run = 'electrostatic'
        kim_m_mode = m_vals(1)
        kim_n_mode = n_vals(1)

        call from_kim_factory_get_kim(kim_type_of_run, kim_instance)
        call kim_instance%init()
        call kim_instance%run()

        ! -----------------------------------------------------------
        ! 5. Extract and interpolate KIM backgrounds onto the
        !    QL-Balance radial grid
        ! -----------------------------------------------------------
        kim_npts = kim_plasma%grid_size
        allocate(kim_r(kim_npts))
        kim_r = kim_plasma%r_grid

        ! Temporary work arrays for interpolation
        allocate(work_old(kim_npts), work_new(nrad))

        ! Safety factor q
        work_old = kim_plasma%q
        call interp_profile(kim_npts, kim_r, work_old, nrad, r, q)

        ! Parallel wavenumber kp
        work_old = kim_plasma%kp
        call interp_profile(kim_npts, kim_r, work_old, nrad, r, kp)

        ! Perpendicular wavenumber ks
        work_old = kim_plasma%ks
        call interp_profile(kim_npts, kim_r, work_old, nrad, r, ks)

        ! ExB drift frequency om_E
        work_old = kim_plasma%om_E
        call interp_profile(kim_npts, kim_r, work_old, nrad, r, om_E)

        ! Equilibrium radial electric field -> dPhi0
        ! Sign convention: dPhi0 = -Er
        work_old = -kim_plasma%Er
        call interp_profile(kim_npts, kim_r, work_old, nrad, r, dPhi0)

        ! Electron density [1/cm^3]
        work_old = kim_plasma%spec(0)%n
        call interp_profile(kim_npts, kim_r, work_old, nrad, r, n)

        ! Electron temperature [eV]
        work_old = kim_plasma%spec(0)%T
        call interp_profile(kim_npts, kim_r, work_old, nrad, r, Te)

        ! Ion temperature [eV]
        work_old = kim_plasma%spec(1)%T
        call interp_profile(kim_npts, kim_r, work_old, nrad, r, Ti)

        ! Electron collision frequency
        work_old = kim_plasma%spec(0)%nu
        call interp_profile(kim_npts, kim_r, work_old, nrad, r, nue)

        ! Ion collision frequency
        work_old = kim_plasma%spec(1)%nu
        call interp_profile(kim_npts, kim_r, work_old, nrad, r, nui)

        ! Background magnetic field B0 (total)
        ! equilibrium_m::B0 is on the plasma r_grid
        work_old = kim_B0
        call interp_profile(kim_npts, kim_r, work_old, nrad, r, B0)

        ! B0 toroidal and poloidal components from equilibrium_m
        call interp_B0_components(kim_npts, kim_r, nrad, r)

        ! Poloidal rotation velocity Vth
        ! KIM does not store Vth; set to zero (can be updated later)
        Vth = 0.0d0

        ! Toroidal rotation Vz
        ! KIM does not store Vz directly in plasma struct;
        ! set to zero (Er already accounts for force balance)
        Vz = 0.0d0

        ! -----------------------------------------------------------
        ! 6. Set up vacuum Br placeholder for form factor
        !    (filled properly in Task 7)
        ! -----------------------------------------------------------
        if (allocated(kim_vac_Br)) deallocate(kim_vac_Br)
        allocate(kim_vac_Br(nrad))
        kim_vac_Br = (1.0d0, 0.0d0)  ! placeholder: unit vacuum Br

        deallocate(kim_r, work_old, work_new)

        write(*, *) "KIM adapter: initialization complete"
        write(*, *) "  Balance grid points: ", nrad
        write(*, *) "  KIM grid points:     ", kim_npts
        write(*, *) "  Number of modes:     ", dim_mn

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

    subroutine kim_run_for_all_modes()
        !! Run KIM solver for all (m,n) modes.
        ! TODO: implement in Task 5
    end subroutine

    subroutine kim_update_profiles()
        !! Update KIM internal profiles from QL-Balance evolved
        !! parameters.
        ! TODO: implement in Task 5
    end subroutine

    subroutine kim_get_wave_fields(i_mn)
        !! Extract wave fields from KIM EBdat into wave_code_data
        !! arrays.
        integer, intent(in) :: i_mn
        ! TODO: implement in Task 6
    end subroutine

    subroutine kim_get_wave_vectors()
        !! Extract wave vectors from KIM plasma into wave_code_data
        !! arrays.
        ! TODO: implement in Task 6
    end subroutine

    subroutine kim_get_background_magnetic_fields()
        !! Extract background B from KIM equilibrium into
        !! wave_code_data arrays.
        ! TODO: implement in Task 6
    end subroutine

    subroutine kim_get_collision_frequencies()
        !! Extract collision frequencies from KIM plasma into
        !! wave_code_data arrays.
        ! TODO: implement in Task 6
    end subroutine

end module kim_wave_code_adapter_m
