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

    !! Per-mode stored field results (nrad, dim_mn)
    !! Filled by kim_run_for_all_modes, read by kim_get_wave_fields
    complex(8), allocatable :: kim_Es_modes(:,:)
    complex(8), allocatable :: kim_Ep_modes(:,:)
    complex(8), allocatable :: kim_Er_modes(:,:)
    complex(8), allocatable :: kim_Et_modes(:,:)
    complex(8), allocatable :: kim_Ez_modes(:,:)
    complex(8), allocatable :: kim_Br_modes(:,:)

    !! Per-mode stored wave vectors (nrad, dim_mn)
    !! kp and ks depend on (m,n) via the equilibrium formulas.
    !! Filled by kim_run_for_all_modes, read by kim_get_wave_vectors
    real(8), allocatable :: kim_kp_modes(:,:)
    real(8), allocatable :: kim_ks_modes(:,:)

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
        !! Run KIM electrostatic solver for each (m,n) mode and
        !! store the resulting fields in per-mode arrays.
        !!
        !! After this call, kim_Es_modes(:,i_mn) etc. hold the
        !! field solutions interpolated onto the QL-Balance radial
        !! grid.  kim_get_wave_fields(i_mn) copies from these
        !! arrays into the wave_code_data module scalars.
        use wave_code_data, only: dim_mn, m_vals, n_vals, &
            dim_r, bal_r => r

        implicit none

        integer :: i_mn, kim_npts, kim_plasma_npts
        real(8), allocatable :: kim_r(:)
        real(8), allocatable :: kim_plasma_r(:), work_real(:)

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

        allocate(kim_Es_modes(dim_r, dim_mn))
        allocate(kim_Ep_modes(dim_r, dim_mn))
        allocate(kim_Er_modes(dim_r, dim_mn))
        allocate(kim_Et_modes(dim_r, dim_mn))
        allocate(kim_Ez_modes(dim_r, dim_mn))
        allocate(kim_Br_modes(dim_r, dim_mn))
        allocate(kim_kp_modes(dim_r, dim_mn))
        allocate(kim_ks_modes(dim_r, dim_mn))

        kim_Es_modes = (0.0d0, 0.0d0)
        kim_Ep_modes = (0.0d0, 0.0d0)
        kim_Er_modes = (0.0d0, 0.0d0)
        kim_Et_modes = (0.0d0, 0.0d0)
        kim_Ez_modes = (0.0d0, 0.0d0)
        kim_Br_modes = (0.0d0, 0.0d0)
        kim_kp_modes = 0.0d0
        kim_ks_modes = 0.0d0

        ! -------------------------------------------------------
        ! 2. Loop over modes: solve and store
        ! -------------------------------------------------------
        do i_mn = 1, dim_mn

            ! Set mode numbers for this solve
            kim_m_mode = m_vals(i_mn)
            kim_n_mode = n_vals(i_mn)

            ! Clean up EBdat from previous solve
            call deallocate_EBdat()

            ! Create a fresh solver, init grids/equilibrium, solve
            call from_kim_factory_get_kim('electrostatic', kim_instance)
            call kim_instance%init()
            call kim_instance%run()

            ! run() calls postprocess_electric_field which computes
            ! Es, Ep, Er, Etheta, Ez from the Poisson-solved Phi.
            ! EBdat%Br is set by set_Br_field inside solve_poisson.

            ! Interpolate KIM fields (on xl_grid%xb) onto
            ! the QL-Balance radial grid and store per-mode.
            kim_npts = size(EBdat%r_grid)
            allocate(kim_r(kim_npts))
            kim_r = EBdat%r_grid

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

            deallocate(kim_r)

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
        use wave_code_data, only: dim_r, &
            wcd_n => n, wcd_Te => Te, wcd_Ti => Ti, &
            wcd_Vz => Vz, wcd_dPhi0 => dPhi0
        use plasma_parameters, only: params_b
        use baseparam_mod, only: ev, rtor
        use grid_mod, only: Ercov

        implicit none

        integer :: k

        do k = 1, dim_r
            wcd_n(k)     = params_b(1, k)
            wcd_Te(k)    = params_b(3, k) / ev   ! erg -> eV
            wcd_Ti(k)    = params_b(4, k) / ev   ! erg -> eV
            wcd_Vz(k)    = params_b(2, k) * rtor  ! rad/s -> cm/s
            wcd_dPhi0(k) = -Ercov(k)
        end do

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

end module kim_wave_code_adapter_m
