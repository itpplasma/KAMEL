module periodic_background_m
    !> Bridges KIM's real background primitives to the vendored forced-periodicity
    !> utility (periodization_m::make_periodic). Provides:
    !>   - kim_aperfuns:  an aperfuns_i-conforming callback that samples the TRUE
    !>                    (aperiodic) primitive profiles of the global `plasma` at
    !>                    an arbitrary radius,
    !>   - sample_periodic_primitives: a thin wrapper that returns the PERIODIZED
    !>                    primitives at a query radius by driving make_periodic.
    !>
    !> LEGACY SAMPLING CONTRACT: the public five-value sampler returns
    !>     funs = [ n_e, Te, Ti, q, Er ].
    !> The full periodic build caches and samples n_s, T_s, and V_parallel,s for
    !> every species.

    use KIM_kinds_m, only: dp
    use periodization_m, only: make_periodic

    implicit none
    private

    public :: kim_aperfuns, kim_geom_aperfuns
    public :: sample_periodic_primitives, build_periodic_plasma
    public :: sample_periodic_species_profiles
    public :: reset_true_background_cache

    integer, parameter, public :: PERIODIC_BACKGROUND_OK = 0
    integer, parameter, public :: PERIODIC_BACKGROUND_INVALID_MAIN_ION = 1
    integer, parameter, public :: PERIODIC_BACKGROUND_NONPHYSICAL_PROFILE = 2

    !> Number of values in the legacy public primitive-sampling interface.
    integer, parameter :: n_legacy_primitives = 5

    !> Number of GEOMETRY functions periodized, and their FIXED order.
    !>     geom = [ B0, ks, kp, om_E ].
    integer, parameter :: n_geom = 4

    !> Cached TRUE (aperiodic) background, captured ONCE from the global plasma on
    !> the FIRST build_periodic_plasma call and reused on every later call.
    !>
    !> WHY A CACHE: build_periodic_plasma REDIRECTS the module-global `plasma`
    !> onto the window on its first call, destroying the true global-grid
    !> primitives and geometry. The convergence harness then calls it AGAIN (once
    !> per n_rg) with `plasma` already pointing at the previous window, so a naive
    !> re-sample of `plasma` would periodize an already-periodized profile. Caching
    !> the true background r-grid + primitives + geometry on the first call, and
    !> having BOTH the primitive callback (kim_aperfuns) and the geometry callback
    !> (kim_geom_aperfuns) read the cache, makes repeated calls exactly correct:
    !> every call periodizes the SAME true background.
    !>
    !> AUTO-INVALIDATION: the cache records the species_m::profiles_generation it
    !> was captured at (cached_generation). set_profiles_from_arrays bumps that
    !> counter on every profile swap (new case / time step / parameter-scan point),
    !> so capture_true_background re-captures whenever the generations differ. This
    !> means a second case in the same process gets a FRESH capture of ITS true
    !> background -- no silent staleness -- without species_m having to know about
    !> this cache (which would create a circular module dependency). No manual
    !> reset is needed after set_profiles_from_arrays.
    !>
    !> Order in cache_prim for S species:
    !>   [ n_0..n_(S-1), T_0..T_(S-1), Vpar_0..Vpar_(S-1), q, Er ].
    !> Order in cache_geom: [ B0, ks, kp, om_E ]     (matches kim_geom_aperfuns).
    integer, save :: cached_generation = -1
    integer, save :: cached_n_species = 0
    real(dp), allocatable, save :: cache_r(:)
    real(dp), allocatable, save :: cache_prim(:, :)   ! (grid_size, 3*S + 2)
    real(dp), allocatable, save :: cache_geom(:, :)   ! (grid_size, n_geom)

contains

    !> Force the cached true background to be re-captured on the next
    !> capture_true_background call. Not needed for correctness -- invalidation is
    !> automatic via species_m::profiles_generation (see the cache header) -- but
    !> kept for callers that want to drop the snapshot explicitly.
    subroutine reset_true_background_cache()
        cached_generation = -1
        cached_n_species = 0
        if (allocated(cache_r))    deallocate(cache_r)
        if (allocated(cache_prim)) deallocate(cache_prim)
        if (allocated(cache_geom)) deallocate(cache_geom)
    end subroutine reset_true_background_cache

    !> Capture the TRUE background (r-grid, primitives, geometry) from the global
    !> plasma into the module cache. Self-guarding: a no-op while the cache matches
    !> the current species_m::profiles_generation, and a (re)capture otherwise.
    !> Must be called while `plasma` still holds the true global-grid background
    !> (i.e. at the top of build_periodic_plasma, before the redirect).
    subroutine capture_true_background()
        use species_m, only: plasma, profiles_generation

        integer :: gs, i, sp, ncols

        if (cached_generation == profiles_generation) return

        if (allocated(cache_r))    deallocate(cache_r)
        if (allocated(cache_prim)) deallocate(cache_prim)
        if (allocated(cache_geom)) deallocate(cache_geom)

        gs = plasma%grid_size
        cached_n_species = plasma%n_species
        ncols = 3*cached_n_species + 2
        allocate(cache_r(gs), cache_prim(gs, ncols), cache_geom(gs, n_geom))

        cache_r = plasma%r_grid(1:gs)
        do i = 1, gs
            do sp = 0, cached_n_species - 1
                cache_prim(i, sp + 1) = plasma%spec(sp)%n(i)
                cache_prim(i, cached_n_species + sp + 1) = plasma%spec(sp)%T(i)
                cache_prim(i, 2*cached_n_species + sp + 1) = plasma%spec(sp)%Vpar(i)
            end do
            cache_prim(i, 3*cached_n_species + 1) = plasma%q(i)
            cache_prim(i, 3*cached_n_species + 2) = plasma%Er(i)

            cache_geom(i, 1) = plasma%B0(i)          ! equilibrium field magnitude
            cache_geom(i, 2) = plasma%ks(i)          ! perpendicular wavenumber
            cache_geom(i, 3) = plasma%kp(i)          ! parallel wavenumber
            cache_geom(i, 4) = plasma%om_E(i)        ! ExB rotation frequency
        end do

        cached_generation = profiles_generation
    end subroutine capture_true_background

    !> 4-point Lagrange interpolation of ALL columns of a cached-background table
    !> at radius x, using the same binsrc + plag_coeff stencil as the rest of KIM,
    !> with the query clamped to the cached r-grid domain (see kim_aperfuns header
    !> note on why clamping is needed and harmless).
    subroutine interp_cache_row(col_table, ncols, x, out)
        real(dp), intent(in) :: col_table(:, :)
        integer, intent(in) :: ncols
        real(dp), intent(in) :: x
        real(dp), intent(out) :: out(ncols)

        integer, parameter :: nlagr = 4, nder = 0
        real(dp) :: coef(0:nder, nlagr), xc
        integer :: gs, ir, ibeg, iend, c

        gs = size(cache_r)
        xc = max(cache_r(1), min(x, cache_r(gs)))
        call binsrc(cache_r, 1, gs, xc, ir)
        ibeg = max(1, ir - nlagr/2)
        iend = ibeg + nlagr - 1
        if (iend > gs) then
            iend = gs
            ibeg = iend - nlagr + 1
        end if
        call plag_coeff(nlagr, nder, xc, cache_r(ibeg:iend), coef)
        do c = 1, ncols
            out(c) = sum(coef(0, :) * col_table(ibeg:iend, c))
        end do
    end subroutine interp_cache_row

    !> Aperiodic sample callback conforming to periodization_m::aperfuns_i.
    !> Evaluates the TRUE primitive profiles at radius x via 4-point Lagrange
    !> interpolation of the CACHED true background, in the FIXED order
    !>     funs = [ n_e(x), Te(x), Ti(x), q(x), Er(x) ].
    !>
    !> Reads the module cache (captured once from the global plasma) rather than
    !> the live `plasma`, so repeated build_periodic_plasma calls (which redirect
    !> `plasma` onto the window) still periodize the SAME true background.
    !>
    !> Clamping: x is clamped to [cache_r(1), cache_r(end)] before interpolation,
    !> returning the endpoint value for out-of-range x. This is required because
    !> make_periodic queries x_inperiod +- L, which can fall outside the raw
    !> profile domain. In the deep as-is / periodic-image regions the far images
    !> are weighted 0 by the localizer, so clamping is harmless there; it only
    !> matters inside the transition region.
    subroutine kim_aperfuns(nfuns, x, funs)
        integer, intent(in) :: nfuns
        real(dp), intent(in) :: x
        real(dp), intent(out) :: funs(nfuns)

        ! Ensure the cache is fresh for the current profile generation. If this
        ! callback is invoked directly (e.g. a unit test exercising kim_aperfuns /
        ! sample_periodic_primitives without going through build_periodic_plasma),
        ! this captures on first use; it also re-captures if the profiles changed.
        ! MUST happen before cache_prim is passed as an argument (passing an
        ! unallocated allocatable to an assumed-shape dummy is illegal).
        ! capture_true_background self-guards on generation, so it is a cheap
        ! no-op in the normal solver path (build_periodic_plasma already captured).
        call capture_true_background()

        block
            real(dp) :: all_funs(3*cached_n_species + 2)

            call interp_cache_row(cache_prim, size(all_funs), x, all_funs)
            funs(1) = all_funs(1)
            funs(2) = all_funs(cached_n_species + 1)
            funs(3) = all_funs(cached_n_species + 2)
            funs(4) = all_funs(3*cached_n_species + 1)
            funs(5) = all_funs(3*cached_n_species + 2)
        end block
    end subroutine kim_aperfuns

    !> Internal callback used by build_periodic_plasma. Unlike the legacy public
    !> sampler, this passes every species profile through make_periodic in the
    !> cache order [n_0..n_(S-1), T_0..T_(S-1), Vpar_0..Vpar_(S-1), q, Er].
    subroutine kim_all_profiles_aperfuns(nfuns, x, funs)
        integer, intent(in) :: nfuns
        real(dp), intent(in) :: x
        real(dp), intent(out) :: funs(nfuns)

        call capture_true_background()
        call interp_cache_row(cache_prim, nfuns, x, funs)
    end subroutine kim_all_profiles_aperfuns

    !> Aperiodic sample callback for the equilibrium GEOMETRY, conforming to
    !> periodization_m::aperfuns_i. Evaluates the TRUE geometry at radius x via
    !> 4-point Lagrange interpolation of the CACHED true background, in the FIXED
    !> order
    !>     funs = [ B0(x), ks(x), kp(x), om_E(x) ].
    !>
    !> These carry the explicit 1/r geometric factors (ks, kp, om_E) and the
    !> force-balance-integrated field (B0), which make_periodic did NOT touch when
    !> only the primitives n,T,q,Er were periodized. Their true radial variation
    !> (including 1/r) drifts ~0.3% across one period, so the raw window integrand
    !> is NOT L-periodic and the periodic-trapezoidal r_g quadrature plateaus.
    !> Periodizing them here (faithfully: the localizer blends the TRUE profile
    !> into a C-infinity-periodic window function, keeping the true variation in
    !> the as-is core) restores integrand periodicity and spectral convergence.
    subroutine kim_geom_aperfuns(nfuns, x, funs)
        integer, intent(in) :: nfuns
        real(dp), intent(in) :: x
        real(dp), intent(out) :: funs(nfuns)

        ! Ensure a fresh cache before passing cache_geom (see kim_aperfuns note).
        call capture_true_background()

        call interp_cache_row(cache_geom, nfuns, x, funs)
    end subroutine kim_geom_aperfuns

    !> Return the PERIODIZED primitives at query radius x, driving make_periodic
    !> with kim_aperfuns. Layout: half_period = dx_asis + dx_tr, period
    !> L = 2*half_period, an as-is core of half-width dx_asis centred on x_mid,
    !> flanked by transition bands of width dx_tr. Output order matches
    !> kim_aperfuns: funs_per = [ n_e, Te, Ti, q, Er ].
    subroutine sample_periodic_primitives(x_mid, dx_asis, dx_tr, x, funs_per)
        real(dp), intent(in) :: x_mid, dx_asis, dx_tr, x
        real(dp), intent(out) :: funs_per(n_legacy_primitives)

        call make_periodic(kim_aperfuns, n_legacy_primitives, x_mid, dx_asis, dx_tr, &
                           x, funs_per)
    end subroutine sample_periodic_primitives

    !> Return independently periodized densities and temperatures for every
    !> configured species, plus q and Er. When requested, Vpar_per returns each
    !> species' field-parallel Maxwellian drift. Arrays use indices 0:S-1.
    subroutine sample_periodic_species_profiles(x_mid, dx_asis, dx_tr, x, &
                                                n_per, T_per, q_per, Er_per, Vpar_per)
        real(dp), intent(in) :: x_mid, dx_asis, dx_tr, x
        real(dp), intent(out) :: n_per(0:), T_per(0:)
        real(dp), intent(out) :: q_per, Er_per
        real(dp), intent(out), optional :: Vpar_per(0:)

        integer :: sp, n_profile_values

        call capture_true_background()
        if (size(n_per) /= cached_n_species .or. &
            size(T_per) /= cached_n_species) then
            error stop 'periodic species sampler received the wrong species count'
        end if
        if (present(Vpar_per)) then
            if (size(Vpar_per) /= cached_n_species) then
                error stop 'periodic flow sampler received the wrong species count'
            end if
        end if

        n_profile_values = 3*cached_n_species + 2
        block
            real(dp) :: all_funs(n_profile_values)

            call make_periodic(kim_all_profiles_aperfuns, n_profile_values, x_mid, &
                               dx_asis, dx_tr, x, all_funs)
            do sp = 0, cached_n_species - 1
                n_per(sp) = all_funs(sp + 1)
                T_per(sp) = all_funs(cached_n_species + sp + 1)
                if (present(Vpar_per)) then
                    Vpar_per(sp) = all_funs(2*cached_n_species + sp + 1)
                end if
            end do
            q_per = all_funs(3*cached_n_species + 1)
            Er_per = all_funs(3*cached_n_species + 2)
        end block
    end subroutine sample_periodic_species_profiles

    !> Return the PERIODIZED geometry [ B0, ks, kp, om_E ] at query radius x,
    !> driving make_periodic with kim_geom_aperfuns. Same window layout as
    !> sample_periodic_primitives.
    subroutine sample_periodic_geometry(x_mid, dx_asis, dx_tr, x, geom_per)
        real(dp), intent(in) :: x_mid, dx_asis, dx_tr, x
        real(dp), intent(out) :: geom_per(n_geom)

        call make_periodic(kim_geom_aperfuns, n_geom, x_mid, dx_asis, dx_tr, x, geom_per)
    end subroutine sample_periodic_geometry

    !> Build a fully-populated PERIODIC `plasma` on an equidistant window grid
    !> centred on rm. The window
    !>     [rm - L/2, rm + L/2],   L = 2*(dx_asis + dx_tr),   n_rg boundary pts
    !> serves as BOTH plasma%r_grid AND the global rg_grid, so the derived-
    !> quantity recompute (calculate_equil -> set_plasma_quantities) runs on the
    !> window with NO cross-grid interpolation (interpolate_plasma_backs becomes
    !> the identity because plasma%r_grid == rg_grid%xb).
    !>
    !> On entry the global `plasma` must hold the TRUE (aperiodic) primitives on
    !> the global grid (as after an electrostatic %init()): kim_aperfuns samples
    !> them while the periodized primitives are laid down on the window.
    !>
    !> On exit plasma%{spec%n,T,dndr,dTdr,q,Er} are the PERIODIZED primitives on
    !> the window and every derived quantity (B0, kp, om_E, vT, nu, omega_c,
    !> rho_L, lambda_D, A1, A2, x1, x2 and the susceptibility functions I00..I21)
    !> is recomputed from them and is therefore L-periodic.
    !>
    !> NOTE: this MUTATES the module-global rg_grid and plasma. It is meant to be
    !> called once per forced-periodicity assembly on an independent process.
    !>
    !> SIDE EFFECT: not filesystem-side-effect-free. set_plasma_quantities
    !> unconditionally calls write_species_backs / write_plasma_backs, which
    !> create output_path//backs/ and write the background profile .dat files
    !> (rho_L, lambda_D, kp, om_E, ...) for the window plasma.
    subroutine build_periodic_plasma(rm, dx_asis, dx_tr, n_rg, info)
        use species_m, only: plasma, reallocate, deallocate_plasma_derived, &
                             calc_plasma_parameter_derivs, calculate_plasma_backs, &
                             interpolate_plasma_backs, compute_rg_cell_centers, &
                             calculate_thermodynamic_forces_and_susc, &
                             write_species_backs, write_species_cc_quantities, &
                             write_plasma_backs
        use equilibrium_m, only: calculate_equil, B0, B0z, B0th, hz, hth, &
                                 equil_grid, u, dpress_prof
        use grid_m, only: rg_grid
        use config_m, only: number_of_ion_species

        real(dp), intent(in) :: rm, dx_asis, dx_tr
        integer, intent(in) :: n_rg
        integer, intent(out), optional :: info

        real(dp) :: L, r_lo, r_hi, u0_seed
        real(dp), allocatable :: r_win(:)
        real(dp), allocatable :: n_win(:, :), T_win(:, :), Vpar_win(:, :)
        real(dp), allocatable :: q_win(:), Er_win(:)
        integer :: j, sigma, npts, n_species, profile_status

        if (present(info)) info = PERIODIC_BACKGROUND_OK

        ! Capture the TRUE background (r-grid, primitives, geometry) from the
        ! global plasma ONCE, before any redirect destroys it; all periodization
        ! (primitives AND geometry) then reads this cache, so repeated calls stay
        ! correct even though `plasma` is redirected onto the window after call 1.
        call capture_true_background()
        n_species = cached_n_species

        ! ---- 1. Construct an equidistant window without mutating globals. -----
        L = 2.0_dp * (dx_asis + dx_tr)
        r_lo = rm - 0.5_dp * L
        r_hi = rm + 0.5_dp * L

        ! ---- 2. Sample and validate PERIODIZED primitives before redirect. ----
        ! kim_aperfuns reads the CACHED true primitives (captured above), so this
        ! is robust to the redirect and to repeated calls. Sample into temporaries;
        ! only a valid profile is allowed to mutate rg_grid or plasma below.
        ! Dynamic order: [n_0..n_(S-1), T_0..T_(S-1),
        !                 Vpar_0..Vpar_(S-1), q, Er].
        allocate(r_win(n_rg), n_win(n_rg, 0:n_species - 1), &
                 T_win(n_rg, 0:n_species - 1), &
                 Vpar_win(n_rg, 0:n_species - 1), q_win(n_rg), Er_win(n_rg))
        do j = 1, n_rg
            r_win(j) = r_lo + L*real(j - 1, dp)/real(n_rg, dp)
            call sample_periodic_species_profiles(rm, dx_asis, dx_tr, r_win(j), &
                n_win(j, :), T_win(j, :), q_win(j), Er_win(j), Vpar_win(j, :))
        end do

        ! Preserve independently periodized impurity densities and use the main
        ! ion (species 1) as the dependent density. This documents one unique
        ! quasineutrality policy and avoids changing impurity concentration
        ! ratios merely because a species has a higher charge state.
        call enforce_dependent_main_ion_density(n_win, T_win, Vpar_win, profile_status)
        if (profile_status /= PERIODIC_BACKGROUND_OK) then
            deallocate(r_win, n_win, T_win, Vpar_win, q_win, Er_win)
            if (present(info)) then
                info = profile_status
                return
            end if
            if (profile_status == PERIODIC_BACKGROUND_INVALID_MAIN_ION) then
                error stop 'periodic background requires a positive-Z main ion'
            end if
            error stop 'periodic background contains a nonphysical species profile'
        end if

        ! Validation succeeded. Redirect the global rg grid to the same
        ! endpoint-exclusive equidistant window used for sampling.
        npts = n_rg
        call rg_grid%grid_init_equidistant(npts, r_lo, r_hi, 'rg')
        call rg_grid%grid_generate_equidistant()
        r_win = rg_grid%xb

        ! Capture the TRUE B0 at the window's left edge (from the cached true
        ! geometry) as the force-balance ODE seed u0 = B0(r_lo)**2. Without this
        ! the window solve restarts from the pressureless vacuum initial condition
        ! and loses the pressure work accumulated over [r_min, r_lo], offsetting
        ! the absolute B0 level (hence omega_c, rho_L) by ~O(0.5%). See
        ! calculate_equil's u0_seed argument.
        u0_seed = interp_cache_B0(r_lo)**2

        ! ---- 3. Redirect the plasma onto the window grid. ---------------------
        ! Clear ALL derived state so the recompute chain reallocates cleanly at
        ! the new (window) size: equilibrium_m arrays, plasma-level derived
        ! arrays (deallocate_plasma_derived), and the thermodynamic/
        ! susceptibility arrays (guarded by `.not. allocated`, so a stale old
        ! size would otherwise be kept).
        if (allocated(B0))          deallocate(B0)
        if (allocated(B0z))         deallocate(B0z)
        if (allocated(B0th))        deallocate(B0th)
        if (allocated(hz))          deallocate(hz)
        if (allocated(hth))         deallocate(hth)
        if (allocated(equil_grid))  deallocate(equil_grid)
        if (allocated(u))           deallocate(u)
        if (allocated(dpress_prof)) deallocate(dpress_prof)
        call deallocate_plasma_derived()
        call deallocate_thermo_and_susc()

        plasma%grid_size = rg_grid%npts_b
        call reallocate(plasma%r_grid, rg_grid%npts_b)
        plasma%r_grid = r_win

        ! Write the periodized primitives onto the redirected plasma. Ion
        ! density follows from n_e by quasineutrality (as elsewhere in KIM).
        call reallocate(plasma%q, rg_grid%npts_b)
        call reallocate(plasma%Er, rg_grid%npts_b)
        call reallocate(plasma%spec(0)%n, rg_grid%npts_b)
        call reallocate(plasma%spec(0)%T, rg_grid%npts_b)
        call reallocate(plasma%spec(0)%Vpar, rg_grid%npts_b)
        plasma%q  = q_win
        plasma%Er = Er_win
        plasma%spec(0)%n = n_win(:, 0)
        plasma%spec(0)%T = T_win(:, 0)
        plasma%spec(0)%Vpar = Vpar_win(:, 0)
        do sigma = 1, number_of_ion_species
            call reallocate(plasma%spec(sigma)%n, rg_grid%npts_b)
            call reallocate(plasma%spec(sigma)%T, rg_grid%npts_b)
            call reallocate(plasma%spec(sigma)%Vpar, rg_grid%npts_b)
            plasma%spec(sigma)%n = n_win(:, sigma)
            plasma%spec(sigma)%T = T_win(:, sigma)
            plasma%spec(sigma)%Vpar = Vpar_win(:, sigma)
        end do

        ! Force calculate_equil / calc_plasma_parameter_derivs to recompute the
        ! profile derivatives (dndr, dTdr, dqdr) on the window instead of reusing
        ! stale global-grid values.
        do sigma = 0, number_of_ion_species
            if (allocated(plasma%spec(sigma)%dndr)) deallocate(plasma%spec(sigma)%dndr)
            if (allocated(plasma%spec(sigma)%dTdr)) deallocate(plasma%spec(sigma)%dTdr)
        end do
        if (allocated(plasma%dqdr)) deallocate(plasma%dqdr)

        deallocate(r_win, n_win, T_win, Vpar_win, q_win, Er_win)

        ! ---- 4./5. Recompute derivatives + all derived quantities. ------------
        ! calculate_equil: recomputes dndr/dTdr/dqdr (via calc_plasma_parameter_
        ! derivs on the equidistant window) and B0, ks, kp, om_E from the
        ! periodized q/Er and btor/R0/m_mode/n_mode. write_out=.false. suppresses
        ! only the equilibrium (B0z/B0th/...) file dump; the profile .dat writes
        ! below (via set_plasma_quantities) still happen -- see SIDE EFFECT note.
        ! u0_seed carries the true B0(r_lo)**2 so the window B0 level matches the
        ! global solve (see the capture above).
        call calculate_equil(.false., u0_seed=u0_seed)

        ! ---- 6. Periodize the equilibrium GEOMETRY on the window. --------------
        ! calculate_equil just wrote B0/ks/kp/om_E on the window with their TRUE
        ! radial variation -- but that variation (explicit 1/r in ks/kp/om_E, plus
        ! the force-balance-integrated B0) is NOT L-periodic, so the raw window
        ! integrand is non-periodic (seam jump ~3e-3 in ks/om_E) and the periodic-
        ! trapezoidal r_g quadrature plateaus. Overwrite B0/ks/kp/om_E with their
        ! FAITHFULLY periodized values (make_periodic blends the TRUE geometry into
        ! a C-infinity-periodic window function; true variation kept in the as-is
        ! core). Done BEFORE calculate_plasma_backs so omega_c/rho_L (from B0) and
        ! the susceptibility (x1 from kp, x2 from om_E) inherit the periodicity.
        call periodize_window_geometry()

        ! Inlined set_plasma_quantities so the Fix-1 periodic-derivative overwrite
        ! can be injected between calc_plasma_parameter_derivs (which writes dndr/
        ! dTdr/dqdr with NON-periodic one-sided endpoint FD) and the A1/A2 forming
        ! in calculate_thermodynamic_forces_and_susc. Steps mirror set_plasma_
        ! quantities (species_mod.f90): calc_plasma_parameter_derivs -> [FIX 1:
        ! periodic central-difference dndr/dTdr/dqdr] -> calculate_plasma_backs
        ! (vT, nu, omega_c, rho_L, lambda_D, z0 -- omega_c/rho_L now from the
        ! PERIODIZED B0) -> interpolate_plasma_backs(rg_grid%xb) (IDENTITY here
        ! because plasma%r_grid == rg_grid%xb; Lagrange is exact at nodes so the
        ! periodized geometry/derivatives survive) -> compute_rg_cell_centers ->
        ! calculate_thermodynamic_forces_and_susc (A1, A2, x1, x2, I00..I21 from
        ! the PERIODIZED kp/om_E) -> write_*_backs (profile .dat dumps).
        call calc_plasma_parameter_derivs

        call make_window_derivs_periodic()

        call calculate_plasma_backs(plasma)
        call interpolate_plasma_backs(plasma, rg_grid%xb)
        call compute_rg_cell_centers(plasma)
        call calculate_thermodynamic_forces_and_susc(plasma)

        do sigma = 0, plasma%n_species - 1
            call write_species_backs(plasma%spec(sigma), plasma%r_grid)
            call write_species_cc_quantities(plasma%spec(sigma), rg_grid%xc)
        end do
        call write_plasma_backs(plasma, plasma%r_grid)

    contains

        !> 4-point Lagrange interpolation of the CACHED true B0 at radius x. Used
        !> for the force-balance ODE seed u0 = B0(r_lo)**2; reads the cache so it
        !> is robust to the redirect and to repeated calls.
        real(dp) function interp_cache_B0(x) result(b0x)
            real(dp), intent(in) :: x
            real(dp) :: geom(n_geom)
            call interp_cache_row(cache_geom, n_geom, x, geom)
            b0x = geom(1)
        end function interp_cache_B0

        !> Overwrite the window B0/ks/kp/om_E (and the equilibrium_m B0 array, kept
        !> consistent for interpolate_plasma_backs / any downstream reader) with
        !> their make_periodic-periodized values, sampled per window node from the
        !> cached true geometry. omega_c and rho_L are re-derived from the
        !> periodized B0 inside calculate_plasma_backs, so periodizing B0 here
        !> periodizes them too.
        subroutine periodize_window_geometry()
            real(dp) :: geom(n_geom)
            integer :: i

            do i = 1, rg_grid%npts_b
                call sample_periodic_geometry(rm, dx_asis, dx_tr, rg_grid%xb(i), geom)
                plasma%B0(i)   = geom(1)
                plasma%ks(i)   = geom(2)
                plasma%kp(i)   = geom(3)
                plasma%om_E(i) = geom(4)
                if (allocated(B0)) B0(i) = geom(1)
            end do
        end subroutine periodize_window_geometry

        !> FIX 1: overwrite the window profile derivatives dndr/dTdr (all species)
        !> and dqdr with PERIODIC central differences on the equidistant window
        !> grid, replacing the non-periodic one-sided endpoint FD that
        !> calc_plasma_parameter_derivs just wrote (without touching that shared
        !> routine, which the global solver uses).
        !>
        !> The window grid rg_grid%xb is equidistant over exactly one period L,
        !> endpoint-exclusive on the right (xb(i) = r_lo + (i-1)*h, h = L/n_rg, so
        !> xb(n_rg) = r_hi - h). Because the periodized primitives are L-periodic,
        !> the point one step left of xb(1) is xb(N) and one step right of xb(N) is
        !> xb(1); the wrap-around central difference
        !>     dfdr(i) = (f(i+1) - f(i-1)) / (2 h),  f(0) := f(N), f(N+1) := f(1)
        !> is exact-to-O(h^2) AND seam-consistent.
        subroutine make_window_derivs_periodic()
            real(dp) :: h
            integer :: sp, N

            N = rg_grid%npts_b
            if (N < 3) return
            h = rg_grid%xb(2) - rg_grid%xb(1)   ! = L / n_rg

            do sp = 0, number_of_ion_species
                call periodic_central_diff(plasma%spec(sp)%n, plasma%spec(sp)%dndr, h, N)
                call periodic_central_diff(plasma%spec(sp)%T, plasma%spec(sp)%dTdr, h, N)
            end do
            call periodic_central_diff(plasma%q, plasma%dqdr, h, N)
        end subroutine make_window_derivs_periodic

        !> Periodic central difference of an L-periodic array f sampled on an
        !> equidistant, endpoint-exclusive one-period grid of N points, spacing h.
        subroutine periodic_central_diff(f, dfdr, h, N)
            real(dp), intent(in) :: f(:)
            real(dp), intent(inout) :: dfdr(:)
            real(dp), intent(in) :: h
            integer, intent(in) :: N
            integer :: i

            dfdr(1) = (f(2) - f(N)) / (2.0_dp * h)
            do i = 2, N - 1
                dfdr(i) = (f(i+1) - f(i-1)) / (2.0_dp * h)
            end do
            dfdr(N) = (f(1) - f(N-1)) / (2.0_dp * h)
        end subroutine periodic_central_diff

        subroutine enforce_dependent_main_ion_density(density, temperature, &
                                                      parallel_velocity, status)
            use, intrinsic :: ieee_arithmetic, only: ieee_is_finite

            real(dp), intent(inout) :: density(:, 0:)
            real(dp), intent(in) :: temperature(:, 0:)
            real(dp), intent(in) :: parallel_velocity(:, 0:)
            integer, intent(out) :: status
            real(dp) :: impurity_charge_density(size(density, 1))
            integer :: sp

            status = PERIODIC_BACKGROUND_INVALID_MAIN_ION
            if (number_of_ion_species < 1 .or. plasma%spec(1)%Zspec <= 0) then
                return
            end if
            do sp = 2, number_of_ion_species
                if (plasma%spec(sp)%Zspec <= 0) return
            end do

            status = PERIODIC_BACKGROUND_NONPHYSICAL_PROFILE
            if (any(.not. ieee_is_finite(density(:, 0))) .or. &
                any(density(:, 0) <= 0.0_dp) .or. &
                any(.not. ieee_is_finite(temperature)) .or. &
                any(temperature <= 0.0_dp) .or. &
                any(.not. ieee_is_finite(parallel_velocity))) then
                return
            end if
            do sp = 2, number_of_ion_species
                if (any(.not. ieee_is_finite(density(:, sp))) .or. &
                    any(density(:, sp) <= 0.0_dp)) return
            end do

            impurity_charge_density = 0.0_dp
            do sp = 2, number_of_ion_species
                impurity_charge_density = impurity_charge_density + &
                    real(plasma%spec(sp)%Zspec, dp)*density(:, sp)
            end do
            density(:, 1) = (density(:, 0) - impurity_charge_density) / &
                            real(plasma%spec(1)%Zspec, dp)

            if (any(.not. ieee_is_finite(density(:, 1))) .or. &
                any(density(:, 1) <= 0.0_dp)) then
                return
            end if
            status = PERIODIC_BACKGROUND_OK
        end subroutine enforce_dependent_main_ion_density

        !> Deallocate the thermodynamic-force and susceptibility arrays that
        !> calculate_thermodynamic_forces_and_susc allocates behind a
        !> `.not. allocated` guard. Without this the arrays keep their previous
        !> (global-grid) size and the recompute silently skips reallocation.
        subroutine deallocate_thermo_and_susc()
            integer :: sp

            do sp = 0, plasma%n_species - 1
                if (allocated(plasma%spec(sp)%symbI)) deallocate(plasma%spec(sp)%symbI)

                if (allocated(plasma%spec(sp)%A1)) deallocate(plasma%spec(sp)%A1)
                if (allocated(plasma%spec(sp)%A2)) deallocate(plasma%spec(sp)%A2)
                if (allocated(plasma%spec(sp)%x1)) deallocate(plasma%spec(sp)%x1)
                if (allocated(plasma%spec(sp)%x2)) deallocate(plasma%spec(sp)%x2)
                if (allocated(plasma%spec(sp)%I00)) deallocate(plasma%spec(sp)%I00)
                if (allocated(plasma%spec(sp)%I01)) deallocate(plasma%spec(sp)%I01)
                if (allocated(plasma%spec(sp)%I20)) deallocate(plasma%spec(sp)%I20)
                if (allocated(plasma%spec(sp)%I21)) deallocate(plasma%spec(sp)%I21)
                if (allocated(plasma%spec(sp)%I22)) deallocate(plasma%spec(sp)%I22)
                if (allocated(plasma%spec(sp)%I02)) deallocate(plasma%spec(sp)%I02)
                if (allocated(plasma%spec(sp)%I11)) deallocate(plasma%spec(sp)%I11)
                if (allocated(plasma%spec(sp)%I13)) deallocate(plasma%spec(sp)%I13)

                if (allocated(plasma%spec(sp)%n_cc)) deallocate(plasma%spec(sp)%n_cc)
                if (allocated(plasma%spec(sp)%dndr_cc)) deallocate(plasma%spec(sp)%dndr_cc)
                if (allocated(plasma%spec(sp)%T_cc)) deallocate(plasma%spec(sp)%T_cc)
                if (allocated(plasma%spec(sp)%dTdr_cc)) deallocate(plasma%spec(sp)%dTdr_cc)
                if (allocated(plasma%spec(sp)%Vpar_cc)) deallocate(plasma%spec(sp)%Vpar_cc)
                if (allocated(plasma%spec(sp)%nu_cc)) deallocate(plasma%spec(sp)%nu_cc)
                if (allocated(plasma%spec(sp)%vT_cc)) deallocate(plasma%spec(sp)%vT_cc)
                if (allocated(plasma%spec(sp)%omega_c_cc)) deallocate(plasma%spec(sp)%omega_c_cc)
                if (allocated(plasma%spec(sp)%lambda_D_cc)) deallocate(plasma%spec(sp)%lambda_D_cc)
                if (allocated(plasma%spec(sp)%rho_L_cc)) deallocate(plasma%spec(sp)%rho_L_cc)
                if (allocated(plasma%spec(sp)%A1_cc)) deallocate(plasma%spec(sp)%A1_cc)
                if (allocated(plasma%spec(sp)%A2_cc)) deallocate(plasma%spec(sp)%A2_cc)
                if (allocated(plasma%spec(sp)%x1_cc)) deallocate(plasma%spec(sp)%x1_cc)
                if (allocated(plasma%spec(sp)%x2_cc)) deallocate(plasma%spec(sp)%x2_cc)
                if (allocated(plasma%spec(sp)%I00_cc)) deallocate(plasma%spec(sp)%I00_cc)
                if (allocated(plasma%spec(sp)%I01_cc)) deallocate(plasma%spec(sp)%I01_cc)
                if (allocated(plasma%spec(sp)%I10_cc)) deallocate(plasma%spec(sp)%I10_cc)
                if (allocated(plasma%spec(sp)%I20_cc)) deallocate(plasma%spec(sp)%I20_cc)
                if (allocated(plasma%spec(sp)%I21_cc)) deallocate(plasma%spec(sp)%I21_cc)
                if (allocated(plasma%spec(sp)%I12_cc)) deallocate(plasma%spec(sp)%I12_cc)
                if (allocated(plasma%spec(sp)%I22_cc)) deallocate(plasma%spec(sp)%I22_cc)
                if (allocated(plasma%spec(sp)%I02_cc)) deallocate(plasma%spec(sp)%I02_cc)
                if (allocated(plasma%spec(sp)%I13_cc)) deallocate(plasma%spec(sp)%I13_cc)
                if (allocated(plasma%spec(sp)%I11_cc)) deallocate(plasma%spec(sp)%I11_cc)
            end do

            if (allocated(plasma%ks_cc)) deallocate(plasma%ks_cc)
            if (allocated(plasma%Er_cc)) deallocate(plasma%Er_cc)
            if (allocated(plasma%om_E_cc)) deallocate(plasma%om_E_cc)
        end subroutine deallocate_thermo_and_susc

    end subroutine build_periodic_plasma

end module periodic_background_m
