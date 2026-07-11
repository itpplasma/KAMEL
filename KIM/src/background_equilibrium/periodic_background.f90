module periodic_background_m
    !> Bridges KIM's real background primitives to the vendored forced-periodicity
    !> utility (periodization_m::make_periodic). Provides:
    !>   - kim_aperfuns:  an aperfuns_i-conforming callback that samples the TRUE
    !>                    (aperiodic) primitive profiles of the global `plasma` at
    !>                    an arbitrary radius,
    !>   - sample_periodic_primitives: a thin wrapper that returns the PERIODIZED
    !>                    primitives at a query radius by driving make_periodic.
    !>
    !> ORDER CONTRACT (relied upon by Phase-2.3): the five sampled functions are
    !> returned in the FIXED order
    !>     funs = [ n_e, Te, Ti, q, Er ].
    !> Do not reorder without updating every consumer.

    use KIM_kinds_m, only: dp
    use periodization_m, only: make_periodic

    implicit none
    private

    public :: kim_aperfuns, sample_periodic_primitives, build_periodic_plasma
    public :: reset_true_background_cache

    !> Number of primitive functions sampled, and their FIXED order.
    integer, parameter :: n_primitives = 5

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
    !> Order in cache_prim: [ n_e, Te, Ti, q, Er ]  (matches kim_aperfuns).
    !> Order in cache_geom: [ B0, ks, kp, om_E ]     (matches kim_geom_aperfuns).
    logical, save :: cache_ready = .false.
    real(dp), allocatable, save :: cache_r(:)
    real(dp), allocatable, save :: cache_prim(:, :)   ! (grid_size, n_primitives)
    real(dp), allocatable, save :: cache_geom(:, :)   ! (grid_size, n_geom)

contains

    !> Clear the cached true background so the NEXT build_periodic_plasma call
    !> re-captures it from the current global plasma. Call this whenever the
    !> global profiles change (e.g. a new time step / a fresh set_profiles_from_
    !> arrays) so the periodic window tracks the updated background.
    subroutine reset_true_background_cache()
        cache_ready = .false.
        if (allocated(cache_r))    deallocate(cache_r)
        if (allocated(cache_prim)) deallocate(cache_prim)
        if (allocated(cache_geom)) deallocate(cache_geom)
    end subroutine reset_true_background_cache

    !> Capture the TRUE background (r-grid, primitives, geometry) from the global
    !> plasma into the module cache, ONCE. No-op if already cached. Must be called
    !> while `plasma` still holds the true global-grid background (i.e. at the top
    !> of build_periodic_plasma, before the redirect).
    subroutine capture_true_background()
        use species_m, only: plasma

        integer :: gs, i

        if (cache_ready) return

        gs = plasma%grid_size
        allocate(cache_r(gs), cache_prim(gs, n_primitives), cache_geom(gs, n_geom))

        cache_r = plasma%r_grid(1:gs)
        do i = 1, gs
            cache_prim(i, 1) = plasma%spec(0)%n(i)   ! electron density
            cache_prim(i, 2) = plasma%spec(0)%T(i)   ! electron temperature
            cache_prim(i, 3) = plasma%spec(1)%T(i)   ! ion temperature
            cache_prim(i, 4) = plasma%q(i)           ! safety factor
            cache_prim(i, 5) = plasma%Er(i)          ! radial electric field

            cache_geom(i, 1) = plasma%B0(i)          ! equilibrium field magnitude
            cache_geom(i, 2) = plasma%ks(i)          ! perpendicular wavenumber
            cache_geom(i, 3) = plasma%kp(i)          ! parallel wavenumber
            cache_geom(i, 4) = plasma%om_E(i)        ! ExB rotation frequency
        end do

        cache_ready = .true.
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

        ! Lazy capture: if this callback is invoked directly (e.g. a unit test
        ! exercising kim_aperfuns / sample_periodic_primitives without going
        ! through build_periodic_plasma), populate the cache from the current
        ! global plasma on first use. MUST happen before cache_prim is passed as
        ! an argument (passing an unallocated allocatable to an assumed-shape
        ! dummy is illegal). build_periodic_plasma captures explicitly at its top
        ! before the redirect, so this is a no-op in the normal solver path.
        if (.not. cache_ready) call capture_true_background()

        call interp_cache_row(cache_prim, nfuns, x, funs)
    end subroutine kim_aperfuns

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

        ! Lazy capture before passing cache_geom (see kim_aperfuns note).
        if (.not. cache_ready) call capture_true_background()

        call interp_cache_row(cache_geom, nfuns, x, funs)
    end subroutine kim_geom_aperfuns

    !> Return the PERIODIZED primitives at query radius x, driving make_periodic
    !> with kim_aperfuns. Layout: half_period = dx_asis + dx_tr, period
    !> L = 2*half_period, an as-is core of half-width dx_asis centred on x_mid,
    !> flanked by transition bands of width dx_tr. Output order matches
    !> kim_aperfuns: funs_per = [ n_e, Te, Ti, q, Er ].
    subroutine sample_periodic_primitives(x_mid, dx_asis, dx_tr, x, funs_per)
        real(dp), intent(in) :: x_mid, dx_asis, dx_tr, x
        real(dp), intent(out) :: funs_per(n_primitives)

        call make_periodic(kim_aperfuns, n_primitives, x_mid, dx_asis, dx_tr, x, funs_per)
    end subroutine sample_periodic_primitives

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
    subroutine build_periodic_plasma(rm, dx_asis, dx_tr, n_rg)
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

        real(dp) :: L, r_lo, r_hi, u0_seed
        real(dp) :: funs(n_primitives)
        real(dp), allocatable :: r_win(:)
        real(dp), allocatable :: n_win(:), Te_win(:), Ti_win(:), q_win(:), Er_win(:)
        integer :: j, sigma, npts

        ! Capture the TRUE background (r-grid, primitives, geometry) from the
        ! global plasma ONCE, before any redirect destroys it; all periodization
        ! (primitives AND geometry) then reads this cache, so repeated calls stay
        ! correct even though `plasma` is redirected onto the window after call 1.
        call capture_true_background()

        ! ---- 1. Equidistant window grid, reused as the global rg_grid. --------
        L = 2.0_dp * (dx_asis + dx_tr)
        r_lo = rm - 0.5_dp * L
        r_hi = rm + 0.5_dp * L

        ! grid_init_equidistant takes npts as intent(inout); pass a variable.
        npts = n_rg
        call rg_grid%grid_init_equidistant(npts, r_lo, r_hi, 'rg')
        call rg_grid%grid_generate_equidistant()

        ! ---- 2. Sample the PERIODIZED primitives on the window FIRST. ----------
        ! kim_aperfuns reads the CACHED true primitives (captured above), so this
        ! is robust to the redirect and to repeated calls. Sample into temporaries;
        ! write them onto the redirected plasma in step 3.
        ! FIXED order from sample_periodic_primitives: [ n_e, Te, Ti, q, Er ].
        allocate(r_win(rg_grid%npts_b), n_win(rg_grid%npts_b), Te_win(rg_grid%npts_b), &
                 Ti_win(rg_grid%npts_b), q_win(rg_grid%npts_b), Er_win(rg_grid%npts_b))
        r_win = rg_grid%xb
        do j = 1, rg_grid%npts_b
            call sample_periodic_primitives(rm, dx_asis, dx_tr, r_win(j), funs)
            n_win(j)  = funs(1)
            Te_win(j) = funs(2)
            Ti_win(j) = funs(3)
            q_win(j)  = funs(4)
            Er_win(j) = funs(5)
        end do

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
        plasma%q  = q_win
        plasma%Er = Er_win
        plasma%spec(0)%n = n_win
        plasma%spec(0)%T = Te_win
        do sigma = 1, number_of_ion_species
            call reallocate(plasma%spec(sigma)%n, rg_grid%npts_b)
            call reallocate(plasma%spec(sigma)%T, rg_grid%npts_b)
            plasma%spec(sigma)%T = Ti_win
        end do
        call set_ion_density_quasineutral()

        ! Force calculate_equil / calc_plasma_parameter_derivs to recompute the
        ! profile derivatives (dndr, dTdr, dqdr) on the window instead of reusing
        ! stale global-grid values.
        do sigma = 0, number_of_ion_species
            if (allocated(plasma%spec(sigma)%dndr)) deallocate(plasma%spec(sigma)%dndr)
            if (allocated(plasma%spec(sigma)%dTdr)) deallocate(plasma%spec(sigma)%dTdr)
        end do
        if (allocated(plasma%dqdr)) deallocate(plasma%dqdr)

        deallocate(r_win, n_win, Te_win, Ti_win, q_win, Er_win)

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

        subroutine set_ion_density_quasineutral()
            integer :: total_Z, sp, jj

            total_Z = 0
            do sp = 1, number_of_ion_species
                total_Z = total_Z + plasma%spec(sp)%Zspec
            end do
            if (total_Z == 0) return
            do sp = 1, number_of_ion_species
                do jj = 1, rg_grid%npts_b
                    plasma%spec(sp)%n(jj) = plasma%spec(0)%n(jj) * plasma%spec(sp)%Zspec / total_Z
                end do
            end do
        end subroutine set_ion_density_quasineutral

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
