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

    !> Number of primitive functions sampled, and their FIXED order.
    integer, parameter :: n_primitives = 5

contains

    !> Aperiodic sample callback conforming to periodization_m::aperfuns_i.
    !> Evaluates the TRUE primitive profiles of the global `plasma` at radius x
    !> via 4-point Lagrange interpolation, in the FIXED order
    !>     funs = [ n_e(x), Te(x), Ti(x), q(x), Er(x) ].
    !>
    !> Clamping: x is clamped to [r_grid(1), r_grid(grid_size)] before
    !> interpolation, returning the endpoint value for out-of-range x. This is
    !> required because make_periodic queries x_inperiod +- L, which can fall
    !> outside the raw profile domain. In the deep as-is / periodic-image regions
    !> the far images are weighted 0 by the localizer, so clamping is harmless
    !> there; it only matters inside the transition region.
    subroutine kim_aperfuns(nfuns, x, funs)
        use species_m, only: plasma

        integer, intent(in) :: nfuns
        real(dp), intent(in) :: x
        real(dp), intent(out) :: funs(nfuns)

        integer, parameter :: nlagr = 4, nder = 0
        real(dp) :: coef(0:nder, nlagr)
        real(dp) :: xc
        integer :: gs, ir, ibeg, iend

        gs = plasma%grid_size

        ! Clamp query point to the raw profile domain (see header note).
        xc = max(plasma%r_grid(1), min(x, plasma%r_grid(gs)))

        ! 4-point Lagrange stencil selection (mirrors calculate_equil.f90).
        call binsrc(plasma%r_grid, 1, gs, xc, ir)
        ibeg = max(1, ir - nlagr/2)
        iend = ibeg + nlagr - 1
        if (iend > gs) then
            iend = gs
            ibeg = iend - nlagr + 1
        end if

        call plag_coeff(nlagr, nder, xc, plasma%r_grid(ibeg:iend), coef)

        ! FIXED order: [ n_e, Te, Ti, q, Er ] (contract with Phase-2.3).
        funs(1) = sum(coef(0, :) * plasma%spec(0)%n(ibeg:iend))   ! electron density
        funs(2) = sum(coef(0, :) * plasma%spec(0)%T(ibeg:iend))   ! electron temperature
        funs(3) = sum(coef(0, :) * plasma%spec(1)%T(ibeg:iend))   ! ion temperature
        funs(4) = sum(coef(0, :) * plasma%q(ibeg:iend))           ! safety factor
        funs(5) = sum(coef(0, :) * plasma%Er(ibeg:iend))          ! radial electric field
    end subroutine kim_aperfuns

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
        use species_m, only: plasma, set_plasma_quantities, reallocate, &
                             deallocate_plasma_derived
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

        ! ---- 1. Equidistant window grid, reused as the global rg_grid. --------
        L = 2.0_dp * (dx_asis + dx_tr)
        r_lo = rm - 0.5_dp * L
        r_hi = rm + 0.5_dp * L

        ! grid_init_equidistant takes npts as intent(inout); pass a variable.
        npts = n_rg
        call rg_grid%grid_init_equidistant(npts, r_lo, r_hi, 'rg')
        call rg_grid%grid_generate_equidistant()

        ! ---- 2. Sample the PERIODIZED primitives on the window FIRST. ----------
        ! kim_aperfuns reads the TRUE (global-grid) primitives of `plasma`, so
        ! this MUST happen before we redirect `plasma` onto the window (which
        ! reallocates and thereby destroys those true profiles). Sample into
        ! temporaries; write them onto the redirected plasma in step 3.
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

        ! Capture the TRUE B0 at the window's left edge (still on the global grid
        ! here) as the force-balance ODE seed u0 = B0(r_lo)**2. Without this the
        ! window solve restarts from the pressureless vacuum initial condition
        ! and loses the pressure work accumulated over [r_min, r_lo], offsetting
        ! the absolute B0 level (hence omega_c, rho_L) by ~O(0.5%). See
        ! calculate_equil's u0_seed argument.
        u0_seed = interp_global_B0(r_lo)**2

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

        ! set_plasma_quantities: calculate_plasma_backs (vT, nu, omega_c, rho_L,
        ! lambda_D, z0), interpolate_plasma_backs(rg_grid%xb) -- the IDENTITY
        ! here because plasma%r_grid == rg_grid%xb -- cell centres,
        ! calculate_thermodynamic_forces_and_susc (A1, A2, x1, x2, I00..I21) on
        ! rg_grid%npts_b == the window, and write_species_backs /
        ! write_plasma_backs (writes background profile .dat files to
        ! output_path//backs/).
        call set_plasma_quantities(plasma)

    contains

        !> 4-point Lagrange interpolation of the global plasma%B0 at radius x,
        !> using the same binsrc + plag_coeff stencil as the rest of KIM. Called
        !> while `plasma` still holds the TRUE global-grid B0 (before redirect).
        real(dp) function interp_global_B0(x) result(b0x)
            real(dp), intent(in) :: x
            integer, parameter :: nlagr = 4, nder = 0
            real(dp) :: coef(0:nder, nlagr)
            integer :: gs, ir, ibeg, iend

            gs = plasma%grid_size
            call binsrc(plasma%r_grid, 1, gs, x, ir)
            ibeg = max(1, ir - nlagr/2)
            iend = ibeg + nlagr - 1
            if (iend > gs) then
                iend = gs
                ibeg = iend - nlagr + 1
            end if
            call plag_coeff(nlagr, nder, x, plasma%r_grid(ibeg:iend), coef)
            b0x = sum(coef(0, :) * plasma%B0(ibeg:iend))
        end function interp_global_B0

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
