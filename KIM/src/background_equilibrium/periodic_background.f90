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

    public :: kim_aperfuns, sample_periodic_primitives

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

end module periodic_background_m
