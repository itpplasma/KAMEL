module kim_diagnostics_m

    use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
    use KIM_kinds_m, only: dp

    implicit none

    private
    public :: integrate_Ipar
    public :: integrate_Ipar_window
    public :: interp_local_complex
    public :: compute_and_write_diagnostics

contains

    pure function integrate_Ipar(npts, r, jpar) result(Ipar)
        ! Integrated parallel current I_par = 2*pi * int(r' * jpar(r') dr')
        ! over the full passed grid (complex trapezoid rule). The caller is
        ! expected to pass the restricted layer domain.
        integer, intent(in) :: npts
        real(dp), intent(in) :: r(npts)
        complex(dp), intent(in) :: jpar(npts)
        complex(dp) :: Ipar

        real(dp), parameter :: pi = 4.0_dp * atan(1.0_dp)
        integer :: i

        Ipar = (0.0_dp, 0.0_dp)
        do i = 1, npts - 1
            Ipar = Ipar + (r(i) * jpar(i) + r(i+1) * jpar(i+1)) &
                          * (r(i+1) - r(i)) / 2.0_dp
        end do
        Ipar = 2.0_dp * pi * Ipar

    end function

    pure function integrate_Ipar_window(npts, r, jpar, r_mid, half_width) result(Ipar)
        ! Integrate over an exact window. Grid values at the two endpoints are
        ! obtained by linear interpolation; integrate_Ipar remains full-grid.
        integer, intent(in) :: npts
        real(dp), intent(in) :: r(npts), r_mid, half_width
        complex(dp), intent(in) :: jpar(npts)
        complex(dp) :: Ipar

        real(dp), parameter :: pi = 4.0_dp * atan(1.0_dp)
        real(dp) :: lower, upper, left, right
        complex(dp) :: jl, jr
        integer :: i

        Ipar = (0.0_dp, 0.0_dp)
        if (npts < 2 .or. half_width <= 0.0_dp) return
        lower = r_mid - half_width
        upper = r_mid + half_width
        if (lower < r(1) .or. upper > r(npts) .or. upper <= lower) return
        do i = 1, npts - 1
            left = max(lower, r(i))
            right = min(upper, r(i + 1))
            if (right <= left) cycle
            jl = linear_complex_value(r(i), jpar(i), r(i + 1), jpar(i + 1), left)
            jr = linear_complex_value(r(i), jpar(i), r(i + 1), jpar(i + 1), right)
            Ipar = Ipar + (left * jl + right * jr) * (right - left) / 2.0_dp
        end do
        Ipar = 2.0_dp * pi * Ipar
    end function integrate_Ipar_window

    pure function linear_complex_value(x1, y1, x2, y2, x) result(y)
        real(dp), intent(in) :: x1, x2, x
        complex(dp), intent(in) :: y1, y2
        complex(dp) :: y

        if (x2 == x1) then
            y = y1
        else
            y = y1 + (y2 - y1) * (x - x1) / (x2 - x1)
        end if
    end function linear_complex_value

    subroutine compute_and_write_diagnostics(require_window)
        ! Evaluate D_ql,e22 at the resonant surface r_res itself (fields
        ! and plasma parameters interpolated to r_res; the dqle22 layer
        ! can be narrower than the grid spacing, so the nearest node is
        ! not good enough), |Br(r_res)|, and the integrated parallel
        ! currents (total and electron-only) from the electromagnetic
        ! solution in EBdat, then hand them to the IO layer. No-op when
        ! both the HDF5 and the flat-file diagnostics outputs are
        ! disabled.
        use config_m, only: hdf5_output, write_diagnostics_dat
        use fields_m, only: EBdat
        use kim_resonances_m, only: r_res, prop
        use kim_qldiff_m, only: calc_dqle22
        use IO_collection_m, only: write_kim_diagnostics

        real(dp) :: vTe, nue, om_E, B0, kpar, dqle22, br_abs_res
        real(dp) :: rm, dx_asis
        complex(dp) :: Ipar, Ipar_e, Es_res, Br_res
        integer :: npts
        logical :: use_window, packed_window
        logical, intent(in), optional :: require_window

        packed_window = .false.
        if (present(require_window)) packed_window = require_window

        if (.not. (hdf5_output .or. write_diagnostics_dat)) return
        ! Equidistant grids never call kim_recnsplit, so trigger the lazy
        ! resonance detection here (same guard as kim_recnsplit).
        if (prop) then
            prop = .false.
            call kim_prepare_resonances
        end if

        if (.not. allocated(EBdat%r_grid) .or. .not. allocated(EBdat%Es) .or. &
            .not. allocated(EBdat%Br) .or. .not. allocated(EBdat%jpar) .or. &
            .not. allocated(EBdat%jpar_e)) return
        if (.not. allocated(EBdat%Phi) .or. .not. allocated(EBdat%jpar_i)) return
        npts = size(EBdat%r_grid)
        if (npts < 4 .or. size(EBdat%Es) /= npts .or. size(EBdat%Br) /= npts .or. &
            size(EBdat%jpar) /= npts .or. size(EBdat%jpar_e) /= npts) return
        if (size(EBdat%Phi) /= npts .or. size(EBdat%jpar_i) /= npts) return

        if (packed_window) then
            ! Periodic and FLR2 outputs are valid only with an explicit,
            ! finite as-is window. Never reinterpret bad metadata as a
            ! legacy full-grid request.
            if (.not. ieee_is_finite(EBdat%r_resonance) .or. &
                .not. ieee_is_finite(EBdat%dx_asis) .or. &
                EBdat%r_resonance <= 0.0_dp .or. EBdat%dx_asis <= 0.0_dp) return
            rm = EBdat%r_resonance
            dx_asis = EBdat%dx_asis
        else
            rm = EBdat%r_resonance
            if (rm <= 0.0_dp) rm = r_res
        end if

        ! kim_prepare_resonances leaves r_res = 0 when |m/n| is outside the
        ! q range. Skip all diagnostics output then: a missing
        ! kim_diagnostics.dat is the scan driver's failure signal,
        ! whereas values evaluated at a bogus radius would be silently
        ! machine-read. Non-fatal on purpose (logger_m's log_error
        ! aborts the run, and the field solution itself is still valid).
        if (.not. ieee_is_finite(rm) .or. rm <= 0.0_dp .or. &
            rm < EBdat%r_grid(1) .or. rm > EBdat%r_grid(npts)) return
        if (.not. all(ieee_is_finite(EBdat%r_grid)) .or. &
            .not. all(ieee_is_finite(real(EBdat%Es, dp))) .or. &
            .not. all(ieee_is_finite(aimag(EBdat%Es))) .or. &
            .not. all(ieee_is_finite(real(EBdat%Br, dp))) .or. &
            .not. all(ieee_is_finite(aimag(EBdat%Br))) .or. &
            .not. all(ieee_is_finite(real(EBdat%jpar, dp))) .or. &
            .not. all(ieee_is_finite(aimag(EBdat%jpar))) .or. &
            .not. all(ieee_is_finite(real(EBdat%jpar_e, dp))) .or. &
            .not. all(ieee_is_finite(aimag(EBdat%jpar_e)))) return
        if (.not. all(ieee_is_finite(real(EBdat%Phi, dp))) .or. &
            .not. all(ieee_is_finite(aimag(EBdat%Phi))) .or. &
            .not. all(ieee_is_finite(real(EBdat%jpar_i, dp))) .or. &
            .not. all(ieee_is_finite(aimag(EBdat%jpar_i)))) return

        ! kpar interpolated at r_res is tiny but generically nonzero;
        ! calc_dqle22's collisional response is regular in the limit
        ! kpar -> 0, so it is used as-is.
        call interp_local_plasma(rm, vTe, nue, om_E, B0, kpar)

        Es_res = interp_local_complex(npts, EBdat%r_grid, EBdat%Es, rm)
        Br_res = interp_local_complex(npts, EBdat%r_grid, EBdat%Br, rm)
        br_abs_res = abs(Br_res)

        dqle22 = calc_dqle22(vTe, nue, om_E, B0, kpar, Es_res, Br_res)
        use_window = packed_window .or. (EBdat%r_resonance > 0.0_dp .and. &
            EBdat%dx_asis > 0.0_dp)
        if (use_window) then
            if (.not. packed_window) dx_asis = EBdat%dx_asis
            if (.not. valid_window(npts, EBdat%r_grid, rm, dx_asis)) return
            Ipar = integrate_Ipar_window(npts, EBdat%r_grid, EBdat%jpar, rm, dx_asis)
            Ipar_e = integrate_Ipar_window(npts, EBdat%r_grid, EBdat%jpar_e, rm, dx_asis)
        else
            Ipar = integrate_Ipar(npts, EBdat%r_grid, EBdat%jpar)
            Ipar_e = integrate_Ipar(npts, EBdat%r_grid, EBdat%jpar_e)
        end if
        if (.not. ieee_is_finite(dqle22) .or. .not. ieee_is_finite(br_abs_res) .or. &
            .not. ieee_is_finite(real(Ipar, dp)) .or. .not. ieee_is_finite(aimag(Ipar)) .or. &
            .not. ieee_is_finite(real(Ipar_e, dp)) .or. .not. ieee_is_finite(aimag(Ipar_e))) return

        call write_kim_diagnostics(dqle22, Ipar, Ipar_e, br_abs_res)

    end subroutine

    logical function valid_window(npts, r, r_mid, half_width)
        integer, intent(in) :: npts
        real(dp), intent(in) :: r(npts), r_mid, half_width

        valid_window = .false.
        if (npts < 2) return
        valid_window = ieee_is_finite(r_mid) .and. ieee_is_finite(half_width) &
                       .and. all(ieee_is_finite(r)) .and. half_width > 0.0_dp &
                       .and. all(r(2:npts) > r(1:npts - 1)) &
                       .and. r_mid - half_width >= r(1) &
                       .and. r_mid + half_width <= r(npts)
    end function valid_window

    subroutine interp_local_plasma(r_pt, vTe, nue, om_E, B0, kpar)
        ! Interpolate the local electron plasma parameters from the
        ! plasma profile grid to r_pt with the 4-point Lagrange stencil
        ! used throughout KIM (binsrc + plag_coeff).
        use species_m, only: plasma

        real(dp), intent(in) :: r_pt
        real(dp), intent(out) :: vTe, nue, om_E, B0, kpar

        integer, parameter :: nlagr = 4, nder = 0
        real(dp) :: coef(0:nder, nlagr)
        integer :: ir, ibeg, iend

        call binsrc(plasma%r_grid, 1, plasma%grid_size, r_pt, ir)
        ibeg = max(1, ir - nlagr / 2)
        iend = ibeg + nlagr - 1
        if (iend > plasma%grid_size) then
            iend = plasma%grid_size
            ibeg = iend - nlagr + 1
        end if

        call plag_coeff(nlagr, nder, r_pt, plasma%r_grid(ibeg:iend), coef)

        vTe = sum(coef(0, :) * plasma%spec(0)%vT(ibeg:iend))
        nue = sum(coef(0, :) * plasma%spec(0)%nu(ibeg:iend))
        om_E = sum(coef(0, :) * plasma%om_E(ibeg:iend))
        B0 = sum(coef(0, :) * plasma%B0(ibeg:iend))
        kpar = sum(coef(0, :) * plasma%kp(ibeg:iend))

    end subroutine

    function interp_local_complex(npts, r_grid, f, r_pt) result(f_pt)
        ! Interpolate the complex profile f from r_grid to r_pt with the
        ! same 4-point Lagrange stencil (binsrc + plag_coeff) as
        ! interp_local_plasma; the real weights act on real and
        ! imaginary parts alike.
        integer, intent(in) :: npts
        real(dp), intent(in) :: r_grid(npts), r_pt
        complex(dp), intent(in) :: f(npts)
        complex(dp) :: f_pt

        integer, parameter :: nlagr = 4, nder = 0
        real(dp) :: coef(0:nder, nlagr)
        integer :: ir, ibeg, iend

        call binsrc(r_grid, 1, npts, r_pt, ir)
        ibeg = max(1, ir - nlagr / 2)
        iend = ibeg + nlagr - 1
        if (iend > npts) then
            iend = npts
            ibeg = iend - nlagr + 1
        end if

        call plag_coeff(nlagr, nder, r_pt, r_grid(ibeg:iend), coef)

        f_pt = sum(coef(0, :) * f(ibeg:iend))

    end function

end module
