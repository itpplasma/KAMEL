module kim_diagnostics_m

    use KIM_kinds_m, only: dp

    implicit none

    private
    public :: integrate_Ipar
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

    subroutine compute_and_write_diagnostics()
        ! Evaluate D_ql,e22 at the resonant surface r_res itself (fields
        ! and plasma parameters interpolated to r_res; the dqle22 layer
        ! can be narrower than the grid spacing, so the nearest node is
        ! not good enough), |Br(r_res)|, and the integrated parallel
        ! currents (total and electron-only) from the electromagnetic
        ! solution in EBdat, then hand them to the IO layer. No-op when
        ! both the HDF5 and the flat-file diagnostics outputs are
        ! disabled.
        use, intrinsic :: iso_fortran_env, only: error_unit
        use config_m, only: hdf5_output, write_diagnostics_dat
        use fields_m, only: EBdat
        use kim_resonances_m, only: r_res, prop
        use kim_qldiff_m, only: calc_dqle22
        use IO_collection_m, only: write_kim_diagnostics

        real(dp) :: vTe, nue, om_E, B0, kpar, dqle22, br_abs_res
        complex(dp) :: Ipar, Ipar_e, Es_res, Br_res
        integer :: npts

        if (.not. (hdf5_output .or. write_diagnostics_dat)) return

        ! Equidistant grids never call recnsplit, so trigger the lazy
        ! resonance detection here (same guard as recnsplit).
        if (prop) then
            prop = .false.
            call prepare_resonances
        end if

        npts = size(EBdat%r_grid)

        ! prepare_resonances leaves r_res = 0 when |m/n| is outside the
        ! q range. Skip all diagnostics output then: a missing
        ! kim_diagnostics.dat is the scan driver's failure signal,
        ! whereas values evaluated at a bogus radius would be silently
        ! machine-read. Non-fatal on purpose (logger_m's log_error
        ! aborts the run, and the field solution itself is still valid).
        if (r_res <= 0.0_dp .or. r_res < EBdat%r_grid(1) &
            .or. r_res > EBdat%r_grid(npts)) then
            write(error_unit, '(A, ES12.5, A)') &
                '[ERROR] kim_diagnostics: no resonant surface inside the ' // &
                'radial grid (r_res = ', r_res, '); diagnostics output skipped'
            return
        end if

        ! kpar interpolated at r_res is tiny but generically nonzero;
        ! calc_dqle22's collisional response is regular in the limit
        ! kpar -> 0, so it is used as-is.
        call interp_local_plasma(r_res, vTe, nue, om_E, B0, kpar)

        Es_res = interp_local_complex(npts, EBdat%r_grid, EBdat%Es, r_res)
        Br_res = interp_local_complex(npts, EBdat%r_grid, EBdat%Br, r_res)
        br_abs_res = abs(Br_res)

        dqle22 = calc_dqle22(vTe, nue, om_E, B0, kpar, Es_res, Br_res)
        Ipar = integrate_Ipar(npts, EBdat%r_grid, EBdat%jpar)
        Ipar_e = integrate_Ipar(npts, EBdat%r_grid, EBdat%jpar_e)

        call write_kim_diagnostics(dqle22, Ipar, Ipar_e, br_abs_res)

    end subroutine

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
