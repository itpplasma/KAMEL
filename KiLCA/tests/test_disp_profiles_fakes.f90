!> Minimal domain parameters for this isolated native-interface test.
module constants
    use iso_c_binding, only: c_double, c_double_complex, c_intptr_t
    implicit none
    integer, parameter :: dp = c_double, dpc = c_double_complex, pp = c_intptr_t
end module constants

module flre_sett
    implicit none
    integer :: nwaves = 2, flre_order = 1
end module flre_sett

module disp_profiles_test_state
    use, intrinsic :: iso_c_binding, only: c_int
    implicit none
    integer(c_int) :: failures = 0, saves = 0
contains
    integer(c_int) function get_disp_test_failures()
        if (saves /= 1) failures = failures + 1
        get_disp_test_failures = failures
    end function
end module disp_profiles_test_state

subroutine calc_dispersion(r, flagback, flagprint, kval, polvec)
    use constants, only: dp, dpc
    use flre_sett, only: nwaves
    real(dp), intent(in) :: r
    character(*), intent(in) :: flagback
    integer, intent(in) :: flagprint
    complex(dpc), intent(out) :: kval(nwaves), polvec(nwaves, nwaves)
    integer :: idx, wave

    ! r holds the 0-based grid index (the test sets x(i) = i-1).
    idx = nint(r)
    kval = (0.0_dp, 0.0_dp)
    polvec = (0.0_dp, 0.0_dp)
    do wave = 1, nwaves
        kval(wave) = cmplx(1000.0_dp * wave + idx, 2000.0_dp * wave + idx, dpc)
    end do
    polvec(1, 1) = cmplx(3000.0_dp + idx, 0.0_dp, dpc)
end subroutine

module kilca_inout_m
    implicit none
contains
    function save_cmplx_matrix_to_one_file_(Nrows, Ncols, Npoints, xgrid, arr, &
                                            full_name) &
        result(ierr)
        use, intrinsic :: iso_c_binding, only: c_int, c_double
        use disp_profiles_test_state, only: failures, saves
        integer(c_int), value :: Nrows, Ncols, Npoints
        real(c_double), intent(in) :: xgrid(0:Npoints - 1)
        real(c_double), intent(in) :: arr(0:2 * Nrows * Ncols * Npoints - 1)
        character(len=*), intent(in) :: full_name
        integer(c_int) :: ierr
        integer(c_int) :: i, dimk, wave

        saves = saves + 1
        if (full_name /= 'out.d') failures = failures + 1
        dimk = 2 * Nrows
        do i = 0, Npoints - 1
            do wave = 1, Nrows
                if (abs(arr(dimk * i + 2 * (wave - 1)) - (1000.0d0 * wave + i)) > 1.0d-12) &
                    failures = failures + 1
                if (abs(arr(dimk * i + 2 * (wave - 1) + 1) - (2000.0d0 * wave + i)) > 1.0d-12) &
                    failures = failures + 1
            end do
        end do
        ierr = 0
    end function

end module kilca_inout_m
