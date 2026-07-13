program test_fp_susceptibility_oracle
    use KIM_kinds_m, only: dp

    implicit none

    integer, parameter :: nmmax = 3
    real(dp), parameter :: map_tol = 5.0e-14_dp
    real(dp), parameter :: oracle_tol = 1.0e-12_dp
    complex(dp) :: symbI(0:nmmax, 0:nmmax), expected
    real(dp) :: kpar, temperature_erg, mass, nu, omega_E, Vpar
    real(dp) :: omega_c, omega, x1_ref, x2_ref, expected_re, expected_im
    real(dp) :: vT, x1, x2, error_scale, scaled_error, max_scaled_error
    integer :: case_id, mphi, k, l, iunit, ios, nrows
    character(len=1024) :: line
    logical :: seen_pair(0:nmmax, 0:nmmax), seen_harmonic(-1:1)

    interface
        subroutine getIfunc(x1, x2, symbI)
            double precision, intent(in) :: x1, x2
            double complex, intent(out) :: symbI(0:3, 0:3)
        end subroutine getIfunc
    end interface

    seen_pair = .false.
    seen_harmonic = .false.
    nrows = 0
    max_scaled_error = 0.0_dp
    open(newunit=iunit, file='fp_susceptibilities.dat', status='old', &
         action='read', iostat=ios)
    if (ios /= 0) error stop 'cannot open fp_susceptibilities.dat'

    do
        read(iunit, '(A)', iostat=ios) line
        if (ios < 0) exit
        if (ios > 0) error stop 'cannot read susceptibility oracle line'
        if (len_trim(line) == 0 .or. line(1:1) == '#') cycle
        read(line, *, iostat=ios) case_id, kpar, temperature_erg, mass, nu, &
            omega_E, Vpar, mphi, omega_c, omega, k, l, x1_ref, x2_ref, &
            expected_re, expected_im
        if (ios /= 0) error stop 'malformed susceptibility oracle row'
        if (k < 0 .or. k > nmmax .or. l < 0 .or. l > nmmax) then
            error stop 'oracle susceptibility index out of range'
        end if
        if (mphi < -1 .or. mphi > 1) error stop 'oracle harmonic out of range'

        vT = sqrt(temperature_erg/mass)
        x1 = kpar*vT/nu
        x2 = -(omega_E + kpar*Vpar + real(mphi, dp)*omega_c - omega)/nu
        call assert_close('x1 mapping', x1, x1_ref, map_tol, case_id, k, l)
        call assert_close('x2 mapping', x2, x2_ref, map_tol, case_id, k, l)

        call getIfunc(x1, x2, symbI)
        expected = cmplx(expected_re, expected_im, kind=dp)
        error_scale = max(1.0_dp, abs(expected))
        scaled_error = abs(symbI(k, l) - expected)/error_scale
        max_scaled_error = max(max_scaled_error, scaled_error)
        if (scaled_error > oracle_tol) then
            print *, 'FAIL: susceptibility oracle mismatch'
            print *, ' case, (k,l), mphi = ', case_id, k, l, mphi
            print *, ' x1, x2 = ', x1, x2
            print *, ' actual   = ', symbI(k, l)
            print *, ' expected = ', expected
            print *, ' scaled error = ', scaled_error
            error stop
        end if

        nrows = nrows + 1
        seen_pair(k, l) = .true.
        seen_harmonic(mphi) = .true.
    end do
    close(iunit)

    if (nrows < 56) error stop 'susceptibility oracle is missing regime coverage'
    call require_pair(0, 0)
    call require_pair(2, 0)
    call require_pair(0, 2)
    call require_pair(0, 1)
    call require_pair(2, 1)
    call require_pair(2, 2)
    call require_pair(1, 1)
    call require_pair(1, 3)
    if (.not. all(seen_harmonic)) error stop 'oracle lacks mphi=-1,0,+1'

    print *, 'PASS: ', nrows, ' independent FP susceptibility oracle rows'
    print *, 'Maximum scaled oracle error: ', max_scaled_error
    print *, 'All tests PASSED'
    stop 0

contains

    subroutine assert_close(label, actual, reference, tolerance, icase, ik, il)
        character(len=*), intent(in) :: label
        real(dp), intent(in) :: actual, reference, tolerance
        integer, intent(in) :: icase, ik, il

        if (abs(actual - reference) > tolerance*max(1.0_dp, abs(reference))) then
            print *, 'FAIL: ', trim(label), ' case/(k,l) = ', icase, ik, il
            print *, ' actual/reference = ', actual, reference
            error stop
        end if
    end subroutine assert_close

    subroutine require_pair(ik, il)
        integer, intent(in) :: ik, il
        if (.not. seen_pair(ik, il)) then
            print *, 'FAIL: missing susceptibility pair ', ik, il
            error stop
        end if
    end subroutine require_pair

end program test_fp_susceptibility_oracle
