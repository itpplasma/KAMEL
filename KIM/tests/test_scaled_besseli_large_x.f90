program test_scaled_besseli_large_x
    ! The large-argument expansion of exp(-x) I_n(x) used by the periodic ion
    ! tensor must agree with fortnum's independent complex Miller evaluation
    ! to rounding for every argument above the 500 threshold.
    use KIM_kinds_m, only: dp
    use quasilinear_flr_m, only: scaled_besseli_large_x
    use fortnum_special_complex_bessel, only: bessel_i_complex_array
    use fortnum_status, only: fortnum_status_t, FORTNUM_OK
    implicit none

    integer, parameter :: nmax = 8
    real(dp), parameter :: arguments(7) = &
        [500.0001_dp, 512.3_dp, 750.0_dp, 1.0e3_dp, 4.7e3_dp, 1.0e5_dp, 1.0e6_dp]
    complex(dp) :: reference(nmax + 1)
    type(fortnum_status_t) :: status
    real(dp) :: got, want, worst
    integer :: i, n

    worst = 0.0_dp
    do i = 1, size(arguments)
        call bessel_i_complex_array(0, nmax + 1, cmplx(arguments(i), 0.0_dp, dp), &
            .true., reference, status)
        if (status%code /= FORTNUM_OK) error stop 'fortnum reference evaluation failed'
        do n = 0, nmax
            got = scaled_besseli_large_x(n, arguments(i))
            want = real(reference(n + 1), dp)
            worst = max(worst, abs(got - want) / abs(want))
            if (abs(got - want) > 1.0e-13_dp * abs(want)) then
                print *, 'FAIL: n =', n, ' x =', arguments(i)
                print *, '  got  =', got
                print *, '  want =', want
                error stop 'large-argument scaled Bessel I disagrees with fortnum'
            end if
        end do
    end do
    print '(a,es10.2)', 'PASS: exp(-x) I_n(x), n <= 8, x > 500; max rel diff ', worst
end program test_scaled_besseli_large_x
