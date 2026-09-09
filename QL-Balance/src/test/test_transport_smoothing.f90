program test_transport_smoothing
    use QLBalance_kinds, only: dp
    use transport_smoothing_m, only: smooth_transport_profile
    implicit none

    integer, parameter :: n = 101
    real(dp) :: input(n), periodic(n), legacy(n), expected(n)

    input = 0.0_dp
    input(49:53) = [1.0_dp, 4.0_dp, 9.0_dp, 4.0_dp, 1.0_dp]
    periodic = input
    call smooth_transport_profile(periodic, .true.)
    if (any(periodic /= input)) then
        error stop 'periodic transport changed its trusted core or compact support'
    end if

    legacy = input
    call smooth_array_gauss(n, 30, input, expected)
    call smooth_transport_profile(legacy, .false.)
    if (any(legacy /= expected)) error stop 'legacy transport smoothing changed'
    if (legacy(48) <= 0.0_dp) error stop 'test must exercise support spreading by legacy filter'
    print *, 'transport smoothing policy tests passed'
end program test_transport_smoothing
