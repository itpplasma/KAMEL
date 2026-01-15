program test_root_finding

    use muller_root_finding
    use KIM_kinds_m, only: dp

    implicit none

    integer :: known, nmore, km, guess, maxits, nfound, ifail
    logical :: realrt
    real(dp) :: ep1, ep2
    complex(dp), allocatable :: rts(:), fnv(:)
    integer :: i, j

    km = 5
    allocate(rts(km), fnv(km))
    rts = (0.0_dp, 0.0_dp)
    fnv = (0.0_dp, 0.0_dp)

    known = 0
    nmore = 2
    guess = 0
    maxits = 1000
    ep1 = 1.0d-6
    ep2 = 1.0d-8
    realrt = .false.

    call roots(calcF, known, nmore, km, realrt, ep1, ep2, &
            guess, maxits, rts, nfound, fnv, ifail)

    if (ifail /= 0) then
        print *, 'Root finding failed with IFAIL = ', ifail
        error stop
    end if

    if (rts(1) /= (1.0_dp, 0.0_dp)) then
        print *, 'Root is wrong.'
        error stop
    end if

    if (rts(2) /= (0.0_dp, 1.0_dp)) then
        print *, 'Root is wrong.'
        error stop
    end if

    print *, 'Roots found: ', nfound

    do i=1, nfound
        print *, 'Root ', i, ' = ', rts(i), ' f = ', fnv(i)
    end do

    contains

    subroutine calcF(kr, f)

        implicit none

        complex(dp), intent(in)  :: kr
        complex(dp), intent(out) :: f

        ! Example function with roots at kr = 1.0 + 0.0i and kr = 0.0 + 1.0i
        f = (kr - (1.0_dp, 0.0_dp)) * (kr - (0.0_dp, 1.0_dp))

    end subroutine

end program