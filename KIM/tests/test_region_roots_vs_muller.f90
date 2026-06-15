program test_region_roots_vs_muller

    ! Equivalence check for the ZEAL replacement (Z1): the fortnum
    ! complex_region_roots finder must return the same dispersion zeros as the
    ! in-tree Muller solver (run_Muller_dispersion uses the same roots() driver).
    ! Both solvers run on one analytic test function with known complex zeros
    ! inside the search rectangle; the root sets must agree to ~1e-6.

    use KIM_kinds_m, only: dp
    use muller_root_finding, only: roots
    use fortnum_roots_complex, only: complex_region_roots
    use fortnum_status, only: fortnum_status_t, FORTNUM_OK

    implicit none

    ! Search box matches the ZEAL/Muller usage: a rectangle in the complex
    ! plane centered near the cluster of dispersion roots.
    complex(dp), parameter :: ll = (-3.0_dp, -3.0_dp)
    complex(dp), parameter :: ur = ( 3.0_dp,  3.0_dp)
    real(dp),    parameter :: tol = 1.0e-6_dp

    ! Exact zeros built into ftest below.
    complex(dp), parameter :: exact(3) = [ &
        (1.5_dp,  0.5_dp), &
        (-2.0_dp, 1.0_dp), &
        (0.25_dp, -1.75_dp)]

    ! fortnum region finder
    complex(dp), allocatable :: rg_roots(:), rg_fvals(:)
    integer,     allocatable :: rg_mult(:)
    integer :: rg_nfound
    type(fortnum_status_t) :: rstatus

    ! Muller deflation finder
    integer :: known, nmore, km, guess, maxits, m_nfound, ifail
    logical :: realrt
    real(dp) :: ep1, ep2
    complex(dp), allocatable :: m_rts(:), m_fnv(:)

    integer :: i

    ! --- fortnum region root finder over [ll, ur] ---
    call complex_region_roots(ftest, ll, ur, rg_roots, rg_fvals, rg_mult, &
        rg_nfound, rstatus, m_max=5)
    if (rstatus%code /= FORTNUM_OK) then
        print *, 'complex_region_roots failed: ', trim(rstatus%msg)
        error stop 1
    end if
    if (rg_nfound /= size(exact)) then
        print *, 'region finder found ', rg_nfound, ' zeros, expected ', size(exact)
        error stop 2
    end if
    do i = 1, rg_nfound
        if (.not. matches_exact(rg_roots(i))) then
            print *, 'region root not in exact set: ', rg_roots(i)
            error stop 3
        end if
    end do

    ! --- Muller deflation finder on the same function ---
    km = 10
    allocate(m_rts(km), m_fnv(km))
    m_rts = (0.0_dp, 0.0_dp)
    m_fnv = (0.0_dp, 0.0_dp)
    known = 0
    nmore = size(exact)
    guess = 0
    maxits = 10000
    ep1 = 0.5e-6_dp
    ep2 = 1.0e-7_dp
    realrt = .false.

    call roots(ftest_scalar, known, nmore, km, realrt, ep1, ep2, &
        guess, maxits, m_rts, m_nfound, m_fnv, ifail)
    if (ifail /= 0) then
        print *, 'Muller roots() failed, ifail = ', ifail
        error stop 4
    end if
    if (m_nfound /= size(exact)) then
        print *, 'Muller found ', m_nfound, ' zeros, expected ', size(exact)
        error stop 5
    end if

    ! --- cross-check: every region root has a Muller partner within tol ---
    do i = 1, rg_nfound
        if (.not. has_partner(rg_roots(i), m_rts(1:m_nfound))) then
            print *, 'region root has no Muller partner within tol: ', rg_roots(i)
            error stop 6
        end if
    end do

    print *, 'Region/Muller equivalence OK: ', rg_nfound, ' shared zeros within', tol
    do i = 1, rg_nfound
        print *, '  zero ', i, ' = ', rg_roots(i), ' mult ', rg_mult(i)
    end do

contains

    ! Analytic function with three simple complex zeros, no poles or branch
    ! cuts in [ll, ur], so the argument-principle finder and Muller agree.
    subroutine ftest(kr, fk, ctx)
        complex(dp), intent(in)            :: kr
        complex(dp), intent(out)           :: fk
        class(*),    intent(in), optional  :: ctx
        fk = (kr - exact(1)) * (kr - exact(2)) * (kr - exact(3))
    end subroutine ftest

    subroutine ftest_scalar(kr, f)
        complex(dp), intent(in)  :: kr
        complex(dp), intent(out) :: f
        call ftest(kr, f)
    end subroutine ftest_scalar

    logical function matches_exact(z)
        complex(dp), intent(in) :: z
        integer :: k
        matches_exact = .false.
        do k = 1, size(exact)
            if (abs(z - exact(k)) < tol) matches_exact = .true.
        end do
    end function matches_exact

    logical function has_partner(z, set)
        complex(dp), intent(in) :: z
        complex(dp), intent(in) :: set(:)
        integer :: k
        has_partner = .false.
        do k = 1, size(set)
            if (abs(z - set(k)) < tol) has_partner = .true.
        end do
    end function has_partner

end program test_region_roots_vs_muller
