program test_quadpack_qagi
    ! Sanity test for QUADPACK DQAGI (infinite intervals).
    ! Integrates exp(-x) from 0 to +inf; exact result = 1.

    use KIM_kinds_m, only: dp

    implicit none

    external :: dqagi
    real(dp) :: result, abserr
    integer :: neval, ier, last
    integer :: limit, lenw
    integer, allocatable :: iwork(:)
    real(dp), allocatable :: work(:)
    real(dp) :: epsabs, epsrel
    real(dp) :: bound
    integer :: inf

    ! Use contained function below; no explicit interface needed here

    epsabs = 1.0d-10
    epsrel = 1.0d-10
    bound = 0.0d0
    inf = 1  ! integrate (bound, +inf)
    limit = 200
    lenw = 4 * limit
    allocate(iwork(limit), work(lenw))

    call dqagi(f_exp_minus_x, bound, inf, epsabs, epsrel, result, abserr, neval, ier, &
               limit, lenw, last, iwork, work)

    if (ier /= 0) then
        print *, 'DQAGI returned error code: ', ier
        stop 1
    end if

    print *, 'Integral exp(-x) from 0..inf: ', result, ' abs err: ', abserr

    if (abs(result - 1.0d0) > 1.0d-8) then
        print *, 'ERROR: result deviates from 1 beyond tolerance'
        stop 1
    else
        print *, 'DQAGI test passed.'
    end if

contains

    function f_exp_minus_x(x) result(val)
        use KIM_kinds_m, only: dp
        real(dp), intent(in) :: x
        real(dp) :: val
        val = exp(-x)
    end function f_exp_minus_x

end program test_quadpack_qagi
