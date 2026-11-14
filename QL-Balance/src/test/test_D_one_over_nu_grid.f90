program test_D_one_over_nu_grid
    !
    ! The 1/nu coefficients are added to the quasilinear ones on the cell
    ! boundaries, so they must exist on all npoib boundaries. In a uniform
    ! plasma every boundary must carry the same value, and rescaling to a
    ! plasma with twice the temperature and three times the density must
    ! change them by the 1/nu factor T/n used in the code.
    !
    use QLBalance_kinds, only: dp
    use grid_mod, only: npoib, npoic, ipbeg, ipend, reint_coef, &
        Donue11, Donue22, Donui11, Donui22, init_Donue11, init_Donui22
    use plasma_parameters, only: params, init_params
    use baseparam_mod, only: btor, rtor, ev

    implicit none

    real(dp), parameter :: tolerance = 1.0e-12_dp
    integer :: ipoi
    logical :: failed

    failed = .false.
    npoic = 4
    npoib = npoic + 1
    btor = -2.75e4_dp
    rtor = 360.0_dp

    ! two-point averages inside, one-sided values at both ends
    allocate (ipbeg(npoib), ipend(npoib), reint_coef(2, npoib))
    do ipoi = 1, npoib
        ipbeg(ipoi) = max(1, min(ipoi - 1, npoic - 1))
        ipend(ipoi) = ipbeg(ipoi) + 1
        reint_coef(:, ipoi) = 0.5_dp
    end do
    reint_coef(:, 1) = [1.0_dp, 0.0_dp]
    reint_coef(:, npoib) = [0.0_dp, 1.0_dp]

    allocate (params(4, npoic), init_params(4, npoic))
    params(1, :) = 1.0e13_dp
    params(2, :) = 0.0_dp
    params(3, :) = 1.8e3_dp*ev
    params(4, :) = 1.0e3_dp*ev
    init_params = params

    call initialize_D_one_over_nu

    call check(size(Donue11) == npoib .and. size(Donui22) == npoib, &
        "coefficients cover all cell boundaries")
    if (failed) stop 1
    call check(all(Donue11 > 0.0_dp) .and. all(Donui22 > 0.0_dp), "coefficients are positive")
    call check(all(abs(Donue11 - Donue11(1)) <= tolerance*Donue11(1)) .and. &
        all(abs(Donui22 - Donui22(1)) <= tolerance*Donui22(1)), &
        "uniform plasma gives the same value on every boundary")

    params(1, :) = 3.0_dp*init_params(1, :)
    params(3, :) = 2.0_dp*init_params(3, :)
    params(4, :) = 2.0_dp*init_params(4, :)
    call rescale_D_one_over_nu_to_new_n_and_T
    call check(all(abs(Donue11 - init_Donue11*2.0_dp/3.0_dp) <= tolerance*init_Donue11) .and. &
        all(abs(Donui22 - init_Donui22*2.0_dp/3.0_dp) <= tolerance*init_Donui22), &
        "rescaling follows T/n on every boundary")

    if (failed) stop 1
    print *, "ALL TESTS PASSED"

contains

    subroutine check(ok, what)
        logical, intent(in) :: ok
        character(*), intent(in) :: what
        if (ok) then
            print *, "PASS: ", what
        else
            print *, "FAIL: ", what
            failed = .true.
        end if
    end subroutine

end program
