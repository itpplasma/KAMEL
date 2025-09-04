subroutine kim_init_for_test

    use KIM_kinds_m, only: dp
    use plasma_parameter, only: r_prof, iprof_length, n_prof, Te_prof, Ti_prof, &
        Er_prof, q_prof, ni_prof
    use equilibrium_m, only: B0z, B0th, B0, hz, hth
    use setup_m, only: btor, R0

    implicit none

    integer :: ierr = 0
    integer :: i
    real(dp) :: h

    call read_config

    iprof_length = 100

    if (.not. allocated(r_prof)) allocate(r_prof(iprof_length), stat=ierr)
    if (ierr /= 0) print *, "array: Allocation request denied"

    h = (70.0d0-3.0d0) / (iprof_length-1)
    r_prof(1) = 3.0d0
    do i=2, iprof_length
        r_prof(i) = r_prof(i-1) + h
    end do
        
    allocate(n_prof(iprof_length), Te_prof(iprof_length), Ti_prof(1, iprof_length), &
        Er_prof(iprof_length), q_prof(iprof_length), ni_prof(1, iprof_length))
    allocate(B0z(iprof_length), B0th(iprof_length), B0(iprof_length), hz(iprof_length), hth(iprof_length))

    n_prof = 2d13
    ni_prof = 2d13
    Te_prof = 1d3
    Ti_prof = 1d3
    Er_prof = 0d0
    q_prof = 1.5d0

    ! calculate equilibrium B field and J
    B0z = btor
    B0th = B0z * r_prof /(q_prof * R0)
    B0 = sqrt(B0th**2d0 + B0z**2d0)

    hz = B0z / B0
    hth = B0th / B0
    ! calculate quantities used for the kernels, e.g. A1, A2, dndr, omega_c,...
    call allocate_backs
    call calculate_backs(.false.)

end subroutine