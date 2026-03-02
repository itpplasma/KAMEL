subroutine prepare_resonances

    use kim_resonances_m
    use config_m, only: hdf5_output
    use setup_m, only: m_mode, n_mode, type_br_field
    use species_m, only: plasma
    use KIM_kinds_m, only: dp

    implicit none

    integer :: j
    real(dp) :: qres,qmin,qmax
    real(dp), dimension(:), allocatable :: q
    integer :: lb, ub

    iunit_res=157

    allocate(q(plasma%grid_size))

    q = abs(plasma%q)
    qmin = minval(q)
    qmax = maxval(q)

    qres = abs(dfloat(m_mode)/dfloat(n_mode))
    if(qres.lt.qmin.or.qres.gt.qmax) then
        write(*,*) "Resonance location not found in q"
        r_res = 0.0d0
        return
    end if

    r_res = qres

    do j= 2, plasma%grid_size
        if(qres .gt. q(j-1) .and. qres .le. q(j)) then
            r_res = (plasma%r_grid(j-1) * (q(j) - qres) + plasma%r_grid(j) * (qres-q(j-1))) / (q(j)-q(j-1))
            exit
        endif
    enddo

    if (type_br_field == 2) then
        r_res = plasma%r_grid(plasma%grid_size)/2
    end if

    write(*,*) 'resonant radius: ',r_res
    deallocate(q)

end subroutine prepare_resonances
