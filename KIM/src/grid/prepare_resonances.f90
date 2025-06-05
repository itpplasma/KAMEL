subroutine prepare_resonances

    use resonances_mod
    use grid, only: gg_width, gg_factor, grid_spacing
    use config, only: hdf5_output, fdebug
    use setup, only: m_mode, n_mode, type_br_field
    use plasma_parameter, only: iprof_length, r_prof, q_prof
    !use h5mod
    !use mpi

    implicit none

    integer :: j
    double precision :: qres,qmin,qmax
    double precision, dimension(:), allocatable :: q
    integer :: lb, ub

    iunit_res=157

    allocate(q(iprof_length))

    q = abs(q_prof)
    qmin = minval(q)
    qmax = maxval(q)

    !width_res = gg_width
    !ampl_res = gg_factor

    qres = abs(dfloat(m_mode)/dfloat(n_mode))
    if(qres.lt.qmin.or.qres.gt.qmax) write(*,*) "Resonance location not found in q"

    r_res = qres

    do j=2,iprof_length
      if(qres .gt. q(j-1) .and. qres .lt. q(j)) then
        r_res = (r_prof(j-1) * (q(j) - qres) + r_prof(j) * (qres-q(j-1))) / (q(j)-q(j-1))
        exit
      endif
    enddo

     if (type_br_field == 2) then
        r_res = r_prof(iprof_length)/2
    end if   

    write(*,*) 'resonant radius: ',r_res
    deallocate(q)

end subroutine prepare_resonances


