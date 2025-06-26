subroutine prepare_resonances

    use resonances_mod
    use grid, only: gg_width, gg_factor, grid_spacing
    use config, only: hdf5_output, fdebug
    use setup, only: m_mode, n_mode, type_br_field
    use species, only: plasma
    use KIM_kinds, only: dp

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

    !width_res = gg_width
    !ampl_res = gg_factor

    qres = abs(dfloat(m_mode)/dfloat(n_mode))
    if(qres.lt.qmin.or.qres.gt.qmax) write(*,*) "Resonance location not found in q"

    r_res = qres

    do j= 2, plasma%grid_size
      if(qres .gt. q(j-1) .and. qres .lt. q(j)) then
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


