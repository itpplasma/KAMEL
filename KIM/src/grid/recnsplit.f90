subroutine recnsplit(r,recnsp)

    use resonances_mod
!  use grid_mod, only: gg_width, gg_factor, gg_r_res;
!
    implicit none;
!
    logical :: prop=.true.
    double precision :: r, recnsp;
!
    if(prop) then
        prop=.false.
        call prepare_resonances
    endif
!
!  recnsp = 1.d0 + gg_factor*exp(-((r-gg_r_res)/gg_width)**2);
    recnsp = 1.d0
  !do k=1,numres
    recnsp = recnsp +  ampl_res * exp(-((r - r_res) / width_res)**2)
  !recnsp = recnsp + ampl_res(k)*exp(-((r-r_res(k))/width_res(k))**2)
  !enddo
!
    return
end subroutine recnsplit