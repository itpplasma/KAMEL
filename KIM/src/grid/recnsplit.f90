subroutine recnsplit(r,recnsp)

    use resonances_mod, only: r_res, width_res, ampl_res

    implicit none;

    logical :: prop=.true.
    double precision :: r, recnsp;

    if(prop) then
        prop=.false.
        call prepare_resonances
    endif

    r_res = 35.0d0

    recnsp = 1.0d0 +  ampl_res * exp(-((r - r_res) / width_res)**2.0d0)
end subroutine recnsplit