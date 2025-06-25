subroutine recnsplit(r,recnsp)

    use resonances_mod, only: r_res, width_res, ampl_res, prop
    use KIM_kinds, only: dp

    implicit none;

    real(dp) :: r, recnsp;

    if(prop) then
        prop=.false.
        call prepare_resonances
    endif

    recnsp = 1.0d0 +  ampl_res * exp(-((r - r_res) / width_res)**2.0d0)

end subroutine recnsplit