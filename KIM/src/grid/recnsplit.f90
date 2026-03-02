subroutine recnsplit(r,recnsp)

    use kim_resonances_m, only: r_res, prop
    use grid_m, only: width_res, ampl_res
    use KIM_kinds_m, only: dp

    implicit none;

    real(dp) :: r, recnsp;

    if(prop) then
        prop=.false.
        call prepare_resonances
    endif

    recnsp = 1.0d0 +  ampl_res * exp(-((r - r_res) / width_res)**2.0d0)

end subroutine recnsplit
