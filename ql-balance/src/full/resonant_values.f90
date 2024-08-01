module resonantValues

    implicit none


    contains

    subroutine interpBrAndDqlAtResonance(i)

        use time_evolution, only: br_abs, dqle22_res_time, br_abs_time, br_abs_antenna_factor, time
        use PolyLagrangeInterpolation
        use grid_mod, only: npoib, r_resonant, rb, dqle22
        use wave_code_data, only: antenna_factor, Br

        implicit none

        integer, intent(in) :: i
        integer :: indResRadius
        
        call binsrc(rb, 1, npoib, r_resonant(1), indResRadius)
        call getIndicesForLagrangeInterp(indResRadius)
        call plag_coeff(nlagr, nder, r_resonant(1), rb(indBeginInterp:indEndInterp), coef)

        br_abs(i) = sum(coef(0, :)*abs(Br(indBeginInterp:indEndInterp)))*sqrt(antenna_factor)
        dqle22_res_time(i) = sum(coef(0, :)*dqle22(indBeginInterp:indEndInterp))

        ! save the time for the improved stopping criterion
        br_abs_time(i) = time
		br_abs_antenna_factor(i) = antenna_factor


        write(*,*) 'Br abs res * C_mn= ', br_abs(i)
        write(*,*) 'Br abs res       = ', br_abs(i)/sqrt(antenna_factor)
        write(*,*) 'Dqle22 res       = ', dqle22_res_time(i)
        write(*,*) 'Antenna factor   = ', antenna_factor
        write(*,*) 'time = ', br_abs_time(i)


    end subroutine


end module