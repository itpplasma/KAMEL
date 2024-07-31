module linear_run


    contains

    subroutine allocate_prev_variables

        use time_evolution
        use recstep_mod, only: timstep_arr, tim_stack
        use grid_mod, only: neqset, npoib

        implicit none

        allocate (yprev(neqset))
        allocate (dqle11_prev(npoib))
        allocate (dqle12_prev(npoib))
        allocate (dqle21_prev(npoib))
        allocate (dqle22_prev(npoib))
        allocate (dqli11_prev(npoib))
        allocate (dqli12_prev(npoib))
        allocate (dqli21_prev(npoib))
        allocate (dqli22_prev(npoib))
        allocate (timstep_arr(neqset), tim_stack(neqset))

    end subroutine

    subroutine init_background_profiles

        use grid_mod, only: npoic, qsaf, params
        use wave_code_data, only: q, n, Vz, Te, Ti
        use baseparam_mod, only: ev, rtor


        implicit none    
        integer :: ipoi                    !initial background profiles:

        do ipoi = 1, npoic
            !safety factor:
            qsaf(ipoi) = 0.5*(q(ipoi) + q(ipoi + 1))
            !electron density :
            params(1, ipoi) = 0.5*(n(ipoi) + n(ipoi + 1))
            !toroidal rotation frequency :
            params(2, ipoi) = 0.5*(Vz(ipoi) + Vz(ipoi + 1))/rtor
            !electron temeperature :
            params(3, ipoi) = 0.5*(Te(ipoi) + Te(ipoi + 1))*ev
            !ion temeperature :
            params(4, ipoi) = 0.5*(Ti(ipoi) + Ti(ipoi + 1))*ev
        end do

    end subroutine

end module