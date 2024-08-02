module linear_run


    contains

    
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

    subroutine finalizeLinearRun

        use h5mod
        use mpi
        use parallelTools, only: ierror
        use control_mod, only: ihdf5IO

        implicit none

        write(*,*) '-> Finalize linear run |'
        if (ihdf5IO .eq. 1) then
            call writeDqle22
        end if
        call MPI_finalize(ierror);
        stop  !! <<----- Stop for linear code usage

    end subroutine

end module