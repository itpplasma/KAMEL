module singleStep

    use balanceBase, only: balance_t

    double precision :: dqle22_res_single, br_abs_res_single

    type, extends(balance_t) :: singleStep_t
        contains
            procedure :: initBalance => initSingleStep
            procedure :: runBalance => runSingleStep
    end type

    contains
    
    subroutine initSingleStep(this)

        use balance_mod, only: balanceInit

        implicit none

        class(SingleStep_t), intent(inout) :: this
        this%runType = "SingleStep"
        print *, "Initialize Single Step run"

        call balanceInit

    end subroutine

    subroutine runSingleStep(this)

        use control_mod, only: irf

        implicit none

        class(SingleStep_t), intent(inout) :: this

        print *, "Running SingleStep"

        irf = 2 ! initialize transport coefficients
        call get_dql
        irf = 1 ! calculate transport coefficients
        call get_dql
        call rescale_transp_coeffs_by_ant_fac

        call interpBrAndDqlAtResonance

    end subroutine

    subroutine interpBrAndDqlAtResonance

        use PolyLagrangeInterpolation
        use wave_code_data, only: Br, antenna_factor
        use grid_mod, only: npoib, rb, r_resonant, dqle22
        
        implicit none

        if (.not. allocated(coef)) allocate (coef(0:nder, nlagr))
        call binsrc(rb, 1, npoib, r_resonant(1), indResRadius)
        call getIndicesForLagrangeInterp(indResRadius)
        call plag_coeff(nlagr, nder, r_resonant(1), rb(indBeginInterp:indEndInterp), coef)
        
        dqle22_res_single = sum(coef(0, :) * dqle22(indBeginInterp:indEndInterp))
        br_abs_res_single = sum(coef(0, :) * abs(Br(indBeginInterp:indEndInterp)))*sqrt(antenna_factor)

        write(*,*) ""
        write(*,*) "Dqle22 res = ", dqle22_res_single
        write(*,*) "|Br| res   = ", br_abs_res_single
		write(*,*) "Antenna factor = ", antenna_factor
        write(*,*) ""

    end subroutine

    subroutine init_background_profiles

        use grid_mod, only: npoic
        use plasma_parameters, only: params, qsaf
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

    subroutine finalizeSingleStepRun

        use h5mod
        use mpi
        use parallelTools, only: ierror
        use control_mod, only: ihdf5IO
        use paramscan_mod, only: writeDqle22

        implicit none

        
        if (ihdf5IO .eq. 1) then
            call writeDqle22
        end if
        call MPI_finalize(ierror);
        stop '-> Finished linear run |'

    end subroutine


    subroutine rescale_transp_coeffs_by_ant_fac

        use grid_mod, only: dqle11, dqle12, dqle21, dqle22, &
                            dqli11, dqli12, dqli21, dqli22
        use wave_code_data, only: antenna_factor

        implicit none

        dqle11 = dqle11*antenna_factor
        dqle12 = dqle12*antenna_factor
        dqle21 = dqle21*antenna_factor
        dqle22 = dqle22*antenna_factor
        dqli11 = dqli11*antenna_factor
        dqli12 = dqli12*antenna_factor
        dqli21 = dqli21*antenna_factor
        dqli22 = dqli22*antenna_factor

    end subroutine

end module