  module time_evolution

    use control_mod
    use parallelTools
    use h5mod
    use balanceBase, only: balance_t

    implicit none

    logical :: flag_run_time_evolution !Added by Philipp Ulbl 12.05.2020
    logical :: br_stopping ! trigger Br stopping criterion
    logical :: discr_reached = .false. ! variable to say if discrepancy to linear regression
    logical :: scratch

    integer :: Nstorage
    integer :: ramp_up_mode !> control ramp up mode of the RMP coil current amplitude
    integer :: save_prof_time_step ! added by Markus Markl 11.03.2021
    integer :: iexit ! used for ramp-up skipping of saving
    integer :: ramp_up_down = 0 !> used in hysteresis mode, tells if ramp-up (0) or ramp-down (1)
    integer :: timeIndex, timescale

    double precision :: tmax
    DOUBLE PRECISION :: br_beta = 0
    DOUBLE PRECISION :: br_predicted

    double precision :: tmax_factor!, antenna_factor
    double precision :: stop_time_step !Added by Philipp Ulbl 13.05.2020
    double precision :: timstep_min
    DOUBLE PRECISION :: t_max_ramp_up = 1e-2 !> 10ms ramp up until antenna_factor_max is reached
    double precision :: timstep
    double precision :: time
    double precision :: t_hysteresis_turn = 0

    double precision, dimension(:), allocatable :: timscal

    DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE :: yprev
    DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE :: dqle11_prev
    DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE :: dqle12_prev
    DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE :: dqle21_prev
    DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE :: dqle22_prev
    DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE :: dqli11_prev
    DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE :: dqli12_prev
    DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE :: dqli21_prev
    DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE :: dqli22_prev

    double precision :: antenna_factor_max
    DOUBLE PRECISION :: antenna_max_stopping

    !needed for interpolation of br abs and stopping criterion
    DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE :: br_abs
    DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE :: br_abs_time
	DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE :: br_abs_antenna_factor
    DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE :: dqle22_res_time

    logical :: firstiterationdone = .false. !Some steps in saving the data to hdf5 file 
    !need to be done only the first time iteration

    integer(HID_T) :: time_dataset_id !> variable to save the time dataset id


    type, extends(balance_t) :: TimeEvolution_t
        contains
            procedure :: initBalance => initTimeEvolution
            procedure :: runBalance => runTimeEvolution
    end type

    contains

    subroutine initTimeEvolution(this)

        implicit none

        class(TimeEvolution_t), intent(inout) :: this
        this%runType = "TimeEvolution"

        call balanceInit

    end subroutine

    subroutine runTimeEvolution(this)

        use parallelTools, only: irank
        use baseparam_mod, only: factolmax, factolred, tol_max
        use recstep_mod, only: tol
        use plasma_parameters, only: params, params_beg, params_begbeg, limitTemperaturesFromBelow
        use restart_mod, only: redostep, inquiry_to_restart, redoTimeStep
        use recstep_mod, only: tim_stack, timstep_arr
        use transp_coeffs_mod, only: rescale_transp_coeffs_by_ant_fac

        implicit none

        integer :: timeIndex
        class(TimeEvolution_t), intent(inout) :: this

        write(*,*) "Running TimeEvolution"

        call calc_geometric_parameter_profiles
        call initialize_get_dql
        call initAntennaFactor
        call genstartsource

        if (irank .eq. 0) then
            call writeKinProfileDataToDisk(0)
        end if

        time = 0.0d0
        tol = tol_max

        call inquiry_to_restart

        timstep_arr = timstep
        tim_stack = timstep_arr

        call get_dql
        call rescale_transp_coeffs_by_ant_fac

        call allocateBrAndDqleForTimeEvolution
        call savePrevTranspCoefficients
        call allocate_timscal_and_params

        params_begbeg = params

        do timeIndex = 1, Nstorage
            print *, "TimeIndex: ", timeIndex
            call saveKinProfilesToYPrev
            redostep = .false.

            call get_dql
            call stopIfTimeStepTooSmall
            call interpBrAndDqlAtResonanceTimeEvol
            call write_br_time_data
            call rescale_transp_coeffs_by_ant_fac

            call writefort9999

            if (.not. redostep) then
                call savePrevTranspCoefficients
                params_begbeg = params
            else 
                call redoTimeStep
            end if

            do ! redo step loop
                params_beg = params
                call evolvestep(timstep, eps)
                call limitTemperaturesFromBelow
                call calcParamsNumAndDenom
                call smoothParamsNumAndDenom
                call determineTimscal
                if (maxval(timscal) .lt. tol * factolmax) exit
                timstep_arr = timstep_arr * factolred
                params = params_beg

                if (irank .eq. 0) then
                    print *, "Redoing step"
                end if
            end do

            call rescaleTimStepArr
            call setTimStep
            call resetTimStepArrWithTimstep
            call writeTimeInfoToDisk
            call relaxPlasmaParameters

            timstep_arr = 0.0d0
            call evolvestep(timstep, eps)
            timstep_arr = timstep
            time = time + timstep

            call messageTimeInfo
            call writeKinProfileAtTimeIndex
            call setFirstIterationTrue
            call checkIfLinearDiscrepancyOfPenRatioReached

            call ramp_coil(timeIndex)
        end do



    end subroutine


    subroutine allocate_prev_variables

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

    subroutine saveKinProfilesToYPrev

        use grid_mod, only: npoi, nbaleqs
        use plasma_parameters, only: params

        implicit none

        integer :: ipoi, ieq, k

        if (debug_mode) write(*,*) 'Debug: yprev loop'
        do ipoi = 1, npoi
            do ieq = 1, nbaleqs
                k = nbaleqs*(ipoi - 1) + ieq
                yprev(k) = params(ieq, ipoi)
            end do
        end do

    end subroutine

    subroutine initAntennaFactor
        !For time evolution mode use antenna_factor as maximum
        !and start with a very small value and ramp this up

        use wave_code_data, only: antenna_factor
        implicit none

        antenna_factor_max = antenna_factor
        if (ramp_up_mode .eq. 4) then 
            antenna_factor = 0d0
        else
            antenna_factor = 1.d-4
        end if

    end subroutine

    subroutine allocateBrAndDqleForTimeEvolution

        implicit none

        allocate(br_abs(Nstorage))
        allocate(br_abs_antenna_factor(Nstorage))
        allocate(br_abs_time(Nstorage))
        allocate(dqle22_res_time(Nstorage))

    end subroutine

    subroutine savePrevTranspCoefficients

        use grid_mod, only: dqle11, dqle12, dqle21, dqle22, &
                            dqli11, dqli12, dqli21, dqli22
        implicit none

        dqle11_prev = dqle11
        dqle12_prev = dqle12
        dqle21_prev = dqle21
        dqle22_prev = dqle22
        dqli11_prev = dqli11
        dqli12_prev = dqli12
        dqli21_prev = dqli21
        dqli22_prev = dqli22
    
    end subroutine

    subroutine allocate_timscal_and_params

        use grid_mod, only: npoi, npoic, nbaleqs, dummy
        use plasma_parameters, only: params_beg, params_num, params_denom, params_begbeg

        implicit none

        allocate(timscal(npoi), dummy(npoic))
        allocate(params_beg(nbaleqs, npoic), params_num(nbaleqs, npoic))
        allocate(params_denom(nbaleqs, npoic))
        allocate(params_begbeg(nbaleqs, npoic))

    end subroutine


    !> @brief subroutine write_br_time_data. Writes radial magnetic field perturbation evaluated at the resonant
    !> surface, the antenna factor, the time and Dqle22 evaluated at the resonant surface for a given 
    !> time step to the hdf5 file
    !> @param[in] i Integer of time step to which the data will be saved. Goes from 1:i.
    !> @param[in] br_abs_time Time value of the time evolution.
    !> @param[in] br_abs_antenna_factor Value of the antenna factor, i.e. the RMP coil current.
    !> @param[in] br_abs Absolute value of the radial magnetic field evaluated at the resonant surface in question.
    !> @param[in] dqle22_res_time Value of Dqle22 evaluated at the resonant surface during the time evolution.
    subroutine write_br_time_data

        use control_mod
        use baseparam_mod
        use h5mod
        use hdf5_tools
        use wave_code_data, only: antenna_factor

        implicit none

	    if (debug_mode) write(*,*) "Debug: writing out br time evolution data"

        if (ihdf5IO .eq. 1) then
       	    CALL h5_init()
            CALL h5_open_rw(path2out, h5_id)

 		    h5_currentgrp = "/"//trim(h5_mode_groupname) //"/br_abs_time"
		    CALL h5_add_double_1(h5_id, trim(h5_currentgrp), br_abs_time(1:timeIndex), &
			    lbound(br_abs_time(1:timeIndex)), ubound(br_abs_time(1:timeIndex)))

 		    h5_currentgrp = "/"//trim(h5_mode_groupname) //"/br_abs_antenna_factor"
		    CALL h5_add_double_1(h5_id, trim(h5_currentgrp), br_abs_antenna_factor(1:timeIndex), &
			    lbound(br_abs_antenna_factor(1:timeIndex)), ubound(br_abs_antenna_factor(1:timeIndex)))

 		    h5_currentgrp = "/"//trim(h5_mode_groupname) //"/br_abs_res"
		    CALL h5_add_double_1(h5_id, trim(h5_currentgrp), br_abs(1:timeIndex), &
			    lbound(br_abs(1:timeIndex)), ubound(br_abs(1:timeIndex)))

 		    h5_currentgrp = "/"//trim(h5_mode_groupname) //"/dqle22_res_time"
		    CALL h5_add_double_1(h5_id, trim(h5_currentgrp), dqle22_res_time(1:timeIndex), &
			    lbound(dqle22_res_time(1:timeIndex)), ubound(dqle22_res_time(1:timeIndex)))


            CALL h5_close(h5_id)
            CALL h5_deinit()

        else
            open (777, file='br_abs_res.dat', position='append')
            write (777, *) timeIndex, time, antenna_factor, br_abs(timeIndex)
            close (777)
        end if
    end subroutine ! write_br_time_data

    !> @brief Ramp up/down the RMP coil current.
    !> @author Markus Markl
    !> @date 13.03.2023
    subroutine ramp_coil(i)

        use paramscan_mod
        use control_mod
        use wave_code_data, only: antenna_factor
        use h5mod
        use parallelTools
        

        implicit none
        integer, intent(in) :: i
        !Ramp up antenna_factor: linear in Icoil or quadratic in D
        !Added by Philipp Ulbl 12.05.2020, (strongly) edited by Markus Markl
        !
        ! Stopping criterion that is always active
        if (antenna_factor .lt. (antenna_factor_max*antenna_max_stopping * 1.5)) then
            if (debug_mode) write(*,*) "Debug: - - ramp up antenna_factor - -"
        
            if (ramp_up_mode .eq. 0) then
                if (debug_mode) write(*,*) "Debug: ramp-up mode linear"
                ! initial ramp up. This is a linear ramp up of the antenna factor in time.
                antenna_factor = time**2 + 1.d-4
            else if (ramp_up_mode .eq. 1) then
                if (debug_mode) write(*,*) "Debug: ramp-up mode faster 1"
                ! First try for faster ramp up. Is (usually) not stable.
                antenna_factor = antenna_factor + antenna_factor_max * (timstep/t_max_ramp_up)
            else if (ramp_up_mode .eq. 2) then
                if (debug_mode) write(*,*) "Debug: ramp-up mode faster 2, stable"
                ! Use this for faster ramp up. Ramps up antenna factor to 100% of max value
                ! in t_max_ramp_up time, but ramps-up further after that.
                if (time .eq. 0) then
                    antenna_factor = 1.d-4
                else
                    antenna_factor = antenna_factor_max * (time/t_max_ramp_up)
                end if
            else if (ramp_up_mode .eq. 3) then ! ramp up instantly to 100% of max antenna factor
                if (debug_mode) write(*,*) "Debug: ramp-up mode instant"
                if (time .eq. 0) then
                    antenna_factor = 1.d-4
                else if (time .ge. 10*t_max_ramp_up) then ! if max time value is reached, stop the code
                    write(*,*) 'stop: reached time max: ', 10*t_max_ramp_up
                    if (suppression_mode .eqv. .false.) then
                        call writeKinProfileDataToDisk(i)
                    end if
                    ! Write the cause of the stopping into the hdf5 file
                    if (ihdf5IO .eq. 1) then
                        CALL h5_init()
                        CALL h5_open_rw(path2out, h5_id)
                        CALL h5_add_string(h5_id, trim(h5_mode_groupname)// &
                                '/stopping_criterion', 'reached time max')
                        CALL h5_close(h5_id)
                        CALL h5_deinit()
                    end if
                    if (paramscan) then
                        if (debug_mode) write(*,*) "Debug: Write br_time _data"
                        if (ihdf5IO .eq. 1) then
                            CALL write_br_time_data!, br_abs_time, br_abs_antenna_factor, br_abs, dqle22_res_time)
                        end if
                        if (ifac_n + ifac_Te + ifac_Ti + ifac_vz .eq. size(fac_n) + &
                                size(fac_Ti) + size(fac_Te) + size(fac_vz)) then
                            CALL MPI_finalize(ierror)
                            stop
                        else
                            ! if it is not the last scan, skip the rest of the
                            ! code and continue with the next loop iteration
                            !call deallocate_wave_code_data()
                            ! exit this time evolution loop specific to a certain set
                            ! of parameters
                            iexit = 1
                        end if
                    else
                        if (debug_mode) write(*,*) "Debug: Write br_time _data"
                        if (ihdf5IO .eq. 1) then
                            CALL write_br_time_data!, br_abs_time, br_abs_antenna_factor, br_abs, dqle22_res_time)
                        end if
                        CALL MPI_finalize(ierror);
                        write(*,*) 'stop'
                        stop
                    end if
 
                else if (antenna_factor .ge. antenna_factor_max) then
                    antenna_factor = antenna_factor_max
                    write(*,*) "- - - - - - - - - - - - - - - - - - - - - - - - - - - -"
                    write(*,*) "Antenna factor reached max value, will not be changed!"
                    write(*,*) "- - - - - - - - - - - - - - - - - - - - - - - - - - - -"
                else
                    !antenna_factor = antenna_factor_max * (time/(t_max_ramp_up*2))
                    antenna_factor = antenna_factor_max!(time/t_max_ramp_up)**2 + 1.d-4
                    !antenna_factor = antenna_factor + antenna_factor_max * (timstep/t_max_ramp_up)
                end if

            else if (ramp_up_mode .eq. 4) then ! no ramp-up at all
                if (debug_mode) write(*,*) "Debug: ramp-up mode none"
                if (time .eq. 0) then
                    antenna_factor = 0d0
                else if (time .ge. 10*t_max_ramp_up) then
                    ! if max time value is reached, stop the code
                    write(*,*) 'stop: time limit reached: ', time
                    if (suppression_mode .eqv. .false.) then
                        call writeKinProfileDataToDisk(i)
                    end if ! suppression mode
                    ! Write the cause of the stopping into the hdf5 file
                    if (ihdf5IO .eq. 1) then
                        CALL h5_init()
                        CALL h5_open_rw(path2out, h5_id)
                        CALL h5_add_string(h5_id, trim(h5_mode_groupname)// &
                                '/stopping_criterion', 'time limit reached')
                        CALL h5_close(h5_id)
                        CALL h5_deinit()
                    end if ! ihdf5IO
                    if (paramscan) then
                        if (debug_mode) write(*,*) "Debug: Write br_time _data"
                        if (ihdf5IO .eq. 1) then 
                            CALL write_br_time_data!, br_abs_time, br_abs_antenna_factor, br_abs, dqle22_res_time)
                        end if
                        if (ifac_n + ifac_Te + ifac_Ti + ifac_vz .eq. size(fac_n) + &
                            size(fac_Ti) + size(fac_Te) + size(fac_vz)) then
                            CALL MPI_finalize(ierror)
                            stop "Finished parameter scan"
                        else
                            ! if it is not the last scan, skip the rest of the
                            ! code and continue with the next loop iteration
                            !call deallocate_wave_code_data()
                            ! exit this time evolution loop specific to a certain set
                            ! of parameters
                            iexit = 1
                        end if
                    else
                        if (debug_mode) write(*,*) "Debug: Write br_time _data"
                        if (ihdf5IO .eq. 1) then
                            CALL write_br_time_data!, br_abs_time, br_abs_antenna_factor, br_abs, dqle22_res_time)
                        end if
                        write(*,*) 'stop'
                        call MPI_finalize(ierror);
                        stop
                    end if
  
                else
                    antenna_factor = 0d0!antenna_factor_max * (time/t_max_ramp_up)
                end if

            else if (ramp_up_mode .eq. 5) then ! hysteresis mode
                if (debug_mode) write(*,*) "Debug: ramp-up mode hysteresis"
                if (debug_mode) write(*,*) "Debug: ramp_up_down = ", ramp_up_down
                if (ramp_up_down .eq. 0) then ! ramp-up
                    if (debug_mode) write(*,*) "Ramp up ^"
                    antenna_factor = time**2 + 1.d-4
                    if (antenna_factor .ge. antenna_factor_max * antenna_max_stopping) then ! switch to ramp-down
                        ramp_up_down = 1
                        t_hysteresis_turn = time
                        if (debug_mode) write(*,*) "Debug: Ramp turn at t = ", time
                    end if
                else if (ramp_up_down .eq. 1) then !ramp down
                    if (debug_mode) write(*,*) "Debug: Ramp down v"
                    ! get linear ramp-up in coil current (scales with sqrt of antenna_factor):
                    antenna_factor = (sqrt(antenna_factor_max * antenna_max_stopping) - (time - t_hysteresis_turn))**2
                    if ((antenna_factor/antenna_factor_max * 100) .le. 1.0) then ! if antenna factor would get below some percentage of the max value
                        write(*,*) 'stop: ramp-up/down finished '
                        if (suppression_mode .eqv. .false.) then
                            call writeKinProfileDataToDisk(i)
                        end if
                        ! Write the cause of the stopping into the hdf5 file
                        if (ihdf5IO .eq. 1) then
                            CALL h5_init()
                            CALL h5_open_rw(path2out, h5_id)
                            CALL h5_add_string(h5_id, trim(h5_mode_groupname)// &
                                '/stopping_criterion', 'ramp-up/down finished')
                            CALL h5_close(h5_id)
                            CALL h5_deinit()
                        end if
                        if (debug_mode) write(*,*) "Debug: Write br_time _data"
                        if (ihdf5IO .eq. 1) then
                            CALL write_br_time_data!, br_abs_time, br_abs_antenna_factor, br_abs, dqle22_res_time)
                        end if
                        if (paramscan) then
                            if (ifac_n + ifac_Te + ifac_Ti + ifac_vz .eq. size(fac_n) + &
                                size(fac_Ti) + size(fac_Te) + size(fac_vz)) then
                                CALL MPI_finalize(ierror)
                                stop
                            else
                                iexit = 1
                            end if ! if last parameter scan
                        else
                            CALL MPI_finalize(ierror)
                            stop
                        end if ! if parameter scan
                    end if ! antenna_facotr less eq 0
                end if ! ramp_up_down eq 1

            else if (ramp_up_mode .eq. 6) then ! fast hysteresis mode
                if (debug_mode) write(*,*) "Debug: ramp-up mode fast hysteresis"
                if (debug_mode) write(*,*) "Debug: ramp_up_down = ", ramp_up_down
                if (ramp_up_down .eq. 0) then ! ramp-up
                    if (debug_mode) write(*,*) "Ramp up ^"
                    !antenna_factor = time**2 + 1.d-4
                    antenna_factor = (antenna_factor_max/t_max_ramp_up * time)**2
                    if (antenna_factor .ge. antenna_factor_max * antenna_max_stopping) then ! switch to ramp-down
                        ramp_up_down = 1
                        t_hysteresis_turn = time
                        if (debug_mode) write(*,*) "Debug: Ramp turn at t = ", time
                    end if
                else if (ramp_up_down .eq. 1) then !ramp down
                    if (debug_mode) write(*,*) "Debug: Ramp down v"
                    ! get linear ramp-up in coil current (scales with sqrt of antenna_factor):
                    antenna_factor = (sqrt(antenna_factor_max * antenna_max_stopping) - (time - t_hysteresis_turn) &
                    * antenna_factor_max/t_max_ramp_up)**2
                    if ((antenna_factor/antenna_factor_max * 100) .le. 0.1) then ! if antenna factor would get below some percentage of the max value
                        write(*,*) 'stop: ramp-up/down finished '
                        if (suppression_mode .eqv. .false.) then
                            call writeKinProfileDataToDisk(i)
                        end if
                        ! Write the cause of the stopping into the hdf5 file
                        if (ihdf5IO .eq. 1) then
                            CALL h5_init()
                            CALL h5_open_rw(path2out, h5_id)
                            CALL h5_add_string(h5_id, trim(h5_mode_groupname)// &
                                '/stopping_criterion', 'ramp-up/down finished')
                            CALL h5_close(h5_id)
                            CALL h5_deinit()
                        end if
                        if (debug_mode) write(*,*) "Debug: Write br_time _data"
                        if (ihdf5IO .eq. 1) then
                            CALL write_br_time_data!, br_abs_time, br_abs_antenna_factor, br_abs, dqle22_res_time)
                        end if
                        if (paramscan) then
                            if (ifac_n + ifac_Te + ifac_Ti + ifac_vz .eq. size(fac_n) + &
                                size(fac_Ti) + size(fac_Te) + size(fac_vz)) then
                                CALL MPI_finalize(ierror)
                                stop
                            else
                                iexit = 1
                            end if ! if last parameter scan
                        else
                            CALL MPI_finalize(ierror)
                            stop
                        end if ! if parameter scan
                    end if ! antenna_facotr less eq 0
                end if ! ramp_up_down eq 1

            end if !ramp_up_mode
                !This can be activated for runs without QL evolution to check steady state behaviour
                !antenna_factor = 1.d-4
                !if(i .gt. 200) then
                !    stop
                !endif
            write(*,*) 'antenna_factor = ', antenna_factor
            write(*,*) "which are ", antenna_factor/antenna_factor_max*100, "% of the max"
            write(*,*) " - - - "
        else
                write(*,*) 'stop: reached antenna_factor_max * ', antenna_max_stopping
                if (suppression_mode .eqv. .false.) then
                    call writeKinProfileDataToDisk(i)
                end if
                ! Write the cause of the stopping into the hdf5 file
                if (ihdf5IO .eq. 1) then
                    CALL h5_init()
                    CALL h5_open_rw(path2out, h5_id)
                    CALL h5_add_string(h5_id, trim(h5_mode_groupname)// &
                        '/stopping_criterion', 'reached antenna_factor_max * antenna_max_stopping')
                    CALL h5_close(h5_id)
                    CALL h5_deinit()
                end if

                if (debug_mode) write(*,*) "Debug: Write br_time_data"
                if (ihdf5IO .eq. 1) then 
                    CALL write_br_time_data!, br_abs_time, br_abs_antenna_factor, br_abs, dqle22_res_time)
                end if

                if (paramscan) then
                    if (ifac_n + ifac_Te + ifac_Ti + ifac_vz .eq. size(fac_n) + &
                        size(fac_Ti) + size(fac_Te) + size(fac_vz)) then
                        CALL MPI_finalize(ierror)
                        stop
                    else
                    ! if it is not the last scan, skip the rest of the
                    ! code and continue with the next loop iteration
                    !call deallocate_wave_code_data()
                    ! exit this time evolution loop specific to a certain set
                    ! of parameters
                        iexit = 1
                    end if
                else
                CALL MPI_finalize(ierror);
                write(*,*) 'stop'
                stop
            end if

            call MPI_finalize(ierror);
            stop
        end if

    end subroutine !ramp_coil


    !@> brief Check if discrepancy between penetration ratio and linear extrapolation exceeds critical value
    !> author> Markus Markl
    !> created 13.03.2023
    subroutine checkIfLinearDiscrepancyOfPenRatioReached

        use paramscan_mod, only: fac_n, fac_Te, fac_Ti, fac_vz, ifac_n, ifac_Te, &
                                ifac_Ti, ifac_vz
        implicit none
        ! Stopping by linear regression of Br_abs_res
        !
        ! Calculate br_beta from br_abs_res = br_beta * br_abs_time, which is
        ! essentially the slope of the curve. It is assumed that the curve
        ! intersects the y axis at 0. The slope is then used to calculate
        ! the "predicted" value for br_abs, which is compared to the actual
        ! value thereof. If the discrepancy becomes too big, the program is
        ! stopped.
        !
        if (timeIndex .gt. 50 .and. .not. discr_reached) then
            ! calculate beta only once
            if (br_beta .eq. 0) then
                ! Calculate slope from the data until this time step. (Simple linear regression algorithm)
                br_beta = sum(br_abs(3:timeIndex)*br_abs_time(3:timeIndex))/sum(br_abs_time(1:timeIndex)**2)
                write(*,*) 'br_beta = ', br_beta
            end if
            ! calculate the from the linear regression predicted value of Br_abs
            br_predicted = br_beta*br_abs_time(timeIndex)
            write(*,*) "Delta = ", abs(br_abs(timeIndex) - br_predicted)
            if (abs(br_abs(timeIndex) - br_predicted) .gt. 0.1) then
                write(*,*) 'discrepancy to linearly predicted value of Br_abs_res > delta'
                if (modulo(timeIndex, save_prof_time_step) .ne. 0) then
                    if (suppression_mode .eqv. .false.) then
                        CALL writeFieldsCurrentsAndTranspCoeffsToH5
                    end if
                end if
                if (suppression_mode .eqv. .false.) then
                    CALL writeKinProfileDataToDisk(timeIndex)
                end if
                if (br_stopping) then
                    ! Write the cause of the stopping into the hdf5 file
                    CALL h5_init()
                    CALL h5_open_rw(path2out, h5_id)
                    CALL h5_add_string(h5_id, trim(h5_mode_groupname)// &
                        '/stopping_criterion', 'discrepancy to linearly predicted value of Br_abs_res > delta')
                    CALL h5_add_double_0(h5_id, trim(h5_mode_groupname)//'/br_beta', br_beta, &
                        'linear slope of br', 'G/s')
                    CALL h5_close(h5_id)
                    CALL h5_deinit()

                    if (paramscan) then
                        if (ifac_n + ifac_Te + ifac_Ti + ifac_vz .eq. size(fac_n) + &
                            size(fac_Ti) + size(fac_Te) + size(fac_vz)) then
                            CALL MPI_finalize(ierror);
                            stop
                        else
                        ! if it is not the last scan, skip the rest of the
                        ! code and continue with the next loop iteration
                        !call deallocate_wave_code_data()
                        write(*,*) 'not last scan; will end time evolution of this parameter choice'
                        iexit = 1
                        end if
                    else
                        if (debug_mode) write(*,*) "Debug: Write br_time _data"
                        CALL write_br_time_data!, br_abs_time, br_abs_antenna_factor, br_abs, dqle22_res_time)
                        CALL MPI_finalize(ierror);
                        write(*,*) 'stop'
                        stop
                    end if
                else
                    ! Write the info into the hdf5 file
                    CALL h5_init()
                    CALL h5_open_rw(path2out, h5_id)
                    CALL h5_add_string(h5_id, trim(h5_mode_groupname)// &
                        '/info', 'discrepancy to linearly predicted value of Br_abs_res > delta')
                    CALL h5_add_double_1(h5_id, trim(h5_mode_groupname)// &
                        '/discrep_time', (/timeIndex*1.d0, time/), (/1/), (/2/))
                    CALL h5_close(h5_id)
                    CALL h5_deinit()
                    discr_reached = .true.
                end if
            end if
        end if

    end subroutine


    !> @brief subroutine writeKinProfileDataToDisk(istep). Writes the profile data to hdf5 files. 
    !> Formerly, this data was written to fort.1xxx ascii files.
    !> This routine was added because of the change that only every
    !>  "save_prof_time_step"th timestep is written. If the program is to be stopped
    !> because a stopping criterion was met, the profiles should be written for that
    !> last time step. Because this occurs more than once, it is more convenient to
    !> summarize this in a subroutine.
    !> @author Markus Markl
    !> @date 12.03.2021
    !> @param[in] istep Current step of the time evolution. Used to name the fort.1000 group in which
    !> the data is written.
    subroutine writeKinProfileDataToDisk(istep)

        use grid_mod
        use plasma_parameters
        use control_mod
        use baseparam_mod
        use h5mod
        use hdf5_tools
        use wave_code_data, only: Vth

        implicit none
        integer, intent(in) :: istep
        integer :: ipoi

        if (ihdf5IO .eq. 1) then
            print *, "Write kinetic profiles"
            do ipoi = 1, npoic
                sqg_bthet_overcavg(ipoi) = 0.5d0*(sqg_bthet_overc(ipoi) &
                                                + sqg_bthet_overc(ipoi + 1))
                Ercovavg(ipoi) = 0.5d0*(Ercov(ipoi) + Ercov(ipoi + 1))
            end do
            ! h5_mode_groupname
            h5_currentgrp = "/"//trim(h5_mode_groupname) &
                            //"/KinProfiles"

            CALL h5_init()
            CALL h5_open_rw(path2out, h5_id)

            ! create datasets
            write (h5_currentgrp, "(A,A,I4,A)") trim(h5_currentgrp), &
                "/", 1000 + istep, "/"

            write (*, *) "h5_currentgrp ", trim(h5_currentgrp)
            write (*, *) "defining KinProfiles/1000 group ", 1000 + istep
            ! define group 1000+istep
            CALL h5_obj_exists(h5_id, trim(h5_currentgrp), h5_exists_log)
            if (.not. h5_exists_log) then
        	    CALL h5_define_group(h5_id, trim(h5_currentgrp), group_id_1)
		    end if



            write (*, *) "group defined"
            ! edited 26.03.2021, Markus Markl
            ! The profiles are now saved in single precision
            CALL h5_add_float_1(h5_id, trim(h5_currentgrp)//"rc", &
                                real(rc), lbound(rc), ubound(rc))

            CALL h5_add_float_1(h5_id, trim(h5_currentgrp)//"n", &
                                real(params(1, :)), lbound(params(1, :)), ubound(params(1, :)))

            CALL h5_add_float_1(h5_id, trim(h5_currentgrp)//"Vz", &
                                real(params(2, :)), lbound(params(2, :)), ubound(params(2, :)))

            CALL h5_add_float_1(h5_id, trim(h5_currentgrp)//"Te", &
                                real(params(3, :)/ev), lbound(params(3, :)), ubound(params(3, :)))

            CALL h5_add_float_1(h5_id, trim(h5_currentgrp)//"Ti", &
                                real(params(4, :)/ev), lbound(params(4, :)), ubound(params(4, :)))

            CALL h5_add_float_1(h5_id, trim(h5_currentgrp)//"Er", &
                                real(Ercovavg), lbound(Ercovavg), ubound(Ercovavg))

            CALL h5_add_float_1(h5_id, trim(h5_currentgrp)//"sqg_btheta_overc", &
                                real(sqg_bthet_overcavg), lbound(sqg_bthet_overcavg), &
                                ubound(sqg_bthet_overcavg))

            CALL h5_add_float_1(h5_id, trim(h5_currentgrp)//"Vth", &
                                real(Vth), lbound(Vth), ubound(Vth))

            CALL h5_close_group(group_id_1)
            CALL h5_close(h5_id)
            CALL h5_deinit()
            write (*, *) "finished writing KinProfiles"

        else
            do ipoi = 1, npoic
                write (1000 + istep, *) rc(ipoi), params(1:2, ipoi) &
                    , params(3, ipoi)/ev &
                    , params(4, ipoi)/ev &
                    , 0.5d0*(Ercov(ipoi) + Ercov(ipoi + 1)) &
                    , 0.5d0*(sqg_bthet_overc(ipoi) + &
                            sqg_bthet_overc(ipoi + 1))
            end do
            close (1000 + istep)
        end if

    end subroutine writeKinProfileDataToDisk

    subroutine writeKinProfileAtTimeIndex

        implicit none
        if (irank .eq. 0) then
            if (modulo(timeIndex, save_prof_time_step) .eq. 0) then
                if (suppression_mode .eqv. .false.) then
                    CALL writeKinProfileDataToDisk(timeIndex)
                end if
            end if
        end if

    end subroutine


    subroutine interpBrAndDqlAtResonanceTimeEvol

        use PolyLagrangeInterpolation
        use grid_mod, only: npoib, r_resonant, rb, dqle22
        use wave_code_data, only: antenna_factor, Br

        implicit none
        
        call binsrc(rb, 1, npoib, r_resonant(1), indResRadius)
        call getIndicesForLagrangeInterp(indResRadius)
        call plag_coeff(nlagr, nder, r_resonant(1), rb(indBeginInterp:indEndInterp), coef)

        br_abs(timeIndex) = sum(coef(0, :)*abs(Br(indBeginInterp:indEndInterp)))*sqrt(antenna_factor)
        dqle22_res_time(timeIndex) = sum(coef(0, :)*dqle22(indBeginInterp:indEndInterp))

        ! save the time for the improved stopping criterion
        br_abs_time(timeIndex) = time
		br_abs_antenna_factor(timeIndex) = antenna_factor

        write(*,*) 'Br abs res * C_mn= ', br_abs(timeIndex)
        write(*,*) 'Br abs res       = ', br_abs(timeIndex)/sqrt(antenna_factor)
        write(*,*) 'Dqle22 res       = ', dqle22_res_time(timeIndex)
        write(*,*) 'Antenna factor   = ', antenna_factor
        write(*,*) 'time = ', br_abs_time(timeIndex)

    end subroutine

    subroutine stopIfTimeStepTooSmall

        use h5mod
        use parallelTools, only: ierror
        use mpi

        implicit none

        character(*), parameter :: reason = 'timestep < stop time step'

        if (timstep .lt. stop_time_step .and. time .gt. 1.0d-3) then
            write(*,*) 'stop: timestep smaller than stop limit'
            if (suppression_mode .eqv. .false.) then
                CALL writeKinProfileDataToDisk(timeIndex)
            end if

            call writeReasonForStopToH5(reason)
            call MPI_finalize(ierror);
            stop
        end if

    end subroutine


    subroutine calcParamsNumAndDenom

        use plasma_parameters, only: params_num, params_denom, params, params_beg

        implicit none

        params_num = (params - params_beg)**2
        params_denom = params**2 + params_beg**2

    end subroutine

    subroutine smoothParamsNumAndDenom
    
        use grid_mod, only: npoi, nbaleqs, mwind, dummy
        use plasma_parameters, only: params_num, params_denom

        implicit none
        
        integer :: ieq

        do ieq = 1, nbaleqs
            call smooth_array_gauss(npoi, mwind, params_num(ieq, :), dummy)
            params_num(ieq, :) = dummy
            call smooth_array_gauss(npoi, mwind, params_denom(ieq, :), dummy)
            params_denom(ieq, :) = dummy
        end do

    end subroutine

    subroutine determineTimscal

        use grid_mod, only: npoi, rc
        use plasma_parameters, only: params_num, params_denom
        use baseparam_mod, only: factolmax
        use recstep_mod, only: tol

        implicit none

        integer :: ipoi

        do ipoi = 1, npoi
            if (rc(ipoi) .lt. 0.95d0*rc(npoi)) then
                timscal(ipoi) = sum(sqrt(params_num(3:4, ipoi)/params_denom(3:4, ipoi)))
            else
                timscal(ipoi) = 1d-30 !0.d0
            end if
        end do

        if (debug_mode) then
            if (irank .eq. 0) then
                write(*,*) "maxval(timscal) = ", maxval(timscal)
                write(*,*) "tol*factolmax = ", tol*factolmax
            end if
        end if

    end subroutine

    subroutine rescaleTimStepArr

        use grid_mod, only: npoi, nbaleqs
        use recstep_mod, only: tim_stack, timstep_arr, tol
        use diag_mod, only: timscal_dql

        implicit none

        integer :: ipoi, ieq, k

        timscal = timscal + timscal_dql
                        
        do ipoi = 1, npoi
            do ieq = 1, nbaleqs
                k = nbaleqs*(ipoi - 1) + ieq
                !timstep_arr(k)=timstep_arr(k)/timscal(ipoi)*tol
                timstep_arr(k) = timstep_arr(k)/max(timscal(ipoi), epsilon(1.d0))*tol
                ! steady state solution:
                !if (ieq .gt. 1 .and. r(ipoi) .gt. rsepar-0.5d0) then
                !    timstep_arr(k) = 0d0
                !end if
            end do
        end do
        
        timstep_arr = timstep_arr*timescale/(timstep_arr + timescale)
        if (scratch) then
            scratch = .false.
            tim_stack = timstep_arr
        end if
        timstep_arr = 2.d0*timstep_arr*tim_stack/(timstep_arr + tim_stack)

    end subroutine

    subroutine setTimStep

        use recstep_mod, only: timstep_arr, tol
        
        implicit none

        timstep = minval(timstep_arr)

        if (.true.) then
            ! limit time step from below:
            timstep = max(timstep, timstep_min)
            ! limit timestep from above:
            !if (ramp_up_mode .ne. 0) timstep = min(timstep,0.1)
            !timstep = min(timstep,0.005)
        else
        ! use for constant time step:
            timstep = 0.5
            write(*,*) "constant time step = ", timstep
        end if

        if (irank .eq. 0) then
            write(*,*) 'timstep', real(timstep), '   timescale', real(timescale), &
                'tolerance', real(tol)
        end if

    end subroutine

    subroutine resetTimStepArrWithTimstep
        
        use recstep_mod, only: tim_stack, timstep_arr

        implicit none
        
        timstep_arr = timstep
        tim_stack = timstep_arr

    end subroutine

    subroutine writeTimeInfoToDisk

        implicit none
        
        if (irank .eq. 0) then
            if (ihdf5IO .eq. 1) then
                call writeTimeInfoToH5
            else
                call writeTimeInfoToTxt
            end if
        end if


    end subroutine

    subroutine writeTimeInfoToTxt

        use diag_mod, only: rate_dql, timscal_dql

        implicit none

        open (4321, file='timstep_evol.dat', position='append')
        write (4321, *) timeIndex, timstep, timscal_dql, timscal(1), rate_dql, time
        close (4321)

    end subroutine

    subroutine writeTimeInfoToH5

        use diag_mod, only: rate_dql, timscal_dql

        implicit none

        h5_currentgrp = trim("/"//trim(h5_mode_groupname)//"/timstep_evol.dat")

        CALL h5_init()
        CALL h5_open_rw(path2out, h5_id)

        if (.not. firstiterationdone) then
            if (diagnostics_output) then
                CALL h5_define_unlimited_matrix(h5_id, trim(h5_currentgrp), &
                                                H5T_NATIVE_DOUBLE, (/6, -1/), time_dataset_id)
            else
                CALL h5_define_unlimited_matrix(h5_id, trim(h5_currentgrp), &
                                                H5T_NATIVE_DOUBLE, (/3, -1/), time_dataset_id)
            end if
        end if
        ! this output works differently, because it is done for each
        ! time step and not all at once. For every time step a row
        ! is appended to the data. This is different to the case
        ! where the whole data is written at once. In the latter
        ! the data for one variable is saved in one row, i.e. there
        ! are length(variable) number of columns.
        if (diagnostics_output) then
            CALL h5_append_double_1(time_dataset_id, (/timeIndex*1.d0, timstep, &
                                                        timscal_dql, timscal(1), rate_dql, time/), timeIndex)
        else
            CALL h5_append_double_1(time_dataset_id, (/timeIndex*1.d0, timstep, &
                                                        time/), timeIndex)
        end if
        CALL h5_close(h5_id)
        CALL h5_deinit()

    end subroutine

    subroutine relaxPlasmaParameters

        use grid_mod, only: npoi, nbaleqs
        use plasma_parameters, only: params
        use baseparam_mod, only: urelax

        implicit none

        integer :: ipoi, ieq, k

        do ipoi = 1, npoi
            do ieq = 1, nbaleqs
                k = nbaleqs*(ipoi - 1) + ieq
                params(ieq, ipoi) = yprev(k)*urelax + params(ieq, ipoi)*(1.d0 - urelax)
            end do
        end do

    end subroutine

    subroutine messageTimeInfo

        implicit none
        if (irank .eq. 0) then
            write(*,*) ' '
            write(*,*) 'i = ', int2(timeIndex), 'time = ', real(time)
            write(*,*) ' '
        end if

    end subroutine

    subroutine setFirstIterationTrue

        implicit none

        if (firstiterationdone .eqv. .false.) firstiterationdone = .true.

    end subroutine


    subroutine determineDqlDiagnostic

        use grid_mod, only: dqle11, dqli11
        use diag_mod

        implicit none

        timscal_dql = maxval(abs(dqle11_prev - dqle11))/maxval(dqle11_prev + dqle11)
        ind_dqle = maxloc(abs(dqle11_prev - dqle11))
        timscal_dqli = maxval(abs(dqli11_prev - dqli11))/maxval(dqli11_prev + dqli11)
        ind_dqli = maxloc(abs(dqli11_prev - dqli11))
        rate_dql = timscal_dql/timstep

    end subroutine

end module