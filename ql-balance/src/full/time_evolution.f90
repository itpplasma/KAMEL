  module time_evolution

    implicit none

    logical :: flag_run_time_evolution !Added by Philipp Ulbl 12.05.2020
    logical :: br_stopping ! trigger Br stopping criterion

    integer :: Nstorage
    integer :: ramp_up_mode !> control ramp up mode of the RMP coil current amplitude
    integer :: save_prof_time_step ! added by Markus Markl 11.03.2021
    integer :: iexit ! used for ramp-up skipping of saving
    integer :: ramp_up_down = 0 !> used in hysteresis mode, tells if ramp-up (0) or ramp-down (1)

    double precision :: tmax_factor!, antenna_factor
    double precision :: stop_time_step !Added by Philipp Ulbl 13.05.2020
    double precision :: timstep_min
    DOUBLE PRECISION :: t_max_ramp_up = 1e-2 !> 10ms ramp up until antenna_factor_max is reached
    double precision :: timstep
    double precision :: time
    double precision :: t_hysteresis_turn = 0
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

    contains

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

    subroutine rescaleTranspCoefficientsByAntennaFac

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


    !> @brief subroutine write_br_time_data. Writes radial magnetic field perturbation evaluated at the resonant
    !> surface, the antenna factor, the time and Dqle22 evaluated at the resonant surface for a given 
    !> time step to the hdf5 file
    !> @param[in] i Integer of time step to which the data will be saved. Goes from 1:i.
    !> @param[in] br_abs_time Time value of the time evolution.
    !> @param[in] br_abs_antenna_factor Value of the antenna factor, i.e. the RMP coil current.
    !> @param[in] br_abs Absolute value of the radial magnetic field evaluated at the resonant surface in question.
    !> @param[in] dqle22_res_time Value of Dqle22 evaluated at the resonant surface during the time evolution.
    subroutine write_br_time_data(i)

        use control_mod
        use baseparam_mod
        use h5mod
        use hdf5_tools
        use wave_code_data, only: antenna_factor

        implicit none
	    integer, intent(in) :: i
        !double precision, dimension(:), intent(in) :: br_abs_time
        !double precision, dimension(:), intent(in) :: br_abs_antenna_factor
        !double precision, dimension(:), intent(in) :: br_abs
        !double precision, dimension(:), intent(in) :: dqle22_res_time
        character(len=1024) :: h5_currentgrp

	    if (debug_mode) write(*,*) "Debug: writing out br time evolution data"

        if (ihdf5IO .eq. 1) then
       	    CALL h5_init()
            CALL h5_open_rw(path2out, h5_id)

 		    h5_currentgrp = "/"//trim(h5_mode_groupname) //"/br_abs_time"
		    CALL h5_add_double_1(h5_id, trim(h5_currentgrp), br_abs_time(1:i), &
			    lbound(br_abs_time(1:i)), ubound(br_abs_time(1:i)))

 		    h5_currentgrp = "/"//trim(h5_mode_groupname) //"/br_abs_antenna_factor"
		    CALL h5_add_double_1(h5_id, trim(h5_currentgrp), br_abs_antenna_factor(1:i), &
			    lbound(br_abs_antenna_factor(1:i)), ubound(br_abs_antenna_factor(1:i)))

 		    h5_currentgrp = "/"//trim(h5_mode_groupname) //"/br_abs_res"
		    CALL h5_add_double_1(h5_id, trim(h5_currentgrp), br_abs(1:i), &
			    lbound(br_abs(1:i)), ubound(br_abs(1:i)))

 		    h5_currentgrp = "/"//trim(h5_mode_groupname) //"/dqle22_res_time"
		    CALL h5_add_double_1(h5_id, trim(h5_currentgrp), dqle22_res_time(1:i), &
			    lbound(dqle22_res_time(1:i)), ubound(dqle22_res_time(1:i)))


            CALL h5_close(h5_id)
            CALL h5_deinit()

        else
            open (777, file='br_abs_res.dat', position='append')
            write (777, *) i, time, antenna_factor, br_abs(i)
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
                        call writefort1000(i)
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
                            CALL write_br_time_data(i)!, br_abs_time, br_abs_antenna_factor, br_abs, dqle22_res_time)
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
                            CALL write_br_time_data(i)!, br_abs_time, br_abs_antenna_factor, br_abs, dqle22_res_time)
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
                        call writefort1000(i)
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
                            CALL write_br_time_data(i)!, br_abs_time, br_abs_antenna_factor, br_abs, dqle22_res_time)
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
                            CALL write_br_time_data(i)!, br_abs_time, br_abs_antenna_factor, br_abs, dqle22_res_time)
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
                            call writefort1000(i)
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
                            CALL write_br_time_data(i)!, br_abs_time, br_abs_antenna_factor, br_abs, dqle22_res_time)
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
                            call writefort1000(i)
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
                            CALL write_br_time_data(i)!, br_abs_time, br_abs_antenna_factor, br_abs, dqle22_res_time)
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
                    call writefort1000(i)
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
                    CALL write_br_time_data(i)!, br_abs_time, br_abs_antenna_factor, br_abs, dqle22_res_time)
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

    !> @brief subroutine writefort1000(istep). Writes the profile data to hdf5 files. 
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
    subroutine writefort1000(istep)

        use grid_mod
        use control_mod
        use baseparam_mod
        use h5mod
        use hdf5_tools
        use wave_code_data, only: Vth

        implicit none
        integer, intent(in) :: istep
        integer :: ipoi
        character(len=1024) :: h5_currentgrp

        if (ihdf5IO .eq. 1) then
            !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
            ! write profiles to hdf5 (former fort.1000+ files)
            !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
            do ipoi = 1, npoic
                sqg_bthet_overcavg(ipoi) = 0.5d0*(sqg_bthet_overc(ipoi) &
                                                + sqg_bthet_overc(ipoi + 1))
                Ercovavg(ipoi) = 0.5d0*(Ercov(ipoi) + Ercov(ipoi + 1))
            end do
            ! h5_mode_groupname
            h5_currentgrp = "/"//trim(h5_mode_groupname) &
                            //"/fort.1000"

            CALL h5_init()
            CALL h5_open_rw(path2out, h5_id)

            ! create datasets
            write (h5_currentgrp, "(A,A,I4,A)") trim(h5_currentgrp), &
                "/", 1000 + istep, "/"

            write (*, *) "h5_currentgrp ", trim(h5_currentgrp)
            write (*, *) "defining fort.1000/1000 group ", 1000 + istep
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
            write (*, *) "finished writing fort.1000"

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

    end subroutine writefort1000


end module