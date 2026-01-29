!> @brief Ramp up/down the RMP coil current.
!> @author Markus Markl
!> @date 13.03.2023
subroutine ramp_coil

    use control_mod, only: debug_mode
    use time_evolution, only: ramp_up_mode

    implicit none

    if (debug_mode) write(*,*) "Debug: - - ramp up antenna_factor - -"

    select case(ramp_up_mode)
        case(0)
            call ramp_up_linear
        case(1)
            call ramp_up_faster_1
        case(2)
            call ramp_up_faster_2 !stable compared to 1
        case(3)
            call ramp_up_instant
        case(4)
            call ramp_up_none
        case(5)
            call ramp_up_hysteresis
        case(6)
            call ramp_up_fast_hysteresis
        case(10)
            call ramp_up_hyst_mod
        case(11)
            call ramp_oscillation
        case(12)
            call ramp_up_down_zero
        case default
            stop 'Error: ramp_up_mode not defined'
    end select

end subroutine


subroutine ramp_up_linear

    use wave_code_data, only: antenna_factor
    use time_evolution, only: time
    use control_mod, only: debug_mode

    implicit none

    if (debug_mode) write(*,*) "Debug: ramp-up mode linear"
    ! initial ramp up. This is a linear ramp up of the antenna factor in time.
    antenna_factor = time**2 + 1.d-4

end subroutine

subroutine ramp_up_faster_1

    use wave_code_data, only: antenna_factor
    use time_evolution, only: timstep, t_max_ramp_up, antenna_factor_max
    use control_mod, only: debug_mode

    implicit none

    if (debug_mode) write(*,*) "Debug: ramp-up mode faster 1"
    ! First try for faster ramp up. Is (usually) not stable.
    antenna_factor = antenna_factor + antenna_factor_max * (timstep/t_max_ramp_up)

end subroutine

subroutine ramp_up_faster_2

    use wave_code_data, only: antenna_factor
    use time_evolution, only: time, t_max_ramp_up, antenna_factor_max
    use control_mod, only: debug_mode

    implicit none

    if (debug_mode) write(*,*) "Debug: ramp-up mode faster 2, stable"
    ! Use this for faster ramp up. Ramps up antenna factor to 100% of max value
    ! in t_max_ramp_up time, but ramps-up further after that.
    if (time .eq. 0) then
        antenna_factor = 1.d-4
    else
        antenna_factor = antenna_factor_max * (time/t_max_ramp_up)**2.0d0
    end if

end subroutine

subroutine ramp_up_instant

    use wave_code_data, only: antenna_factor
    use time_evolution, only: time, antenna_factor_max
    use control_mod, only: debug_mode

    implicit none

    if (debug_mode) write(*,*) "Debug: ramp-up mode instant"

    call stop_if_t_max_reached

    if (time .eq. 0) then
        antenna_factor = 1.d-4
    else
        antenna_factor = antenna_factor_max
    end if

end subroutine

subroutine ramp_up_none

    use wave_code_data, only: antenna_factor
    use control_mod, only: debug_mode

    implicit none

    if (debug_mode) write(*,*) "Debug: ramp-up mode none"
    call stop_if_t_max_reached
    antenna_factor = 0d0

end subroutine

subroutine ramp_up_hysteresis

    use time_evolution, only: time, ramp_up_down, t_hysteresis_turn, antenna_factor_max, antenna_max_stopping, &
        write_kin_prof_data_to_disk, write_br_dqle22_time_data
    use wave_code_data, only: antenna_factor
    use control_mod, only: debug_mode
    use h5mod

    implicit none

    integer :: ierror

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
                call write_kin_prof_data_to_disk
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
                CALL write_br_dqle22_time_data!, br_abs_time, br_abs_antenna_factor, br_abs, dqle22_res_time)
            end if

            CALL MPI_finalize(ierror)
            stop
        end if
    end if ! ramp_up_down eq 1

end subroutine

subroutine ramp_up_fast_hysteresis

    use time_evolution, only: time, t_max_ramp_up, ramp_up_down, t_hysteresis_turn, antenna_factor_max, &
        antenna_max_stopping, write_kin_prof_data_to_disk, write_br_dqle22_time_data
    use wave_code_data, only: antenna_factor
    use control_mod, only: debug_mode
    use h5mod

    implicit none

    integer :: ierror

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
                call write_kin_prof_data_to_disk
            end if
            ! Write the cause of the stopping into the hdf5 file
            if (ihdf5IO .eq. 1) then
                CALL h5_init()
                CALL h5_open_rw(path2out, h5_id)
                CALL h5_add_string(h5_id, trim(h5_mode_groupname)// &
                    '/stopping_criterion', 'ramp-up/down finished')
                CALL h5_close(h5_id)
                CALL h5_deinit()

                CALL write_br_dqle22_time_data!, br_abs_time, br_abs_antenna_factor, br_abs, dqle22_res_time)
            end if

            CALL MPI_finalize(ierror)
            stop
        end if ! antenna_facotr less eq 0
    end if ! ramp_up_down eq 1

end subroutine

subroutine stop_if_antenna_fac_max_reached

    use control_mod, only: suppression_mode, ihdf5IO
    use time_evolution, only: antenna_max_stopping, antenna_factor_max, write_br_dqle22_time_data, &
        write_kin_prof_data_to_disk
    use wave_code_data, only: antenna_factor
    use h5mod

    implicit none

    integer :: ierror

    if (antenna_factor .gt. (antenna_factor_max * antenna_max_stopping)) then
        write(*,*) 'stop: reached antenna_factor_max * ', antenna_max_stopping

        if (suppression_mode .eqv. .false.) then
            call write_kin_prof_data_to_disk
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
            CALL write_br_dqle22_time_data!, br_abs_time, br_abs_antenna_factor, br_abs, dqle22_res_time)
        end if

        call MPI_finalize(ierror);
        stop
    end if

end subroutine

subroutine ramp_up_hyst_mod
    ! 1: ramp up linearlly,
    ! 2: flat top
    ! 3: modulation of flat top with sine function

    use time_evolution, only: bif_criterion, time_ind, time
    use wave_code_data, only: antenna_factor
    use time_evolution_stellarator, only: hyst_mod_stage, hyst_mod_amp_fac, hyst_mod_freq, hyst_mod_phase, t_flattop_begin, &
        bif_crit_flattop, t_flattop_end, ant_fac_flattop, delta_t_flattop
    use baseparam_mod, only: pi

    implicit none

    if (bif_criterion(time_ind) .lt. bif_crit_flattop .and. hyst_mod_stage .eq. 0) then
        call ramp_up_faster_2
        print *, " ramp up phase 0: linear"
    else if(bif_criterion(time_ind) .ge. bif_crit_flattop .and. hyst_mod_stage .eq. 0) then
        hyst_mod_stage = 1
        t_flattop_begin = time
        print *, " begin of ramp up phase 1: flat top"
        print *, " time = ", time
    else if (hyst_mod_stage .eq. 1 .and. abs(time - t_flattop_begin) .gt. delta_t_flattop) then
        ! flat top
        hyst_mod_stage = 2
        ant_fac_flattop = antenna_factor
        t_flattop_end = time
        print *, " end of ramp up phase 1: flat top"
        print *, " time = ", time
    else if (hyst_mod_stage .eq. 2) then
        print *, " ramp up phase 2: oscillation"
        antenna_factor = ant_fac_flattop * (1.0d0 + hyst_mod_amp_fac * sin(2.0d0 * pi * (time - t_flattop_end) &
            * hyst_mod_freq + hyst_mod_phase))
        if (time .ge. t_flattop_end + 2.0d0/hyst_mod_freq) then
            call stop_evolution
        end if
    end if

end subroutine

subroutine ramp_oscillation

    use time_evolution, only: time_ind, time
    use wave_code_data, only: antenna_factor
    use time_evolution_stellarator, only: hyst_mod_freq, hyst_mod_phase, hyst_mod_amp_fac
    use time_evolution, only: antenna_factor_max, t_max_ramp_up
    use baseparam_mod, only: pi

    implicit none

    antenna_factor = antenna_factor_max * (hyst_mod_amp_fac * sin(2.0d0 * pi * time * hyst_mod_freq + hyst_mod_phase))**2

    if (time .ge. t_max_ramp_up) then
        call stop_evolution
    end if

end subroutine

subroutine ramp_up_down_zero

    use time_evolution, only: time_ind, time
    use wave_code_data, only: antenna_factor
    use time_evolution_stellarator, only: hyst_mod_freq, hyst_mod_phase, hyst_mod_amp_fac
    use time_evolution, only: antenna_factor_max, t_max_ramp_up
    use baseparam_mod, only: pi

    implicit none

    if (time .lt. t_max_ramp_up/2.0d0) then
        antenna_factor = antenna_factor_max * (hyst_mod_amp_fac * sin(2.0d0 * pi * time * hyst_mod_freq + hyst_mod_phase))**2
    else
        antenna_factor = 0.0d0
    end if

    if (time .ge. t_max_ramp_up) then
        call stop_evolution
    end if

end subroutine

subroutine stop_if_t_max_reached

    use time_evolution, only: time, t_max_ramp_up, write_kin_prof_data_to_disk, write_br_dqle22_time_data
    use control_mod, only: suppression_mode, ihdf5IO
    use h5mod

    implicit none

    integer :: ierror

    if (time .ge. t_max_ramp_up) then ! if max time value is reached, stop the code
        write(*,*) 'stop: reached time max: ', t_max_ramp_up
        if (suppression_mode .eqv. .false.) then
            call write_kin_prof_data_to_disk
        end if
        ! Write the cause of the stopping into the hdf5 file
        if (ihdf5IO .eq. 1) then
            CALL h5_init()
            CALL h5_open_rw(path2out, h5_id)
            CALL h5_add_string(h5_id, trim(h5_mode_groupname)// &
                    '/stopping_criterion', 'reached time max')
            CALL h5_close(h5_id)
            CALL h5_deinit()

            CALL write_br_dqle22_time_data!, br_abs_time, br_abs_antenna_factor, br_abs, dqle22_res_time)
        end if
        if (debug_mode) write(*,*) "Debug: Write br_time _data"

        CALL MPI_finalize(ierror)
        stop
    end if

end subroutine

subroutine stop_evolution

    use time_evolution, only: time, t_max_ramp_up, write_kin_prof_data_to_disk, write_br_dqle22_time_data
    use control_mod, only: suppression_mode, ihdf5IO
    use h5mod

    implicit none

    integer :: ierror

    if (suppression_mode .eqv. .false.) then
        call write_kin_prof_data_to_disk
    end if
    ! Write the cause of the stopping into the hdf5 file
    CALL write_br_dqle22_time_data!, br_abs_time, br_abs_antenna_factor, br_abs, dqle22_res_time)

    CALL MPI_finalize(ierror)
    stop

end subroutine


subroutine check_linear_discr_pen_ratio

    use time_evolution, only: br_beta, time_ind, br_abs_time, save_prof_time_step, &
        br_stopping, discr_reached, br_abs, br_predicted, write_br_dqle22_time_data, write_kin_prof_data_to_disk
    use control_mod, only: suppression_mode
    use mpi
    use h5mod, only: write_reason_for_stop_to_h5

    implicit none

    integer :: ierror

    if (time_ind .gt. 50 .and. .not. discr_reached) then
        ! calculate beta only once
        if (br_beta .eq. 0) then
            ! Calculate slope from the data until this time step. (Simple linear regression algorithm)
            br_beta = sum(br_abs(3:time_ind)*br_abs_time(3:time_ind))/sum(br_abs_time(1:time_ind)**2)
            write(*,*) 'br_beta = ', br_beta
        end if
        ! calculate the from the linear regression predicted value of Br_abs
        br_predicted = br_beta*br_abs_time(time_ind)
        write(*,*) "Delta = ", abs(br_abs(time_ind) - br_predicted)
        if (abs(br_abs(time_ind) - br_predicted) .gt. 0.1) then
            write(*,*) 'discrepancy to linearly predicted value of Br_abs_res > delta'
            if (modulo(time_ind, save_prof_time_step) .ne. 0) then
                if (suppression_mode .eqv. .false.) then
                    CALL write_fields_currs_transp_coefs_to_h5
                end if
            end if
            if (suppression_mode .eqv. .false.) then
                CALL write_kin_prof_data_to_disk
            end if
            if (br_stopping) then
                ! Write the cause of the stopping into the hdf5 file
                call write_reason_for_stop_to_h5("discrepancy to " //&
                    "linearly predicted value of Br_abs_res > delta")
                CALL write_br_dqle22_time_data!, br_abs_time, br_abs_antenna_factor, br_abs, dqle22_res_time)
                CALL MPI_finalize(ierror);
                stop "Finished time evolution: br_stopping"

            else
                call write_br_discrepancy_reached_info
                discr_reached = .true.
            end if
        end if
    end if

end subroutine

subroutine write_br_discrepancy_reached_info

    use h5mod
    use time_evolution, only: time_ind, time

    implicit none

    CALL h5_init()
    CALL h5_open_rw(path2out, h5_id)
    CALL h5_add_string(h5_id, trim(h5_mode_groupname)// &
        '/info', 'discrepancy to linearly predicted value of Br_abs_res > delta')
    CALL h5_add_double_1(h5_id, trim(h5_mode_groupname)// &
        '/discrep_time', (/time_ind*1.d0, time/), (/1/), (/2/))
    CALL h5_close(h5_id)
    CALL h5_deinit()

end subroutine
