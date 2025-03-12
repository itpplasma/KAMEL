module time_evolution_stellarator

    use control_mod
    use parallelTools
    use h5mod
    use balance_base, only: balance_t
    use time_evolution, only: Nstorage, ramp_up_mode, save_prof_time_step, iexit, ramp_up_down, time_ind, &
        antenna_max_stopping, timstep_min, tmax_factor, t_max_ramp_up, br_formfactor, br_vac_res, &
        antenna_factor_max, br_beta, br_predicted, scratch, br_abs, br_abs_time, br_abs_antenna_factor, &
        dqle22_res_time, dae22_res_time, bif_criterion, firstiterationdone, dqle11_prev, dqle12_prev, &
        dqle21_prev, dqle22_prev, dqli11_prev, dqli12_prev, dqli21_prev, dqli22_prev, yprev, stop_time_step, &
        write_kin_prof_data_to_disk, write_br_dqle22_time_data
    use QLBalance_kinds, only: dp

    implicit none

    logical :: set_momentum_source_to_zero, set_Q_neo_to_zero
    integer :: update_transport_coefficients
    logical :: reduce_time_step, turn_off_heat_sources

    real(dp) :: tmax, timescale

    real(dp) :: timstep
    real(dp) :: time

    real(dp), dimension(:), allocatable :: timscal


    integer(HID_T) :: time_dataset_id !> variable to save the time dataset id

    type, extends(balance_t) :: time_evolution_stellarator_t
        contains
            procedure :: init_balance => initTimeEvolution
            procedure :: run_balance => runTimeEvolution
    end type

    contains

    subroutine initTimeEvolution(this)

        use recstep_mod, only: tol
        use transp_coeffs_mod, only: rescale_transp_coeffs_by_ant_fac
        use grid_mod, only: mwind, rmax, rmin, set_boundary_condition, npoib, rb
        use baseparam_mod, only: dperp, tol_max
        use QLbalance_diag, only: write_diag, write_diag_b
        use QLBalance_hdf5_tools, only: h5overwrite
        use h5mod, only: mode_m, mode_n
        use control_mod, only: gyro_current_study, write_gyro_current, debug_mode, &
                        ihdf5IO
        use parallelTools, only: initMPI, irank
        use wave_code_data, only: m_vals, n_vals
        use plasma_parameters, only: write_initial_parameters, alloc_hold_parameters, &
                                params, params_begbeg, init_background_profiles
        implicit none

        class(time_evolution_stellarator_t), intent(inout) :: this
        this%runType = "TimeEvolution"
        
        if (irank .eq. 0) then

            call read_stell_config
            call print_stell_config

            iexit = 0 ! 0 - dont skip, 1 - skip, 2 - stop
            mwind = 10
            write_diag = .false.
            write_diag_b = .false.

            ! if h5overwrite = true, existing data will be deleted
            ! before new one is written
            ! This is contained in hdf5_tools module
            h5overwrite = .true.
    
            if (gyro_current_study .ne. 0) then
                write_gyro_current = .true.
            else
                write_gyro_current = .false.
            end if

            !call read_config
            timescale = (rmax - rmin)**2 / dperp
            tmax = timescale * tmax_factor
            timstep = tmax / Nstorage
            time = 0.0d0
            tol = tol_max

            call gengrid
            call set_boundary_condition

            CALL initialize_wave_code_interface(npoib, rb);

            mode_m = m_vals(1)
            mode_n = n_vals(1)
            if (ihdf5IO .eq. 1) then
                CALL create_group_structure_timeevol
            end if
            if (debug_mode) write(*,*) 'Debug: mode_m = ', mode_m, 'mode_n = ', mode_n

            call allocate_prev_variables
            call init_background_profiles
            CALL write_initial_parameters
        end if

        call alloc_Br_Dqle_for_timeevol
        call calc_geometric_parameter_profiles
        call initialize_get_dql
        call initialize_D_one_over_nu
        call initialize_antenna_factor
        call det_balance_eqs_source_terms_stell

        if (irank .eq. 0) then
            if (.not. suppression_mode) call write_kin_prof_data_to_disk
        end if

        call allocate_timscal_and_params
        timstep = timstep*tol
        scratch = .true.

        call reset_timstep_arr_w_timstep

        call get_dql
        call rescale_transp_coeffs_by_ant_fac
        call hold_prev_transp_coeffs

        params_begbeg = params

    end subroutine

    subroutine runTimeEvolution(this)

        use parallelTools, only: irank
        use baseparam_mod, only: factolmax, factolred
        use recstep_mod, only: tol
        use plasma_parameters, only: params, params_beg, params_begbeg, limit_temps_from_below
        use restart_mod, only: redostep
        use transp_coeffs_mod, only: rescale_transp_coeffs_by_ant_fac
        use recstep_mod, only: timstep_arr

        implicit none

        class(time_evolution_stellarator_t), intent(inout) :: this
        integer :: iredo = 0

        do time_ind = 1, Nstorage

            call copy_kin_profs_to_yprev
            redostep = .false.

            if (update_transport_coefficients .eq. 1) then
                call get_D_one_over_nu
            else if (update_transport_coefficients .eq. 2) then
                call rescale_D_one_over_nu_to_new_n_and_T
            end if

            call get_dql
            call rescale_transp_coeffs_by_ant_fac
            call interp_Br_Dql_at_resonance_timeevol
            call determine_Dql_diagnostic

            call write_br_dqle22_time_data
            call message_Br_Dqle_values

            if (diagnostics_output)then
                call writefort9999_stellarator
            end if

            if (.true.) then
                call hold_prev_transp_coeffs
                params_begbeg = params
            else 
                call redoTimeStep
            end if

            do ! evolve in time, reduce time step if too large
                iredo = iredo + 1
                params_beg = params

                print *, ""
                if (debug_mode) write(*,*) "Debug: Timstep before evolvestep is ", timstep, " eps = " , eps
                call evolvestep_stell(timstep, eps)

                call limit_temps_from_below

                call calc_params_num_and_denom
                call smooth_params_num_and_denom
                call determine_timscal

                if (maxval(timscal) .lt. tol * factolmax) then
                    exit
                end if

                timstep_arr = timstep_arr * factolred
                params = params_beg

                if (irank .eq. 0) then
                    if (debug_mode) then
                        print *, "Redoing step: Maxval(timscal) is not lesser than tol * factolmax"
                        print *, "Maxval(timscal) = ", maxval(timscal)
                        print *, "tol = ", tol
                        print *, "factolmax = ", factolmax
                        print *, "tol * factolmax = ", tol * factolmax
                        print *, "timstep_arr(1) = ", timstep_arr(1)
                        print *, "timstep_arr(100) = ", timstep_arr(100)
                        print *, ""
                    end if
                end if
                if (iredo > 300) then
                    stop "Redoing step: Maxval(timscal) is not lesser than tol * factolmax after 100 redos"
                end if
            end do

            call rescale_time_step_array
            call set_time_step
            call stop_if_time_step_too_small

            call reset_timstep_arr_w_timstep
            call write_time_info

            call relax_plasma_parameters

            timstep_arr = 0.0d0
            call evolvestep_stell(timstep, eps)
            timstep_arr = timstep
            time = time + timstep

            if (debug_mode) call messageTimeInfo
            if (.not. suppression_mode) call write_kin_profile_at_time_index
            call set_first_iteration_true
            call check_linear_discr_pen_ratio
            call stop_if_antenna_fac_max_reached

            call ramp_coil
        end do

    end subroutine

    subroutine read_stell_config

        implicit none

        NAMELIST /BALANCE_STELL/ set_momentum_source_to_zero, set_Q_neo_to_zero, &
            update_transport_coefficients, reduce_time_step, turn_off_heat_sources

        ! read the parameters from namelist file
        open (22, file='stell_conf.nml');
        read (22, NML=BALANCE_STELL)
        close (22);

    end subroutine

    subroutine print_stell_config

        implicit none

        print *, "============================================================================="
        print *, "    Stellarator mode configuration"
        print *, "        set_momentum_source_to_zero           = ", set_momentum_source_to_zero
        print *, "        set_Q_neo_to_zero                     = ", set_Q_neo_to_zero
        print *, "        update_transport_coefficients         = ", update_transport_coefficients
        print *, "        reduce_time_step                      = ", reduce_time_step
        print *, "        turn_off_heat_sources                 = ", turn_off_heat_sources
        print *, "============================================================================="

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

    subroutine copy_kin_profs_to_yprev

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

    subroutine initialize_antenna_factor
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

    subroutine alloc_Br_Dqle_for_timeevol

        implicit none

        allocate(br_abs(Nstorage))
        allocate(br_formfactor(Nstorage))
        allocate(br_vac_res(Nstorage))
        allocate(br_abs_antenna_factor(Nstorage))
        allocate(br_abs_time(Nstorage))
        allocate(dqle22_res_time(Nstorage))
        allocate(dae22_res_time(Nstorage))
        allocate(bif_criterion(Nstorage))

    end subroutine

    subroutine hold_prev_transp_coeffs

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


    subroutine write_kin_profile_at_time_index

        implicit none

        if (irank .eq. 0) then
            if (debug_mode) write(*,*) "Debug: Write kinetic profiles at time index: ", time_ind
            if (modulo(time_ind, save_prof_time_step) .eq. 0) then
                CALL write_kin_prof_data_to_disk
            end if
        end if

    end subroutine


    subroutine interp_Br_Dql_at_resonance_timeevol

        use PolyLagrangeInterpolation
        use grid_mod, only: npoib, r_resonant, rb, dqle22, dae22
        use wave_code_data, only: antenna_factor, Br

        implicit none

        integer :: indResRadius, ind_begin_interp, ind_end_interp
        
        call binsrc(rb, 1, npoib, r_resonant(1), indResRadius)
        call get_ind_Lagr_interp(indResRadius, ind_begin_interp, ind_end_interp)
        call plag_coeff(nlagr, nder, r_resonant(1), rb(ind_begin_interp:ind_end_interp), coef)

        br_abs(time_ind) = sum(coef(0, :)*abs(Br(ind_begin_interp:ind_end_interp)))*sqrt(antenna_factor)
        dqle22_res_time(time_ind) = sum(coef(0, :)*dqle22(ind_begin_interp:ind_end_interp))
        dae22_res_time(time_ind) = sum(coef(0, :)*dae22(ind_begin_interp:ind_end_interp))
        bif_criterion(time_ind) = dqle22_res_time(time_ind)/dae22_res_time(time_ind)

        if (bif_criterion(time_ind) .ge. 1.0d0) then
            write(*,*) "!!! Bifurcation criterion met !!!"
            write(*,*) "bif_criterion = ", bif_criterion(time_ind)
        end if

        br_abs_time(time_ind) = time
        br_abs_antenna_factor(time_ind) = antenna_factor

    end subroutine


    subroutine message_Br_Dqle_values

        use wave_code_data, only: antenna_factor

        implicit none

        write(*,*) " "
        write(*,*) "+ + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + +"
        WRITE(*,'(A9,F10.4,A16,I4)') '  time = ', br_abs_time(time_ind), " s, time_ind = ", time_ind
        WRITE(*,'(A23,E10.5,A11,F6.2,A12)') '    Antenna_factor   = ', antenna_factor, " which are ", &
            antenna_factor/antenna_factor_max*100, "% of the max"
        WRITE(*,'(A23,F10.5,A2)') '    Br abs res * C_mn= ', br_abs(time_ind), " G"
        WRITE(*,'(A23,F10.5,A2)') '    Br abs res       = ', br_abs(time_ind)/SQRT(antenna_factor), " G"
        WRITE(*,'(A23,F10.3,A7)') '    Dqle22 res       = ', dqle22_res_time(time_ind), " cm^2/s"
        WRITE(*,'(A23,F10.5)')    '    bif crit         = ', bif_criterion(time_ind)
        write(*,*) '   Form factor      = ', abs(br_formfactor(time_ind))
        write(*,*) "+ + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + +"
        write(*,*) " "

    end subroutine

    subroutine stop_if_time_step_too_small

        use h5mod
        use parallelTools, only: ierror
        use mpi

        implicit none

        character(*), parameter :: reason = 'timestep < stop time step'

        if (timstep .lt. stop_time_step .and. time .gt. 1.0d-3) then
            write(*,*) 'stop: timestep smaller than stop limit'
            if (suppression_mode .eqv. .false.) then
                CALL write_kin_prof_data_to_disk
            end if

            call write_reason_for_stop_to_h5(reason)
            call MPI_finalize(ierror);
            stop
        end if

    end subroutine


    subroutine calc_params_num_and_denom

        use plasma_parameters, only: params_num, params_denom, params, params_beg

        implicit none

        params_num = (params - params_beg)**2
        params_denom = params**2 + params_beg**2

    end subroutine

    subroutine smooth_params_num_and_denom
    
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

    subroutine determine_timscal

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

    subroutine rescale_time_step_array

        use grid_mod, only: npoi, nbaleqs
        use recstep_mod, only: tim_stack, timstep_arr, tol
        use QLbalance_diag, only: timscal_dql

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

    subroutine set_time_step

        use recstep_mod, only: timstep_arr, tol
        use time_evolution, only: set_constant_time_step, constant_time_step
        
        implicit none

        timstep = minval(timstep_arr)

        if (.not. set_constant_time_step) then
            ! limit time step from below:
            timstep = max(timstep, timstep_min)
            if (reduce_time_step) then
                timstep = timstep * 0.3
            end if
            ! limit timestep from above:
            !if (ramp_up_mode .ne. 0) timstep = min(timstep,0.1)
            !timstep = min(timstep,0.005)
        else
        ! use for constant time step:
            timstep = constant_time_step
            write(*,*) "constant time step = ", timstep
        end if

        if (irank .eq. 0) then
            if (debug_mode) write(*,*) 'Debug: timstep', real(timstep), '   timescale', real(timescale), &
                'tolerance', real(tol)
        end if

    end subroutine

    subroutine reset_timstep_arr_w_timstep
        
        use recstep_mod, only: tim_stack, timstep_arr

        implicit none
        
        timstep_arr = 0.0d0
        timstep_arr = timstep
        tim_stack = timstep_arr

    end subroutine

    subroutine write_time_info

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

        use QLbalance_diag, only: rate_dql, timscal_dql

        implicit none

        open (4321, file='timstep_evol.dat', position='append')
        write (4321, *) time_ind, timstep, timscal_dql, timscal(1), rate_dql, time
        close (4321)

    end subroutine

    subroutine writeTimeInfoToH5

        use QLbalance_diag, only: rate_dql, timscal_dql

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
            CALL h5_append_double_1(time_dataset_id, (/time_ind*1.d0, timstep, &
                                                        timscal_dql, timscal(1), rate_dql, time/), time_ind)
        else
            CALL h5_append_double_1(time_dataset_id, (/time_ind*1.d0, timstep, &
                                                        time/), time_ind)
        end if
        CALL h5_close(h5_id)
        CALL h5_deinit()

    end subroutine

    subroutine relax_plasma_parameters

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
            write(*,*) 'Debug: i = ', int2(time_ind), 'time = ', real(time)
            write(*,*) ' '
        end if

    end subroutine

    subroutine set_first_iteration_true

        implicit none

        if (firstiterationdone .eqv. .false.) firstiterationdone = .true.

    end subroutine


    subroutine determine_Dql_diagnostic

        use grid_mod, only: dqle11, dqli11
        use QLbalance_diag

        implicit none

        timscal_dql = maxval(abs(dqle11_prev - dqle11))/maxval(dqle11_prev + dqle11)
        ind_dqle = maxloc(abs(dqle11_prev - dqle11))
        timscal_dqli = maxval(abs(dqli11_prev - dqli11))/maxval(dqli11_prev + dqli11)
        ind_dqli = maxloc(abs(dqli11_prev - dqli11))
        rate_dql = timscal_dql/timstep

    end subroutine

    subroutine create_group_structure_timeevol

        use control_mod
        use wave_code_data, only: m_vals, n_vals
        use resonances_mod, only: numres
        use h5mod

        implicit none

        if (debug_mode) write (*, *) "Debug: Creating group structure for TimeEvol"

        if (numres .eq. 1) then
            write (h5_mode_groupname, "(A,I1,A,I1)") "f_", m_vals(1), "_", n_vals(1)
        else
            write (h5_mode_groupname, "(A,I1,A,I1)") "multi_mode"
        end if

        CALL h5_init()
        CALL h5_open_rw(path2out, h5_id)

        if (.not. suppression_mode) then
            if (debug_mode) write(*,*) "Debug: h5_mode_groupname ", trim(h5_mode_groupname)
            CALL h5_create_parent_groups(h5_id, trim(h5_mode_groupname) //'/')
            CALL h5_create_parent_groups(h5_id, trim(h5_mode_groupname)//"/KinProfiles/")

            call create_group_if_not_existent(trim(h5_mode_groupname)//"/LinearProfiles/")
            call create_group_if_not_existent("/init_params")
        else
            if (debug_mode) write (*,*) "Debug: h5_mode_groupname: ", trim(h5_mode_groupname)
            call create_group_if_not_existent(trim(h5_mode_groupname))
            call create_group_if_not_existent("/init_params")
        end if

        CALL h5_close(h5_id)
        CALL h5_deinit()

        if (debug_mode) write (*, *) "Debug: finished creating group structure for TimeEvol"
    end subroutine

    subroutine redoTimeStep

        use parallelTools, only: irank
        use recstep_mod, only: timstep_arr
        use grid_mod, only: npoic, rc, Ercov
        use plasma_parameters, only: params, params_begbeg
        use baseparam_mod, only: eV, factolmax
        use restart_mod

        implicit none

        integer :: ipoi

        if (irank .eq. 0) then
            print *, 'redo step with old DQL'
        end if

        call hold_prev_transp_coeffs
        iunit_redo = 137

        if (irank .eq. 0) then
            open (iunit_redo, file='params_redostep.after')
            do ipoi = 1, npoic
                write (iunit_redo, *) rc(ipoi), params(1:2, ipoi) &
                    , params(3, ipoi)/ev &
                    , params(4, ipoi)/ev &
                    , 0.5d0*(Ercov(ipoi) + Ercov(ipoi + 1))
            end do
            close (iunit_redo)
        end if
        params = params_begbeg
        if (irank .eq. 0) then
            open (iunit_redo, file='params_redostep.before')
            do ipoi = 1, npoic
                write (iunit_redo, *) rc(ipoi), params(1:2, ipoi) &
                    , params(3, ipoi)/ev &
                    , params(4, ipoi)/ev &
                    , 0.5d0*(Ercov(ipoi) + Ercov(ipoi + 1))
            end do
            close (iunit_redo)
        end if

        timstep = timstep/factolmax
        timstep_arr = timstep

    end subroutine

end module