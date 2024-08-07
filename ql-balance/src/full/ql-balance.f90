!> @file
!> This is the main program file.

!> @details This program runs the balance code that solves the balance equations described in Heyn et al. NF2014.
!> The solution of the balance code includes envoking KiLCA. Note that KiLCA will read the profiles still from
!> ASCII files, in contrast to the balance code, which can read them from a HDF5 file.
program ql_balance

    use grid_mod
    use baseparam_mod
    use control_mod
    use h5mod 
    use wave_code_data
    use recstep_mod, only: nstack, tol, tim_stack, y_stack, timstep_arr
    use diag_mod
    use resonances_mod, only: numres
    use hdf5_tools
    use paramscan_mod
    use mpi
    use time_evolution
    use linear_run
    use parallelTools
    use restart_mod
    use PolyLagrangeInterpolation
    use plasma_parameters

    implicit none

    integer :: np_num

    logical :: discr_reached = .false. ! variable to say if discrepancy to linear regression
    ! of Br is reached, only used if br_stopping = .false.
    integer :: ipoi, i, ieq, l, k
    integer :: ioddeven
    double precision :: evoltime, timescale, tmax

    DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE :: ych_one, ych_tot
    double precision, dimension(:), allocatable :: timscal
    
    integer :: lb, ub
    
    DOUBLE PRECISION :: br_beta = 0
    DOUBLE PRECISION :: br_predicted
    
    integer(HID_T) :: time_dataset_id !> variable to save the time dataset id

    call read_config

    iexit = 0 ! 0 - don't skip, 1 - skip, 2 - stop
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

    call MPI_Init(ierror);
    call MPI_Comm_size(MPI_COMM_WORLD, np_num, ierror);
    call MPI_Comm_rank(MPI_COMM_WORLD, irank, ierror);

    if (irank .eq. 0) then
        write(*,*) ' '
        write(*,*) '******************************'
        write(*,*) 'number of processes:', np_num
        write(*,*) '              irank:', irank
        write(*,*) '******************************'
    end if


    timescale = (rmax - rmin)**2/dperp
    tmax = timescale*tmax_factor
    timstep = tmax/Nstorage
    if (irank .eq. 0) then
        write(*,*) "timstep = ", timstep
    end if

    call gengrid

    ! boundary condition
    if (iboutype .eq. 1) then
        npoi = npoic - 1
    else
        npoi = npoic
    end if

    CALL initialize_wave_code_interface(npoib, rb);
    CALL initialize_parameter_scan_vars

    mode_m = m_vals(1)
    mode_n = n_vals(1)
    if (debug_mode) write(*,*) 'Debug: mode_m = ', mode_m, 'mode_n = ', mode_n

    if (ihdf5IO .eq. 1) then
        CALL creategroupstructure
    end if

    call allocate_prev_variables

    call init_background_profiles

    if (irank .eq. 0) then
        CALL write_init_profiles
        call alloc_hold_parameters
    end if

    ! parameter scan loops that span over (nearly) the rest of the code
    do ifac_n = 1, size(fac_n)
        do ifac_Te = 1, size(fac_Te)
            do ifac_Ti = 1, size(fac_Ti)
                do ifac_vz = 1, size(fac_vz)
                    if (paramscan) then
                        ! change parameter scan string used for navigating the hdf5 file, only if
                        ! the suppression_mode is not activated
                        if (.not. suppression_mode) then
                            write (parscan_str, "(A,F0.3,A,F0.3,A,F0.3,A,F0.3,A)") "n", fac_n(ifac_n), &
                                "Te", fac_Te(ifac_Te), "Ti", fac_Ti(ifac_Ti), "vz", fac_vz(ifac_vz) &
                                , "/"
                        else
                            parscan_str = ""
                        end if
                    else
                        write(*,*) " "
                        write(*,*) "- No parameter scan -"
                        write(*,*) " "
                        ! leave it empty if no parameter scan
                        parscan_str = ""
                    end if

                    write (h5_mode_groupname, "(A)") trim(parscan_str)
                    write(*,*) "h5_mode_groupname: ", trim(h5_mode_groupname)

                    ! allocate variables in first parameter scan loop iteration
                    !if (ifac_n + ifac_Ti + ifac_Te + ifac_vz .eq. 4) then
                    !    call allocate_prev_variables
                    !end if
				
                    ! if more than one RMP mode is used, use different group name
                    if (numres .eq. 1) then
                        write (h5_mode_groupname, "(A,A,I1,A,I1)") trim(h5_mode_groupname), &
                            "f_", m_vals(1), "_", n_vals(1)
                    else
                        write (h5_mode_groupname, "(A,A,I1,A,I1)") trim(h5_mode_groupname), &
                            "multi_mode"
                    end if

                    write(*,*) "h5_mode_groupname after f_m_n: ", trim(h5_mode_groupname)

                    call rescale_profiles
                    call geomparprof

                    irf = 2 ! initialize dql variables and set to zero !
                    call get_dql

                    if (flag_run_time_evolution) then
                        call initAntennaFactor
                    end if

                    call rescaleTranspCoefficientsByAntennaFac
                    
                    irf = 1

                    call genstartsource

					write(*,*) "h5_mode_groupname before writeKinProfileDataToDisk: ", trim(h5_mode_groupname)
                    if (irank .eq. 0) then
                        if (suppression_mode .eqv. .false.) then
                            CALL writeKinProfileDataToDisk(0) ! write the profiles to hdf5 file
                        end if
                    end if

                    time = 0.d0
                    tol = tol_max

                    call InquiryToRestart

                    timstep_arr = timstep
                    tim_stack = timstep_arr

                    iunit_diag = 5000

                    call get_dql ! also writes out diffusion coefficients and other data
                    call rescaleTranspCoefficientsByAntennaFac

                    if (flag_run_time_evolution) then
                        if (ifac_n + ifac_Ti + ifac_Te + ifac_vz .eq. 4) then
                            call allocateBrAndDqleForTimeEvolution
                        end if
                    end if

                    call savePrevTranspCoefficients

                    if (ifac_n + ifac_Ti + ifac_Te + ifac_vz .eq. 4) then
                        allocate(timscal(npoi), dummy(npoic))
                        allocate(params_beg(nbaleqs, npoic), params_num(nbaleqs, npoic))
                        allocate(params_denom(nbaleqs, npoic))
                        allocate(params_begbeg(nbaleqs, npoic))
                    end if

                    params_begbeg = params

                    if (.not. flag_run_time_evolution) then
                        ! linear run
                        ! if velocity scan, determine Er_res for v_ExB velocity at resonant surface
                        call interpBrAndDqlAtResonanceParamScan

                        if (paramscan) then
                            ! if the last parameter scan is done, write data and stop the code
                            if (ifac_n + ifac_Te + ifac_Ti + ifac_vz .eq. size(fac_n) + size(fac_Ti) + &
                                size(fac_Te) + size(fac_vz)) then
                                if (debug_mode) write(*,*) "Debug: Last parameter done. Finalize MPI"

                                call writeBrAndDqlAtResonanceToH5
                                CALL deallocate_wave_code_data()
                                CALL MPI_finalize(ierror);
                                stop
                            else
                                ! if it is not the last scan, skip the rest of the code and continue
                                ! with the next loop iteration
                                CYCLE
                            end if
                        else
                            !Stop if mode is not time evolution
                            call finalizeLinearRun
                        end if
                    end if

                    iunit_diag = 137
                    iunit_diag_b = 8138
                    ioddeven = 1
                    
                    ! time evolution

                    if (ihdf5IO .eq. 0) then
                        ! sweep files that are only appended to
                        open(4321, file='timstep_evol.dat', status='replace')
                        close(4321)
                        open(777, file='br_abs_res.dat', status='replace')
                        close(777)
                    end if

                    do i = 1, Nstorage ! loop over time steps
                        timeStep = i
                        write (*, *) "Time Step = ", timeStep
                        
                        call saveKinProfilesToYPrev

                        redostep = .false.

                        if (irank .eq. 0) then
                            iunit_diag = 5000 + i
                            if (write_diag_b) then
                                if (ioddeven/2*2 .eq. ioddeven) then
                                    open (iunit_diag_b, file='params_b_redostep.even')
                                else
                                    open (iunit_diag_b, file='params_b_redostep.odd')
                                end if
                            end if
                        end if

                        ! in get_dql the fort.5000 data is written. The argument is used
                        ! to restrict the writing of the data. Only every "save_prof_time_step"th
                        ! step the data is written
                        call get_dql
                        call stopIfTimeStepTooSmall
                        call interpBrAndDqlAtResonanceTimeEvol
						call write_br_time_data
                        call rescaleTranspCoefficientsByAntennaFac

                        if (write_diag_b) close (iunit_diag_b)

                        call writefort9999
                        
                        if (.not. redostep) then
                            call savePrevTranspCoefficients
                            params_begbeg = params
                            ioddeven = ioddeven + 1
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
                                write(*,*) 'redo step'
                            end if
                        end do ! end of redo step loop
                        
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
                        !
                        timstep_arr = timstep_arr*timescale/(timstep_arr + timescale)
                        if (scratch) then
                            scratch = .false.
                            tim_stack = timstep_arr
                        end if
                        timstep_arr = 2.d0*timstep_arr*tim_stack/(timstep_arr + tim_stack)
                        !timstep = minval(timstep_arr, MASK = r .gt. rsepar - 0.5d0)
                        timstep = minval(timstep_arr)

                        !This can be used to limit timestep
                        !Added by Philipp Ulbl, June 2020
                        !if(time .gt. 1.2 .and. timstep .gt. 1.d-4) then
                        !    timstep = 1.d-4
                        !endif
                        ! the rimstep_min variable was added to the namelist, the file
                        ! timstep_min.inp is redundant now^A&M2B5* d3v F MB(e)A(v)
                        !open(5432,file='timstep_min.inp')
                        !read (5432,*) timstep_min
                        !close(5432)
                        ! non-constant time step:
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
                        
                        timstep_arr = timstep
                        !
                        tim_stack = timstep_arr
                        !
                        
                        if (irank .eq. 0) then
                            write(*,*) 'timstep', real(timstep), '   timescale', real(timescale), &
                                'tolerance', real(tol)
                        end if
                        !
                        if (irank .eq. 0) then
                            if (ihdf5IO .eq. 1) then
                                ! write timstep_evol data to hdf5 file
                                h5_currentgrp = trim("/"//trim(h5_mode_groupname)// &
                                                     "/timstep_evol.dat")

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
                                    CALL h5_append_double_1(time_dataset_id, (/i*1.d0, timstep, &
                                                                               timscal_dql, timscal(1), rate_dql, time/), i)
                                else
                                    CALL h5_append_double_1(time_dataset_id, (/i*1.d0, timstep, &
                                                                               time/), i)
                                end if
                                CALL h5_close(h5_id)
                                CALL h5_deinit()
                            else
                                open (4321, file='timstep_evol.dat', position='append')
                                write (4321, *) i, timstep, timscal_dql, timscal(1), rate_dql, time
                                close (4321)
                            end if
                        end if
                        

                        do ipoi = 1, npoi
                            do ieq = 1, nbaleqs
                                k = nbaleqs*(ipoi - 1) + ieq
                                params(ieq, ipoi) = yprev(k)*urelax + params(ieq, ipoi)*(1.d0 - urelax)
                            end do
                        end do
                        timstep_arr = 0.d0
                        call evolvestep(timstep, eps)
                        timstep_arr = timstep
                        !
                        time = time + timstep
                        !
                        !
                        if (irank .eq. 0) then
                            write(*,*) 'i = ', int2(i), 'time = ', real(time)
                            write(*,*) ' '
                        end if
                        !
                        !for debugging:
                        if (irank .eq. 0) then
                            if (modulo(i, save_prof_time_step) .eq. 0) then
                                if (suppression_mode .eqv. .false.) then
                                    CALL writeKinProfileDataToDisk(i)
                                end if
                            end if
                        end if
!
                        if (firstiterationdone .eqv. .false.) firstiterationdone = .true.

                        ! check if linear discepancy in penetration ratio is reached
                        CALL linear_discrepancy_pen_ratio
 
                        if (iexit .eq. 1) then
                            iexit = 0
                            EXIT
                        end if                       

                        ! ramp-up RMP coil current
                        CALL ramp_coil(i)

                        if (iexit .eq. 1) then
                            iexit = 0
                            EXIT
                        end if
                    end do ! end of time evol loop



                    if (debug_mode) write(*,*) 'Debug: deallocate data for next parameter scan'
                    call deallocate_wave_code_data();
!deallocate(yprev)
                    deallocate (coef)
!deallocate(tim_stack)
!deallocate(timscal,params_beg,params_num,params_denom,dummy)

! end of loops of the parameter scan
                end do
            end do
        end do
    end do

    write (*, *) 'Programm is finalized without stopping criterion met';
    call MPI_finalize(ierror)

contains

!> @brief subroutine write_init_profiles. Write initial profiles to hdf5 or ascii.
!> @author Markus Markl
!> @date 05.10.2022
subroutine write_init_profiles

    use plasma_parameters, only: params, qsaf
    use control_mod, only: debug_mode, ihdf5IO
    use h5mod, only: h5_exists_log, h5_id, path2out
    use wave_code_data, only: r

    implicit none

    if (debug_mode) write(*,*) "Debug: writing initial background profiles"
    if (ihdf5IO .eq. 1) then
        CALL h5_init()
        ! open hdf5 file
        CALL h5_open_rw(path2out, h5_id)
        CALL h5_obj_exists(h5_id, "/init_params/n", h5_exists_log)
        if (.not. h5_exists_log) then
            CALL h5_add_double_1(h5_id, "/init_params/n", &
                params(1, :), lbound(params(1, :)), ubound(params(1, :)))
            CALL h5_add_double_1(h5_id, "/init_params/Vz", &
                params(2, :), lbound(params(2, :)), ubound(params(2, :)))
            CALL h5_add_double_1(h5_id, "/init_params/Te", &
                params(3, :)/ev, lbound(params(3, :)), ubound(params(3, :)))
            CALL h5_add_double_1(h5_id, "/init_params/Ti", &
                params(4, :)/ev, lbound(params(4, :)), ubound(params(4, :)))
            CALL h5_add_double_1(h5_id, "/init_params/qsaf", &
                qsaf(:), lbound(qsaf(:)), ubound(qsaf(:)))
            CALL h5_add_double_1(h5_id, "/init_params/r", &
                r, lbound(r), ubound(r))
        else
            if (debug_mode) write(*,*) "Debug: they are already there -> skiping"
        end if

        CALL h5_close(h5_id)
        CALL h5_deinit()
        if (debug_mode) write(*,*) "Debug: finished writing initial background profiles"
                                !stop ! for test purposes

    else
        open (123, form='unformatted', file='init_params.dat')
        write (123) params
        close (123)
    end if

end subroutine ! write_initial_profiles



!@> brief Check if discrepancy between penetration ratio and linear extrapolation exceeds critical value
!> author> Markus Markl
!> created 13.03.2023
subroutine linear_discrepancy_pen_ratio

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
    if (i .gt. 50 .and. .not. discr_reached) then
        ! calculate beta only once
        if (br_beta .eq. 0) then
            ! Calculate slope from the data until this time step. (Simple linear regression algorithm)
            br_beta = sum(br_abs(3:i)*br_abs_time(3:i))/sum(br_abs_time(1:i)**2)
            write(*,*) 'br_beta = ', br_beta
        end if
        ! calculate the from the linear regression predicted value of Br_abs
        br_predicted = br_beta*br_abs_time(i)
        write(*,*) "Delta = ", abs(br_abs(i) - br_predicted)
        if (abs(br_abs(i) - br_predicted) .gt. 0.1) then
            write(*,*) 'discrepancy to linearly predicted value of Br_abs_res > delta'
            if (modulo(i, save_prof_time_step) .ne. 0) then
                if (suppression_mode .eqv. .false.) then
                    CALL writeFieldsCurrentsAndTranspCoeffsToH5
                end if
            end if
            if (suppression_mode .eqv. .false.) then
                CALL writeKinProfileDataToDisk(i)
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
                    '/discrep_time', (/i*1.d0, time/), (/1/), (/2/))
                CALL h5_close(h5_id)
                CALL h5_deinit()
                discr_reached = .true.
            end if
        end if
    end if

end subroutine


end program ql_balance
