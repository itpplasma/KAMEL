!> @file
!> This is the main program file.

!> @details This program runs the balance code that solves the balance equations described in Heyn et al. NF2014.
!> The solution of the balance code includes envoking KiLCA. Note that KiLCA will read the profiles still from
!> ASCII files, in contrast to the balance code, which can read them from a HDF5 file.
program ql_balance

    use grid_mod
    use baseparam_mod
    use control_mod
    use h5mod ! added by Markus Markl, 25.02.2021
    use wave_code_data
    use recstep_mod, only: nstack, tol, tim_stack, y_stack, timstep_arr
    !DIAG:
    use diag_mod, only: write_diag, iunit_diag, write_diag_b, iunit_diag_b
    !END DIAG
    use resonances_mod, only: numres
    use hdf5_tools
    use paramscan_mod
    use mpi
    use time_evolution
    use linear_run
    use parallelTools
    use restart_mod
    use PolyLagrangeInterpolation

    implicit none

    integer :: ierror, np_num
    integer :: iexit ! used for ramp-up skipping of saving

    logical :: firstiterationdone !Added by Markus Markl 25.02.2021. Some steps
    ! in saving the data to hdf5 file need to be done only the first time iteration
    logical :: discr_reached = .false. ! variable to say if discrepancy to linear regression
    ! of Br is reached, only used if br_stopping = .false.
    integer :: ipoi, i, ieq, l, k
    integer :: iunit_redo, ioddeven
    double precision :: evoltime, timescale, tmax
    double precision :: time
    double precision :: timscal_dql, rate_dql, timscal_dqli

    double precision :: urelax = 0.5e0 !0.5d0  !0.9d0
    double precision :: tol_max = 3.d-2 !3.d-4 !3.d-3 !3.d-2
    double precision :: factolmax = 3.d0 ! keep
    double precision :: factolred = 0.5d0 ! keep

    DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE :: ych_one, ych_tot
    double precision, dimension(:), allocatable :: timscal
    double precision, dimension(:), allocatable :: dummy
    
    double precision, dimension(:, :), allocatable :: params_beg, params_begbeg
    double precision, dimension(:, :), allocatable :: params_num, params_denom
    
    integer, dimension(1) :: ind_dqle, ind_dqli
    integer :: lb, ub
    
    DOUBLE PRECISION :: br_beta = 0
    DOUBLE PRECISION :: br_predicted

    integer :: ramp_up_down = 0 !> used in hysteresis mode, tells if ramp-up (0) or ramp-down (1)
    integer ::  ibrabsres, ibeg, iend
    
    character(len=1024) :: h5_currentgrp !> current hdf5 group string
    integer(HID_T) :: time_dataset_id !> variable to save the time dataset id
    double precision :: t_hysteresis_turn = 0


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

                    irf = 2
                    call get_dql(0)

                    if (flag_run_time_evolution) then
                        call initAntennaFactor
                    end if

                    dqle11 = dqle11*antenna_factor
                    dqle12 = dqle12*antenna_factor
                    dqle21 = dqle21*antenna_factor
                    dqle22 = dqle22*antenna_factor
                    dqli11 = dqli11*antenna_factor
                    dqli12 = dqli12*antenna_factor
                    dqli21 = dqli21*antenna_factor
                    dqli22 = dqli22*antenna_factor
                    irf = 1
!
                    call genstartsource

					write(*,*) "h5_mode_groupname before writefort1000: ", trim(h5_mode_groupname)
                    if (irank .eq. 0) then
                        if (suppression_mode .eqv. .false.) then
                            CALL writefort1000(0) ! write the profiles to hdf5 file
                        end if
                    end if

                    if (ifac_n + ifac_Ti + ifac_Te + ifac_vz .eq. 4) then
                        allocate(timscal(npoi), dummy(npoic))
                        allocate(params_beg(nbaleqs, npoic), params_num(nbaleqs, npoic))
                        allocate(params_denom(nbaleqs, npoic))
                        allocate(params_begbeg(nbaleqs, npoic))
                    end if
                    time = 0.d0
                    tol = tol_max

                    call InquiryToRestart

                    timstep_arr = timstep
                    tim_stack = timstep_arr

                    write(*,*) 'start balance, irank = ', irank
                    iunit_diag = 5000

                    call get_dql(0) ! also writes out diffusion coefficients and other data
                    call rescaleTranspCoefficientsByAntennaFac
!
                    if (flag_run_time_evolution) then
                        if (ifac_n + ifac_Ti + ifac_Te + ifac_vz .eq. 4) then
                            call allocateBrAndDqleForTimeEvolution
                        end if
                    end if

                    call savePrevTranspCoefficients

                    params_begbeg = params
                    if (debug_mode) write(*,*) 'Debug: dql ready'


                    if (.not. allocated(coef)) allocate (coef(0:nder, nlagr))
                    !binsearch, get index of r_resonant
                    call binsrc(rb, 1, npoib, r_resonant(1), ibrabsres)

                    ibeg = max(1, ibrabsres - nlagr/2)
                    iend = ibeg + nlagr - 1
                    if (iend .gt. npoib) then
                        iend = npoib
                        ibeg = iend - nlagr + 1
                    end if

                    if (.not. flag_run_time_evolution) then
                        ! linear run
                        call plag_coeff(nlagr, nder, r_resonant(1), rb(ibeg:iend), coef)
                        dqle22_res(ifac_n, ifac_Te, ifac_Ti, ifac_vz) = sum(coef(0, :) * dqle22(ibeg:iend))
                        br_abs_res_parscan(ifac_n, ifac_Te, ifac_Ti, ifac_vz) = sum(coef(0, :) &
                            * abs(Br(ibeg:iend)))*sqrt(antenna_factor)
                        ! if velocity scan, determine Er_res for v_ExB velocity at resonant surface
                        if (size(fac_vz) .ne. 1) then
                            if (debug_mode) write (*, *) "Debug: determine Er_res"
                            do ipoi = 1, npoic
                                Ercovavg(ipoi) = 0.5d0*(Ercov(ipoi) + Ercov(ipoi + 1))
                            end do
                            Er_res(ifac_n, ifac_Te, ifac_Ti, ifac_vz) = sum(coef(0, :)*Ercovavg(ibeg:iend))
                            if (debug_mode) write (*, *) "Debug: Er_res = ", Er_res(ifac_n, ifac_Te, ifac_Ti, ifac_vz)

                        end if
                        write(*,*) "dqle22 res = ", dqle22_res(ifac_n, ifac_Te, ifac_Ti, ifac_vz)
                        write(*,*) "Br abs res = ", br_abs_res_parscan(ifac_n, ifac_Te, ifac_Ti, ifac_vz)
						write(*,*) "Antenna factor = ", antenna_factor

                        if (paramscan) then
                            ! if the last parameter scan is done, write data and stop the code
                            if (ifac_n + ifac_Te + ifac_Ti + ifac_vz .eq. size(fac_n) + size(fac_Ti) + &
                                size(fac_Te) + size(fac_vz)) then
                                if (debug_mode) write(*,*) "Debug: Last parameter done. Finalize MPI"
                                write (h5_mode_groupname, "(A,I1,A,I1)") "f_", m_vals(1), "_", n_vals(1)


                                ! write the diffusion coefficient
                                !write (h5_mode_groupname, "(A,I1,A,I1)") "f_", &
                                !mode_m, "_", mode_n

                                if (debug_mode) write(*,*) "Debug: Write out results"

                                CALL h5_init()
                                CALL h5_open_rw(path2out, h5_id)
                                CALL h5_obj_exists(h5_id, trim(h5_mode_groupname), &
                                    h5_exists_log)
                                if (.not. h5_exists_log) then
                                    CALL h5_define_group(h5_id, &
                                        trim(h5_mode_groupname), group_id_2)
                                    CALL h5_close_group(group_id_2)
                                end if

                                CALL h5_add_double_1(h5_id, trim(h5_mode_groupname)//'/dqle22_res', &
                                                     reshape(dqle22_res, (/size(dqle22_res)/)), &
                                                     lbound(reshape(dqle22_res, (/size(dqle22_res)/))), &
                                                     ubound(reshape(dqle22_res, (/size(dqle22_res)/))))
                                CALL h5_add_double_1(h5_id, trim(h5_mode_groupname)//'/br_abs_res', &
                                                     reshape(br_abs_res_parscan, (/size(br_abs_res_parscan)/)), &
                                                     lbound(reshape(br_abs_res_parscan, (/size(br_abs_res_parscan)/))), &
                                                     ubound(reshape(br_abs_res_parscan, (/size(br_abs_res_parscan)/))))

                                if (size(fac_vz) .ne. 1) then
                                    CALL h5_add_double_1(h5_id, trim(h5_mode_groupname)//'/Er_res', &
                                                         reshape(Er_res, (/size(Er_res)/)), &
                                                         lbound(reshape(Er_res, (/size(Er_res)/))), &
                                                         ubound(reshape(Er_res, (/size(Er_res)/))))
                                end if

                                CALL h5_close(h5_id)
                                CALL h5_deinit()

                                CALL deallocate_wave_code_data()

                                CALL MPI_finalize(ierror);
                                stop
                            else
                                ! if it is not the last scan, skip the rest of the code and continue
                                ! with the next loop iteration
                                !call deallocate_wave_code_data();
                                CYCLE
                            end if
                        else
                            !Stop if mode is not time evolution
                            !Added by Philipp Ulbl 12.05.2020
                            write(*,*) 'stop: linear code only'
                            !write (h5_mode_groupname, "(A,I1,A,I1)") "f_", &
                            !    mode_m, "_", mode_n
                            if (debug_mode) write(*,*) "Debug: Write out results"

                            if (ihdf5IO .eq. 1) then
                                CALL h5_init()
                                CALL h5_open_rw(path2out, h5_id)
                                CALL h5_obj_exists(h5_id, trim(h5_mode_groupname), &
                                    h5_exists_log)
                                if (.not. h5_exists_log) then
                                    CALL h5_define_group(h5_id, &
                                        trim(h5_mode_groupname), group_id_2)
                                    CALL h5_close_group(group_id_2)
                                end if
!
                                CALL h5_add_double_1(h5_id, trim(h5_mode_groupname)//'/dqle22_res', &
                                                 reshape(dqle22_res, (/size(dqle22_res)/)), &
                                                 lbound(reshape(dqle22_res, (/size(dqle22_res)/))), &
                                                 ubound(reshape(dqle22_res, (/size(dqle22_res)/))))
                                CALL h5_add_double_1(h5_id, trim(h5_mode_groupname)//'/dqle22', &
                                                     dqle22, &
                                                     lbound(dqle22), &
                                                     ubound(dqle22))
 
                                CALL h5_close(h5_id)
                                CALL h5_deinit()
                            end if
                            call MPI_finalize(ierror);
                            stop  !! <<----- Stop for linear code usage
                        end if
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
                    end if

                    iunit_diag = 137
                    iunit_diag_b = 8138
                    ioddeven = 1
! #########################################################################################
! Time evolution
!

                    !if (flag_run_time_evolution) then
                    !    if (ramp_up_mode .eq. 5) then
                    !        ramp_up_down = 0 ! ramp up
                    !        t_hysteresis_turn = 0
                    !    end if
                    !end if
                    firstiterationdone = .false. ! if first iteration is done, some variables are already allocated

                    if (ihdf5IO .eq. 0) then
                        ! sweep files that are only appended to
                        open(4321, file='timstep_evol.dat', status='replace')
                        close(4321)
                        open(777, file='br_abs_res.dat', status='replace')
                        close(777)
                    end if

                    do i = 1, Nstorage ! loop over time steps
                        step_counter = i
                        write (*, *) "i = ", i
                        !
                        if (debug_mode) write(*,*) 'Debug: yprev loop'
                        do ipoi = 1, npoi
                            do ieq = 1, nbaleqs
                                k = nbaleqs*(ipoi - 1) + ieq
                                yprev(k) = params(ieq, ipoi)
                            end do
                        end do
                        !
                        !
                        dostep = .true.
                        !
                        !do while(dostep)

                        if (irank .eq. 0) then
                            !
                            !DIAG:
                            iunit_diag = 5000 + i
                            if (write_diag_b) then
                                if (ioddeven/2*2 .eq. ioddeven) then
                                    open (iunit_diag_b, file='params_b_redostep.even')
                                else
                                    open (iunit_diag_b, file='params_b_redostep.odd')
                                end if
                            end if
                            !END DIAG
                            !
                        end if

                        ! in get_dql the fort.5000 data is written. The argument is used
                        ! to restrict the writing of the data. Only every "save_prof_time_step"th
                        ! step the data is written
                        if (debug_mode) write(*,*) "Debug: get_dql(", i, ")"
                        call get_dql(i)
                        if (debug_mode) write(*,*) "Debug: coming out get_dql(", i, ")"

                        ! stopping criterion
                        !Stop if timestep becomes too small
                        !Added by Philipp Ulbl 13.05.2020
                        if (timstep .lt. stop_time_step .and. time .gt. 1.0d-3) then
                            write(*,*) 'stop: timestep smaller than stop limit'
                            if (suppression_mode .eqv. .false.) then
                                CALL writefort1000(i)
                            end if
                           if (ihdf5IO .eq. 1) then
                                CALL h5_init()
                                CALL h5_open_rw(path2out, h5_id)
                                CALL h5_add_string(h5_id, trim(h5_mode_groupname)// &
                                '/stopping_criterion', 'timestep < stop time step')
                                CALL h5_close(h5_id)
                                CALL h5_deinit()
                            end if

                            call MPI_finalize(ierror);
                            stop
                        end if

                        !calculate Br abs at the resonant surface for stopping criterion
                        !Added by Philipp Ulbl 04.06.2020

                        !binsearch - is also already done before time evolution
                        call binsrc(rb, 1, npoib, r_resonant(1), ibrabsres)
                        !
                        ibeg = max(1, ibrabsres - nlagr/2)
                        iend = ibeg + nlagr - 1
                        if (iend .gt. npoib) then
                            iend = npoib
                            ibeg = iend - nlagr + 1
                        end if

                        !lagrange interpolation with order 4 only for function (0)
                        call plag_coeff(nlagr, nder, r_resonant(1), rb(ibeg:iend), coef)
                        if (debug_mode) then
                            write(*,*) "Debug: nlagr = ", nlagr
                            write(*,*) "Debug: nder = ", nder
                            write(*,*) "Debug: r_resonant = ", r_resonant
                            write(*,*) "Debug: rb(ibeg:iend) = ", rb(ibeg:iend)
                            write(*,*) "Debug: coef = ", coef
                        end if
                        br_abs(i) = sum(coef(0, :)*abs(Br(ibeg:iend)))*sqrt(antenna_factor)
                        dqle22_res_time(i) = sum(coef(0, :)*dqle22(ibeg:iend))

                        !output on console and save to file

                        ! save the time for the improved stopping criterion
                        br_abs_time(i) = time
						br_abs_antenna_factor(i) = antenna_factor


                        ! save the data
                        write(*,*) 'Br abs res * C_mn= ', br_abs(i)
                        write(*,*) 'Br abs res       = ', br_abs(i)/sqrt(antenna_factor)
                        write(*,*) 'Dqle22 res       = ', dqle22_res_time(i)
                        write(*,*) 'Antenna factor   = ', antenna_factor
                        write(*,*) 'time = ', br_abs_time(i)

                        ! most important data to be saved

						CALL write_br_time_data(i, br_abs_time, br_abs_antenna_factor, br_abs, dqle22_res_time)

                        ! check if linear discepancy in penetration ratio is reached
                        !CALL linear_discrepancy_pen_ratio
 
                        !if (iexit .eq. 1) then
                        !    iexit = 0
                        !    EXIT
                        !end if                       

                        ! ramp-up RMP coil current
                        !CALL ramp_coil

                        !if (iexit .eq. 1) then
                        !    iexit = 0
                        !    EXIT
                        !end if
                        !Old ramp up by Martin: ramp up to 1s then ramp down to 2s
                        !to check hysteresis behaviour. afterwards small increases
                        !if (time.lt.1.d0) then
                        !    antenna_factor = time**2
                        !elseif (time.gt.1.d0.AND.time.lt.2.d0 ) then
                        !    antenna_factor = (2.d0 - time)**2
                        !endif
                        !antenna_factor = antenna_factor + 1.d-4

                        dqle11 = dqle11*antenna_factor
                        dqle12 = dqle12*antenna_factor
                        dqle21 = dqle21*antenna_factor
                        dqle22 = dqle22*antenna_factor
                        dqli11 = dqli11*antenna_factor
                        dqli12 = dqli12*antenna_factor
                        dqli21 = dqli21*antenna_factor
                        dqli22 = dqli22*antenna_factor

                        !
                        !DIAG:
                        if (write_diag_b) close (iunit_diag_b)
                        !END DIAG
                        timscal_dql = maxval(abs(dqle11_prev - dqle11))/maxval(dqle11_prev + dqle11)
                        ind_dqle = maxloc(abs(dqle11_prev - dqle11))
                        timscal_dqli = maxval(abs(dqli11_prev - dqli11))/maxval(dqli11_prev + dqli11)
                        ind_dqli = maxloc(abs(dqli11_prev - dqli11))
                        rate_dql = timscal_dql/timstep
                        ! write fort.9999 data if diagnostics output is done
                        if (diagnostics_output) then
                            if (irank .eq. 0) then
                                print *, 'timscal_dqle = ', sngl(timscal_dql) &
                                    , 'timscal_dqli = ', sngl(timscal_dqli)
                                print *, 'maximum dqle at r = ', rc(ind_dqle(1)) &
                                    , 'maximum dqli at r = ', rc(ind_dqli(1))
                                ! Edited by Markus Markl, 26.02.2021
                                if (ihdf5IO .eq. 1) then
                                    ! write fort.9999 data to hdf5 file
                                    h5_currentgrp = trim("/"//trim(h5_mode_groupname) &
                                                         //"/fort.9999")
                                    CALL h5_init()
                                    CALL h5_open_rw(path2out, h5_id)
                                    CALL h5_obj_exists(h5_id, trim(h5_currentgrp), h5_exists_log)
                                    if (h5_exists_log) then
                                        CALL h5_delete(h5_id, trim(h5_currentgrp))
                                    end if

                                    CALL h5_define_unlimited_matrix(h5_id, trim(h5_currentgrp), &
                                                                    H5T_NATIVE_DOUBLE, (/-1, 3/), dataset_id)
                                    CALL h5_append_double_1(dataset_id, rb, 1)
                                    CALL h5_append_double_1(dataset_id, abs(dqle11_prev - dqle11), 2)
                                    CALL h5_append_double_1(dataset_id, abs(dqli11_prev - dqli11), 3)

                                    CALL h5_close(h5_id)
                                    CALL h5_deinit()

                                else
                                    do ipoi = 1, npoib
                                        write (9999, *) rb(ipoi), abs(dqle11_prev(ipoi) - dqle11(ipoi)), &
                                            abs(dqli11_prev(ipoi) - dqli11(ipoi))
                                    end do
                                    close (9999)
                                end if
                            end if
                        end if
                        !    timscal_dql=timscal_dql+timscal_dqli
                        !    if(timscal_dql.lt.tol*factolmax) then
                        if (.true.) then
                            dostep = .false.
                            dqle11_prev = dqle11
                            dqle12_prev = dqle12
                            dqle21_prev = dqle21
                            dqle22_prev = dqle22
                            dqli11_prev = dqli11
                            dqli12_prev = dqli12
                            dqli21_prev = dqli21
                            dqli22_prev = dqli22
                            params_begbeg = params
                            ioddeven = ioddeven + 1
                        else
                            if (irank .eq. 0) then
                                print *, 'redo step with old DQL'
                            end if
                            dqle11 = dqle11_prev
                            dqle12 = dqle12_prev
                            dqle21 = dqle21_prev
                            dqle22 = dqle22_prev
                            dqli11 = dqli11_prev
                            dqli12 = dqli12_prev
                            dqli21 = dqli21_prev
                            dqli22 = dqli22_prev
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
                        end if
                        !
                        do ! redo step loop
                            !
                            params_beg = params
                            !
                            write(*,*) "timstep before evolvestep is ", timstep
                            write(*,*) "eps before evolvestep is ", eps
                            call evolvestep(timstep, eps)
                            write(*,*) "timstep after evolvestep is ", timstep
                            !
                            !limits ion and electron temperatures from below (by 10 eV in this example).
                            !Added by Philipp Ulbl on 09.06.2020

                            !write(*,*) "Limit Te and Ti from below"
                            do ipoi = 1, npoic
                                if (.true.) then
                                    params(3, ipoi) = max(params(3, ipoi), temperature_limit*ev)
                                    params(4, ipoi) = max(params(4, ipoi), temperature_limit*ev)
                                else 
                                    !> Quick fix of steady state solution. Keep boundary inside the separatrix.
                                    if (r(ipoi) > rsepar-0.5d0) then
                                        params(3, ipoi) = hold_Te(ipoi)
                                        params(4, ipoi) = hold_Ti(ipoi)
                                    end if
                                end if
                            end do
                            !
                            params_num = (params - params_beg)**2
                            !write(1234,*) params_num
                            params_denom = params**2 + params_beg**2
                            !write(1235, *) params_denom
                            !
                            do ieq = 1, nbaleqs
                                !
                                call smooth_array_gauss(npoi, mwind, params_num(ieq, :), dummy)
                                !
                                params_num(ieq, :) = dummy
                                !
                                call smooth_array_gauss(npoi, mwind, params_denom(ieq, :), dummy)
                                !
                                params_denom(ieq, :) = dummy
                            end do
                            !
                            do ipoi = 1, npoi
                                if (rc(ipoi) .lt. 0.95d0*rc(npoi)) then
                                    timscal(ipoi) = sum(sqrt(params_num(3:4, ipoi)/params_denom(3:4, ipoi)))
                                    !write(*,*) "timscal(", ipoi,")= ", timscal(ipoi)
                                else
                                    timscal(ipoi) = 1d-30 !0.d0
                                end if
                            end do
                            !
                            if (irank .eq. 0) then
                                write(*,*) "maxval(timscal) = ", maxval(timscal)
                                write(*,*) "tol*factolmax = ", tol*factolmax
                            end if
                            if (maxval(timscal) .lt. tol*factolmax) exit
                            !
                            timstep_arr = timstep_arr*factolred
                            params = params_beg
                            if (irank .eq. 0) then
                                write(*,*) 'redo step'
                            end if
                        end do ! end of redo step loop
                        !
                        timscal = timscal + timscal_dql
                        !
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
                        !
                        !enddo ! do while(dostep)
                        !
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
                                    CALL writefort1000(i)
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
                        CALL ramp_coil

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

!> @brief subroutine creategroupstructure. Creates the group structure in the hdf5 file.
!> @author  Markus Markl
!> @date 12.03.2021
subroutine creategroupstructure

    use h5mod
    use paramscan_mod
    use control_mod
    use wave_code_data, only: m_vals, n_vals
    use resonances_mod, only: numres

    implicit none
    !logical :: suppression_mode = .true.

    write (*, *) "Creating group structure"
    ! if the profiles should be written out, i.e. suppression_mode = false, then an extended group
    ! structure is created to save them
    if (.not. suppression_mode) then
        CALL h5_init()
        CALL h5_open_rw(path2out, h5_id)
        do ifac_n = 1, size(fac_n)
            do ifac_Te = 1, size(fac_Te)
                do ifac_Ti = 1, size(fac_Ti)
                    do ifac_vz = 1, size(fac_vz)
                        write (*, *) ifac_n, ifac_Te, ifac_Ti, ifac_vz
                        if (paramscan) then
                            ! change parameter scan string used for the
                            !group structure in hdf5 file
                            write (parscan_str, "(A,F0.3,A,F0.3,A,F0.3,A,F0.3)") &
                                "n", fac_n(ifac_n), "Te", fac_Te(ifac_Te), &
                                "Ti", fac_Ti(ifac_Ti), "vz", fac_vz(ifac_vz)
                            
                            if (numres .eq. 1) then
                                write (h5_mode_groupname, "(A,A,A,I1,A,I1)") &
                                    trim(parscan_str), "/", "f_", m_vals(1), &
                                    "_", n_vals(1)
                            else
                                write (h5_mode_groupname, "(A,A,A,I1,A,I1)") &
                                    trim(parscan_str), "/", "multi_mode"
                            end if

                        else
                            ! leave it empty if no parameter scan
                            parscan_str = ""
                            ! if more than one RMP mode is used, use different group name
                            if (numres .eq. 1) then
                                write (h5_mode_groupname, "(A,I1,A,I1)") &
                                    "f_", m_vals(1), "_", n_vals(1)
                            else
                                write (h5_mode_groupname, "(A,I1,A,I1)") &
                                    "multi_mode"
                            end if
                        end if
                        ! create the groups that are furthest down: fort.1000,
                        ! fort.5000 and init_params
                        write(*,*) "h5_mode_groupname ", trim(h5_mode_groupname)
                        CALL h5_create_parent_groups(h5_id, trim(h5_mode_groupname) &
                                                     //'/')

                        if (suppression_mode .eqv. .false.) then
                            CALL h5_create_parent_groups(h5_id, &
                                                         trim(h5_mode_groupname)//"/fort.1000/")
                            CALL h5_define_group(h5_id, &
                                                 trim(h5_mode_groupname)//"/fort.5000/", group_id_1)
                            CALL h5_close_group(group_id_1)
                        end if
                        CALL h5_obj_exists(h5_id, "/init_params", &
                                           h5_exists_log)
                        if (.not. h5_exists_log) then
                            CALL h5_define_group(h5_id, &
                                                 "/init_params", group_id_2)
                            CALL h5_close_group(group_id_2)
                        end if
                    end do
                end do
            end do
        end do
        CALL h5_close(h5_id)
        CALL h5_deinit()

        ! reset loop variables, since they are also used in main code
        ifac_n = 1
        ifac_Te = 1
        ifac_Ti = 1
        ifac_vz = 1
        ! reset h5_mode_groupname string
        if (paramscan) then
            ! change parameter scan string used for the
            !group structure in hdf5 file
            write (parscan_str, "(A,F0.3,A,F0.3,A,F0.3,A,F0.3,A)") &
                "n", fac_n(ifac_n), "Te", fac_Te(ifac_Te), &
                "Ti", fac_Ti(ifac_Ti), "vz", fac_vz(ifac_vz), "/"
        else
            ! leave it empty if no parameter scan
            parscan_str = ""
        end if
    else
        ! if suppression_mode is true, only a simple group structure is created
        ! i.e. /f_m_n and /init_params
        ! if more than one RMP mode is used, use different group name
        if (numres .eq. 1) then
            write (h5_mode_groupname, "(A,I1,A,I1)") &
                "f_", m_vals(1), "_", n_vals(1)
        else
            write (h5_mode_groupname, "(A,I1,A,I1)") &
                "multi_mode"
        end if

        if (debug_mode) write (*,*) "Debug: h5_mode_groupname: ", trim(h5_mode_groupname)
        CALL h5_init()
        CALL h5_open_rw(path2out, h5_id)
        CALL h5_define_group(h5_id, trim(h5_mode_groupname), group_id_2)
        CALL h5_close_group(group_id_2)
        CALL h5_obj_exists(h5_id, "/init_params", &
                           h5_exists_log)
        if (.not. h5_exists_log) then
            CALL h5_define_group(h5_id, &
                                 "/init_params", group_id_2)
            CALL h5_close_group(group_id_2)
        end if

        CALL h5_close(h5_id)
        CALL h5_deinit()
    end if

    write (*, *) "finished creating group structure"
end subroutine


!> @brief subroutine write_br_time_data. Writes radial magnetic field perturbation evaluated at the resonant
!> surface, the antenna factor, the time and Dqle22 evaluated at the resonant surface for a given 
!> time step to the hdf5 file
!> @param[in] i Integer of time step to which the data will be saved. Goes from 1:i.
!> @param[in] br_abs_time Time value of the time evolution.
!> @param[in] br_abs_antenna_factor Value of the antenna factor, i.e. the RMP coil current.
!> @param[in] br_abs Absolute value of the radial magnetic field evaluated at the resonant surface in question.
!> @param[in] dqle22_res_time Value of Dqle22 evaluated at the resonant surface during the time evolution.
subroutine write_br_time_data(i, br_abs_time, br_abs_antenna_factor, br_abs, dqle22_res_time)

    use control_mod
    use baseparam_mod
    use h5mod
    use hdf5_tools

    implicit none
	integer, intent(in) :: i
    double precision, dimension(:), intent(in) :: br_abs_time
    double precision, dimension(:), intent(in) :: br_abs_antenna_factor
    double precision, dimension(:), intent(in) :: br_abs
    double precision, dimension(:), intent(in) :: dqle22_res_time
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


!> @brief subroutine write_init_profiles. Write initial profiles to hdf5 or ascii.
!> @author Markus Markl
!> @date 05.10.2022
subroutine write_init_profiles

    use grid_mod, only: params, qsaf
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


!> @brief subroutine rescale_profiles. Rescales kinetic profiles (n,Vz,Te,Ti).
!> @author Markus Markl
!> @date 05.10.2022
subroutine rescale_profiles

    use grid_mod, only: params
    use control_mod, only: debug_mode
    use paramscan_mod
    implicit none

    double precision, dimension(:), allocatable :: ErVzfac ! factor to rescale Er

    if (debug_mode) write(*,*) "Debug: coming into rescaling profiles"

    params(1, :) = hold_n * fac_n(ifac_n)
    params(2, :) = hold_vz * fac_vz(ifac_vz)
    params(3, :) = hold_Te * fac_Te(ifac_Te)
    params(4, :) = hold_Ti * fac_Ti(ifac_Ti)

    if (fac_vz(ifac_vz) .ne. 1.d0) then
        if (debug_mode) write(*,*) "Debug: fac_vz not equal 1. need to rescale Er as well"
        CALL h5_init()
        CALL h5_open_rw(path2out, h5_id)
        CALL h5_get_bounds_1(h5_id, '/factors/ErVzfac', lb, ub)
        allocate (ErVzfac(ub))
        CALL h5_get_double_1(h5_id, '/factors/ErVzfac', ErVzfac)
        CALL h5_close(h5_id)
        CALL h5_deinit()
        write (*, *) "rescale Er"
        idPhi0 = hold_dphi0 + ErVzfac*params(2, :)*(fac_vz(ifac_vz) - 1.d0)
        deallocate (ErVzfac)
    end if

    write(*,*) "Parameter scan, current factors: "
    write(*,*) "fac_n = ", fac_n(ifac_n), "   ", ifac_n, " of ", size(fac_n)
    write(*,*) "fac_Ti = ", fac_Ti(ifac_Ti), "   ", ifac_Ti, " of ", size(fac_Ti)
    write(*,*) "fac_Te = ", fac_Te(ifac_Te), "   ", ifac_Te, " of ", size(fac_Te)
    write(*,*) "fac_vz = ", fac_vz(ifac_vz), "   ", ifac_vz, " of ", size(fac_vz)


    if (debug_mode) write(*,*) "Debug: going out of rescaling profiles"
end subroutine !rescale_profiles

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
                    CALL writefort5000(i)
                end if
            end if
            if (suppression_mode .eqv. .false.) then
                CALL writefort1000(i)
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
                    CALL write_br_time_data(i, br_abs_time, br_abs_antenna_factor, br_abs, dqle22_res_time)
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

!> @brief Ramp up/down the RMP coil current.
!> @author Markus Markl
!> @date 13.03.2023
subroutine ramp_coil

    use paramscan_mod
    use control_mod
    implicit none
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
                        CALL write_br_time_data(i, br_abs_time, br_abs_antenna_factor, br_abs, dqle22_res_time)
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
                        CALL write_br_time_data(i, br_abs_time, br_abs_antenna_factor, br_abs, dqle22_res_time)
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
                        CALL write_br_time_data(i, br_abs_time, br_abs_antenna_factor, br_abs, dqle22_res_time)
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
                        CALL write_br_time_data(i, br_abs_time, br_abs_antenna_factor, br_abs, dqle22_res_time)
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
                        CALL write_br_time_data(i, br_abs_time, br_abs_antenna_factor, br_abs, dqle22_res_time)
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
                        CALL write_br_time_data(i, br_abs_time, br_abs_antenna_factor, br_abs, dqle22_res_time)
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
                CALL write_br_time_data(i, br_abs_time, br_abs_antenna_factor, br_abs, dqle22_res_time)
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

end program ql_balance
