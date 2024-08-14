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
    use singleStep
    use parallelTools
    use restart_mod
    use PolyLagrangeInterpolation
    use plasma_parameters
    use balance_mod

    implicit none

    integer :: ipoi, i, ieq, l, k
    integer :: ioddeven
    character(100) :: typeOfRun = "ParameterScan"
    
    class(balance_t), allocatable :: balanceInstance

    call fromBalanceFactoryGetBalance(typeOfRun, balanceInstance)
    
    !call balanceInstance%runBalance()

    !call initialize_balance_code
    call balanceInstance%initBalance()
    call balanceInstance%runBalance()

    call MPI_finalize(ierror)
    stop "to test SingleStep run"

    ! parameter scan loops that span over (nearly) the rest of the code
    do ifac_n = 1, size(fac_n)
        do ifac_Te = 1, size(fac_Te)
            do ifac_Ti = 1, size(fac_Ti)
                do ifac_vz = 1, size(fac_vz)
                    !if (paramscan) then
                        ! change parameter scan string used for navigating the hdf5 file, only if
                        ! the suppression_mode is not activated
                    !    if (.not. suppression_mode) then
                    !        write (parscan_str, "(A,F0.3,A,F0.3,A,F0.3,A,F0.3,A)") "n", fac_n(ifac_n), &
                    !            "Te", fac_Te(ifac_Te), "Ti", fac_Ti(ifac_Ti), "vz", fac_vz(ifac_vz) &
                    !            , "/"
                    !    else
                    !        parscan_str = ""
                    !    end if
                    !else
                    !    write(*,*) " "
                    !    write(*,*) "- No parameter scan -"
                    !    write(*,*) " "
                        ! leave it empty if no parameter scan
                    !    parscan_str = ""
                    !end if

                    !write (h5_mode_groupname, "(A)") trim(parscan_str)
                    !write(*,*) "h5_mode_groupname: ", trim(h5_mode_groupname)

				
                    ! if more than one RMP mode is used, use different group name
                    !if (numres .eq. 1) then
                    !    write (h5_mode_groupname, "(A,A,I1,A,I1)") trim(h5_mode_groupname), &
                    !        "f_", m_vals(1), "_", n_vals(1)
                    !else
                    !    write (h5_mode_groupname, "(A,A,I1,A,I1)") trim(h5_mode_groupname), &
                    !        "multi_mode"
                    !end if

                    !write(*,*) "h5_mode_groupname after f_m_n: ", trim(h5_mode_groupname)

                    call rescale_profiles
                    call calc_geometric_parameter_profiles

                    irf = 2 ! initialize dql variables and set to zero !
                    call get_dql

                    if (flag_run_time_evolution) then
                        call initialize_antenna_factor
                    end if

                    call rescale_transp_coeffs_by_ant_fac
                    
                    irf = 1

                    call det_balance_eqs_source_terms ! calcs source term in balance equations, needed for time evolution

					write(*,*) "h5_mode_groupname before writeKinProfileDataToDisk: ", trim(h5_mode_groupname)
                    if (irank .eq. 0) then
                        if (suppression_mode .eqv. .false.) then
                            CALL writeKinProfileDataToDisk(0) ! write the profiles to hdf5 file
                        end if
                    end if

                    time = 0.d0
                    tol = tol_max

                    call inquiry_to_restart

                    timstep_arr = timstep
                    tim_stack = timstep_arr

                    iunit_diag = 5000

                    call get_dql ! also writes out diffusion coefficients and other data
                    call rescale_transp_coeffs_by_ant_fac

                    if (flag_run_time_evolution) then
                        if (ifac_n + ifac_Ti + ifac_Te + ifac_vz .eq. 4) then
                            call alloc_Br_Dqle_for_timeevol
                        end if
                    end if

                    call hold_prev_transp_coeffs

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
                        call interpolate_Br_Dql_at_res_parscan

                        if (paramscan) then
                            ! if the last parameter scan is done, write data and stop the code
                            if (ifac_n + ifac_Te + ifac_Ti + ifac_vz .eq. size(fac_n) + size(fac_Ti) + &
                                size(fac_Te) + size(fac_vz)) then
                                if (debug_mode) write(*,*) "Debug: Last parameter done. Finalize MPI"

                                call write_Br_Dql_at_res_to_hdf5
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
                            call finalizeSingleStepRun
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
                        timeIndex = i
                        write (*, *) "Time Index = ", timeIndex
                        
                        call copy_kin_profs_to_yprev

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
						call write_br_dqle22_time_data
                        call rescale_transp_coeffs_by_ant_fac

                        if (write_diag_b) close (iunit_diag_b)

                        call writefort9999
                        
                        if (.not. redostep) then
                            call hold_prev_transp_coeffs
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
                        
                        call rescaleTimStepArr
                        call setTimStep
                        call reset_timstep_arr_w_timstep
                        call writeTimeInfoToDisk
                        call relaxPlasmaParameters

                        timstep_arr = 0.d0
                        call evolvestep(timstep, eps)
                        timstep_arr = timstep
                        time = time + timstep

                        call messageTimeInfo
                        call writeKinProfileAtTimeIndex
                        call setFirstIterationTrue
                        CALL checkIfLinearDiscrepancyOfPenRatioReached
 
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
                    !deallocate (coef)

                end do
            end do
        end do
    end do

    write (*, *) 'Programm is finalized without stopping criterion met';
    call MPI_finalize(ierror)

end program ql_balance
