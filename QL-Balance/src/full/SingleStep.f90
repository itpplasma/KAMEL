module singleStep

    use balanceBase, only: balance_t

    double precision :: dqle22_res_single, br_abs_res_single

    type, extends(balance_t) :: singleStep_t
        contains
            procedure :: init_balance => initSingleStep
            procedure :: run_balance => runSingleStep
    end type

    contains
    
    subroutine initSingleStep(this)

        use grid_mod, only: mwind, rmax, rmin, setBoundaryCondition, npoib, rb
        use baseparam_mod, only: dperp
        use diag_mod, only: write_diag, write_diag_b
        use hdf5_tools, only: h5overwrite
        use h5mod, only: mode_m, mode_n
        use control_mod, only: gyro_current_study, write_gyro_current, debug_mode, &
                          ihdf5IO
        use parallelTools, only: irank
        use wave_code_data, only: m_vals, n_vals
        use plasma_parameters, only: write_initial_parameters, alloc_hold_parameters, &
                                init_background_profiles

        implicit none

        class(SingleStep_t), intent(inout) :: this
        this%runType = "SingleStep"
        

        if (irank .eq. 0) then
            print *, "Initialize Single Step run"
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

            call gengrid
            call setBoundaryCondition
            CALL initialize_wave_code_interface(npoib, rb);

            mode_m = m_vals(1)
            mode_n = n_vals(1)
            if (debug_mode) write(*,*) 'Debug: mode_m = ', mode_m, 'mode_n = ', mode_n
            !call allocate_prev_variables
            if (ihdf5IO .eq. 1) then
                call create_group_structure_singlestep
            end if
            call init_background_profiles
            CALL write_initial_parameters
            !call alloc_hold_parameters
        end if

    end subroutine

    subroutine runSingleStep(this)

        use control_mod, only: irf
        use transp_coeffs_mod, only: rescale_transp_coeffs_by_ant_fac
        use parallelTools, only: irank

        implicit none

        class(SingleStep_t), intent(inout) :: this

        if (irank .eq. 0) then
            print *, "Running SingleStep"
        end if

        call initialize_get_dql
        call get_dql
        call rescale_transp_coeffs_by_ant_fac

        call interpBrAndDqlAtResonance
        call finalizeSingleStepRun

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

        write(*,*) " "
        write(*,*) "=== === === Results of Single Step: === === ==="
        write(*,*) "    Dqle22 res  = ", dqle22_res_single
        write(*,*) "    |Br| res    = ", br_abs_res_single
        write(*,*) " "
		write(*,*) "    With the Antenna factor = ", antenna_factor
        write(*,*) "=== === === === === === === === === === === ==="
        write(*,*) " "

    end subroutine


    subroutine finalizeSingleStepRun

        use h5mod
        use mpi
        use parallelTools, only: ierror
        use control_mod, only: ihdf5IO

        implicit none

        if (ihdf5IO .eq. 1) then
            call writeDqle22_SingleStep
        end if
        call MPI_finalize(ierror);
        stop '-> Finished linear run |'

    end subroutine

    subroutine writeDqle22_SingleStep

        use grid_mod, only: dqle22, rb, r_resonant
        use h5mod

        implicit none

        !print *, "In write dqle22"
        
        CALL h5_init()
        CALL h5_open_rw(path2out, h5_id)
        !CALL h5_obj_exists(h5_id, trim(h5_mode_groupname), &
        !    h5_exists_log)
        !if (.not. h5_exists_log) then
        !    CALL h5_define_group(h5_id, trim(h5_mode_groupname), group_id_2)
        !    CALL h5_close_group(group_id_2)
        !end if

        call create_group_if_not_existent(h5_mode_groupname)

        CALL h5_add_double_0(h5_id, trim(h5_mode_groupname)//'/dqle22_res', dqle22_res_single)
        CALL h5_add_double_0(h5_id, trim(h5_mode_groupname)//'/Br_abs_res', br_abs_res_single)
        CALL h5_add_double_1(h5_id, trim(h5_mode_groupname)//'/dqle22', &
                                dqle22, lbound(dqle22), ubound(dqle22))
        CALL h5_add_double_1(h5_id, trim(h5_mode_groupname)//'/r_eff', &
                                rb, lbound(rb), ubound(rb))
        CALL h5_add_double_0(h5_id, trim(h5_mode_groupname)//'/r_res', r_resonant(1))
 
        CALL h5_close(h5_id)
        CALL h5_deinit()

    end subroutine

    subroutine create_group_structure_singlestep

        use control_mod
        use wave_code_data, only: m_vals, n_vals
        use resonances_mod, only: numres
        use h5mod

        implicit none

        if (debug_mode) then
            print *, "Creating group structure for Single Step"
        end if

        if (numres .eq. 1) then
            if (m_vals(1) <= 9) then
                write (h5_mode_groupname, "(A,I1,A,I1)") "f_", m_vals(1), "_", n_vals(1)
            else
                write (h5_mode_groupname, "(A,I2,A,I1)") "f_", m_vals(1), "_", n_vals(1)
            end if
        else
            write (h5_mode_groupname, "(A)") "multi_mode"
        end if

        CALL h5_init()
        CALL h5_open_rw(path2out, h5_id)

        if (.not. suppression_mode) then
            if (debug_mode) then
                write(*,*) "h5_mode_groupname ", trim(h5_mode_groupname)
            end if
            CALL h5_create_parent_groups(h5_id, trim(h5_mode_groupname)//'/')
            CALL h5_create_parent_groups(h5_id, trim(h5_mode_groupname)//"/KinProfiles/")
            CALL h5_define_group(h5_id, trim(h5_mode_groupname)//"/LinearProfiles/", group_id_1)
            CALL h5_close_group(group_id_1)
            CALL h5_obj_exists(h5_id, "/init_params", h5_exists_log)
            if (.not. h5_exists_log) then
                CALL h5_define_group(h5_id, "/init_params", group_id_2)
                CALL h5_close_group(group_id_2)
            end if
        else
            if (debug_mode) then 
                write (*,*) "Debug: h5_mode_groupname: ", trim(h5_mode_groupname)
            end if
            CALL h5_define_group(h5_id, trim(h5_mode_groupname), group_id_2)
            CALL h5_close_group(group_id_2)
            CALL h5_obj_exists(h5_id, "/init_params", h5_exists_log)
            if (.not. h5_exists_log) then
                CALL h5_define_group(h5_id, "/init_params", group_id_2)
                CALL h5_close_group(group_id_2)
            end if
        end if

        CALL h5_close(h5_id)
        CALL h5_deinit()

        if (debug_mode) write (*, *) "Debug: finished creating group structure for Single Step"
    end subroutine



end module