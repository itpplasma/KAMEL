module singleStep

    use balance_base, only: balance_t
    use QLBalance_kinds, only: dp

    implicit none

    real(dp) :: dqle22_res_single, br_abs_res_single

    type, extends(balance_t) :: singleStep_t
        contains
            procedure :: init_balance => initSingleStep
            procedure :: run_balance => runSingleStep
    end type

    contains
    
    subroutine initSingleStep(this)

        use grid_mod, only: mwind, set_boundary_condition, npoib, rb
        use KAMEL_hdf5_tools, only: h5overwrite
        use h5mod, only: mode_m, mode_n
        use control_mod, only: gyro_current_study, write_gyro_current, debug_mode, &
                        ihdf5IO
        use parallelTools, only: irank
        use wave_code_data, only: m_vals, n_vals
        use plasma_parameters, only: write_initial_parameters, alloc_hold_parameters, &
                                init_background_profiles
        use QLbalance_diag, only: write_diag, write_diag_b

        implicit none

        class(SingleStep_t), intent(inout) :: this
        this%runType = "SingleStep"
        

        if (irank .eq. 0) then
            print *, "Initialize Single Step run"
            mwind = 10
            write_diag = .false.
            write_diag_b = .false.
    
            if (gyro_current_study .ne. 0) then
                write_gyro_current = .true.
            else
                write_gyro_current = .false.
            end if

            call gengrid
            call set_boundary_condition
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
            call calc_geometric_parameter_profiles
            !call alloc_hold_parameters
        end if

    end subroutine

    subroutine runSingleStep(this)

        use transp_coeffs_mod, only: rescale_transp_coeffs_by_ant_fac
        use parallelTools, only: irank
        use control_mod, only: debug_mode

        implicit none

        class(SingleStep_t), intent(inout) :: this

        if (irank .eq. 0) then
            print *, "Running SingleStep"
        end if

        call initialize_get_dql

        if (debug_mode) write(*,*) "Debug: before get_dql"
        call get_dql
        if (debug_mode) write(*,*) "Debug: after get_dql"
        call rescale_transp_coeffs_by_ant_fac

        call interp_Br_Dql_at_res
        call finalize_singlestep_run

    end subroutine

    subroutine interp_Br_Dql_at_res

        use PolyLagrangeInterpolation, only: coef, plag_coeff, binsrc, get_ind_Lagr_interp, nder, nlagr
        use wave_code_data, only: Br, antenna_factor
        use grid_mod, only: npoib, rb, r_resonant, dqle22
        
        implicit none

        integer :: indResRadius, ind_begin_interp, ind_end_interp

        ! TODO: Correct for multi mode
        if (.not. allocated(coef)) allocate (coef(0:nder, nlagr))
        call binsrc(rb, 1, npoib, r_resonant(1), indResRadius)
        call get_ind_Lagr_interp(indResRadius, ind_begin_interp, ind_end_interp)
        call plag_coeff(nlagr, nder, r_resonant(1), rb(ind_begin_interp:ind_end_interp), coef)
        
        dqle22_res_single = sum(coef(0, :) * dqle22(ind_begin_interp:ind_end_interp))
        br_abs_res_single = sum(coef(0, :) * abs(Br(ind_begin_interp:ind_end_interp)))*sqrt(antenna_factor)

        write(*,*) " "
        write(*,*) "=== === === Results of Single Step: === === ==="
        write(*,*) "    Dqle22 res  = ", dqle22_res_single
        write(*,*) "    |Br| res    = ", br_abs_res_single
        write(*,*) " "
        write(*,*) "    With the Antenna factor = ", antenna_factor
        write(*,*) "=== === === === === === === === === === === ==="
        write(*,*) " "

    end subroutine


    subroutine finalize_singlestep_run

        use h5mod
        use mpi
        use parallelTools, only: ierror
        use control_mod, only: ihdf5IO

        implicit none

        if (ihdf5IO .eq. 1) then
            call write_Dqle22_SingleStep
        end if
        call MPI_finalize(ierror);
        stop '-> Finished linear run |'

    end subroutine

    subroutine write_Dqle22_SingleStep

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