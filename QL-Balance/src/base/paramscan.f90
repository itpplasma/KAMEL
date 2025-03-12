
module paramscan_mod

    use control_mod
    use balance_base, only: balance_t
    use QLBalance_kinds, only: dp

    implicit none

    integer :: ifac_n, ifac_Te, ifac_Ti, ifac_vz ! counter for the do loops
    integer :: numoffac                          ! total number of factors
    real(dp), dimension(:), allocatable :: fac_n, fac_Te, &
        fac_Ti, fac_vz
    character(len=1024) :: parscan_str
    real(dp) :: viscosity_factor
    real(dp), DIMENSION(:,:,:,:), ALLOCATABLE :: Er_res
    real(dp), DIMENSION(:, :, :, :), ALLOCATABLE :: br_abs_res_parscan
    real(dp), DIMENSION(:, :, :, :), ALLOCATABLE :: dqle22_res
    
    type, extends(balance_t) :: ParameterScan_t
        contains
            procedure :: init_balance => initParameterScan
            procedure :: run_balance => runParameterScan
    end type

    contains

    subroutine initParameterScan(this)

        use grid_mod, only: mwind, set_boundary_condition, npoib, rb
        use QLbalance_diag, only: write_diag, write_diag_b
        use QLBalance_hdf5_tools, only: h5overwrite
        use h5mod, only: mode_m, mode_n
        use control_mod, only: gyro_current_study, write_gyro_current, debug_mode, &
                        ihdf5IO
        use parallelTools, only: irank
        use wave_code_data, only: m_vals, n_vals
        use plasma_parameters, only: write_initial_parameters, alloc_hold_parameters, &
                                init_background_profiles

        implicit none

        class(ParameterScan_t), intent(inout) :: this
        this%runType = "ParameterScan"

        paramscan = .true.

        if (irank .eq. 0) then
            !iexit = 0 ! 0 - dont skip, 1 - skip, 2 - stop
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
            call set_boundary_condition
            CALL initialize_wave_code_interface(npoib, rb);
            CALL initialize_parameter_scan_vars

            mode_m = m_vals(1)
            mode_n = n_vals(1)
            if (debug_mode) write(*,*) 'Debug: mode_m = ', mode_m, 'mode_n = ', mode_n
            !call allocate_prev_variables

            if (ihdf5IO .eq. 1) then
                CALL create_group_structure_paramscan
            end if
            call init_background_profiles
            CALL write_initial_parameters

            !call alloc_hold_parameters

        end if

    end subroutine

    subroutine runParameterScan(this)

        use transp_coeffs_mod, only: rescale_transp_coeffs_by_ant_fac
        use parallelTools, only: irank
        use plasma_parameters, only: alloc_hold_parameters
            
        implicit none

        class(ParameterScan_t), intent(inout) :: this

        if (irank .eq. 0) then
            print *, ""
            print *, "  Running ParameterScan   "
            print *, ""
        end if

        call alloc_hold_parameters

        do ifac_n = 1, size(fac_n)
            do ifac_Te = 1, size(fac_Te)
                do ifac_Ti = 1, size(fac_Ti)
                    do ifac_vz = 1, size(fac_vz)
                        print *, "Before prepare h5 group_name"
                        call prepare_h5_group_name
                        print *, "Before rescale"
                        call rescale_profiles
                        print *, "Before calc geom"
                        call calc_geometric_parameter_profiles
                        print *, "Before initialize get dql"
                        call initialize_get_dql
                        print *, "Before get dql"
                        call get_dql
                        print *, "Before rescale transp"
                        call rescale_transp_coeffs_by_ant_fac
                        print *, "Before interpolate Br Dql"
                        call interpolate_Br_Dql_at_res_parscan
                    end do
                end do
            end do
        end do

        call finalize_parscan_run

    end subroutine

    subroutine prepare_h5_group_name

        use resonances_mod, only: add_mode_group_to_h5_mode_groupname
        use h5mod, only: h5_mode_groupname

        implicit none

        call add_mode_group_to_h5_mode_groupname
        call determine_h5_mode_groupname_of_scan
        print *, "Prepare_h5_group_name : h5_mode_groupname: ", trim(h5_mode_groupname)
        if (debug_mode) then
            print *, "Prepare_h5_group_name : h5_mode_groupname: ", trim(h5_mode_groupname)
        end if

    end subroutine

    subroutine determine_h5_mode_groupname_of_scan

        use h5mod

        implicit none

        if (.not. suppression_mode) then
            write (parscan_str, "(A,F0.3,A,F0.3,A,F0.3,A,F0.3,A)") "n", fac_n(ifac_n), &
                "Te", fac_Te(ifac_Te), "Ti", fac_Ti(ifac_Ti), "vz", fac_vz(ifac_vz) &
                , "/"
        else
            parscan_str = ""
        end if
        write (h5_mode_groupname, "(A,A,A)") trim(h5_mode_groupname), '/', trim(parscan_str)

    end subroutine

    subroutine finalize_parscan_run

        use parallelTools

        implicit none

        call write_Br_Dql_at_res_to_hdf5
        call deallocate_wave_code_data
        call MPI_finalize(ierror)
        stop '-> Finished parameter scan run |'

    end subroutine


    !> @brief subroutine initialize_parameter_scan_vars. Read factors for parameter scan and allocate variables. Still needed if no parameter scan is done.
    !> @author Markus Markl
    !> @date 05.10.2022
    subroutine initialize_parameter_scan_vars

        implicit none

        if (paramscan) then
            write(*,*) "Parameter scan: fetch factors for parameter scan"
            CALL getfactors
            write(*,*) "Got out of getfactors"
            if (size(fac_vz) .ne. 1) allocate(Er_res(size(fac_n), size(fac_Te), size(fac_Ti), size(fac_vz)))
            write(*,*) "allocated Er_res"
        else
            allocate(fac_n(1))
            allocate(fac_Ti(1))
            allocate(fac_Te(1))
            allocate(fac_vz(1))
            fac_n = (/1.d0/)
            fac_Ti = (/1.d0/)
            fac_Te = (/1.d0/)
            fac_vz = (/1.d0/)
        end if

        allocate(dqle22_res(size(fac_n), size(fac_Te), size(fac_Ti), size(fac_vz)))
        allocate(br_abs_res_parscan(size(fac_n), size(fac_Te), size(fac_Ti), size(fac_vz)))
        write(*,*) "Finished initialize parameter scan vars"

    end subroutine ! initialize_parameter_scan_vars


    subroutine getfactors

        use h5mod

        implicit none

        integer :: lb, ub

        CALL h5_init()
        CALL h5_open_rw(path2out, h5_id)
        CALL h5_get_bounds_1(h5_id, "/factors/fac_n", lb, ub)
        write (*, *) "lower bound fac_n", lb, " upper bound ", ub
        allocate(fac_n(ub))
        CALL h5_get_bounds_1(h5_id, "/factors/fac_Te", lb, ub)
        write (*, *) "lower bound fac_Te", lb, " upper bound ", ub
        allocate(fac_Te(ub))
        CALL h5_get_bounds_1(h5_id, "/factors/fac_Ti", lb, ub)
        write (*, *) "lower bound fac_Ti", lb, " upper bound ", ub
        allocate(fac_Ti(ub))
        CALL h5_get_bounds_1(h5_id, "/factors/fac_vz", lb, ub)
        write (*, *) "lower bound fac_vz", lb, " upper bound ", ub
        allocate(fac_vz(ub))

        CALL h5_get_double_1(h5_id, "/factors/fac_n", fac_n)
        CALL h5_get_double_1(h5_id, "/factors/fac_Te", fac_Te)
        CALL h5_get_double_1(h5_id, "/factors/fac_Ti", fac_Ti)
        CALL h5_get_double_1(h5_id, "/factors/fac_vz", fac_vz)

        CALL h5_close(h5_id)
        CALL h5_deinit()

    end subroutine


    !> @brief subroutine rescale_profiles. Rescales kinetic profiles (n,Vz,Te,Ti).
    !> @author Markus Markl
    !> @date 05.10.2022
    subroutine rescale_profiles

        use plasma_parameters, only: params, hold_n, hold_vz, hold_Te, hold_Ti, hold_dphi0
        use control_mod, only: debug_mode
        use h5mod
        use wave_code_data, only: idPhi0

        implicit none

        double precision, dimension(:), allocatable :: ErVzfac ! factor to rescale Er
        integer :: lowerBound, upperBound

        if (debug_mode) write(*,*) "Debug: coming into rescaling profiles"

        print *, "hold times fac"
        params(1, :) = hold_n * fac_n(ifac_n)
        params(2, :) = hold_vz * fac_vz(ifac_vz)
        params(3, :) = hold_Te * fac_Te(ifac_Te)
        params(4, :) = hold_Ti * fac_Ti(ifac_Ti)
        print *, "After hold times fac"

        if (fac_vz(ifac_vz) .ne. 1.d0) then
            if (debug_mode) write(*,*) "Debug: fac_vz not equal 1. need to rescale Er as well"
            CALL h5_init()
            CALL h5_open_rw(path2out, h5_id)
            CALL h5_get_bounds_1(h5_id, '/factors/ErVzfac', lowerBound, upperBound)
            allocate(ErVzfac(upperBound))
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


    subroutine interpolate_Br_Dql_at_res_parscan

        use PolyLagrangeInterpolation
        use grid_mod, only: npoib, rb, r_resonant, dqle22, npoic, Ercovavg, Ercov
        use wave_code_data, only: Br, antenna_factor

        implicit none

        integer :: indResRadius, ind_begin_interp, ind_end_interp
        integer :: ipoi

        if (.not. allocated(coef)) allocate (coef(0:nder, nlagr))
        call binsrc(rb, 1, npoib, r_resonant(1), indResRadius)
        call get_ind_Lagr_interp(indResRadius, ind_begin_interp, ind_end_interp)
        call plag_coeff(nlagr, nder, r_resonant(1), rb(ind_begin_interp:ind_end_interp), coef)
        
        dqle22_res(ifac_n, ifac_Te, ifac_Ti, ifac_vz) = sum(coef(0, :) * dqle22(ind_begin_interp:ind_end_interp))
        br_abs_res_parscan(ifac_n, ifac_Te, ifac_Ti, ifac_vz) = sum(coef(0, :) &
            * abs(Br(ind_begin_interp:ind_end_interp)))*sqrt(antenna_factor)

        write(*,*) ""
        write(*,*) "Dqle22 res = ", dqle22_res(ifac_n, ifac_Te, ifac_Ti, ifac_vz)
        write(*,*) "|Br| res   = ", br_abs_res_parscan(ifac_n, ifac_Te, ifac_Ti, ifac_vz)
        write(*,*) "Antenna factor = ", antenna_factor
        write(*,*) ""

        if (size(fac_vz) .ne. 1) then
            do ipoi = 1, npoic
                Ercovavg(ipoi) = 0.5d0*(Ercov(ipoi) + Ercov(ipoi + 1))
            end do
            Er_res(ifac_n, ifac_Te, ifac_Ti, ifac_vz) = sum(coef(0, :)*Ercovavg(ind_begin_interp:ind_end_interp))
            if (debug_mode) write (*, *) "Er_res = ", Er_res(ifac_n, ifac_Te, ifac_Ti, ifac_vz)
        end if

    end subroutine 

    !> @brief subroutine create_group_structure_paramscan. Creates the group structure in the hdf5 file.
    !> @author  Markus Markl
    !> @date 12.03.2021
    subroutine create_group_structure_paramscan

        use control_mod
        use wave_code_data, only: m_vals, n_vals
        use resonances_mod, only: numres
        use h5mod

        implicit none
        !logical :: suppression_mode = .true.

        if (debug_mode) write (*, *) "Debug: Creating group structure"
        ! if the profiles should be written out, i.e. suppression_mode = false, then an extended group
        ! structure is created to save them
        CALL h5_init()
        CALL h5_open_rw(path2out, h5_id)
        CALL h5_obj_exists(h5_id, "/init_params", h5_exists_log)
        if (.not. h5_exists_log) then
            CALL h5_define_group(h5_id, "/init_params", group_id_2)
            CALL h5_close_group(group_id_2)
        end if

        if (.not. suppression_mode) then
            do ifac_n = 1, size(fac_n)
                do ifac_Te = 1, size(fac_Te)
                    do ifac_Ti = 1, size(fac_Ti)
                        do ifac_vz = 1, size(fac_vz)
                            write (*, *) ifac_n, ifac_Te, ifac_Ti, ifac_vz
                            write (parscan_str, "(A,F0.3,A,F0.3,A,F0.3,A,F0.3)") &
                                "n", fac_n(ifac_n), "Te", fac_Te(ifac_Te), &
                                "Ti", fac_Ti(ifac_Ti), "vz", fac_vz(ifac_vz)
                            
                            if (numres .eq. 1) then
                                if (m_vals(1) < 10) then
                                    write (h5_mode_groupname, "(A,I1,A,I1,A,A)") &
                                        "f_", m_vals(1), "_", n_vals(1), '/', &
                                        trim(parscan_str)
                                else
                                    write (h5_mode_groupname, "(A,I2,A,I1,A,A)") &
                                        "f_", m_vals(1), "_", n_vals(1), '/', trim(parscan_str)
                                end if
                            else
                                write (h5_mode_groupname, "(A,A,A,I1,A,I1)") &
                                    "multi_mode/", trim(parscan_str)
                            end if

                            ! create the groups that are furthest down: fort.1000,
                            ! fort.5000 and init_params
                            if (debug_mode) write(*,*) "Debug: h5_mode_groupname ", trim(h5_mode_groupname)
                            CALL h5_create_parent_groups(h5_id, trim(h5_mode_groupname) //'/')
                            CALL h5_create_parent_groups(h5_id, trim(h5_mode_groupname)//"/KinProfiles/")
                            CALL h5_create_parent_groups(h5_id, trim(h5_mode_groupname)//"/LinearProfiles/")
                        end do
                    end do
                end do
            end do

            ! reset loop variables, since they are also used in main code
            ifac_n = 1
            ifac_Te = 1
            ifac_Ti = 1
            ifac_vz = 1
            write (parscan_str, "(A,F0.3,A,F0.3,A,F0.3,A,F0.3,A)") &
                    "n", fac_n(ifac_n), "Te", fac_Te(ifac_Te), &
                    "Ti", fac_Ti(ifac_Ti), "vz", fac_vz(ifac_vz), "/"
        else
            if (numres .eq. 1) then
                if (m_vals(1) < 10) then
                    write (h5_mode_groupname, "(A,I1,A,I1)") "f_", m_vals(1), "_", n_vals(1)
                else
                    write (h5_mode_groupname, "(A,I2,A,I1)") "f_", m_vals(1), "_", n_vals(1)
                end if
            else
                write (h5_mode_groupname, "(A)") "multi_mode"
            end if

            if (debug_mode) write (*,*) "Debug: h5_mode_groupname: ", trim(h5_mode_groupname)
            CALL h5_define_group(h5_id, trim(h5_mode_groupname), group_id_2)
            CALL h5_close_group(group_id_2)

        end if
        CALL h5_close(h5_id)
        CALL h5_deinit()

        if (debug_mode) write (*, *) "Debug: finished creating group structure"
    end subroutine

    subroutine write_Br_Dql_at_res_to_hdf5

        !use paramscan_mod, only: dqle22_res, br_abs_res_parscan, Er_res, &
        !    fac_vz
        use wave_code_data, only: m_vals, n_vals
        use resonances_mod, only: numres
        use h5mod

        implicit none

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

    end subroutine


    subroutine writeDqle22

        use grid_mod, only: dqle22
        use h5mod

        implicit none

        !print *, "In write dqle22"
        
        CALL h5_init()
        CALL h5_open_rw(path2out, h5_id)
        CALL h5_obj_exists(h5_id, trim(h5_mode_groupname), &
            h5_exists_log)
        if (.not. h5_exists_log) then
            CALL h5_define_group(h5_id, trim(h5_mode_groupname), group_id_2)
            CALL h5_close_group(group_id_2)
        end if

        CALL h5_add_double_1(h5_id, trim(h5_mode_groupname)//'/dqle22_res', &
                            reshape(dqle22_res, (/size(dqle22_res)/)), &
                            lbound(reshape(dqle22_res, (/size(dqle22_res)/))), &
                            ubound(reshape(dqle22_res, (/size(dqle22_res)/))))
        CALL h5_add_double_1(h5_id, trim(h5_mode_groupname)//'/dqle22', &
                                dqle22, lbound(dqle22), ubound(dqle22))
 
        CALL h5_close(h5_id)
        CALL h5_deinit()

    end subroutine

end module