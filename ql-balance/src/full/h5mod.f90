
module h5mod

    use hdf5_tools

    implicit none

    integer(HID_T) :: h5_id, group_id_1, group_id_2, group_id_3, dataset_id
    logical :: h5_exists_log
    integer :: mode_n, mode_m
    character(len=1024) :: path2inp ! path to hdf5 file with input data
    character(len=1024) :: path2time ! path to hdf5 file from which time evolved profiles are read
    character(len=1024) :: path2out ! path to hdf5 file where output is written
    character(len=1024) :: h5_mode_groupname

    contains
    !> @brief subroutine creategroupstructure. Creates the group structure in the hdf5 file.
    !> @author  Markus Markl
    !> @date 12.03.2021
    subroutine creategroupstructure

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



end module h5mod