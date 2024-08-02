
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

    subroutine writeBrAndDqlAtResonanceToH5

        use paramscan_mod, only: dqle22_res, br_abs_res_parscan, Er_res, &
            fac_vz
        use wave_code_data, only: m_vals, n_vals

        implicit none
        
        write (h5_mode_groupname, "(A,I1,A,I1)") "f_", m_vals(1), "_", n_vals(1)

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
        use paramscan_mod, only: dqle22_res

        implicit none
        
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

    subroutine writeReasonForStopToH5(reason)

        implicit none

        character(*), intent(in) :: reason

        CALL h5_init()
        CALL h5_open_rw(path2out, h5_id)
        CALL h5_add_string(h5_id, trim(h5_mode_groupname)// &
        '/stopping_criterion', reason)
        CALL h5_close(h5_id)
        CALL h5_deinit()

    end subroutine

    ! Added by Markus Markl, 18.03.2021
    ! Used to save the fort5000 data.
    subroutine writefort5000
        use grid_mod, only: nbaleqs, neqset, iboutype, npoic, npoib &
                            , params, ddr_params, deriv_coef &
                            , ipbeg, ipend, rb, params_b, reint_coef &
                            , rc, sqg_bthet_overc, Ercov &
                            , ddr_params_nl, y, mwind &
                            , dqle11, dqle12, dqle21, dqle22 &
                            , dqli11, dqli12, dqli21, dqli22, d11_misalign, Es_pert_flux

        use baseparam_mod, only: Z_i, e_charge, am, p_mass, c, btor, e_mass, ev, rtor
        use control_mod, only: ihdf5IO, diagnostics_output, misalign_diffusion
        use wave_code_data
        use diag_mod, only: write_diag, iunit_diag, write_diag_b, iunit_diag_b, i_mn_loop
        use time_evolution, only: timeStep
    
        implicit none

        integer :: ipoi, ieq, i, npoi!,i_mn,ierr,mwind_save
        character(len=1024) :: tempch

        write(*,*) "writing fort.5000"

        if (ihdf5IO .eq. 1) then

            CALL h5_init()
            CALL h5_open_rw(path2out, h5_id)
            tempch = "/"//trim(h5_mode_groupname)//"/fort.5000"

            write (tempch, "(A,A,I4,A)") trim(tempch), "/", timeStep, "/"

            CALL h5_define_group(h5_id, trim(tempch), group_id_1)
            CALL h5_close_group(group_id_1)

            CALL h5_add_double_1(h5_id, trim(tempch)//"r", &
                                r, lbound(r), ubound(r))
            CALL h5_add_double_1(h5_id, trim(tempch)//"Br_abs", &
                                abs(Br), lbound(Br), ubound(Br))
            CALL h5_add_double_1(h5_id, trim(tempch)//"Re_Br", &
                                real(Br), lbound(Br), ubound(Br))
            CALL h5_add_double_1(h5_id, trim(tempch)//"Im_Br", &
                                dimag(Br), lbound(Br), ubound(Br))
 
            CALL h5_add_double_1(h5_id, trim(tempch)//"Jpe_abs", &
                                abs(Jpe), lbound(Jpe), ubound(Jpe))
            CALL h5_add_double_1(h5_id, trim(tempch)//"Jpi_abs", &
                                abs(Jpi), lbound(Jpi), ubound(Jpi))
            CALL h5_add_double_1(h5_id, trim(tempch)//"dqle22", &
                                    dqle22, lbound(dqle22), ubound(dqle22))

            if (misalign_diffusion .eqv. .true.) then
                write(*,*) " "
                write(*,*) "Writing misalignment diffusion to hdf5"
                !CALL h5_init()
                !CALL h5_open_rw(path2out, h5_id)
                CALL h5_add_double_1(h5_id, trim(tempch)//"D11_MA", d11_misalign, lbound(d11_misalign), ubound(d11_misalign))
                CALL h5_add_double_1(h5_id, trim(tempch)//"Es_pert_flux_real", real(Es_pert_flux), & 
                    lbound(real(Es_pert_flux)), ubound(dreal(Es_pert_flux)))
                CALL h5_add_double_1(h5_id, trim(tempch)//"Es_pert_flux_imag", dimag(Es_pert_flux), &
                    lbound(dimag(Es_pert_flux)), ubound(dimag(Es_pert_flux)))

                CALL h5_add_double_1(h5_id, trim(tempch)//"Br_real", &
                    real(Br), lbound(real(Br)), ubound(real(Br)))
                CALL h5_add_double_1(h5_id, trim(tempch)//"Es_real", &
                    real(Es), lbound(real(Es)), ubound(real(Es)))

                CALL h5_add_double_1(h5_id, trim(tempch)//"Br_imag", &
                    dimag(Br), lbound(dimag(Br)), ubound(dimag(Br)))
                CALL h5_add_double_1(h5_id, trim(tempch)//"Es_imag", &
                    dimag(Es), lbound(dimag(Es)), ubound(dimag(Es)))
                !CALL h5_close(h5_id)
                !CALL h5_deinit()
                write(*,*) "Carry on"
                write(*,*) " "
            end if ! hdf5test .eq. 1
 


            ! write the whole content only if diagnostics_output is true
            if (diagnostics_output) then
                !CALL h5_add_double_1(h5_id, trim(tempch)//"dqle22", &
                !                     dqle22, lbound(dqle22), ubound(dqle22))
                CALL h5_add_double_1(h5_id, trim(tempch)//"dqle11", &
                                    dqle11, lbound(dqle11), ubound(dqle11))
                CALL h5_add_double_1(h5_id, trim(tempch)//"dqle12", &
                                    dqle12, lbound(dqle12), ubound(dqle12))
                CALL h5_add_double_1(h5_id, trim(tempch)//"dqli11", &
                                    dqli11, lbound(dqli11), ubound(dqli11))
                CALL h5_add_double_1(h5_id, trim(tempch)//"dqli12", &
                                    dqli12, lbound(dqli12), ubound(dqli12))
                CALL h5_add_double_1(h5_id, trim(tempch)//"dqli22", &
                                    dqli22, lbound(dqli22), ubound(dqli22))
                CALL h5_add_double_1(h5_id, trim(tempch)//"Br-ckpEs_om_E", &
                                    abs(Br - c*kp*Es/om_E), lbound(Br), ubound(Br))
                CALL h5_add_double_1(h5_id, trim(tempch)//"Br-cksEp_om_E", &
                                    abs(Br - c*ks*Ep/om_E), lbound(Br), ubound(Br))
                CALL h5_add_double_1(h5_id, trim(tempch)//"Jpe_abs", &
                                    abs(Jpe), lbound(Jpe), ubound(Jpe))
                CALL h5_add_double_1(h5_id, trim(tempch)//"Jpi_abs", &
                                    abs(Jpi), lbound(Jpi), ubound(Jpi))
                CALL h5_add_double_1(h5_id, trim(tempch)//"JpeJpi_abs", &
                                    abs(Jpe + Jpi), lbound(Jpe), ubound(Jpe))
            end if
            CALL h5_close(h5_id)
            CALL h5_deinit()
        else
            do ipoi = 1, npoib
    !      write(iunit_diag,*) r(ipoi),dqle11(ipoi),dqli11(ipoi) &
    !                            ,real(Br(ipoi)),imag(Br(ipoi))  &
    !                            ,real(Ep(ipoi)),imag(Ep(ipoi))  &
    !                            ,real(Er(ipoi)),imag(Er(ipoi))  &
    !                            ,om_E(ipoi) &
    !                            ,abs(c*ks(ipoi)*Ep(ipoi) - om_E(ipoi)*Br(ipoi))**2
                write (iunit_diag, *) r(ipoi), dqle11(ipoi), dqle12(ipoi) &
                    , dqle22(ipoi), dqli11(ipoi) &
                    , dqli12(ipoi), dqli22(ipoi) &
                    , abs(Br(ipoi)) &
                    , abs(Br(ipoi) - c*kp(ipoi)*Es(ipoi)/om_E(ipoi)) &
                    , abs(Br(ipoi) - c*ks(ipoi)*Ep(ipoi)/om_E(ipoi)) &
                    , abs(Jpe(ipoi)), abs(Jpi(ipoi)) &
                    , abs(Jpe(ipoi) + Jpi(ipoi))
            end do
            close (iunit_diag)
        end if
        print *, "Finished writing fort.5000 data"
    end subroutine writefort5000

    


end module h5mod