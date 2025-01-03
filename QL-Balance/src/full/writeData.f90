subroutine write_fields_currs_transp_coefs_to_h5
    use grid_mod, only: nbaleqs, neqset, iboutype, npoic, npoib &
                        , deriv_coef &
                        , ipbeg, ipend, rb, reint_coef &
                        , rc, sqg_bthet_overc, Ercov &
                        , y, mwind &
                        , dqle11, dqle12, dqle21, dqle22 &
                        , dqli11, dqli12, dqli21, dqli22, d11_misalign, Es_pert_flux
    use plasma_parameters, only: params, ddr_params, params_b, ddr_params_nl
    use baseparam_mod, only: Z_i, e_charge, am, p_mass, c, btor, e_mass, ev, rtor
    use control_mod, only: ihdf5IO, diagnostics_output, misalign_diffusion
    use wave_code_data
    use diag_mod, only: write_diag, iunit_diag, write_diag_b, iunit_diag_b, i_mn_loop
    use time_evolution, only: timeIndex
    use h5mod
    
    implicit none

    integer :: ipoi, ieq, i, npoi!,i_mn,ierr,mwind_save
    character(len=1024) :: tempch

    if (ihdf5IO .eq. 1) then

        CALL h5_init()
        CALL h5_open_rw(path2out, h5_id)
        if (debug_mode) print *, "Debug: ", trim(h5_mode_groupname)
        tempch = "/"//trim(h5_mode_groupname)//"/LinearProfiles"

        write (tempch, "(A,A,I4,A)") trim(tempch), "/", timeIndex, "/"

        call create_group_if_not_existent(tempch)

        CALL h5_add_double_1(h5_id, trim(tempch)//"r", &
                            r, lbound(r), ubound(r))
        CALL h5_add_double_1(h5_id, trim(tempch)//"Br_abs", &
                            abs(Br), lbound(Br), ubound(Br))
        CALL h5_add_double_1(h5_id, trim(tempch)//"Br_Re", &
                            real(Br), lbound(Br), ubound(Br))
        CALL h5_add_double_1(h5_id, trim(tempch)//"Br_Im", &
                            dimag(Br), lbound(Br), ubound(Br))
 
        CALL h5_add_double_1(h5_id, trim(tempch)//"Jpe_abs", &
                            abs(Jpe), lbound(Jpe), ubound(Jpe))
        CALL h5_add_double_1(h5_id, trim(tempch)//"Jpi_abs", &
                            abs(Jpi), lbound(Jpi), ubound(Jpi))
        CALL h5_add_double_1(h5_id, trim(tempch)//"dqle22", &
                                dqle22, lbound(dqle22), ubound(dqle22))

        if (misalign_diffusion .eqv. .true.) then
            call write_misalignment_data_to_hdf5(tempch)
        end if 
        if (diagnostics_output) then
            call write_dql_Br_Jp_profiles_to_hdf5(tempch)
        end if

        CALL h5_close(h5_id)
        CALL h5_deinit()
    else
        do ipoi = 1, npoib
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
end subroutine write_fields_currs_transp_coefs_to_h5


subroutine write_misalignment_data_to_hdf5(tempch)

    use h5mod
    use control_mod, only: debug_mode
    use grid_mod
    use wave_code_data
    
    implicit none

    character(*), intent(in) :: tempch

    if (debug_mode) write(*,*) "Writing misalignment diffusion to hdf5"
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
    if (debug_mode) write(*,*) "Carry on from writing misalignment diffusion to hdf5"

end subroutine

subroutine write_dql_Br_Jp_profiles_to_hdf5(tempch)

    use grid_mod
    use wave_code_data
    use h5mod
    use baseparam_mod, only: c

    implicit none

    character(*), intent(in) :: tempch

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

end subroutine

subroutine writefort9999
        
    use grid_mod, only: dqle11, dqli11, rb, rc, npoib
    use diag_mod, only: timscal_dql, timscal_dqli, ind_dqle, ind_dqli 
    use time_evolution, only: dqle11_prev, dqli11_prev, determine_Dql_diagnostic
    use ParallelTools, only: irank
    use h5mod

    implicit none

    integer :: ipoi
        
    if (irank .eq. 0) then

        call determine_Dql_diagnostic

        if (debug_mode) print *, 'Debug: timscal_dqle = ', sngl(timscal_dql) &
            , 'timscal_dqli = ', sngl(timscal_dqli)
        if (debug_mode) print *, 'Debug: maximum dqle at r = ', rc(ind_dqle(1)) &
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

end subroutine

