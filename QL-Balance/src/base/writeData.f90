module writeData_m
    use QLBalance_kinds, only: dp

    implicit none
    private

    public :: write_fields_currs_transp_coefs_to_h5
    public :: writefort9999
    public :: writefort9999_stellarator
    public :: write_D_one_over_nu_to_h5

contains

subroutine write_fields_currs_transp_coefs_to_h5(time_ind)
    use grid_mod, only: npoib, dqle11, dqle12, dqle22, dqli11, dqli12, dqli22, &
                        T_EM_phi_e, T_EM_phi_i
    use baseparam_mod, only: c
    use control_mod, only: ihdf5IO, data_verbosity, misalign_diffusion
    use logger_m, only: log_debug
    use wave_code_data
    use QLbalance_diag, only: iunit_diag
    use h5mod

    implicit none

    integer, intent(in) :: time_ind
    integer :: ipoi
    character(len=1024) :: tempch

    if (ihdf5IO .eq. 1) then

        CALL h5_init()
        CALL h5_open_rw(path2out, h5_id)
        call log_debug(trim(h5_mode_groupname))
        tempch = "/"//trim(h5_mode_groupname)//"/LinearProfiles"

        write (tempch, '(A,"/",I0,"/")') trim(tempch), time_ind

        call create_group_if_not_existent(tempch)

        CALL h5_add_double_1(h5_id, trim(tempch)//"r", &
                            r, lbound(r), ubound(r))
        CALL h5_add_double_1(h5_id, trim(tempch)//"Br_abs", &
                            abs(Br), lbound(Br), ubound(Br))
        CALL h5_add_double_1(h5_id, trim(tempch)//"Br_Re", &
                            real(Br), lbound(Br), ubound(Br))
        CALL h5_add_double_1(h5_id, trim(tempch)//"Br_Im", &
                            dimag(Br), lbound(Br), ubound(Br))
        CALL h5_add_double_1(h5_id, trim(tempch)//"Es_abs", &
                            abs(Es), lbound(Es), ubound(Es))

        CALL h5_add_double_1(h5_id, trim(tempch)//"Jpe_abs", &
                            abs(Jpe), lbound(Jpe), ubound(Jpe))
        CALL h5_add_double_1(h5_id, trim(tempch)//"Jpe_Re", &
                            real(Jpe), lbound(Jpe), ubound(Jpe))
        CALL h5_add_double_1(h5_id, trim(tempch)//"Jpe_Im", &
                            dimag(Jpe), lbound(Jpe), ubound(Jpe))
        CALL h5_add_double_1(h5_id, trim(tempch)//"Jpi_abs", &
                            abs(Jpi), lbound(Jpi), ubound(Jpi))
        CALL h5_add_double_1(h5_id, trim(tempch)//"Jpi_Re", &
                            real(Jpi), lbound(Jpi), ubound(Jpi))
        CALL h5_add_double_1(h5_id, trim(tempch)//"Jpi_Im", &
                            dimag(Jpi), lbound(Jpi), ubound(Jpi))
        CALL h5_add_double_1(h5_id, trim(tempch)//"dqle22", &
                                dqle22, lbound(dqle22), ubound(dqle22))

        if (misalign_diffusion .eqv. .true.) then
            call write_misalignment_data_to_hdf5(tempch)
        end if
        if (data_verbosity >= 2) then
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
    use logger_m, only: log_debug
    use grid_mod
    use wave_code_data

    implicit none

    character(*), intent(in) :: tempch

    call log_debug("Writing misalignment diffusion to hdf5")
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
    call log_debug("Carry on from writing misalignment diffusion to hdf5")

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
    call h5_add_double_1(h5_id, trim(tempch) // "T_EM_phi_e", T_EM_phi_e, lbound(T_EM_phi_e), &
                         ubound(T_EM_phi_e))
    call h5_add_double_1(h5_id, trim(tempch) // "T_EM_phi_i", T_EM_phi_i, lbound(T_EM_phi_i), &
                         ubound(T_EM_phi_i))
    call h5_add_double_1(h5_id, trim(tempch) // "torque_ntv", torque_ntv, &
                         lbound(torque_ntv), ubound(torque_ntv), &
                         comment="Total NTV torque calculated by NEO-RT")
    CALL h5_add_double_1(h5_id, trim(tempch)//"Br-ckpEs_om_E", &
                        abs(Br - c*kp*Es/om_E), lbound(Br), ubound(Br))
    CALL h5_add_double_1(h5_id, trim(tempch)//"Br-cksEp_om_E", &
                        abs(Br - c*ks*Ep/om_E), lbound(Br), ubound(Br))
    !CALL h5_add_double_1(h5_id, trim(tempch)//"Jpe_abs", &
                        !abs(Jpe), lbound(Jpe), ubound(Jpe))
    !CALL h5_add_double_1(h5_id, trim(tempch)//"Jpi_abs", &
                        !abs(Jpi), lbound(Jpi), ubound(Jpi))
    !CALL h5_add_double_1(h5_id, trim(tempch)//"JpeJpi_abs", &
                        !abs(Jpe + Jpi), lbound(Jpe), ubound(Jpe))

end subroutine

subroutine writefort9999(dqle11_prev, dqli11_prev)

    use grid_mod, only: dqle11, dqli11, rb, rc, npoib
    use QLbalance_diag, only: timscal_dql, timscal_dqli, ind_dqle, ind_dqli
    use h5mod
    use logger_m, only: log_debug

    implicit none

    real(dp), dimension(:), intent(in) :: dqle11_prev
    real(dp), dimension(:), intent(in) :: dqli11_prev
    integer :: ipoi
    character(256) :: buf

    write(buf, '(A,ES12.4,A,ES12.4)') 'timscal_dqle = ', &
        sngl(timscal_dql), ' timscal_dqli = ', sngl(timscal_dqli)
    call log_debug(trim(buf))
    write(buf, '(A,ES12.4,A,ES12.4)') 'maximum dqle at r = ', &
        rc(ind_dqle(1)), ' maximum dqli at r = ', rc(ind_dqli(1))
    call log_debug(trim(buf))
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

end subroutine

subroutine writefort9999_stellarator(dqle11_prev, dqli11_prev)

    use grid_mod, only: dqle11, dqli11, rb, rc, npoib
    use QLbalance_diag, only: timscal_dql, timscal_dqli, ind_dqle, ind_dqli
    use h5mod
    use logger_m, only: log_debug

    implicit none

    real(dp), dimension(:), intent(in) :: dqle11_prev
    real(dp), dimension(:), intent(in) :: dqli11_prev
    integer :: ipoi
    character(256) :: buf

    write(buf, '(A,ES12.4,A,ES12.4)') 'timscal_dqle = ', &
        sngl(timscal_dql), ' timscal_dqli = ', sngl(timscal_dqli)
    call log_debug(trim(buf))
    write(buf, '(A,ES12.4,A,ES12.4)') 'maximum dqle at r = ', &
        rc(ind_dqle(1)), ' maximum dqli at r = ', rc(ind_dqli(1))
    call log_debug(trim(buf))
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

end subroutine

subroutine write_D_one_over_nu_to_h5(time_ind)

    use grid_mod, only: Donue11, Donue12, Donue21, Donue22, &
        Donui11, Donui12, Donui21, Donui22
    use h5mod

    implicit none

    integer, intent(in) :: time_ind
    character(256) :: tempch
    tempch = "/"//trim(h5_mode_groupname)//"/LinearProfiles"
    write (tempch, '(A,"/",I0,"/")') trim(tempch), time_ind

    CALL h5_init()
    CALL h5_open_rw(path2out, h5_id)
    CALL h5_add_double_1(h5_id, trim(tempch)//"Donue11", &
                        Donue11, lbound(Donue11), ubound(Donue11))
    CALL h5_add_double_1(h5_id, trim(tempch)//"Donue12", &
                        Donue12, lbound(Donue12), ubound(Donue12))
    CALL h5_add_double_1(h5_id, trim(tempch)//"Donue21", &
                        Donue21, lbound(Donue21), ubound(Donue21))
    CALL h5_add_double_1(h5_id, trim(tempch)//"Donue22", &
                        Donue22, lbound(Donue22), ubound(Donue22))

    CALL h5_add_double_1(h5_id, trim(tempch)//"Donui11", &
                        Donui11, lbound(Donui11), ubound(Donui11))
    CALL h5_add_double_1(h5_id, trim(tempch)//"Donui12", &
                        Donui12, lbound(Donui12), ubound(Donui12))
    CALL h5_add_double_1(h5_id, trim(tempch)//"Donui21", &
                        Donui21, lbound(Donui21), ubound(Donui21))
    CALL h5_add_double_1(h5_id, trim(tempch)//"Donui22", &
                        Donui22, lbound(Donui22), ubound(Donui22))
    CALL h5_close(h5_id)
    CALL h5_deinit()
end subroutine

end module writeData_m
