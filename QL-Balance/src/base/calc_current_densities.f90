subroutine calc_parallel_current_directly
    ! this subroutine calculates the electron parallel current (eq. (60) in Heyn et. al 2014)

    use grid_mod, only: npoib, rb, Ercov
    use plasma_parameters, only: params_b, ddr_params_nl
    use baseparam_mod, only: e_charge, p_mass, c, e_mass, ev
    use control_mod, only: ihdf5IO, data_verbosity, write_gyro_current, &
        gyro_current_study
    use logger_m, only: log_debug
    use h5mod
    use wave_code_data
    use QLBalance_kinds, only: dp

    implicit none

    integer :: ipoi, i, iunit, mnmax
    real(dp), dimension(:), allocatable :: x1, x2, vT, A1, A2
    complex(dp), dimension(:), allocatable :: curr_e_par
    complex(dp), dimension(:, :, :), allocatable :: symbI
    real(dp), dimension(1) :: coll_fac = (/1/)!(/10.0, 100.0, 1.0e3, 1.0e4, 1.0e5/)
    integer :: study_i_omE
    integer :: study_j_nue

    character(len=1024) :: tempch

    iunit = 731
    mnmax = 3
    allocate (x1(npoib), x2(npoib), A1(npoib), A2(npoib), symbI(0:mnmax, 0:mnmax, npoib))
    allocate (curr_e_par(npoib), vT(npoib))

    vT = sqrt(params_b(3, :)/e_mass)


    if (gyro_current_study .eq. 0) then
        x1 = kp*vT/nue
        x2 = -om_E/nue

        do i = 1, npoib
            call getIfunc(x1(i), x2(i), symbI(:, :, i))
        end do


        A2 = ddr_params_nl(3, :)/params_b(3, :)
        A1 = ddr_params_nl(1, :)/params_b(1, :) + e_charge*Ercov/params_b(3, :) - 1.5d0*A2

        curr_e_par = e_charge*params_b(1, :)*vT/(nue*B0) &
                 *(c*Es*((A1 + A2)*symbI(1, 0, :) + 0.5d0*A2*symbI(2, 1, :)) &
                    + vT*Br*((A1 + A2)*symbI(1, 1, :) + 0.5d0*A2*symbI(3, 1, :)))

        call h5_init()
        call h5_open_rw(path2out, h5_id)
        tempch = "/"//trim(h5_mode_groupname)//"/par_current_e/"
        call log_debug("In group: "//trim(tempch))

        call h5_obj_exists(h5_id, trim(tempch), h5_exists_log)
        if (.not. h5_exists_log) then
            call h5_define_group(h5_id, trim(tempch), group_id_1)
            call h5_close_group(group_id_1)
        end if


        ! TODO: fix this, makes problem when writing
        !call h5_add_double_1(h5_id, trim(tempch)//"rb", &
                !rb, lbound(rb), ubound(rb))
        !call h5_add_double_1(h5_id, trim(tempch)//"x2", &
                !x2, lbound(x2), ubound(x2))

        if (write_gyro_current) then
            call log_debug("writing par_current_e.dat")
            ! Write out gyro current which is different to KiLCA current Jpe.
            ! The gyro current is calculated from (60) in Heyn et. al 2014
            !call h5_add_double_1(h5_id, trim(tempch)//"par_current_e_real", &
            !   real(curr_e_par), lbound(real(curr_e_par)), ubound(real(curr_e_par)))
            !all h5_add_double_1(h5_id, trim(tempch)//"par_current_e_imag", &
            !   dimag(curr_e_par), lbound(dimag(curr_e_par)), ubound(dimag(curr_e_par)))
            !all h5_add_double_1(h5_id, trim(tempch)//"Jpe_real", &
            !   real(Jpe), lbound(real(Jpe)), ubound(real(Jpe)))
            !all h5_add_double_1(h5_id, trim(tempch)//"Jpe_imag", &
            !   dimag(Jpe), lbound(dimag(Jpe)), ubound(dimag(Jpe)))
            !all h5_add_double_1(h5_id, trim(tempch)//"Jse_real", &
            !   real(Jse), lbound(real(Jse)), ubound(real(Jse)))
            !all h5_add_double_1(h5_id, trim(tempch)//"Jse_imag", &
            !   dimag(Jse), lbound(dimag(Jse)), ubound(dimag(Jse)))
            !all h5_add_double_1(h5_id, trim(tempch)//"Jre_real", &
            !   real(Jpe), lbound(real(Jpe)), ubound(real(Jpe)))
            !all h5_add_double_1(h5_id, trim(tempch)//"Jre_imag", &
            !   dimag(Jre), lbound(dimag(Jre)), ubound(dimag(Jre)))
            call h5_add_double_1(h5_id, trim(tempch)//"kp", &
                kp, lbound(kp), ubound(kp))
            call h5_add_double_1(h5_id, trim(tempch)//"ks", &
                ks, lbound(ks), ubound(ks))
            call h5_add_double_1(h5_id, trim(tempch)//"x1", &
                x1, lbound(x1), ubound(x1))
            call h5_add_double_1(h5_id, trim(tempch)//"nue", &
                nue, lbound(nue), ubound(nue))
            call h5_add_double_1(h5_id, trim(tempch)//"om_E", &
                om_E, lbound(om_E), ubound(om_E))
            call h5_add_double_1(h5_id, trim(tempch)//"vT", &
                vT, lbound(vT), ubound(vT))

            call h5_add_double_1(h5_id, trim(tempch)//"I00_re", &
                real(symbI(0, 0, :)), lbound(symbI(0, 0, :)), ubound(symbI(0, 0, :)))
            call h5_add_double_1(h5_id, trim(tempch)//"I00_im", &
                dimag(symbI(0, 0, :)), lbound(symbI(0, 0, :)), ubound(symbI(0, 0, :)))

            call h5_add_double_1(h5_id, trim(tempch)//"I20_re", &
                real(symbI(2, 0, :)), lbound(symbI(2, 0, :)), ubound(symbI(2, 0, :)))
            call h5_add_double_1(h5_id, trim(tempch)//"I20_im", &
                dimag(symbI(2, 0, :)), lbound(symbI(2, 0, :)), ubound(symbI(2, 0, :)))

            call h5_add_double_1(h5_id, trim(tempch)//"I22_re", &
                real(symbI(2, 2, :)), lbound(symbI(2, 2, :)), ubound(symbI(2, 2, :)))
            call h5_add_double_1(h5_id, trim(tempch)//"I22_im", &
                dimag(symbI(2, 2, :)), lbound(symbI(2, 2, :)), ubound(symbI(2, 2, :)))

        end if ! write_gyro_current
        call h5_close(h5_id)
        call h5_deinit()

        if (data_verbosity >= 2) then
            if (ihdf5IO .eq. 1) then
                call log_debug("writing par_current_e.dat")
                call h5_init()
                call h5_open_rw(path2out, h5_id)

                ! par_current_e data
                tempch = "/"//trim(h5_mode_groupname)//"/par_current_e.dat"

                call h5_obj_exists(h5_id, trim(tempch), h5_exists_log)
                if (h5_exists_log) then
                    call h5_delete(h5_id, trim(tempch))
                end if

                call h5_define_unlimited_matrix(h5_id, trim(tempch), &
                                            H5T_NATIVE_DOUBLE, (/-1, 5/), dataset_id)
                call h5_append_double_1(dataset_id, rb, 1)
                call h5_append_double_1(dataset_id, real(curr_e_par), 2)
                call h5_append_double_1(dataset_id, dimag(curr_e_par), 3)
                call h5_append_double_1(dataset_id, real(Jpe), 4)
                call h5_append_double_1(dataset_id, dimag(Jpe), 5)

                ! cond_e data
                tempch = "/"//trim(h5_mode_groupname)//"/cond_e.dat"

                call h5_obj_exists(h5_id, trim(tempch), h5_exists_log)
                if (h5_exists_log) then
                    call h5_delete(h5_id, trim(tempch))
                end if

                call h5_define_unlimited_matrix(h5_id, trim(tempch), &
                                            H5T_NATIVE_DOUBLE, (/-1, 9/), dataset_id)
                call h5_append_double_1(dataset_id, rb, 1)
                call h5_append_double_1(dataset_id, real(symbI(1, 0, :)), 2)
                call h5_append_double_1(dataset_id, dimag(symbI(1, 0, :)), 3)
                call h5_append_double_1(dataset_id, real(symbI(1, 1, :)), 4)
                call h5_append_double_1(dataset_id, dimag(symbI(1, 1, :)), 5)
                call h5_append_double_1(dataset_id, real(symbI(2, 1, :)), 6)
                call h5_append_double_1(dataset_id, dimag(symbI(2, 1, :)), 7)
                call h5_append_double_1(dataset_id, real(symbI(3, 1, :)), 8)
                call h5_append_double_1(dataset_id, dimag(symbI(3, 1, :)), 9)

                call h5_close(h5_id)
                call h5_deinit()
                else ! ihdf5IO .eq. 1
                    open (iunit, file='par_current_e.dat')
                    open (10000, file='cond_e.dat')
                    do ipoi = 1, npoib
                        write (iunit, *) rb(ipoi) &
                            , real(curr_e_par(ipoi)), dimag(curr_e_par(ipoi)) &
                            , real(Jpe(ipoi)), dimag(Jpe(ipoi))
                        write (10000, *) rb(ipoi) &
                            , real(symbI(1, 0, ipoi)), dimag(symbI(1, 0, ipoi)) &
                            , real(symbI(1, 1, ipoi)), dimag(symbI(1, 1, ipoi)) &
                            , real(symbI(2, 1, ipoi)), dimag(symbI(2, 1, ipoi)) &
                            , real(symbI(3, 1, ipoi)), dimag(symbI(3, 1, ipoi))

                    end do
                    close (10000)
                    close (iunit)
                end if
            end if

    elseif (gyro_current_study .eq. 1) then ! used to scan over nue and omega_E, deprecated
        write(*,*) " - - - - - - - - - "
        write(*,*) "Gyro current study"

        CALL h5_init()
        CALL h5_open_rw(path2out, h5_id)
        tempch = "/"//trim(h5_mode_groupname)//"/gyro_current_study/"
        call log_debug("In group: "//trim(tempch))

        CALL h5_define_group(h5_id, trim(tempch), group_id_1)
        CALL h5_close_group(group_id_1)

        ! Quantities that are not affected by parameter study:
        CALL h5_add_double_1(h5_id, trim(tempch)//"rb", &
            rb, lbound(rb), ubound(rb))
        CALL h5_add_double_1(h5_id, trim(tempch)//"nue0", &
            nue, lbound(nue), ubound(nue))
        CALL h5_add_double_1(h5_id, trim(tempch)//"om_E", &
            om_E, lbound(om_E), ubound(om_E))

        CALL h5_add_double_1(h5_id, trim(tempch)//"B0", &
            B0, lbound(B0), ubound(B0))

        CALL h5_add_double_1(h5_id, trim(tempch)//"Br_real", &
            real(Br), lbound(real(Br)), ubound(real(Br)))
        CALL h5_add_double_1(h5_id, trim(tempch)//"Es_real", &
            real(Es), lbound(real(Es)), ubound(real(Es)))

        CALL h5_add_double_1(h5_id, trim(tempch)//"Br_imag", &
            dimag(Br), lbound(dimag(Br)), ubound(dimag(Br)))
        CALL h5_add_double_1(h5_id, trim(tempch)//"Es_imag", &
            dimag(Es), lbound(dimag(Es)), ubound(dimag(Es)))
        ! - - - - - - - - - - - - - - - - - - - - - - - - -

        A2 = ddr_params_nl(3, :)/params_b(3, :)
        A1 = ddr_params_nl(1, :)/params_b(1, :) + e_charge*Ercov/params_b(3, :) - 1.5d0*A2

        CALL h5_add_double_1(h5_id, trim(tempch)//"A1", &
            A1, lbound(A1), ubound(A1))
        CALL h5_add_double_1(h5_id, trim(tempch)//"A2", &
            A2, lbound(A2), ubound(A2))

        ! vary om_E/nue
        do study_i_omE = 1, 1
            tempch = "/"//trim(h5_mode_groupname)//"/gyro_current_study/"
            write(tempch, "(A,i4.4,A)") trim(tempch), study_i_omE, "/"

            CALL h5_define_group(h5_id, trim(tempch), group_id_1)
            CALL h5_close_group(group_id_1)

            do study_j_nue = 1, size(coll_fac)

                om_E = om_E * study_i_omE
                nue = nue * coll_fac(study_j_nue)
                x2 = -om_E/(nue)! * study_i * 0.5
                x1 = kp*vT/(nue)! * study_i * 0.5

                do i = 1, npoib
                    call getIfunc(x1(i), x2(i), symbI(:, :, i))
                end do
                tempch = "/"//trim(h5_mode_groupname)//"/gyro_current_study/"
                write(tempch, "(A,i4.4,A,i4.4,A)") trim(tempch), study_i_omE, "/", &
                    study_j_nue, "/"
                call log_debug("In group: "//trim(tempch))

                CALL h5_define_group(h5_id, trim(tempch), group_id_1)
                CALL h5_close_group(group_id_1)
                ! write I functions
                ! real part
                !CALL h5_add_double_1(h5_id, trim(tempch)//"I10_real", &
                !    real(symbI(1,0,:)), lbound(real(symbI(1,0,:))), &
                !    ubound(real(symbI(1,0,:))))
                !CALL h5_add_double_1(h5_id, trim(tempch)//"I11_real", &
                !    real(symbI(1,1,:)), lbound(real(symbI(1,1,:))), &
                !    ubound(real(symbI(1,1,:))))
                !CALL h5_add_double_1(h5_id, trim(tempch)//"I21_real", &
                !    real(symbI(2,1,:)), lbound(real(symbI(2,1,:))), &
                !    ubound(real(symbI(2,1,:))))
                !CALL h5_add_double_1(h5_id, trim(tempch)//"I31_real", &
                !    real(symbI(3,1,:)), lbound(real(symbI(3,1,:))), &
                !    ubound(real(symbI(3,1,:))))
                ! imag part

                !CALL h5_add_double_1(h5_id, trim(tempch)//"I10_imag", &
                !    dimag(symbI(1,0,:)), lbound(dimag(symbI(1,0,:))),&
                !    ubound(dimag(symbI(1,0,:))))
                !CALL h5_add_double_1(h5_id, trim(tempch)//"I11_imag", &
                !    dimag(symbI(1,1,:)), lbound(dimag(symbI(1,1,:))), &
                !    ubound(dimag(symbI(1,1,:))))
                !CALL h5_add_double_1(h5_id, trim(tempch)//"I21_imag", &
                !    dimag(symbI(2,1,:)), lbound(dimag(symbI(2,1,:))), &
                !    ubound(dimag(symbI(2,1,:))))
                !CALL h5_add_double_1(h5_id, trim(tempch)//"I31_imag", &
                !    dimag(symbI(3,1,:)), lbound(dimag(symbI(3,1,:))), &
                !    ubound(dimag(symbI(3,1,:))))

                curr_e_par = e_charge*params_b(1, :)*vT/(nue*B0) &
                 *(c*Es*((A1 + A2)*symbI(1, 0, :) + 0.5d0*A2*symbI(2, 1, :)) &
                    + vT*Br*((A1 + A2)*symbI(1, 1, :) + 0.5d0*A2*symbI(3, 1, :)))

                call get_current_densities_from_wave_code(flre_cd_ptr(1), dim_r, r, &
                                                    m_vals(1), n_vals(1), Jri, Jsi, Jpi, Jre, Jse, Jpe)

                CALL h5_add_double_1(h5_id, trim(tempch)//"je_real", &
                    real(curr_e_par), lbound(real(curr_e_par)), ubound(real(curr_e_par)))
                CALL h5_add_double_1(h5_id, trim(tempch)//"je_imag", &
                    dimag(curr_e_par), lbound(dimag(curr_e_par)), ubound(dimag(curr_e_par)))

                CALL h5_add_double_1(h5_id, trim(tempch)//"Je_real", &
                    real(Jpe), lbound(real(Jpe)), ubound(real(Jpe)))
                CALL h5_add_double_1(h5_id, trim(tempch)//"Je_imag", &
                    dimag(Jpe), lbound(dimag(Jpe)), ubound(dimag(Jpe)))
            end do ! study_j_nue
        end do ! study_i_omE

        CALL h5_close(h5_id)
        CALL h5_deinit()

        write(*,*) " - - - - - - - - - "

    elseif (gyro_current_study .eq. 2) then
        write(*,*) " - - - - - - - - - "
        write(*,*) "Writing currents"
        x1 = kp*vT/nue
        x2 = -om_E/nue


        CALL h5_init()
        CALL h5_open_rw(path2out, h5_id)
        tempch = "/"//trim(h5_mode_groupname)//"/currents/"
        call log_debug("In group: "//trim(tempch))

        CALL h5_define_group(h5_id, trim(tempch), group_id_1)
        CALL h5_close_group(group_id_1)

        ! Quantities that are not affected by parameter study:
        CALL h5_add_double_1(h5_id, trim(tempch)//"rb", &
            rb, lbound(rb), ubound(rb))
        CALL h5_add_double_1(h5_id, trim(tempch)//"nue0", &
            nue, lbound(nue), ubound(nue))
        CALL h5_add_double_1(h5_id, trim(tempch)//"om_E", &
            om_E, lbound(om_E), ubound(om_E))

        CALL h5_add_double_1(h5_id, trim(tempch)//"B0", &
            B0, lbound(B0), ubound(B0))

        CALL h5_add_double_1(h5_id, trim(tempch)//"Br_real", &
            real(Br), lbound(real(Br)), ubound(real(Br)))
        CALL h5_add_double_1(h5_id, trim(tempch)//"Es_real", &
            real(Es), lbound(real(Es)), ubound(real(Es)))

        CALL h5_add_double_1(h5_id, trim(tempch)//"Br_imag", &
            dimag(Br), lbound(dimag(Br)), ubound(dimag(Br)))
        CALL h5_add_double_1(h5_id, trim(tempch)//"Es_imag", &
            dimag(Es), lbound(dimag(Es)), ubound(dimag(Es)))
        ! - - - - - - - - - - - - - - - - - - - - - - - - -

        A2 = ddr_params_nl(3, :)/params_b(3, :)
        A1 = ddr_params_nl(1, :)/params_b(1, :) + e_charge*Ercov/params_b(3, :) - 1.5d0*A2

        CALL h5_add_double_1(h5_id, trim(tempch)//"A1", &
            A1, lbound(A1), ubound(A1))
        CALL h5_add_double_1(h5_id, trim(tempch)//"A2", &
            A2, lbound(A2), ubound(A2))

        do i = 1, npoib
            call getIfunc(x1(i), x2(i), symbI(:, :, i))
        end do
        curr_e_par = e_charge*params_b(1, :)*vT/(nue*B0) &
                    *(c*Es*((A1 + A2)*symbI(1, 0, :) + 0.5d0*A2*symbI(2, 1, :)) &
                    + vT*Br*((A1 + A2)*symbI(1, 1, :) + 0.5d0*A2*symbI(3, 1, :)))

        call get_current_densities_from_wave_code(flre_cd_ptr(1), dim_r, r, &
                                    m_vals(1), n_vals(1), Jri, Jsi, Jpi, Jre, Jse, Jpe)

        CALL h5_add_double_1(h5_id, trim(tempch)//"je_real", & ! drift kinetic current (gyrocurrent)
            real(curr_e_par), lbound(real(curr_e_par)), ubound(real(curr_e_par)))
        CALL h5_add_double_1(h5_id, trim(tempch)//"je_imag", &
            dimag(curr_e_par), lbound(dimag(curr_e_par)), ubound(dimag(curr_e_par)))

        CALL h5_add_double_1(h5_id, trim(tempch)//"Je_real", & ! kinetic current from KiLCA
            real(Jpe), lbound(real(Jpe)), ubound(real(Jpe)))
        CALL h5_add_double_1(h5_id, trim(tempch)//"Je_imag", &
            dimag(Jpe), lbound(dimag(Jpe)), ubound(dimag(Jpe)))

        CALL h5_close(h5_id)
        CALL h5_deinit()

        write(*,*) " - - - - - - - - - "

    end if

    deallocate (x1, x2, symbI, curr_e_par, vT)

end subroutine calc_parallel_current_directly


subroutine calc_ion_parallel_current_directly

    use grid_mod, only: npoib, rb, Ercov
    use plasma_parameters, only: params_b, ddr_params_nl
    use baseparam_mod, only: Z_i, e_charge, p_mass, c, e_mass, ev
    use wave_code_data
    use control_mod, only: ihdf5IO, data_verbosity, write_gyro_current
    use logger_m, only: log_debug
    use h5mod
    use QLBalance_kinds, only: dp

    implicit none

    integer :: ipoi, i, iunit, mnmax
    real(dp) :: ei_charge
    real(dp), dimension(:), allocatable :: x1, x2, vT
    complex(dp), dimension(:), allocatable :: curr_i_par
    complex(dp), dimension(:, :, :), allocatable :: symbI

    character(len=1024) :: tempch

    iunit = 731
    mnmax = 3
    allocate (x1(npoib), x2(npoib), symbI(0:mnmax, 0:mnmax, npoib))
    allocate (curr_i_par(npoib), vT(npoib))

    ei_charge = -Z_i*e_charge

    vT = sqrt(params_b(4, :)/p_mass)

    x1 = kp*vT/nui
    x2 = -om_E/nui

    do i = 1, npoib
        call getIfunc(x1(i), x2(i), symbI(:, :, i))
    end do

    ! Here x1 and x2 are used for A_1 and A_2:
    x2 = ddr_params_nl(4, :)/params_b(4, :)
    x1 = ddr_params_nl(1, :)/params_b(1, :) + ei_charge*Ercov/params_b(4, :) - 1.5d0*x2

    curr_i_par = ei_charge*params_b(1, :)*vT/(nui*B0) &
                 *(c*Es*((x1 + x2)*symbI(1, 0, :) + 0.5d0*x2*symbI(2, 1, :)) &
                    + vT*Br*((x1 + x2)*symbI(1, 1, :) + 0.5d0*x2*symbI(3, 1, :)))

    if (write_gyro_current) then
        call log_debug("writing par_current_i.dat")
        ! Write out gyro current which is different to KiLCA current Jpe.
        ! The gyro current is calculated from (60) in Heyn et. al 2014
        call h5_init()
        call h5_open_rw(path2out, h5_id)
        tempch = "/"//trim(h5_mode_groupname)//"/par_current_i/"
        call log_debug("In group: "//trim(tempch))

        call h5_define_group(h5_id, trim(tempch), group_id_1)
        call h5_close_group(group_id_1)

        !all h5_add_double_1(h5_id, trim(tempch)//"rb", &
        !   rb, lbound(rb), ubound(rb))
        !all h5_add_double_1(h5_id, trim(tempch)//"par_current_i_real", &
        !   real(curr_i_par), lbound(real(curr_i_par)), ubound(real(curr_i_par)))
        !all h5_add_double_1(h5_id, trim(tempch)//"par_current_i_imag", &
        !   dimag(curr_i_par), lbound(dimag(curr_i_par)), ubound(dimag(curr_i_par)))
        !all h5_add_double_1(h5_id, trim(tempch)//"Jpi_real", &
        !   real(Jpi), lbound(real(Jpi)), ubound(real(Jpi)))
        !all h5_add_double_1(h5_id, trim(tempch)//"Jpi_imag", &
        !   dimag(Jpi), lbound(dimag(Jpi)), ubound(dimag(Jpi)))
        !all h5_add_double_1(h5_id, trim(tempch)//"Jsi_real", &
        !   real(Jsi), lbound(real(Jsi)), ubound(real(Jsi)))
        !all h5_add_double_1(h5_id, trim(tempch)//"Jsi_imag", &
        !   dimag(Jsi), lbound(dimag(Jsi)), ubound(dimag(Jsi)))
        !all h5_add_double_1(h5_id, trim(tempch)//"Jri_real", &
        !   real(Jpi), lbound(real(Jpi)), ubound(real(Jpi)))
        !all h5_add_double_1(h5_id, trim(tempch)//"Jri_imag", &
        !   dimag(Jri), lbound(dimag(Jri)), ubound(dimag(Jri)))
        call h5_add_double_1(h5_id, trim(tempch)//"kp", &
            kp, lbound(kp), ubound(kp))
        call h5_add_double_1(h5_id, trim(tempch)//"ks", &
            ks, lbound(ks), ubound(ks))
        call h5_add_double_1(h5_id, trim(tempch)//"x1", &
            x1, lbound(x1), ubound(x1))
        call h5_add_double_1(h5_id, trim(tempch)//"x2", &
            x2, lbound(x2), ubound(x2))
        call h5_add_double_1(h5_id, trim(tempch)//"nui", &
            nui, lbound(nui), ubound(nui))

        call h5_close(h5_id)
        call h5_deinit()

    end if ! write_gyro_current

    if (data_verbosity >= 2) then
        if (ihdf5IO .eq. 1) then
            call log_debug("writing par_current_i.dat")
            call h5_init()
            call h5_open_rw(path2out, h5_id)
            tempch = "/"//trim(h5_mode_groupname)//"/par_current_i.dat"

            call h5_obj_exists(h5_id, trim(tempch), h5_exists_log)
            if (h5_exists_log) then
                call h5_delete(h5_id, trim(tempch))
            end if

            call h5_define_unlimited_matrix(h5_id, trim(tempch), &
                                            H5T_NATIVE_DOUBLE, (/-1, 5/), dataset_id)
            call h5_append_double_1(dataset_id, rb, 1)
            call h5_append_double_1(dataset_id, real(curr_i_par), 2)
            call h5_append_double_1(dataset_id, dimag(curr_i_par), 3)
            call h5_append_double_1(dataset_id, real(Jpi), 4)
            call h5_append_double_1(dataset_id, dimag(Jpi), 5)

            call h5_close(h5_id)
            call h5_deinit()

        else
                open (iunit, file='par_current_i.dat')
                do ipoi = 1, npoib
                    write (iunit, *) rb(ipoi), real(curr_i_par(ipoi)), dimag(curr_i_par(ipoi)) &
                        , real(Jpi(ipoi)), dimag(Jpi(ipoi))
                end do
                close (iunit)
            end if
        end if

    deallocate (x1, x2, symbI, curr_i_par, vT)

end subroutine calc_ion_parallel_current_directly
