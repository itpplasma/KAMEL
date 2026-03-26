subroutine det_balance_eqs_source_terms

    ! calculates source terms in the balance equations. Is determined by assuming steady state.

    use grid_mod, only : y, dery, dery_equisource, nbaleqs, set_boundary_condition, npoi
    use plasma_parameters, only: params
    use control_mod, only: data_verbosity
    use logger_m, only: log_debug
    use h5mod
    use matrix_mod
    use QLBalance_kinds, only: dp
    use rhs_balance_m, only: rhs_balance, initialize_rhs

    implicit none

    integer :: ipoi, ieq, i, k
    real(dp) :: x

    call log_debug("Generating starting source")

    call set_boundary_condition

    do ipoi=1,npoi
        do ieq=1,nbaleqs
            i=nbaleqs*(ipoi-1)+ieq
            y(i)=params(ieq,ipoi)
        enddo
    enddo

    call log_debug("Before initialize_rhs")
    call initialize_rhs(y,dery)

    dery_equisource=0.d0
    call log_debug("Before rhs_balance")
    call rhs_balance(x,y,dery)

    do k = 1,nz
        dery_equisource(irow(k))=dery_equisource(irow(k))-amat(k)*y(icol(k))
    end do

    dery_equisource=dery_equisource-rhsvec

    if (data_verbosity >= 2) then
        call write_balance_eqs_source_terms
    end if

end subroutine det_balance_eqs_source_terms

subroutine write_balance_eqs_source_terms

    use grid_mod, only : dery_equisource, npoi
    use h5mod
    use control_mod, only: ihdf5IO
    use logger_m, only: log_debug

    implicit none

    character(len=1024) :: tempch
    integer :: ipoi

    call log_debug("Writing equisource")

    if (ihdf5IO .eq. 1) then
        tempch = "/"//trim(h5_mode_groupname)//"/equisource.dat"

        CALL h5_init()
        CALL h5_open_rw(path2out, h5_id)
        call create_group_if_not_existent("/"//trim(h5_mode_groupname))

        CALL h5_obj_exists(h5_id, trim(tempch), h5_exists_log)
        if (h5_exists_log) then
            CALL h5_delete(h5_id, trim(tempch))
        end if

        CALL h5_define_unlimited_matrix(h5_id, trim(tempch), &
            H5T_NATIVE_DOUBLE, (/4,-1/), dataset_id)
        do ipoi = 1, npoi
            CALL h5_append_double_1(dataset_id, &
                dery_equisource(4*(ipoi-1)+1:4*ipoi), ipoi)
        enddo
        CALL h5_close(h5_id)
        CALL h5_deinit()
    else
        open(666,file='equisource.dat')
        do ipoi=1,npoi
            write(666,*) dery_equisource(4*(ipoi-1)+1:4*ipoi)
        enddo
        close (666)
    end if

end subroutine
