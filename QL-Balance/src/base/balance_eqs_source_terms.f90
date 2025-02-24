subroutine det_balance_eqs_source_terms
    
    ! calculates source terms in the balance equations. Is determined by assuming steady state.

    use grid_mod, only : y,dery,dery_equisource &
                    , nbaleqs,neqset,iboutype,npoic
    use plasma_parameters, only: params

    use control_mod, only: iwrite, ihdf5IO, diagnostics_output, debug_mode, irf
    use h5mod
    use matrix_mod

    implicit none

    integer :: ipoi,ieq,i,npoi,icount,k
    double precision :: x
    character(len=1024) :: tempch

    if (debug_mode) write(*,*) "Debug: Generating starting source"

    if(iboutype.eq.1) then
        npoi=npoic-1
    else
        npoi=npoic
    endif

    do ipoi=1,npoi
        do ieq=1,nbaleqs
            i=nbaleqs*(ipoi-1)+ieq
            y(i)=params(ieq,ipoi)
        enddo
    enddo

    if (debug_mode) print *, "Debug: Before initialize_rhs"
    call initialize_rhs(y,dery)

    dery_equisource=0.d0
    if (debug_mode) print *, "Debug: Before rhs_balance"
    call rhs_balance(x,y,dery)

    do k=1,nz
        dery_equisource(irow(k))=dery_equisource(irow(k))-amat(k)*y(icol(k))
    end do

    dery_equisource=dery_equisource-rhsvec

    if (diagnostics_output) then
        if (debug_mode) write(*,*) "Debug: Writing equisource"
        if (ihdf5IO .eq. 1) then
            tempch = "/"//trim(h5_mode_groupname)//"/equisource.dat"

            CALL h5_init()
            CALL h5_open_rw(path2out, h5_id)

            CALL h5_obj_exists(h5_id, "/"//trim(h5_mode_groupname), &
                h5_exists_log)
            if (.not. h5_exists_log) then
                CALL h5_define_group(h5_id, trim(h5_mode_groupname), &
                    group_id_1)
            end if

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
                do ieq=1,nbaleqs
                    i=nbaleqs*(ipoi-1)+ieq
                    params(ieq,ipoi)=y(i)
                enddo
                write(666,*) dery_equisource(4*(ipoi-1)+1:4*ipoi)
            enddo
            close (666)
        end if
    end if

end subroutine det_balance_eqs_source_terms
