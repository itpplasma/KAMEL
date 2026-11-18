subroutine det_balance_eqs_source_terms

    ! Determines the equilibrium sources of the balance equations by assuming that
    ! the current state is stationary without them. They are stored in conserved
    ! form, as sources of particles and of pressure [Markl2023 (22), (24), (25)],
    ! so that the dilution by the particle source enters the temperature equations
    ! (see rhs_balance_m).

    use grid_mod, only : y, dery, dery_equisource, equisource_reference, nbaleqs, &
                         set_boundary_condition, npoi
    use plasma_parameters, only: params
    use control_mod, only: data_verbosity
    use logger_m, only: log_debug
    use h5mod
    use matrix_mod
    use QLBalance_kinds, only: dp
    use rhs_balance_m, only: rhs_balance, initialize_rhs, convert_rates_to_equilibrium_sources

    implicit none

    integer :: ipoi, ieq, i, k
    real(dp) :: x
    real(dp), dimension(4) :: sources

    call log_debug("Generating starting source")

    call set_boundary_condition

    do ipoi=1,npoi
        do ieq=1,nbaleqs
            i=nbaleqs*(ipoi-1)+ieq
            y(i)=params(ieq,ipoi)
        enddo
    enddo

    equisource_reference=y

    call log_debug("Before initialize_rhs")
    call initialize_rhs(y,dery)

    dery_equisource=0.d0
    call log_debug("Before rhs_balance")
    call rhs_balance(x,y,dery)

    do k = 1,nz
        dery_equisource(irow(k))=dery_equisource(irow(k))-amat(k)*y(icol(k))
    end do

    dery_equisource=dery_equisource-rhsvec

    do ipoi=1,npoi
        i=nbaleqs*(ipoi-1)
        call convert_rates_to_equilibrium_sources(dery_equisource(i+1:i+nbaleqs), &
                                                  params(:,ipoi), sources)
        dery_equisource(i+1:i+nbaleqs)=sources
    enddo

    if (data_verbosity >= 2) then
        call write_balance_eqs_source_terms
    end if

end subroutine det_balance_eqs_source_terms

subroutine write_balance_eqs_source_terms

    use grid_mod, only : dery_equisource, npoi
    use h5mod
    use control_mod, only: ihdf5IO
    use logger_m, only: log_debug
    use QLBalance_kinds, only: dp

    implicit none

    character(len=1024) :: tempch
    real(dp), dimension(:, :), allocatable :: diagnostic_data
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

        allocate (diagnostic_data(4, npoi))
        do ipoi = 1, npoi
            diagnostic_data(:, ipoi) = dery_equisource(4 * (ipoi - 1) + 1:4 * ipoi)
        enddo
        call h5_add(h5_id, trim(tempch), diagnostic_data, lbound(diagnostic_data), &
                    ubound(diagnostic_data))
        deallocate (diagnostic_data)
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
