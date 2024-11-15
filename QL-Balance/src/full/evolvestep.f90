  subroutine evolvestep(timstep,eps)

    use grid_mod, only : nbaleqs,neqset,iboutype,npoic,y,dery
    use plasma_parameters, only: params
    USE sparse_mod, ONLY : sparse_talk,sparse_solve_method,sparse_solve, &
        column_full2pointer,remap_rc,sparse_solver_test
    use matrix_mod
    use recstep_mod, only : timstep_arr
    use control_mod, only : debug_mode

    implicit none

    external :: rhs_balance !, rhs_func, gslint
    integer :: ipoi,ieq,i,k,npoi,iopt,nz_sp,nz_sq,nrow,ncol
    double precision :: timstep,x1,x2,eps,time_start!,time_factorization,time_solver
    INTEGER,          DIMENSION(:),   ALLOCATABLE :: ipcol,irow_sp,icol_sp
    DOUBLE PRECISION, DIMENSION(:),   ALLOCATABLE :: amat_sp,bvec_sp

    x1=0.d0
    sparse_talk = .true.
    !  x2=timstep

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

    call initialize_rhs(y,dery)
    call rhs_balance(x1,y,dery)

    call write_amat

    nz_sp=nz+nsize
    nrow=nsize
    ncol=nsize

    allocate(ipcol(nsize))
    allocate(amat_sp(nz_sp),irow_sp(nz_sp),icol_sp(nz_sp),bvec_sp(nsize))
    irow_sp(1:nz)=irow
    icol_sp(1:nz)=icol
    !  amat_sp(1:nz)=-timstep*amat
    amat_sp(1:nz)=-timstep_arr(irow)*amat!

    k=nz
    do i=1,nsize
        k=k+1
        irow_sp(k)=i 
        icol_sp(k)=i 
        amat_sp(k)=1.d0
    enddo

    !  bvec_sp=y+timstep*dery
    bvec_sp=y+timstep_arr*dery

    call  remap_rc(nz_sp,nz_sq,irow_sp,icol_sp,amat_sp)

    nz_sp=nz_sq

    CALL column_full2pointer(icol_sp(1:nz_sp),ipcol)

    if (.false.) then
        iopt=0
        !write(*,*) " timstep_arr(1) in evolvestep = ", timstep_arr(1)
        !write(*,*) " amat_sp(1) = ", amat_sp(1)
        !call cpu_time(time_start)
        call write_bvec_sp_to_txt
        
        !CALL sparse_solve(nrow,ncol,nz_sp,irow_sp(1:nz_sp),ipcol,amat_sp(1:nz_sp),       &
        CALL sparse_solve(nrow, ncol, nz_sp, irow_sp, ipcol, amat_sp, bvec_sp, iopt)
                          !bvec_sp,iopt)
    else
        iopt=1
        CALL sparse_solve(nrow,ncol,nz_sp,irow_sp(1:nz_sp),ipcol,amat_sp(1:nz_sp),       &
                      bvec_sp,iopt)
    
        !  call cpu_time(time_factorization)
        !  print *,'factorization completed ',time_factorization - time_start,' sec'
        print *, "After factorization"
        !
        iopt=2
        !
        ! Solution of inhomogeneus equation (account of sources):
        !
        CALL sparse_solve(nrow,ncol,nz_sp,irow_sp(1:nz_sp),ipcol,amat_sp(1:nz_sp),       &
                            bvec_sp,iopt)

        !iopt=3
        !CALL sparse_solve(nrow,ncol,nz_sp,irow_sp(1:nz),ipcol,amat_sp(1:nz_sp),bvec_sp,iopt)

    end if
    
    y=bvec_sp

    deallocate(ipcol)
    deallocate(amat_sp,irow_sp,icol_sp,bvec_sp)

    do ipoi=1,npoi
        do ieq=1,nbaleqs
            i=nbaleqs*(ipoi-1)+ieq
            params(ieq,ipoi)=y(i)
        enddo
    enddo

    contains

    subroutine write_bvec_sp_to_txt

      implicit none
      
      open(666,file='bvec_sp.txt')
      write(666,*) bvec_sp
      close (666)

    end subroutine

    subroutine write_amat

        implicit none

        open(666, file='amat.txt')
        write(666,*) amat
        close(666)

    end subroutine

end subroutine evolvestep





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
  
    irf = 0
    print *, "Before initialize_rhs"
    call initialize_rhs(y,dery)

    dery_equisource=0.d0
    print *, "Before rhs_balance"
    call rhs_balance(x,y,dery)
    irf = 1

    do k=1,nz
      dery_equisource(irow(k))=dery_equisource(irow(k))-amat(k)*y(icol(k))
    end do

    dery_equisource=dery_equisource-rhsvec

    if (diagnostics_output) then
      write(*,*) "Writing equisource"
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


    if (iwrite .eq. 1) then
        write(666,*) dery_equisource
        close (666)
    end if

end subroutine det_balance_eqs_source_terms
