subroutine evolvestep(timstep,eps)

    use grid_mod, only : nbaleqs,iboutype,npoic,y,dery
    use plasma_parameters, only: params
    USE sparse_mod
    use matrix_mod, only: amat
    use matrix_mod, only: nz, nsize, irow, icol, amat
    use recstep_mod, only : timstep_arr

    implicit none

    external :: rhs_balance !, rhs_func, gslint
    integer :: ipoi,ieq,i,k,npoi,iopt,nz_sp,nz_sq,nrow,ncol
    double precision :: timstep,x1,eps
    INTEGER,          DIMENSION(:),   ALLOCATABLE :: ipcol,irow_sp,icol_sp
    DOUBLE PRECISION, DIMENSION(:),   ALLOCATABLE :: amat_sp,bvec_sp

    x1=0.d0
    sparse_solve_method = 3

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

    nz_sp=nz+nsize
    nrow=nsize
    ncol=nsize

    allocate(ipcol(nsize))
    allocate(amat_sp(nz_sp),irow_sp(nz_sp),icol_sp(nz_sp),bvec_sp(nsize))
    irow_sp(1:nz)=irow
    icol_sp(1:nz)=icol
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

    call remap_rc(nz_sp,nz_sq,irow_sp,icol_sp,amat_sp)

    nz_sp=nz_sq

    CALL column_full2pointer(icol_sp(1:nz_sp),ipcol)

    iopt=0
    CALL sparse_solve(nrow, ncol, nz_sp, irow_sp, ipcol, amat_sp, bvec_sp, iopt)

    ! for debugging:
    !iopt=1
    !CALL sparse_solve(nrow,ncol,nz_sp,irow_sp(1:nz_sp),ipcol,amat_sp(1:nz_sp),       &
                    !bvec_sp,iopt)
    !iopt=2
    !CALL sparse_solve(nrow,ncol,nz_sp,irow_sp(1:nz_sp),ipcol,amat_sp(1:nz_sp),       &
                        !bvec_sp,iopt)

    y=bvec_sp

    deallocate(ipcol)
    deallocate(amat_sp,irow_sp,icol_sp,bvec_sp)

    do ipoi=1,npoi
        do ieq=1,nbaleqs
            i=nbaleqs*(ipoi-1)+ieq
            params(ieq,ipoi)=y(i)
        enddo
    enddo

end subroutine evolvestep